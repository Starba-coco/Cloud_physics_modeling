subroutine adiabatic_process(z, T, p, rho, w, dt)
    use constants, only: g, cp
    implicit none
    real(8), intent(inout) :: z, T, p
    real(8), intent(in)    :: dt, w
    real(8) :: dz, dp, rho

    dz = w * dt
    z  = z + dz
    T = T - (g / cp) * dz
    call cal_rhoa(p, T, rho)
    dp = -rho * g * dz
    p  = p + dp

end subroutine adiabatic_process


subroutine cal_rhoa(p, T, rho)
    use constants, only: R_dry
    implicit none
    real(8), intent(in)  :: p, T
    real(8), intent(out) :: rho

    rho = p / (R_dry * T)                                        ! 밀도 업데이트

end subroutine cal_rhoa

subroutine cal_bin(r, r_center, log_r, log_rmin, delta_logr, &
                   r_drop, r_center_drop, log_r_drop, log_rmin_drop, delta_logr_drop)
    use constants
    implicit none
    ! Aerosol bin
    real(8), intent(out) :: r(nbin+1), r_center(nbin), log_r(nbin+1)
    real(8), intent(in)  :: log_rmin, delta_logr
    ! drop_bin
    real(8), intent(out) :: r_drop(nbin_drop+1), r_center_drop(nbin_drop), log_r_drop(nbin_drop+1)
    real(8), intent(in)  :: log_rmin_drop, delta_logr_drop
    integer              :: i
    
    ! Aerosol bin 경계 계산
    do i = 1, nbin + 1
        log_r(i) = log_rmin + (i - 1) * delta_logr
        r(i)     = exp(log_r(i))
    end do
    
    ! Aerosol bin 중심 계산 (기하 평균)
    do i = 1, nbin
        r_center(i) = sqrt(r(i) * r(i + 1))
    end do
    
    ! drop_bin 경계 계산
    do i = 1, nbin_drop + 1
        log_r_drop(i) = log_rmin_drop + (i - 1) * delta_logr_drop
        r_drop(i)     = exp(log_r_drop(i))
    end do
    
    ! drop_bin 중심 계산 (기하 평균)
    do i = 1, nbin_drop
        r_center_drop(i) = sqrt(r_drop(i) * r_drop(i + 1))
    end do
end subroutine cal_bin

subroutine lognormal(r, pdf_value)
    use constants
    implicit none
    real(8), intent(in)  :: r
    real(8), intent(out) :: pdf_value
    real(8)              :: term1, term2, pdf_i, total_pdf

    total_pdf = 0.0d0

    ! 모드 1
    term1     = - (log(r / rm_val_1))**2 / (2.0d0 * (log(sig_val_1))**2)
    term2     = r * log(sig_val_1) * sqrt(2.0d0 * pi)
    pdf_i     = (1.0d0 / term2) * exp(term1)
    pdf_i     = N0_val_1 * pdf_i
    total_pdf = total_pdf + pdf_i

    ! ! 모드 2
    ! term1     = - (log(r / rm_val_2))**2 / (2.0d0 * (log(sig_val_2))**2)
    ! term2     = r * log(sig_val_2) * sqrt(2.0d0 * pi)
    ! pdf_i     = (1.0d0 / term2) * exp(term1)
    ! pdf_i     = N0_val_2 * pdf_i
    ! total_pdf = total_pdf + pdf_i

    ! ! 모드 3
    ! term1     = - (log(r / rm_val_3))**2 / (2.0d0 * (log(sig_val_3))**2)
    ! term2     = r * log(sig_val_3) * sqrt(2.0d0 * pi)
    ! pdf_i     = (1.0d0 / term2) * exp(term1)
    ! pdf_i     = N0_val_3 * pdf_i
    ! total_pdf = total_pdf + pdf_i

    pdf_value = total_pdf

end subroutine lognormal

subroutine set_aerosol(aerosol_type, Ms, rho_s, i_vant)
    use constants
    implicit none
    character(len=*), intent(in)  :: aerosol_type
    real(8),          intent(out) :: Ms, rho_s
    integer(8),       intent(out) :: i_vant

    select case (trim(adjustl(aerosol_type)))
        case ('NaCl')
            Ms     = 58.44d-3    ! kg/mol
            rho_s  = 2160.0d0    ! kg/m³
            i_vant = 2
        case ('(NH4)2SO4')
            Ms     = 132.14d-3   ! kg/mol
            rho_s  = 1770.0d0    ! kg/m³
            i_vant = 3
        case ('KCl')
            Ms     = 74.55d-3    ! kg/mol
            rho_s  = 1980.0d0    ! kg/m³
            i_vant = 2
        ! 필요한 경우 더 많은 에어로졸 추가
        case default
            print *, 'Error: Unknown Aerosol Type "', trim(adjustl(aerosol_type)), '". Using default Aerosol (NaCl).'
            Ms     = 58.44d-3    ! kg/mol
            rho_s  = 2160.0d0    ! kg/m³
            i_vant = 2
    end select
end subroutine set_aerosol

subroutine write_aerosol_distribution(nbin, r_center, n_bin, delta_logr, time)
    implicit none
    integer, intent(in)    :: nbin
    real(8), intent(in)    :: r_center(nbin), n_bin(nbin), delta_logr, time
    integer                :: i, dist_unit
    real(8)                :: dlogD, diameter_nm, y_value
    character(len=256)     :: filename

    ! 파일 경로와 이름 생성
    write(filename, '("/home/dmafytn89/study/Cloud_physics_modeling/Week6/result/distribution_", I0, "s.txt")') int(time)

    ! 분포 데이터를 출력하기 위한 파일 설정
    dist_unit = 20
    open(dist_unit, file=filename, status='unknown', action='write')

    ! 헤더
    write(dist_unit, '(A)') 'Diameter(nm)    dN/dlogD * 1e-3 (cm^-3/nm)'

    ! dlogD 계산 (상수로 계산)
    dlogD = delta_logr / log(10.0d0)  ! ln(10)으로 나누어 log10 기반으로 변환

    ! 분포 데이터를 파일에 출력
    do i = 1, nbin
        ! 입자 지름 계산 (nm)
        diameter_nm = 2.0d0 * r_center(i) * 1.0d9

        y_value = (n_bin(i) / dlogD) * 1.0d-9  ! #/kg 단위를 #/cm^3로 변환

        write(dist_unit, '(E12.5, 3X, E12.5)') diameter_nm, y_value
    end do

    close(dist_unit)
end subroutine write_aerosol_distribution

subroutine write_drop_distribution(nbin_drop, r_center_drop, n_bin_drop, delta_logr_drop, time)
    implicit none
    integer, intent(in)    :: nbin_drop
    real(8), intent(in)    :: r_center_drop(nbin_drop), n_bin_drop(nbin_drop), delta_logr_drop, time
    integer                :: i, dist_unit_drop
    real(8)                :: dlogD, diameter_nm, y_value
    character(len=256)     :: filename

    ! 파일 경로와 이름 생성
    write(filename, '("/home/dmafytn89/study/Cloud_physics_modeling/Week6/result/drop_distribution_", I0, "s.txt")') int(time)

    ! 분포 데이터를 출력하기 위한 파일 설정
    dist_unit_drop = 30
    open(dist_unit_drop, file=filename, status='unknown', action='write')

    ! 헤더
    write(dist_unit_drop, '(A)') 'Diameter(nm)    dN/dlogD * 1e-3 (cm^-3/nm)'

    ! dlogD 계산 (상수로 계산)
    dlogD = delta_logr_drop / log(10.0d0)  ! ln(10)으로 나누어 log10 기반으로 변환

    ! 분포 데이터를 파일에 출력
    do i = 1, nbin_drop
        ! 입자 지름 계산 (nm)
        diameter_nm = 2.0d0 * r_center_drop(i) * 1.0d9

        y_value = (n_bin_drop(i) / dlogD) * 1.0d-9  ! #/kg 단위를 #/cm^3로 변환

        write(dist_unit_drop, '(E20.10, 3X, E20.10)') diameter_nm, y_value
    end do

    close(dist_unit_drop)
end subroutine write_drop_distribution

subroutine update_saturation(p, T, qv, S)
    use constants
    implicit none
    real(8), intent(in)    :: p, T, qv
    real(8), intent(out)   :: S
    real(8)                :: e, e_s, T_Celsius

    ! 수증기 압력 계산
    e         = p * qv / (epsilon + qv * (1.0d0 - epsilon))

    ! 포화 수증기 압력 계산
    T_Celsius = T - 273.15d0
    e_s       = 611.2d0 * exp((17.67d0 * T_Celsius) / (T_Celsius + 243.5d0))

    ! 상대습도 계산
    S         = (e / e_s) - 1.0d0                                      

end subroutine update_saturation

subroutine activation(p, T, S, nbin, nbin_drop, r, n_bin, n_bin_drop, &
                      Ms, rho_s, i_vant)
    use constants, only: sigma_v, Rv, rho_w, Mw, Lv, cp, pi
    implicit none

    ! 변수 선언
    integer, intent(in)    :: nbin, nbin_drop, i_vant
    real(8), intent(inout) :: p, T, S
    real(8), intent(in)    :: r(nbin+1)
    real(8), intent(inout) :: n_bin(nbin), n_bin_drop(nbin_drop)
    real(8), intent(in)    :: Ms, rho_s

    ! 로컬 변수
    real(8) :: a, b, ln_r0, ln_r1, ln_r2, activate_ratio, n_activated
    integer :: i

    ! 콜러 상수 계산 (온도 의존)
    a = (2.0d0 * sigma_v)      / (Rv * rho_w * T)
    b = i_vant * ((rho_s * Mw) / (Ms * rho_w))

    ln_r0 = log((a / 3.0d0) * ((4.0d0 / (b * (S**2))) ** (1.0d0 / 3.0d0)))

    ! 각 bin에 대해 임계 반경과 임계 과포화도 계산
    do i = 1, nbin
        ln_r1 = log(r(i))
        ln_r2 = log(r(i+1))
    
        ! 활성화 여부 판단
        if (S .ge. 0.0d0) then
            if (ln_r0 > ln_r2) then  ! 활성화되지 않음
                activate_ratio = 0.0d0
                cycle
            else if (ln_r0 > ln_r1) then ! 일부만 활성화
                activate_ratio = (ln_r2 - ln_r0) / (ln_r2 - ln_r1)
            else ! ln_r0 ≤ ln_r1 --> 모두 활성화
                activate_ratio = 1.0d0
            end if
    
            ! 활성화된 입자의 수 계산
            n_activated = min(n_bin(i), n_bin(i) * activate_ratio)

            if (n_activated <= 0.0d0) then
                exit
            end if

            ! 활성화된 입자를 drop bin으로 이동
            n_bin_drop(1) = n_bin_drop(1) + n_activated
            n_bin(i)      = n_bin(i)      - n_activated
            
        end if
    end do
end subroutine activation

subroutine condensation(T, p, S, dt, m, m_new, rho, r_new, n_bin_drop, r_center_drop, qv)
    use constants, only: rho_w, Rv, Lv, nbin_drop, pi, cp
    implicit none

    real(8), intent(in)    :: p, dt, rho
    real(8), intent(in)    :: n_bin_drop(nbin_drop), r_center_drop(nbin_drop)
    real(8), intent(inout) :: T, S, qv
    real(8), intent(inout) :: m(nbin_drop), m_new(nbin_drop), r_new(nbin_drop)
    real(8)                :: T_Celsius, e_s, Ka, Dv, Fd, Fk, dqc
    real(8), allocatable   :: dm(:)
    integer                :: ibin

    ! 배열 할당
    allocate(dm(nbin_drop))

    ! 상수 값 계산
    ! dqc       = 0.0d0
    T_Celsius = T - 273.15d0
    e_s       = 611.2d0 * exp((17.67d0 * T_Celsius) / (T_Celsius + 243.5d0))

    Ka = (4.1868d0 * 1.0d-3) * (5.69d0 + 0.017d0 * T_Celsius)                ! W/m·K
    Dv = (2.1100d0 * 1.0d-5) * ((T / 273.15d0) ** 1.94d0) * (101325.0d0 / p) ! m^2/s
    ! print '(A, F10.7)', "Ka :", Ka
    ! print '(A, F10.7)', "Dv :", Dv

    Fk = ((Lv / (Rv * T)) - 1.0d0) * (Lv / (Ka * T))
    Fd = (Rv * T) / (Dv * e_s)
    
    do ibin = 1, nbin_drop
        m(ibin)     = rho_w * (4.0d0 / 3.0d0) * pi * r_center_drop(ibin)**3
        r_new(ibin) = r_center_drop(ibin) + (S / (r_center_drop(ibin) * rho_w * (Fd + Fk))) * dt
        m_new(ibin) = rho_w * (4.0d0 / 3.0d0) * pi * r_new(ibin)**3
        dm(ibin)    = m_new(ibin) - m(ibin)
    end do

    ! 전체 수증기 질량 변화량 계산
    dqc = dot_product(n_bin_drop, dm)

    ! 수증기 혼합비 업데이트 (kg/kg)
    qv = qv - dqc

    ! 온도 업데이트 (K)
    T = T + (Lv * dqc) / (cp * rho)

    ! 포화도 업데이트
    call update_saturation(p, T, qv, S)

    ! 배열 메모리 해제
    deallocate(dm)
end subroutine condensation

subroutine redistribution(m, m_new, n_bin_drop, nbin_drop, r_new, r_drop, n_bin_new)
    implicit none

    ! 변수 선언
    integer                :: i, j          
    integer, intent(in)    :: nbin_drop 
    real(8), intent(in)    :: m(nbin_drop), m_new(nbin_drop), r_new(nbin_drop), r_drop(nbin_drop+1)
    real(8), intent(inout) :: n_bin_drop(nbin_drop)
    real(8), intent(out)   :: n_bin_new(nbin_drop)          
    real(8)                :: position        

    n_bin_new = 0.0d0
    
    ! Redistribution
    do i = 1, nbin_drop
        do j = 1, nbin_drop
            if ((r_new(i) >= r_drop(j)) .and. (r_new(i) < r_drop(j + 1))) then

                position  = (m_new(i) - m(j)) / (m(j+1) - m(j))

                ! n_bin_drop 값 업데이트
                n_bin_new(j)     = n_bin_new(j)     + n_bin_drop(i) * (1.0d0 - position)
                n_bin_new(j + 1) = n_bin_new(j + 1) + n_bin_drop(i) * position

                exit  ! 해당 bin을 찾았으므로 루프 종료
            end if
        end do
    end do

    n_bin_drop = n_bin_new

end subroutine redistribution

subroutine terminal_velocity(r_center_drop, rho, T, P, Vt)
    use constants, only: R_dry, g, rho_w, nbin_drop
    implicit none

    ! 입력 변수
    real(8), intent(in)    :: r_center_drop(nbin_drop)  ! 빗방울의 반지름 [m]
    real(8), intent(in)    :: T                    ! 온도 [K]
    real(8), intent(in)    :: P                    ! 압력 [Pa]
    real(8), intent(in)    :: rho                  ! 공기의 밀도
    real(8), intent(out)   :: Vt(nbin_drop)        ! 종단속도 [cm/s]

    ! 지역 변수
    real(8), allocatable  :: d0(:)
    real(8) :: l, C1, Csc, eta0
    real(8) :: b0, b1, b2, b3, b4, b5, b6
    real(8) :: C2, Da, X, Y, Re, Cl, C3, sigma, Bo, Np
    integer :: i
    integer :: file_unit

    allocate(d0(nbin_drop))

    file_unit = 40
    open(unit=file_unit, file='vt_output.txt', status='unknown', action='write')
    write(file_unit, '(A)') 'Radius(m)        Vt(cm/s)'

    eta0 = 1.818d-5 
    l    = 6.620d-8 * 0.018d0 * (101325d0 / P) * sqrt(T / 293.15d0)
    C1   = (rho_w - rho) * g / (18.0d0 * eta0)
    C2   = 4.0d0 * rho * (rho_w - rho) * g / (3.0d0 * (eta0 ** 2.0d0))

    do i = 1, nbin_drop

        d0(i) = r_center_drop(i)
        ! print *, d0(i)
        Csc = 1.0d0 + (2.51d0 * l / d0(i))

        ! 0.5 um <= d0 < 19 um
        if (0.5d-6 <= d0(i) .and. d0(i) < 19.0d-6) then
            Vt(i) = C1 * Csc * (d0(i)**2)

        ! 19 um <= d0 < 1070 um
        else if (d0(i) < 1.07d-3) then
            b0 = -0.318657d1
            b1 =  0.992696d0
            b2 = -0.153193d-2
            b3 = -0.987059d-3
            b4 = -0.578878d-3
            b5 =  0.855176d-4
            b6 = -0.327815d-5

            Da    = C2 * (d0(i)**3)
            X     = log(Da)
            Y     = b0 + b1 * X + b2 * X**2 + b3 * X**3 + b4 * X**4 + b5 * X**5 + b6 * X**6
            Re    = Csc * exp(Y)
            Vt(i) = eta0 * Re / (rho * d0(i))

        ! 1070 um <= d0 < 7 mm
        else if (d0(i) < 7.0d-3) then
            b0 = -0.500015d1
            b1 =  0.523778d1
            b2 = -0.204914d1
            b3 =  0.475294d0
            b4 = -0.542819d-1
            b5 =  0.238449d-2

            Cl    = -1.55d-4
            sigma = Cl * T + 0.118d0
            C3    = 4.0d0 * (rho_w - rho) * g / (3.0d0 * sigma)
            Bo    = C3 * d0(i)**2.0d0
            Np    = ((sigma**3.0d0) * (rho**2.0d0)) / ((eta0**4.0d0) * (rho_w - rho) * g)
            X     = log(Bo * (Np**(1.0d0 / 6.0d0)))
            Y     = b0 + b1 * X + b2 * X**2 + b3 * X**3 + b4 * X**4 + b5 * X**5
            Re    = Np**(1.0d0 / 6.0d0) * exp(Y)
            Vt(i) = eta0 * Re / (rho * d0(i))
        else
            Vt(i) = 9.12499177288172d0
        end if

        print *, Vt(i)

        write(file_unit, '(E15.7, 3X, E15.7)') d0(i), Vt(i)
    end do 

    close(file_unit)

    deallocate(d0)
end subroutine terminal_velocity

subroutine collision(r_center_drop, n_bin_drop, dt...)
    use constants, only: nbin_drop
    implicit none

    real(8), intent(in) :: dt
    real(8), intent(in) :: r_center_drop(:), n_bin_drop(:)
    do i = 1, nbin_drop-1
        do j = i+1, nbin_drop
            dN = k(i, j) * n_bin_drop(i) * n_bin_drop(j) * dt
            dN = min(n(i), dN)
            dN = min(n(j), dN)
            n(i) = n(i) - dN
            n(j) = n(j) - dN
            call redistribution(m, m_new, n_bin_drop, nbin_drop, r_new, r_drop, n_bin_new)
        end do
    end do


end subroutine collision
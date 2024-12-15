subroutine adiabatic_process(z, T, p, rho, w, dt)
    use constants, only: g, cp
    implicit none
    real(8), intent(inout) :: z, T, p, rho
    real(8), intent(in)    :: dt, w
    real(8) :: dz, dp

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

    rho = p / (R_dry * T)  ! 밀도 업데이트

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
        ! write(*,*) r_center_drop(i)
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

subroutine activation(T, S, r, n_bin, n_bin_drop, Ms, rho_s, i_vant)
    use constants, only: sigma_v, Rv, rho_w, Mw, Lv, cp, pi, nbin, nbin_drop
    implicit none

    ! 변수 선언
    integer, intent(in)    :: i_vant
    real(8), intent(in)    :: T, S, Ms, rho_s
    real(8), intent(in)    :: r(nbin+1)
    real(8), intent(inout) :: n_bin(nbin), n_bin_drop(nbin_drop)

    ! 로컬 변수
    real(8) :: a, b, ln_r0, ln_r1, ln_r2, activate_ratio, n_activated
    integer :: i

    ! kohler 상수 계산
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
    real(8), intent(out)   :: m(nbin_drop), m_new(nbin_drop), r_new(nbin_drop)
    real(8)                :: T_Celsius, e_s, Ka, Dv, Fd, Fk, dqc
    real(8), allocatable   :: dm(:)
    integer                :: i

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
    
    do i = 1, nbin_drop
        m(i)     = rho_w * (4.0d0 / 3.0d0) * pi * r_center_drop(i)**3
        r_new(i) = r_center_drop(i) + (S / (r_center_drop(i) * rho_w * (Fd + Fk))) * dt
        m_new(i) = rho_w * (4.0d0 / 3.0d0) * pi * r_new(i)**3
        dm(i)    = m_new(i) - m(i)
    end do
    ! print *, dm
    ! print *, "=================================="

    ! 전체 수증기 질량 변화량 계산
    dqc = dot_product(n_bin_drop, dm)
    ! print *, dqc
    ! 수증기 혼합비 업데이트 (kg/kg)
    qv = qv - dqc

    ! 온도 업데이트 (K)
    T = T + (Lv * dqc) / (cp * rho)

    ! =====call rhoa 가 들어가야 하는지???============
    
    ! 포화도 업데이트
    call update_saturation(p, T, qv, S)

    ! 배열 메모리 해제
    deallocate(dm)
end subroutine condensation

subroutine redistribution(m, m_new, n_bin_drop, r_new, r_drop)
    use constants, only: nbin_drop
    implicit none

    ! 변수 선언
    integer                :: i, j
    real(8), intent(in)    :: m(nbin_drop), m_new(nbin_drop), r_new(nbin_drop), r_drop(nbin_drop+1)
    real(8), intent(inout) :: n_bin_drop(nbin_drop)
    real(8)                :: position, n_bin_new(nbin_drop)      

    n_bin_new = 0.0d0
    
    ! Redistribution
    do i = 1, nbin_drop
        do j = 1, nbin_drop
            ! if ((r_new(i) >= r_drop(j)) .and. (r_new(i) < r_drop(j + 1))) then
            if (((m_new(i) >= m(j))) .and. (m_new(i) < m(j + 1))) then

                position  = (m_new(i) - m(j)) / (m(j + 1) - m(j))
                ! n_bin_drop 값 업데이트
                n_bin_new(j)     = n_bin_new(j)     + n_bin_drop(i) * (1.0d0 - position)
                n_bin_new(j + 1) = n_bin_new(j + 1) + n_bin_drop(i) * position

                exit  ! 해당 bin을 찾았으므로 루프 종료
            end if
        end do
    end do

    n_bin_drop = n_bin_new

end subroutine redistribution

subroutine redistribution_for_collision(mass_new, m_drop, dN, dn_drop)
    use constants, only: nbin_drop
    implicit none    

    real(8), intent(in)    :: mass_new     
    real(8), intent(in)    :: m_drop(nbin_drop) 
    real(8), intent(in)    :: dN
    real(8), intent(inout) :: dn_drop(nbin_drop)

    integer :: i
    real(8) :: position
    real(8), dimension(nbin_drop) :: dn_new

    ! if (dN == 0.0d0) return
    ! if (mass_new == 0.0d0) return

    dn_new = 0.0d0

    if (mass_new < m_drop(1)) then
        dn_new(1) = dn_new(1) + dN

    else if (mass_new >= m_drop(nbin_drop)) then
        dn_new(nbin_drop) = dn_new(nbin_drop) + dN

    else
        do i = 1, nbin_drop - 1
            if (mass_new >= m_drop(i) .and. mass_new < m_drop(i+1)) then
                position    = (m_drop(i+1) - mass_new) / (m_drop(i+1) - m_drop(i))
                dn_new(i)   = dn_new(i)   + dN * position
                dn_new(i+1) = dn_new(i+1) + dN * (1.0d0 - position)
                exit
            end if
        end do
    end if

    ! 새로 형성된 입자 수를 기존 분포에 반영
    dn_drop = dn_new
    ! write(*,*) 'redistribution good'
end subroutine redistribution_for_collision

subroutine terminal_velocity(r_center_drop, rho, T, p, Vt)
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

        write(file_unit, '(E15.7, 3X, E15.7)') d0(i), Vt(i)
    end do 

    close(file_unit)

    deallocate(d0)
end subroutine terminal_velocity

subroutine effic(r_center_drop, ec)
    ! collision efficiencies of Hall kernel
    use constants, only: nbin_drop
    implicit none
    real(8), intent(in),  dimension(nbin_drop) :: r_center_drop
    real(8), intent(out), dimension(nbin_drop, nbin_drop) :: ec
    real(8), dimension(11, 21) :: ecoll
    real(8) :: rat(21), r0(11)
    integer :: i, j, k, ir, iq
    real(8) :: rq, p, q
    real(8) :: r_j_um, r_i_um

    ! Initialize r0 and rat
    r0  = (/300.0d0, 200.0d0, 150.0d0, 100.0d0, 70.0d0, 60.0d0, 50.0d0, 40.0d0, 30.0d0, 20.0d0, 10.0d0/)
    rat = (/0.0d0, 0.05d0, 0.1d0, 0.15d0, 0.2d0, 0.25d0, 0.3d0, 0.35d0, 0.4d0, 0.45d0, 0.5d0, 0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0, 0.9d0, 0.95d0, 1.0d0/)

    ! Initialize ecoll
    ecoll = reshape([ &
    0.0d0, 0.97d0,   1.0d0,    1.0d0,    1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   &
    0.0d0, 0.87d0,   0.96d0,   0.98d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   &
    0.0d0, 0.77d0,   0.93d0,   0.97d0,   0.97d0,  1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   &
    0.0d0, 0.5d0,    0.79d0,   0.91d0,   0.95d0,  0.95d0,  1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   &
    0.0d0, 0.2d0,    0.58d0,   0.75d0,   0.84d0,  0.88d0,  0.9d0,   0.92d0,  0.94d0,  0.95d0,  0.95d0,  0.95d0,  0.95d0,  0.95d0,  0.95d0,  0.97d0,  1.0d0,   1.02d0,  1.04d0,  2.3d0,   4.0d0,   &
    0.0d0, 0.05d0,   0.43d0,   0.64d0,   0.77d0,  0.84d0,  0.87d0,  0.89d0,  0.91d0,  0.91d0,  0.91d0,  0.91d0,  0.91d0,  0.91d0,  0.92d0,  0.93d0,  0.95d0,  1.0d0,   1.03d0,  1.7d0,   3.0d0,   &
    0.0d0, 0.005d0,  0.4d0,    0.6d0,    0.7d0,   0.78d0,  0.83d0,  0.86d0,  0.88d0,  0.9d0,   0.9d0,   0.9d0,   0.9d0,   0.89d0,  0.88d0,  0.88d0,  0.89d0,  0.92d0,  1.01d0,  1.3d0,   2.3d0,   &
    0.0d0, 0.001d0,  0.07d0,   0.28d0,   0.5d0,   0.62d0,  0.68d0,  0.74d0,  0.78d0,  0.8d0,   0.8d0,   0.8d0,   0.78d0,  0.77d0,  0.76d0,  0.77d0,  0.77d0,  0.78d0,  0.79d0,  0.95d0,  1.4d0,   &
    0.0d0, 0.0001d0, 0.002d0,  0.02d0,   0.04d0,  0.085d0, 0.17d0,  0.27d0,  0.4d0,   0.5d0,   0.55d0,  0.58d0,  0.59d0,  0.58d0,  0.54d0,  0.51d0,  0.49d0,  0.47d0,  0.45d0,  0.47d0,  0.52d0,  &
    0.0d0, 0.0001d0, 0.0001d0, 0.005d0,  0.016d0, 0.022d0, 0.03d0,  0.043d0, 0.052d0, 0.064d0, 0.072d0, 0.079d0, 0.082d0, 0.08d0,  0.076d0, 0.067d0, 0.057d0, 0.048d0, 0.04d0,  0.033d0, 0.027d0, &
    0.0d0, 0.0001d0, 0.0001d0, 0.0001d0, 0.014d0, 0.017d0, 0.019d0, 0.022d0, 0.027d0, 0.03d0,  0.033d0, 0.035d0, 0.037d0, 0.038d0, 0.038d0, 0.037d0, 0.036d0, 0.035d0, 0.032d0, 0.029d0, 0.027d0], shape=[11, 21])

! Compute collision efficiency
    do j = 1, nbin_drop
        r_j_um = r_center_drop(j) * 1.0d6  ! m -> µm

        do i = 1, nbin_drop
            ! m 단위의 r_center_drop를 µm로 변환    
            r_i_um = r_center_drop(i) * 1.0d6  ! m -> µm
            ! write(*,*) r_j_um, r_i_um
            ! Check if radius exceeds 300 µm
            if (r_j_um > 300.0d0 .or. r_i_um > 300.0d0) then
                ec(j, i) = 1.0
                cycle
            end if
    
            ! Find indices for r(j_um)
            ! r0 배열은 11개 구간점 (300, 200, 150, ... µm)
            ! r_j_um이 어느 r0(k) 사이에 있는지 찾아 ir 결정
            ir = 1
            do k = 2, 11
                if (r_j_um <= r0(k)) then
                    ir = k
                    exit
                end if
            end do
    
            ! Find indices for rq
            ! rq = r_i_um / r_j_um (반경비)
            ! rat 배열은 21개 구간점 (0.0, 0.05, 0.1, ..., 1.0)
            ! rq가 어느 rat(k) 사이에 있는지 찾아 iq 결정
            rq = r_i_um / r_j_um
            iq = 1
            do k = 2, 21
                if (rq <= rat(k)) then
                    iq = k
                    exit
                end if
            end do
    
            ! Interpolation
            ! ir과 iq 구간 내라면 bilinear interpolation 수행
            if (ir < 11 .and. iq < 21) then
                ! p: r_j_um이 r0(ir-1)와 r0(ir) 사이에서의 상대 위치
                p = (r_j_um - r0(ir-1)) / (r0(ir) - r0(ir-1))
                ! q: rq가 rat(iq-1)와 rat(iq) 사이에서의 상대 위치
                q = (rq - rat(iq-1)) / (rat(iq) - rat(iq-1))
                ! write(*,*) p, q
                ! ecoll(ir,iq)는 표준 효율 테이블
                ! (1-p)*(1-q), p*(1-q), q*(1-p), p*q를 가중치로 하여 2차원 보간
                ec(j, i) = (1-p)*(1-q)*ecoll(ir-1, iq-1) + p*(1-q)*ecoll(ir, iq-1) + &
                           q*(1-p)*ecoll(ir-1, iq)       + p*q*ecoll(ir, iq)
            else
                ! 구간 밖이면 효율 0
                ec(j, i) = 0.0
            end if
        end do
    end do
    ! print *, ec
end subroutine effic

subroutine collision_kernel(r_center_drop, Vt, ec, k)
    use constants, only: nbin_drop, pi
    implicit none
    
    real(8), dimension(nbin_drop),            intent(in) :: Vt
    real(8), dimension(nbin_drop),            intent(in) :: r_center_drop
    real(8), dimension(nbin_drop, nbin_drop), intent(in) :: ec
    real(8), dimension(nbin_drop, nbin_drop), intent(out):: k
    integer :: i, j
    
    k = 0.0d0

    do i = 1, nbin_drop - 1
        do j = i + 1, nbin_drop
            k(i,j) = pi * (r_center_drop(i) + r_center_drop(j))**2 * abs(Vt(i) - Vt(j))
            ! k(i,j) = k(i,j) * ec(i,j)
        end do
    end do

end subroutine collision_kernel

subroutine collision(dt, rho, r_center_drop, n_bin_drop, k, total_mass)
    use constants, only: nbin_drop, pi, rho_w
    implicit none

    real(8), intent(in) :: dt, rho
    real(8), intent(in),    dimension(nbin_drop)            :: r_center_drop
    real(8), intent(in),    dimension(nbin_drop, nbin_drop) :: k
    real(8), intent(inout), dimension(nbin_drop)            :: n_bin_drop
    real(8), intent(out),   dimension(nbin_drop)            :: total_mass
    real(8), dimension(nbin_drop) :: dn_new

    real(8), dimension(nbin_drop) :: m_drop
    real(8) :: dN, mm_new, total_mass_before, total_mass_after
    integer :: i, j
    integer :: unit_coll

    ! #/kg -> #/m3 변환
    n_bin_drop = n_bin_drop * rho
    total_mass = 0.0d0

    do i = 1, nbin_drop
        m_drop(i) = rho_w * (4.0d0 / 3.0d0) * pi * r_center_drop(i)**3
    end do 
    ! if (sum(n_bin_drop) == 0.0d0) return

    do i = 1, nbin_drop
        if (n_bin_drop(i) == 0.0) cycle

        do j = i, nbin_drop
            if (n_bin_drop(j) == 0.0) cycle
            ! #/kg -> #/m3 변환
            ! n_bin_drop = n_bin_drop * rho
            ! mm_new     = 0.0d0 
            dn_new = 0.0d0
            dN = k(i, j) * n_bin_drop(i) * n_bin_drop(j) * dt
            dN = min(n_bin_drop(i), dN)
            dN = min(n_bin_drop(j), dN)
            ! print *, dN
            ! if (dN == 0) cycle

            ! 충돌로 i,j bin에서 dN만큼 감소
            n_bin_drop(i) = n_bin_drop(i) - dN
            n_bin_drop(j) = n_bin_drop(j) - dN
        
            ! 충돌해서 만들어진 새 물방울 질량 mm_new
            mm_new = m_drop(i) + m_drop(j)

            ! write(*,*) mm_new
            ! 새로 형성된 물방울을 적절한 bin에 재분배
            call redistribution_for_collision(mm_new, m_drop, dN, dn_new)
            n_bin_drop = n_bin_drop + dn_new
            ! write(*,*) n_bin_drop
        end do
    end do

    n_bin_drop = n_bin_drop / rho

    total_mass = n_bin_drop * m_drop
    ! print *, total_mass
    print *, sum(n_bin_drop)
end subroutine collision    

subroutine adiabatic_process(z, T, p, w, dt)
    use constants, only: g, cp
    implicit none
    real(8), intent(inout) :: z, T, p
    real(8), intent(in)    :: w, dt
    real(8) :: dz, dp, rho

    dz = w * dt
    z  = z + dz

    call cal_rhoa(p, T, rho)

    T         = T - (g / cp) * dz   ! 온도 업데이트
    dp        = -rho * g * dz       ! 압력 업데이트
    p         = p + dp

end subroutine adiabatic_process

subroutine cal_rhoa(p, T, rho)
    use constants, only: R_dry
    implicit none
    real(8), intent(in)  :: p, T
    real(8), intent(out) :: rho

    rho = p / (R_dry * T)           ! 밀도 업데이트

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
    
    ! Aerosol bin 중심 계산
    do i = 1, nbin
        r_center(i) = sqrt(r(i) * r(i + 1))
    end do
    
    ! drop_bin 경계 계산
    do i = 1, nbin_drop + 1
        log_r_drop(i) = log_rmin_drop + (i - 1) * delta_logr_drop
        r_drop(i)     = exp(log_r_drop(i))
    end do
    
    ! drop_bin 중심 계산
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
    integer              :: i

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

subroutine write_distribution(nbin, r_center, n_bin, delta_logr)
    implicit none
    integer, intent(in)    :: nbin
    real(8), intent(in)    :: r_center(nbin), n_bin(nbin), delta_logr
    integer                :: i, dist_unit
    real(8)                :: dlogD, diameter_nm, y_value

    ! 분포 데이터를 출력하기 위한 파일 설정
    dist_unit = 20
    open(dist_unit, file='distribution.txt', status='unknown', action='write')

    ! 헤더
    write(dist_unit, '(A)') 'Diameter(nm)    dN/dlogD * 1e-3 (cm^-3/nm)'

    ! dlogD 계산
    dlogD = delta_logr / log(10.0d0)  ! ln(10)으로 나누어 log10 기반으로 변환

    ! 분포 데이터를 파일에 출력
    do i = 1, nbin
        ! 입자 지름 계산 (nm)
        diameter_nm = 2.0d0 * r_center(i) * 1.0d9

        y_value = (n_bin(i) / dlogD) * 1.0d-9

        write(dist_unit, '(E12.5, 3X, E12.5)') diameter_nm, y_value
    end do

    close(dist_unit)
end subroutine write_distribution

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

subroutine activation(p, T, qv, S, nbin, nbin_drop, r, n_bin, n_bin_drop, &
                      r_center, r_center_drop, Ms, rho_s, i_vant)
    use constants, only: sigma_v, Rv, rho_w, Mw, Lv, cp, pi
    implicit none

    ! 변수 선언
    integer, intent(in)    :: nbin, nbin_drop, i_vant
    real(8), intent(inout) :: p, T, qv, S
    real(8), intent(in)    :: r(nbin+1), r_center(nbin), r_center_drop(nbin_drop)
    real(8), intent(inout) :: n_bin(nbin), n_bin_drop(nbin_drop)
    real(8), intent(in)    :: Ms, rho_s

    ! 로컬 변수
    real(8) :: a, b, ln_r0, ln_r1, ln_r2, activate_ratio, n_activated, dqc
    integer :: i

    ! 콜러 상수 계산 (온도 의존)
    a = (2.0d0 * sigma_v)      / (Rv * rho_w * T)
    b = i_vant * ((rho_s * Mw) / (Ms * rho_w))

    ln_r0 = log( (a / 3.0d0) * ( (4.0d0 / (b * (S**2))) ** (1.0d0 / 3.0d0) ) )

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

            ! 활성화된 입자의 수
            n_activated   = n_bin(i) * activate_ratio
            n_bin_drop(1) = n_bin_drop(1) + n_activated
            n_bin(i)      = n_bin(i)      - n_activated

            dqc = n_activated * rho_w * (4.0d0 / 3.0d0) * pi * r_center_drop(1)**3
            qv  = qv - dqc
            T   = T + (Lv * dqc / cp)

            call update_saturation(p, T, qv, S)

        end if
    end do

end subroutine activation

subroutine update_time_and_w(time, dt, i_time, w_val, w, end_flag)
    implicit none
    ! 변수 선언
    real(8), intent(inout) :: time
    real(8), intent(in)    :: dt
    real(8), intent(in)    :: i_time(2), w_val(2)
    real(8), intent(out)   :: w
    logical, intent(out)   :: end_flag

    ! 시간 업데이트
    time = time + dt

    ! 초기화
    w = 0.0d0

    ! 상승 속도 업데이트
    if (time <= i_time(1)) then
        w = w_val(1)
    else if (time <= i_time(2)) then
        w = w_val(2)
    else
        end_flag = .true.
    end if

end subroutine update_time_and_w

subroutine condensation(T, p, S, dt, m, r, n_particles, rho, delta_qv, delta_T)
    use constants, only: pi, Rv, Lv, rho_w, cp
    implicit none

    ! 입력 변수
    real(8), intent(in)    :: T           ! 온도 (K)
    real(8), intent(in)    :: p           ! 압력 (Pa)
    real(8), intent(in)    :: S           ! 과포화도
    real(8), intent(in)    :: dt          ! 시간 간격 (s)
    real(8), intent(in)    :: n_particles ! 해당 bin의 물방울 수 (#/kg)
    real(8), intent(in)    :: rho         ! 공기 밀도 (kg/m^3)

    ! 입출력 변수
    real(8), intent(inout) :: m           ! 물방울 질량 (kg)
    real(8), intent(inout) :: r           ! 물방울 반경 (m)
    real(8), intent(inout) :: delta_qv    ! 수증기 혼합비 변화량 (kg/kg)
    real(8), intent(inout) :: delta_T     ! 온도 변화량 (K)

    ! 지역 변수
    real(8) :: T_Celsius  ! 섭씨 온도 (°C)
    real(8) :: e_s        ! 포화 수증기 압력 (Pa)
    real(8) :: Ka         ! 열전도도 (W/m·K)
    real(8) :: Dv         ! 수증기 확산 계수 (m^2/s)
    real(8) :: Fd         ! 확산 제한 인자
    real(8) :: Fk         ! 열전도 제한 인자
    real(8) :: dm         ! 개별 물방울의 질량 변화량 (kg)
    real(8) :: delta_m    ! 해당 bin에서의 총 질량 변화량 (kg)

    ! 포화 수증기 압력 계산
    T_Celsius = T - 273.15d0
    e_s       = 611.2d0 * exp( (17.67d0 * T_Celsius) / (T_Celsius + 243.5d0) )

    ! 열전도도 Ka와 확산 계수 Dv 계산
    Ka = (4.1868d0 * 1.0d-3) * (5.69d0 + 0.017d0 * T_Celsius)    ! W/m·K
    Dv = (2.1100d0 * 1.0d-5) * ( (T / 273.15d0) ** 1.94d0 ) * (101325.0d0 / p)   ! m^2/s

    ! 확산 및 열전도 제한 인자 계산
    Fd = (rho_w * Rv * T) / (e_s * Dv)
    Fk = ((Lv * rho_w) / (Ka * T)) - ((Rv * T) / Lv) + 1.0d0

    ! 개별 물방울의 질량 변화량 계산
    dm = ((4.0d0 * pi * r * S) / (Fd + Fk)) * dt

    ! 물방울 질량과 반경 업데이트
    m = m + dm
    r = ( (3.0d0 * m) / (4.0d0 * pi * rho_w) ) ** (1.0d0 / 3.0d0)

    ! 해당 bin에서의 총 질량 변화량 계산
    delta_m = n_particles * dm

    ! 수증기 혼합비와 온도 변화량 누적
    delta_qv = delta_qv - delta_m / rho
    delta_T  = delta_T + (Lv * delta_m) / cp

end subroutine condensation

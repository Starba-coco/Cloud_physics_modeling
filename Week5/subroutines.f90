subroutine adiabatic_process(z, T, p, rho, w, dt, q, S)
    use constants
    implicit none
    real(8), intent(inout) :: z, T, p, rho, q
    real(8), intent(in)    :: w, dt
    real(8), intent(out)   :: S
    real(8) :: dz, dp, e_s, e, T_Celsius

    dz = w * dt
    z  = z + dz

    T         = T - (g / cp) * dz                                      ! 온도 업데이트
    dp        = -rho * g * dz                                          ! 압력 업데이트
    p         = p + dp
    rho       = p / (R_dry * T)                                        ! 밀도 업데이트
    e         = p * q / (epsilon + q * (1.0d0 - epsilon))              ! 수증기 압력 계산
    T_Celsius = T - 273.15d0
    e_s       = 6.112d0 * exp((17.67d0 * T_Celsius) / (T_Celsius + 243.5d0)) ! 포화 수증기 압력 계산
    e_s       = e_s * 100.0d0                                          ! hPa -> Pa 변환
    S         = (e / e_s) - 1.0d0                                      ! 상대습도 계산

end subroutine adiabatic_process

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
    integer              :: i

    total_pdf = 0.0d0

    ! 모드 1
    term1     = - (log(r / rm_val_1))**2 / (2.0d0 * (log(sig_val_1))**2)
    term2     = r * log(sig_val_1) * sqrt(2.0d0 * pi)
    pdf_i     = (1.0d0 / term2) * exp(term1)
    pdf_i     = N0_val_1 * pdf_i
    total_pdf = total_pdf + pdf_i

    ! 모드 2
    term1     = - (log(r / rm_val_2))**2 / (2.0d0 * (log(sig_val_2))**2)
    term2     = r * log(sig_val_2) * sqrt(2.0d0 * pi)
    pdf_i     = (1.0d0 / term2) * exp(term1)
    pdf_i     = N0_val_2 * pdf_i
    total_pdf = total_pdf + pdf_i

    ! 모드 3
    term1     = - (log(r / rm_val_3))**2 / (2.0d0 * (log(sig_val_3))**2)
    term2     = r * log(sig_val_3) * sqrt(2.0d0 * pi)
    pdf_i     = (1.0d0 / term2) * exp(term1)
    pdf_i     = N0_val_3 * pdf_i
    total_pdf = total_pdf + pdf_i

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

    ! dlogD 계산 (상수로 계산)
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
    real(8)                :: e_s, qs

    ! 포화 수증기 압력 계산
    e_s = 610.94d0 * exp(17.625d0 * (T - 273.15d0) / (T - 30.11d0))

    ! 포화 혼합비 계산
    qs = epsilon * e_s / (p - e_s)

    ! 과포화도 계산
    S = qv / qs - 1.0d0
end subroutine update_saturation
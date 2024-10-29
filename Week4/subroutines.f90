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
    rho       = p / (R_gas * T)                                        ! 밀도 업데이트
    e         = p * q / (epsilon + q * (1.0d0 - epsilon))              ! 수증기 압력 계산
    T_Celsius = T - 273.15d0
    e_s       = 6.112d0 * exp((17.67d0 * T_Celsius) / (T_Celsius + 243.5d0)) ! 포화 수증기 압력 계산
    e_s       = e_s * 100.0d0                                          ! hPa -> Pa 변환
    S         = (e / e_s) - 1.0d0                                      ! 상대습도 계산

end subroutine adiabatic_process

subroutine lognormal(r, pdf_value)
    use constants
    implicit none
    real(8), intent(in)  :: r
    real(8), intent(out) :: pdf_value
    real(8)              :: term1, term2

    term1     = - (log(r / rm))**2 / (2.0d0 * (sigma**2))
    term2     = (r * sigma * sqrt(2.0d0 * pi))
    pdf_value = (1.0d0 / term2) * exp(term1)  ! 로그 정규 분포 계산

end subroutine lognormal

subroutine set_aerosol(aerosol_type, Ms, rho_s, i_vant)
    use constants
    implicit none
    character(len=*), intent(in)  :: aerosol_type
    real(8),          intent(out) :: Ms, rho_s
    integer(8),       intent(out) :: i_vant

    select case (trim(aerosol_type))
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
            print *, 'Unknown Aerosol. Proceeding with default Aerosol(NaCl)'
            Ms     = 58.44d-3    ! kg/mol
            rho_s  = 2160.0d0    ! kg/m³
            i_vant = 2
    end select
end subroutine set_aerosol
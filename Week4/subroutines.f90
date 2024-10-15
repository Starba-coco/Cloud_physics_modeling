subroutine adiabatic_process(z, T, p, rho, w, dt, q, RH, S)
    use constants
    implicit none
    real, intent(inout) :: z, T, p, rho, q
    real, intent(in)    :: w, dt
    real, intent(out)   :: RH, S
    real :: dz, dp, e_s, e, T_Celsius

    dz = w * dt
    z  = z + dz

    ! 온도 업데이트
    T = T - (g / cp) * dz

    ! 압력 업데이트
    dp = -rho * g * dz
    p  = p + dp

    ! 밀도 업데이트
    rho = p / (R_gas * T)

    ! 수증기 압력 계산
    e = p * q / (epsilon + q * (1.0 - epsilon))

    ! 포화 수증기 압력 계산 (Magnus 공식)
    T_Celsius = T - 273.15
    e_s = 6.112 * exp((17.67 * T_Celsius) / (T_Celsius + 243.5))
    e_s = e_s * 100.0   ! hPa -> Pa 변환

    ! 상대습도 계산
    RH = (e / e_s) * 100

    ! 과포화도 계산
    S = RH - 100

end subroutine adiabatic_process

subroutine lognormal(r, pdf_value)
    use constants
    implicit none
    real(8), intent(in)  :: r
    real(8), intent(out) :: pdf_value
    real(8)              :: term1, term2

    ! 로그 정규 분포 계산
    term1     = - (log(r / rm))**2 / (2.0d0 * sigma**2)
    term2     = (r * sigma * sqrt(2.0d0 * pi))
    pdf_value = (1.0d0 / term2) * exp(term1)

end subroutine lognormal
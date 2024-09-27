subroutine adiabatic_process(z, T, p, rho, w, dt)
    use constants
    implicit none
    real, intent(inout) :: z, T, p, rho
    real, intent(in)    :: w, dt

    ! 위치 업데이트
    z = z + w * dt

    ! 압력 업데이트
    p = p0 * exp(-g * (z - z0) / (R * T))

    ! 밀도 업데이트
    rho = p / (R * T)

    ! 온도 업데이트
    T = T - (g / cp) * w * dt
end subroutine adiabatic_process
subroutine adiabatic_process(z, T, p, rho, w, dt)
    use constants
    implicit none
    real, intent(inout) :: z, T, p, rho
    real, intent(in)    :: w, dt

    z   = z + w * dt

    p   = p0 * exp(-g * (z - z0) / (R * T))

    rho = p / (R * T)

    T   = T - (g / cp) * w * dt

end subroutine adiabatic_process
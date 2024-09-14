program adiabatic_box_model
    use constants
    implicit none
  
    real    :: dt, total_time, time, z, T, p, rho
    integer :: i
  
    z    = z0
    T    = T0
    p    = p0
    dt   = 1
    time = 0.0
  
    do i = 1, 100
       time = time + dt
  
       z = z + w * dt
  
       p = p0 * exp(-g * (z - z0) / (R * T))
  
       rho = p / (R * T)
  
       T = T - (g / cp) * w * dt
  
       print *, 'Time:', time, 'Temperature:', T
    end do
  
end program adiabatic_box_model

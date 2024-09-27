program adiabatic_box_model
    use constants
    implicit none

    real    :: dt, total_time, time, z, T, p, rho, w
    integer :: unit_num

    z      = z0
    T      = T0
    p      = p0
    dt     = 1.0
    time   = 0.0

    ! 파일 열기
    unit_num = 10
    call open_file(unit_num, 'results.txt')

    do
        time = time + dt

        if (time .gt. i_time) exit

        w = 1.0 + 0.01 * time

        if (w > 5.0) w = 5.0

        call adiabatic_process(z, T, p, rho, w, dt)

        print *, 'Time:', time, 'Temperature:', T

        call write_results(unit_num, time, T, p, rho, z)
    end do

    close(unit=unit_num)
end program adiabatic_box_model
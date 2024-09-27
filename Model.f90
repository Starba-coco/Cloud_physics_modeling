program adiabatic_box_model
    use constants
    implicit none

    real    :: dt, total_time, time, z, T, p, rho, w
    integer :: unit_num, num_w_vals, i
    integer, parameter :: max_size = 10
    real    :: w_vals(max_size), w_time(max_size)
    
    namelist /input_params/ w_time, w_vals

    ! 배열 초기화
    w_time = 0.0
    w_vals = 0.0

    z    = z0
    T    = T0
    p    = p0
    dt   = 1.0
    time = 0.0

    ! 파일 열기
    unit_num = 10
    call open_file(unit_num, 'results.txt')

    ! namelist 파일에서 값 읽기
    open(unit=20, file='./input.nml', status='old')
    read(20, input_params)
    close(20)

    do
        time = time + dt

        if (time > i_time) exit

        ! 시간에 따른 w 값 설정
        if (time <= w_time(1)) then
            w = w_vals(1)
        else if (time <= w_time(2)) then
            w = w_vals(2)
        else if (time <= w_time(3)) then
            w = w_vals(3)
        else
            w = 5.0
        end if

        call adiabatic_process(z, T, p, rho, w, dt)

        print *, 'Time:', time, 'Temperature:', T

        call write_results(unit_num, time, T, p, rho, z)

    end do

    close(unit=unit_num)
end program adiabatic_box_model
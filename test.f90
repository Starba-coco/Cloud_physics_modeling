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

        ! 상태 업데이트 서브루틴 호출
        call update_state(z, T, p, rho, w, dt)

        ! 화면에 출력
        print *, 'Time:', time, 'Temperature:', T

        ! 결과를 파일에 저장하는 서브루틴 호출
        call write_results(unit_num, time, T, p, rho, z)
    end do

    ! 파일 닫기
    close(unit=unit_num)
end program adiabatic_box_model
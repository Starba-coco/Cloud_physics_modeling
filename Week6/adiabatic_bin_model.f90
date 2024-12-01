program adiabatic_bin_model
    use constants
    implicit none

    ! ============================================================
    ! 변수 선언
    integer              :: i, i_vant

    ! 시간 및 공간 변수
    real(8)              :: dt, time, z, i_time(2), w_val(2)
    real(8)              :: save_start_time, save_end_time, save_interval

    ! 물리적 상수 및 초기 상태
    real(8)              :: Ms, rho_s, S, w, T, p, rho, qv

    ! 분포 함수 및 에어로졸 특성 변수
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:), &
                             r_drop(:), r_center_drop(:), n_bin_drop(:), log_r_drop(:), &
                             m(:), dm(:), m_new(:), r_new(:), n_bin_new(:)
    real(8)              :: pdf_value, dr, total_aerosols

    ! 출력 파일 관련 변수
    character(len=20)    :: aerosol_type

    ! 빈(bin) 경계 계산 변수
    real(8)              :: log_rmin, log_rmax, delta_logr
    real(8)              :: log_rmin_drop, log_rmax_drop, delta_logr_drop
    real(8)              :: rmin, rmax, rmin_drop, rmax_drop
    ! ============================================================

    ! namelist에서 읽을 변수 선언
    namelist /input_params/ aerosol_type, w_val, T, z, p, i_time, &
                            rmin, rmax, qv, rmin_drop, rmax_drop
    ! ============================================================
    ! 변수 초기화
    dt               = 1.0d0
    time             = 0.0d0
    aerosol_type     = 'NaCl'   ! 기본값 설정
    ! ============================================================

    ! namelist 파일에서 변수 읽기
    open(unit=15, file='input.nml', status='old', action='read')
    read(15, nml=input_params)
    close(15)

    ! 초기 공기 밀도 계산
    call cal_rhoa(p, T, rho)

    ! 배열 할당
    allocate(m(nbin_drop))
    allocate(m_new(nbin_drop))
    allocate(dm(nbin_drop))
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(n_bin_new(nbin_drop))
    allocate(log_r(nbin+1))
    allocate(r_drop(nbin_drop+1))
    allocate(r_new(nbin_drop))
    allocate(r_center_drop(nbin_drop))
    allocate(n_bin_drop(nbin_drop))
    allocate(log_r_drop(nbin_drop+1))

    ! bin 경계 계산을 위한 로그 값 계산
    log_rmin   = log(rmin)
    log_rmax   = log(rmax)
    delta_logr = (log_rmax - log_rmin) / nbin

    ! drop_bin 경계 계산을 위한 로그 값 계산
    log_rmin_drop   = log(rmin_drop)
    log_rmax_drop   = log(rmax_drop)
    delta_logr_drop = (log_rmax_drop - log_rmin_drop) / nbin_drop

    ! 에어로졸 설정
    call set_aerosol(aerosol_type, Ms, rho_s, i_vant)

    ! Aerosol과 Droplet의 Bin model
    call cal_bin(r, r_center, log_r, log_rmin, delta_logr, &
                 r_drop, r_center_drop, log_r_drop, log_rmin_drop, delta_logr_drop)

    ! 각 bin마다 포함되는 입자 수 계산 (Aerosol)
    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr       = r(i + 1) - r(i)
        n_bin(i) = pdf_value * dr
    end do

    ! from #/m3 to #/kg
    n_bin = n_bin / rho
    total_aerosols = sum(n_bin)

    ! 결과 출력 파일 설정
    open(10, file='output.txt', status='unknown', action='write')

    ! 헤더 출력
    write(10, '(A7,3X,A6,3X,A10,3X,A10,3X,A8,3X,A12,3X,A12,3X,A12)') &
         'Time(s)', 'w(m/s)', 'T(K)', 'p(Pa)', 'z(m)', 'RH(%)', 'Nd(#/kg)', 'qv(g/kg)'

    ! 단열 상승 과정
    do
        time = time + dt

        ! i_time에 따라 w 설정
        if (time <= i_time(1)) then
            w = w_val(1)
        else
            w = w_val(2)
        end if
        if (time > i_time(2)) exit

        ! 단열 과정
        call adiabatic_process(z, T, p, rho, w, dt)

        call cal_rhoa(p, T, rho)

        ! 상대 습도 및 과포화도 업데이트
        call update_saturation(p, T, qv, S)

        call activation(p, T, S, nbin, nbin_drop, r, n_bin, n_bin_drop, &
                        Ms, rho_s, i_vant)

        call condensation(T, p, S, dt, m, m_new, rho, r_new, n_bin_drop, r_center_drop, qv)

        call redistribution(m, m_new, n_bin_drop, nbin_drop, r_new, r_drop, n_bin_new)

        write(10, '(F7.2,3X,F6.2,3X,F10.2,3X,F10.2,3X,F8.2,3X,E12.4,3X,E12.4,3X,F12.6)') &
              time, w, T, p, z, (S+1.0d0) * 100.0d0, sum(n_bin_drop), qv

        if (time >= 205.0d0 .and. time <= 220.0d0) then
            call write_distribution(nbin, r_center, n_bin, delta_logr, time)
        end if
    
        ! 루프 종료 조건
        if (time > 220.0d0) exit
    end do
    close(10)

    ! 메모리 해제
    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)
    deallocate(r_drop)
    deallocate(r_center_drop)
    deallocate(n_bin_drop)
    deallocate(log_r_drop)
    deallocate(m)
    deallocate(m_new)
    deallocate(dm)
    deallocate(r_new)
    deallocate(n_bin_new)
end program adiabatic_bin_model

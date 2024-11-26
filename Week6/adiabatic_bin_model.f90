program adiabatic_bin_model
    use constants
    implicit none

    ! ============================================================
    ! 변수 선언
    integer              :: i, i_vant

    ! 시간 및 공간 변수
    real(8)              :: dt, time, z, i_time(2), w_val(2)

    ! 물리적 상수 및 초기 상태
    real(8)              :: a, b, Ms, rho_s, q, S, w, T, p, rho

    ! 분포 함수 및 에어로졸 특성 변수
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:), r_drop(:), r_center_drop(:), n_bin_drop(:), log_r_drop(:), m(:), radius_drop(:), dm(:), m_new(:), r_new(:), n_bin_new(:)
    real(8)              :: rs, rc, Sc, pdf_value, dr, total_aerosols, activated_drops
    real(8)              :: delta_qv, n_activated
    real(8)              :: delta_volume, delta_mass_per_particle, delta_mass_bin

    ! 출력 파일 관련 변수
    character(len=20)    :: aerosol_type
    integer              :: dist_unit, nt_steps

    ! 에어로졸 직경 관련 변수
    real(8)              :: diameter_nm, dlogD, y_value

    ! 추가된 변수 선언
    real(8)              :: log_rmin, log_rmax, delta_logr, log_rmin_drop, log_rmax_drop, delta_logr_drop
    real(8)              :: rmin, rmax, rmin_drop, rmax_drop, total_mass, mean_radius
    real(8)              :: qv, dqc, activate_ratio, ln_r1, ln_r2, ln_rc, ln_r0
    ! ============================================================

    ! namelist에서 읽을 변수 선언
    namelist /input_params/ aerosol_type, w_val, T, z, p, i_time, rmin, rmax, qv, rmin_drop, rmax_drop
    ! ============================================================
    ! 변수 초기화
    dt           = 1.0d0
    time         = 0.0d0
    aerosol_type = 'NaCl'   ! 기본값 설정

    ! ============================================================

    ! namelist 파일에서 변수 읽기
    open(unit=15, file='input.nml', status='old', action='read')
    read(15, nml=input_params)
    close(15)

    ! calculate the inital air density to convert units from #/m3 to #/kg
    call cal_rhoa(p, T, rho)

    ! 배열 할당
    allocate(m(nbin_drop))
    allocate(m_new(nbin_drop))
    
    allocate(radius_drop(nbin_drop))
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(n_bin_new(nbin_drop))
    allocate(log_r(nbin+1))
    allocate(r_drop(nbin_drop+1))
    allocate(r_new(nbin_drop+1))
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
    call cal_bin(r, r_center, log_r, log_rmin, delta_logr, &
                 r_drop, r_center_drop, log_r_drop, log_rmin_drop, delta_logr_drop)

    ! 각 bin마다 포함되는 입자 수 계산

    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr          = r(i + 1) - r(i)
        n_bin(i)    = pdf_value * dr
    end do

    ! from #/m3 to #/kg
    n_bin = n_bin / rho
    total_aerosols = sum(n_bin)

    call write_distribution(nbin, r_center, n_bin, delta_logr)

    ! 결과 출력 파일 설정
    open(10, file='output.txt', status='unknown', action='write')
    ! open(11, file='test.txt',   status='unknown', action='write')

    ! 헤더
    write(10, '(A)') 'Time(s)   w(m/s)         T(K)        p(Pa)       z(m)    RH(%)  Nd(#/kg), qv(g/kg)'
    ! write(11) 'Time(s)   w(m/s)         T(K)        p(Pa)       z(m)    RH(%)  Nd(#/kg), qv(g/kg)'

    ! 단열 상승 과정
    do
        time = time + dt
    
        ! i_time에 따라 w 설정
        if (time <= i_time(1)) then
            w = 1.0
        else
            w = 0.0
        end if
        if (time > i_time(2)) exit
        ! print *, "1", w
        ! 단열 과정
        call adiabatic_process(z, T, p, rho, w, dt)
        ! print *, "2", w, z
        call cal_rhoa(p, T, rho)
    
        ! 상대 습도 및 과포화도 업데이트
        call update_saturation(p, T, qv, S)

        call activation(p, T, S, nbin, nbin_drop, r, n_bin, n_bin_drop, &
                        Ms, rho_s, i_vant)

        call condensation(T, p, S, dt, m, m_new, rho, r_new, n_bin_drop, r_center_drop, qv)

        call redistribution(m, m_new, n_bin_drop, nbin_drop, r_new, r_drop, n_bin_new)
        
        write(10, '(F7.2,3X,F6.2,3X,F10.2,3X,F10.2,3X,F8.2,3X,E12.4,3X,E12.4,3X,F12.6)') &
              time, w, T, p, z, (S+1.0) * 100.0d0, sum(n_bin_drop), qv
    end do

    close(10)

    ! 메모리 해제
    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)
    ! deallocate(r_drop(nbin_drop+1))
    ! deallocate(r_center_drop(nbin_drop))
    ! deallocate(n_bin_drop(nbin_drop))
    ! deallocate(log_r_drop(nbin_drop+1))

end program adiabatic_bin_model
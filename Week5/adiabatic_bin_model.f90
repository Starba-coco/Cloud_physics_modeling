program adiabatic_bin_model
    use constants
    implicit none

    ! ============================================================
    ! 변수 선언
    integer              :: i, i_vant

    ! 시간 및 공간 변수
    real(8)              :: dt, time, z, i_time(2), w_val(2)
    logical              :: end_flag

    ! 물리적 상수 및 초기 상태
    real(8)              :: a, b, Ms, rho_s, q, S, w, T, p, rho

    ! 분포 함수 및 에어로졸 특성 변수
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:)
    real(8), allocatable :: r_drop(:), m_drop(:), r_center_drop(:), n_bin_drop(:), log_r_drop(:)
    real(8)              :: rs, rc, Sc, pdf_value, dr, total_aerosols, activated_drops
    real(8)              :: delta_qv, n_activated, delta_T
    real(8)              :: delta_volume, delta_mass_per_particle, delta_mass_bin, total_liquid_mass

    ! 출력 파일 관련 변수
    character(len=20)    :: aerosol_type
    integer              :: dist_unit

    ! 에어로졸 직경 관련 변수
    real(8)              :: diameter_nm, dlogD, y_value

    ! 추가된 변수 선언
    real(8)              :: log_rmin, log_rmax, delta_logr, log_rmin_drop, log_rmax_drop, delta_logr_drop
    real(8)              :: rmin, rmax, rmin_drop, rmax_drop
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

    ! 초기 공기 밀도 계산 (단위 변환을 위해)
    call cal_rhoa(p, T, rho)

    ! bin_model 계산을 위한 배열 할당
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(log_r(nbin+1))
    allocate(m_drop(nbin_drop))
    allocate(r_drop(nbin_drop))
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

    ! 단위 변환: #/m3에서 #/kg로
    n_bin = n_bin / rho
    total_aerosols = sum(n_bin)

    ! 활성화된 물방울 수 초기화
    n_bin_drop = 0.0d0

    ! 결과 출력 파일 설정
    open(10, file='output.txt', status='unknown', action='write')

    ! 헤더
    write(10, '(A7,3X,A6,3X,A10,3X,A10,3X,A8,3X,A12,3X,A12,3X,A12)') &
              'Time(s)', 'w(m/s)', 'T(K)', 'p(Pa)', 'z(m)', 'RH(%)', 'Nd(#/kg)', 'qv(g/kg)'

    ! 단열 상승 과정
    do
        call update_time_and_w(time, dt, i_time, w_val, w, end_flag)
        if (end_flag) exit

        ! 단열 과정
        call adiabatic_process(z, T, p, rho, w, dt)
        call cal_rhoa(p, T, rho)

        ! 상대 습도 및 과포화도 업데이트
        call update_saturation(p, T, qv, S)

        ! Activation
        call activation(p, T, qv, S, nbin, nbin_drop, r, n_bin, n_bin_drop, &
                        r_center, r_center_drop, Ms, rho_s, i_vant)

        ! 활성화된 물방울의 질량과 반경 초기화
        do i = 1, nbin_drop
            if (n_bin_drop(i) > 0.0d0) then
                r_drop(i) = r_center_drop(i)
                m_drop(i) = (4.0d0 / 3.0d0) * pi * r_drop(i)**3 * rho_w
            end if
        end do

        ! 응결 과정
        delta_qv = 0.0d0
        delta_T  = 0.0d0
        total_liquid_mass = 0.0d0

        do i = 1, nbin_drop
            if (n_bin_drop(i) > 0.0d0) then
                call condensation(T, p, S, dt, m_drop(i), r_drop(i), n_bin_drop(i), rho, delta_qv, delta_T)
                total_liquid_mass = total_liquid_mass + m_drop(i) * n_bin_drop(i)
            end if
        end do
        print *, time, total_liquid_mass
        ! 수증기 혼합비와 온도 업데이트
        qv = qv + delta_qv
        T  = T + delta_T

        ! 과포화도 재계산
        call update_saturation(p, T, qv, S)

        ! 결과 출력 (파일로 쓰기)
        write(10, '(F7.2,3X,F6.2,3X,F10.2,3X,F10.2,3X,F8.2,3X,E12.4,3X,E12.4,3X,F12.6)') &
              time, w, T, p, z, (S + 1.0d0) * 100.0d0, sum(n_bin_drop), qv
    end do

    close(10)

    ! 메모리 해제
    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)
    deallocate(m_drop)
    deallocate(r_drop)
    deallocate(r_center_drop)
    deallocate(n_bin_drop)
    deallocate(log_r_drop)

end program adiabatic_bin_model

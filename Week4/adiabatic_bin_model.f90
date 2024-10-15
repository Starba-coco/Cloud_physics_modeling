program adiabatic_bin_model
    use constants
    implicit none

    ! 변수 선언
    ! 공통 변수
    integer              :: i
    real(8)              :: delta_logr, log_rmin, log_rmax
    real(8)              :: pdf_value, dr, total_drops
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:)

    ! 단열 상승 관련 변수
    real                 :: dt, total_time, time, z, T, p, rho, w
    integer              :: num_w_vals
    integer, parameter   :: max_size = 10
    real                 :: w_vals(max_size), w_time(max_size)
    
    ! 추가 변수 선언
    real                 :: q, RH, S
    real(8)              :: a, b, rs, rc, Sc, activated_drops

    ! 배열 선언
    real(8), allocatable :: m_s(:), r_c(:), S_c(:)
    allocate(m_s(nbin))
    allocate(r_c(nbin))
    allocate(S_c(nbin))

    namelist /input_params/ w_time, w_vals

    ! 초기화
    dt     = 1.0
    time   = 0.0
    total_time = i_time
    z      = z0
    T      = T0
    p      = p0
    rho    = p / (R_gas * T)
    q      = qv0   ! 초기 수증기 혼합비 설정

    ! 상승 속도 관련 초기화 (필요에 따라 수정 가능)
    w_time = 0.0
    w_vals = w

    ! bin_model 계산을 위한 배열 할당
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(log_r(nbin+1))

    open(unit=20, file='./input.nml', status='old')
    read(20, input_params)
    close(20)

    ! 로그 스케일로 bin 경계 계산
    log_rmin   = log(rmin)
    log_rmax   = log(rmax)
    delta_logr = (log_rmax - log_rmin)/nbin

    do i = 1, nbin+1
        log_r(i) = log_rmin + (i-1)*delta_logr
        r(i)     = exp(log_r(i))
    end do

    ! bin 중심 계산 (기하 평균)
    do i = 1, nbin
        r_center(i) = sqrt(r(i)*r(i+1))
    end do

    ! 각 bin마다 포함되는 drop의 개수 계산
    total_drops = 0.0
    
    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr          = r(i+1) - r(i)
        n_bin(i)    = N0 * pdf_value * dr
        total_drops = total_drops + n_bin(i)
    end do

    ! 용질 질량 계산
    do i = 1, nbin
        m_s(i) = (4.0 / 3.0) * pi * (r_center(i))**3 * rho_s
    end do

    ! 단열 상승 과정 시뮬레이션
    do
        time = time + dt
        if (time > total_time) exit

        ! 상승 속도 결정 (필요에 따라 수정 가능)
        if (time <= w_time(1)) then
            w = w_vals(1)
        else if (time <= w_time(2)) then
            w = w_vals(2)
        else if (time <= w_time(3)) then
            w = w_vals(3)
        else
            w = w
        end if

        ! 단열 과정 계산
        call adiabatic_process(z, T, p, rho, w, dt, q, RH, S)

        ! 코흘러 상수 계산 (온도 의존)
        a = (2.0 * sigma_v) / (Rv * rho_w * T)
        b = (rho_s * Mw) / (Ms * rho_w)

        ! 활성화된 입자 수 초기화
        activated_drops = 0.0

        ! 각 bin에 대해 임계 반경과 임계 과포화도 계산 및 활성화 여부 판단
        do i = 1, nbin
            rs = r_center(i)

            ! 임계 반경 계산
            rc = ( (a * rs**3) / (3.0 * b) ) ** 0.25

            ! 임계 과포화도 계산
            Sc = exp( (a / rc) - (b * rc**3) / rs**3 )

            ! 현재 과포화도와 비교하여 활성화 여부 판단
            if ( (S * 100) >= (Sc - 1.0) ) then
                activated_drops = activated_drops + n_bin(i)
            end if
        end do

        ! 결과 출력
        print *, 'Time:', time, ' w:', w, ' T:', T, ' p:', p, ' z:', z, ' RH:', RH, ' Activated Drops:', activated_drops

    end do

    ! 배열 해제
    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)
    deallocate(m_s)
    deallocate(r_c)
    deallocate(S_c)

end program adiabatic_bin_model

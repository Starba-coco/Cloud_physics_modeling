program adiabatic_bin_model
    use constants
    implicit none

    integer              :: i
    real(8)              :: delta_logr, log_rmin, log_rmax
    real(8)              :: pdf_value, dr, total_drops
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:)
    real(8)              :: dt, time, z, T, p, rho, w
    integer              :: num_w_vals, i_vant
    integer, parameter   :: max_size = 10
    real(8)              :: w_vals(max_size), w_time(max_size)
    real(8)              :: q, S
    real(8)              :: a, b, rs, rc, Sc, activated_drops, Ms, rho_s
    real(8), allocatable :: m_s(:)
    character(len=20)    :: aerosol_type

    ! 출력 파일 관련 변수
    character(len=100)   :: output_file

    namelist /input_params/ aerosol_type, w_vals

    ! 초기화
    dt           = 1.0d0
    time         = 0.0d0
    z            = z0
    T            = T0
    p            = p0
    rho          = p0 / (R_gas * T)
    q            = qv0
    aerosol_type = 'NaCl'

    ! bin_model 계산을 위한 배열 할당
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(log_r(nbin+1))
    allocate(m_s(nbin))

    open(unit=20, file='./input.nml', status='old')
    read(20, nml=input_params)
    close(20)

    ! 에어로졸 특성 설정 서브루틴 호출
    call set_aerosol(aerosol_type, Ms, rho_s, i_vant)

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
    total_drops = 0.0d0

    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr          = r(i+1) - r(i)
        n_bin(i)    = N0 * pdf_value * dr
        total_drops = total_drops + n_bin(i)
    end do

    ! 용질 질량 계산
    do i = 1, nbin
        m_s(i) = (4.0d0 / 3.0d0) * pi * (r_center(i))**3 * rho_s
    end do

    ! 출력 파일 열기
    open(10, file='output.txt', status='unknown', action='write')

    ! 헤더 작성
    write(10, '(A)') 'Time(s)   w(m/s)         T(K)        p(Pa)       z(m)    RH(%)  Activated Drops'

    ! 단열 상승 과정 시뮬레이션
    do
        time = time + dt
        if (time > i_time) exit

        w = w_vals(1)

        ! 단열 과정 계산
        call adiabatic_process(z, T, p, rho, w, dt, q, S)

        ! 코흘러 상수 계산 (온도 의존)
        a = (2.0d0 * sigma_v) / (Rv * rho_w * T)
        b = (rho_s * Mw)      / (Ms * rho_w)

        activated_drops = 0.0d0

        ! 각 bin에 대해 임계 반경과 임계 과포화도 계산 및 활성화 여부 판단
        do i = 1, nbin
            rs = r_center(i)

            ! 임계 반경 계산
            rc = ( (a * rs**3) / (3.0d0 * b) ) ** 0.25d0
            ! 임계 과포화도 계산
            Sc = exp( (a / rc) - (b * rc**3) / rs**3 )

            ! 현재 과포화도와 비교하여 활성화 여부 판단
            if ( S >= (Sc - 1.0d0) ) then
                activated_drops = activated_drops + n_bin(i)
            end if
        end do

        ! 결과 출력 (파일로 쓰기)
        write(10, '(F7.2, 3X, F6.2, 3X, F10.2, 3X, F10.2, 3X, F8.2, 3X, F6.2, 5X, E12.4)') &
        time, w, T, p, z, (S+1.0d0)*100.0d0, activated_drops

    end do

    close(10)

    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)
    deallocate(m_s)

end program adiabatic_bin_model

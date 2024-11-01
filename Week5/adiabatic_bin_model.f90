program adiabatic_bin_model
    use constants
    implicit none

    ! ============================================================
    ! 상수 및 파라미터 선언
    integer              :: i, i_vant
    
    ! 시간 및 공간 변수
    real(8)              :: dt, time, z
    
    ! 물리적 상수 및 초기 상태
    real(8)              :: delta_logr, log_rmin, log_rmax
    real(8)              :: a, b, Ms, rho_s, q, S, w, w_vals, T, p, rho
    
    ! 분포 함수 및 에어로졸 특성 변수
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:)
    real(8)              :: rs, rc, Sc, pdf_value, dr, total_drops, activated_drops
    
    ! 출력 파일 관련 변수
    character(len=20)    :: aerosol_type
    integer              :: dist_unit
    
    ! 에어로졸 직경 관련 변수
    real(8)              :: diameter_nm, dlogD, y_value
    ! ============================================================
    namelist /input_params/ aerosol_type, w
    ! ============================================================
    ! 변수 초기화
    dt           = 1.0d0
    time         = 0.0d0
    z            = z0
    T            = T0
    p            = p0
    rho          = p0 / (R_gas * T)
    q            = qv0
    aerosol_type = 'NaCl'
    ! ============================================================

    open(unit=15, file='input.nml', status='old', action='read')
    read(15, nml=input_params)
    close(15)

    ! bin_model 계산을 위한 배열 할당
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(log_r(nbin+1))
    
    ! 에어로졸 설정
    call set_aerosol(aerosol_type, Ms, rho_s, i_vant)
    
    ! bin 경계 계산
    log_rmin   = log(rmin)
    log_rmax   = log(rmax)
    delta_logr = (log_rmax - log_rmin) / nbin
    
    do i = 1, nbin + 1
        log_r(i) = log_rmin + (i - 1) * delta_logr
        r(i)     = exp(log_r(i))
    end do
    
    ! bin 중심 계산 (기하 평균)
    do i = 1, nbin
        r_center(i) = sqrt(r(i) * r(i + 1))
    end do
    
    ! 각 bin마다 포함되는 입자 수 계산
    total_drops = 0.0d0
    
    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr          = r(i + 1) - r(i)
        n_bin(i)    = pdf_value * dr
        total_drops = total_drops + n_bin(i)
    end do
    
    ! 분포 데이터를 출력하기 위한 파일 설정
    dist_unit = 20
    open(dist_unit, file='distribution.txt', status='unknown', action='write')
    
    ! 헤더
    write(dist_unit, '(A)') 'Diameter(nm)    dN/dlogD * 1e-3 (cm^-3/nm)'
    
    ! dlogD 계산 (상수로 계산)
    dlogD = delta_logr / log(10.0d0)  ! ln(10)으로 나누어 log10 기반으로 변환
    
    ! 분포 데이터를 파일에 출력
    do i = 1, nbin
        ! 입자 지름 계산 (nm)
        diameter_nm = 2.0d0 * r_center(i) * 1.0d9
    
        y_value = (n_bin(i) / dlogD) * 1.0d-9
    
        write(dist_unit, '(E12.5, 3X, E12.5)') diameter_nm, y_value
    end do
    
    close(dist_unit)
    
    ! 결과 출력 파일 설정
    open(10, file='output.txt', status='unknown', action='write')
    
    ! 헤더
    write(10, '(A)') 'Time(s)   w(m/s)         T(K)        p(Pa)       z(m)    RH(%)  Activated Drops'
    
    ! 단열 상승 과정
    do
        time = time + dt
        if (time > i_time) exit
    
        ! 단열 상승
        call adiabatic_process(z, T, p, rho, w, dt, q, S)
    
        ! kohler constants
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
            if (S >= (Sc - 1.0d0)) then
                activated_drops = activated_drops + n_bin(i)
            end if
        end do
    
        ! 결과 출력 (파일로 쓰기)
        write(10, '(F7.2, 3X, F6.2, 3X, F10.2, 3X, F10.2, 3X, F8.2, 3X, F6.2, 5X, E12.4)') &
               time, w, T, p, z, (S + 1.0d0) * 100.0d0, activated_drops
    end do
    
    close(10)
    
    ! 메모리 해제
    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)
    

end program adiabatic_bin_model

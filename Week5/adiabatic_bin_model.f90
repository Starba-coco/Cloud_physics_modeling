program adiabatic_bin_model
    use constants
    implicit none

    ! ============================================================
    ! 변수 선언
    integer              :: i, i_vant

    ! 시간 및 공간 변수
    real(8)              :: dt, time, z, i_time

    ! 물리적 상수 및 초기 상태
    real(8)              :: a, b, Ms, rho_s, q, S, w, T, p, rho

    ! 분포 함수 및 에어로졸 특성 변수
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:)
    real(8)              :: rs, rc, Sc, pdf_value, dr, total_drops, activated_drops
    real(8)              :: delta_qv, n_activated
    real(8)              :: delta_volume, delta_mass_per_particle, delta_mass_bin

    ! 출력 파일 관련 변수
    character(len=20)    :: aerosol_type
    integer              :: dist_unit

    ! 에어로졸 직경 관련 변수
    real(8)              :: diameter_nm, dlogD, y_value

    ! 추가된 변수 선언
    real(8)              :: log_rmin, log_rmax, delta_logr
    real(8)              :: rmin, rmax
    real(8)              :: qv
    ! ============================================================

    ! namelist에서 읽을 변수 선언
    namelist /input_params/ aerosol_type, w, T, z, p, i_time, rmin, rmax, qv
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

    ! namelist에서 읽은 변수들을 사용하여 초기화
    rho = p / (R_dry * T)
    q   = qv

    ! bin_model 계산을 위한 배열 할당
    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(log_r(nbin+1))

    ! bin 경계 계산을 위한 로그 값 계산
    log_rmin   = log(rmin)
    log_rmax   = log(rmax)
    delta_logr = (log_rmax - log_rmin) / nbin

    ! 에어로졸 설정
    call set_aerosol(aerosol_type, Ms, rho_s, i_vant)
    call cal_bin(r, r_center, log_r, log_rmin, delta_logr)

    ! 각 bin마다 포함되는 입자 수 계산
    total_drops = 0.0d0

    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr          = r(i + 1) - r(i)
        n_bin(i)    = pdf_value * dr
        total_drops = total_drops + n_bin(i)
    end do

    call write_distribution(nbin, r_center, n_bin, delta_logr)

    ! 결과 출력 파일 설정
    open(10, file='output.txt', status='unknown', action='write')

    ! 헤더
    write(10, '(A)') 'Time(s)   w(m/s)         T(K)        p(Pa)       z(m)    RH(%)  Activated_Drops, q'

    ! 단열 상승 과정
    do
        time = time + dt
        if (time > i_time) exit

        ! 단열 과정
        call adiabatic_process(z, T, p, rho, w, dt, q, S)

        ! kohler 상수 계산 (온도 의존)
        a = (2.0d0 * sigma_v)      / (Rv * rho_w * T)
        b = i_vant * ((rho_s * Mw) / (Ms * rho_w))

        activated_drops = 0.0d0
        delta_qv        = 0.0d0
        
        ! 각 bin에 대해 임계 반경과 임계 과포화도 계산
        do i = 1, nbin
            rs = r_center(i)
        
            ! 임계 과포화도 계산
            Sc = sqrt((4.0d0 * (a**3)) / (27.0d0 * b * rs**3))
        
            ! 임계 반경 계산
            ! rc = ( (a * rs**3) / (3.0d0 * b) ) ** 0.25d0
            rc = ((a / 3.0d0) / (4.0d0 * b * (Sc**2))) ** (1.0d0 / 3.0d0)

            ! 활성화 여부 판단
            if (S >= Sc) then
                ! 활성화된 입자의 수 (개수/부피, #/m^3)
                n_activated = n_bin(i)
                activated_drops = activated_drops + n_activated
        
                ! 응결된 물의 부피 계산 (m^3)
                delta_volume = ((4.0d0 / 3.0d0) * pi * (rs**3)) * rho_w * n_activated
        
                ! 총 응결된 물의 질량에 합산 (kg/kg)
                delta_qv = delta_qv + delta_mass_bin / rho
            end if
        end do
        
        ! 수증기 혼합비 업데이트 (kg/kg)
        q = q - delta_qv

        ! q가 음수가 되지 않도록 체크
        if (q < 0.0d0) q = 0.0d0

        ! 상대 습도 및 과포화도 업데이트
        call update_saturation(p, T, q, S)

        ! 결과 출력 (파일로 쓰기)
        write(10, '(F7.2,3X,F6.2,3X,F10.2,3X,F10.2,3X,F8.2,3X,E12.4,3X,E12.4,3X,F12.6)') &
        time, w, T, p, z, (S + 1.0d0) * 100.0d0, activated_drops, q
    end do

    close(10)

    ! 메모리 해제
    deallocate(r)
    deallocate(r_center)
    deallocate(n_bin)
    deallocate(log_r)

end program adiabatic_bin_model

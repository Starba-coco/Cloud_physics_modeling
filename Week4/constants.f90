module constants
    implicit none

    ! 물리 상수
    real,    parameter :: g       = 9.80665             ! 중력 가속도        (m/s^2)
    real,    parameter :: R_gas   = 287.058             ! 건조 공기 기체 상수 (J/kg K)
    real,    parameter :: cp      = 1005.0              ! 공기의 비열        (J/kg K)
    real,    parameter :: T0      = 300.0               ! 초기 온도          (K)
    real,    parameter :: z0      = 0.0                 ! 초기 고도          (m)
    real,    parameter :: p0      = 102400.0            ! 초기 압력          (Pa)
    real,    parameter :: i_time  = 3600                ! 적분 시간          (s)
    real,    parameter :: rmin    = 1.00E-08            ! Drop의 최소 크기   (m)
    real,    parameter :: rmax    = 1.00E-04            ! Drop의 최대 크기   (m)
    real,    parameter :: N0      = 1.00E9              ! 수농도             (#/m^3)
    real,    parameter :: sigma   = 1.2                 ! 로그 정규 분포의 표준편차
    real,    parameter :: rm      = 1.00E-06            ! 기하 평균 반경      (m)
    real(8), parameter :: pi      = 4.0d0 * atan(1.0d0) ! 파이

    real,    parameter :: qv0     = 0.01093             ! 초기 수증기 혼합비  (kg/kg)
    integer, parameter :: nbin    = 200                 ! bin의 개수         (#)

    ! 수증기 및 공기 몰 질량
    real,    parameter :: Mv      = 18.01528d-3         ! 수증기 몰 질량      (kg/mol)
    real,    parameter :: Md      = 28.9644d-3          ! 건조 공기 몰 질량   (kg/mol)
    real,    parameter :: epsilon = Mv / Md             ! 0.622 근처
    real,    parameter :: Rv      = 461.5               ! 수증기 기체 상수    (J/kg·K)

    ! 물의 특성
    real,    parameter :: Mw      = 18.01528d-3         ! 물의 몰 질량   (kg/mol)
    real,    parameter :: rho_w   = 1000.0d0            ! 물의 밀도      (kg/m^3)
    real,    parameter :: sigma_v = 0.0728              ! 물의 표면 장력 (J/m^2)

    ! Aerosol (염화나트륨, NaCl)
    real,    parameter :: Ms      = 58.44d-3            ! NaCl의 몰 질량     (kg/mol)
    real,    parameter :: rho_s   = 2160.0d0            ! NaCl의 밀도        (kg/m^3)
    integer, parameter :: i_vant  = 2                   ! 반 호프 인자       (i = 2)

end module constants

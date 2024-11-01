module constants
    implicit none

    ! 물리 상수
    real(8),    parameter :: g       = 9.80665d0              ! 중력 가속도        (m/s^2)
    real(8),    parameter :: R_gas   = 287.058d0              ! 건조 공기 기체 상수 (J/kg·K)
    real(8),    parameter :: cp      = 1005.0d0               ! 공기의 비열        (J/kg·K)
    real(8),    parameter :: T0      = 300.0d0                ! 초기 온도          (K)
    real(8),    parameter :: z0      = 0.0d0                  ! 초기 고도          (m)
    real(8),    parameter :: p0      = 102400.0d0             ! 초기 압력          (Pa)
    real(8),    parameter :: i_time  = 3600.0d0               ! 적분 시간          (s)
    real(8),    parameter :: rmin    = 1.00D-09               ! Drop의 최소 크기   (m) - 수정됨
    real(8),    parameter :: rmax    = 1.00D-04               ! Drop의 최대 크기   (m)
    real(8),    parameter :: pi      = 4.0d0 * atan(1.0d0)    ! 파이

    real(8),    parameter :: qv0     = 0.01093d0              ! 초기 수증기 혼합비  (kg/kg)
    integer(8), parameter :: nbin    = 300                    ! bin의 개수 증가 (#) - 수정됨

    ! 수증기 및 공기의 특성
    real(8),    parameter :: Mv      = 18.01528d-3            ! 수증기 몰 질량      (kg/mol)
    real(8),    parameter :: Md      = 28.9644d-3             ! 건조 공기 몰 질량   (kg/mol)
    real(8),    parameter :: epsilon = Mv / Md                ! 약 0.622
    real(8),    parameter :: Rv      = 461.5d0                ! 수증기 기체 상수    (J/kg·K)

    ! 물의 특성
    real(8),    parameter :: Mw      = 18.01528d-3            ! 물의 몰 질량   (kg/mol)
    real(8),    parameter :: rho_w   = 1000.0d0               ! 물의 밀도      (kg/m^3)
    real(8),    parameter :: sigma_v = 0.0728d0               ! 물의 표면 장력 (J/m^2)

    ! 로그 정규 분포 파라미터 (세 개의 모드)
    integer,    parameter :: n_modes = 3                      ! 모드의 수

    ! 모드 1
    real(8),    parameter :: N0_val_1  = 8.5E9                ! 수농도 (#/m^3)
    real(8),    parameter :: rm_val_1  = 3.65D-9              ! 기하 평균 반경 (m)
    real(8),    parameter :: sig_val_1 = 1.60d0               ! 기하 표준 편차 (무차원)

    ! 모드 2
    real(8),    parameter :: N0_val_2  = 5.8E9
    real(8),    parameter :: rm_val_2  = 22.5D-9
    real(8),    parameter :: sig_val_2 = 1.65d0

    ! 모드 3
    real(8),    parameter :: N0_val_3  = 0.95E9
    real(8),    parameter :: rm_val_3  = 77D-9
    real(8),    parameter :: sig_val_3 = 1.60d0

end module constants

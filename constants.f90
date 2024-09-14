module constants
    implicit none

    real, parameter :: g  = 9.80665    ! 중력 가속도    (m/s^2)
    real, parameter :: R  = 2.87       ! 공기 기체 상수 (J/kg K)
    real, parameter :: cp = 1005.0     ! 공기의 비열    (J/kg K)
    real, parameter :: T0 = 300.0      ! 초기 온도      (K)
    real, parameter :: z0 = 0.0        ! 초기 고도      (m)
    real, parameter :: w  = 5.0        ! 상승 속도      (m/s)
    real, parameter :: p0 = 1024.0     ! 초기 압력      (hPa)

end module constants
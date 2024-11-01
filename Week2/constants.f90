module constants
    implicit none

    real,    parameter :: g      = 9.80665       ! 중력 가속도    (m/s^2)
    real,    parameter :: R      = 287.058       ! 공기 기체 상수 (J/kg K)
    real,    parameter :: cp     = 1005.0        ! 공기의 비열    (J/kg K)
    real,    parameter :: T0     = 300.0         ! 초기 온도      (K)
    real,    parameter :: z0     = 0.0           ! 초기 고도      (m)
    ! real,    parameter :: w      = 5.0           ! 상승 속도      (m/s)
    real,    parameter :: p0     = 102400.0      ! 초기 압력      (Pa)
    real,    parameter :: i_time = 1600          ! 적분 시간      (s)
    real,    parameter :: rmin   = 1.00E-06      ! Drop의 최소 크기 (m)
    real,    parameter :: rmax   = 1.00E-02      ! Drop의 최대 크기 (m)
    real,    parameter :: N0     = 1.00E8        ! 수농도           (#/m^3)
    real,    parameter :: sigma  = 1.2           ! 표준편차
    real,    parameter :: rm     = 5.00E-04
    real,    parameter :: AS     = 0.13214       ! (NH4)2SO4      (kg/mol)
    real,    parameter :: Mv     = 
    real(8), parameter :: pi     = 4.0*atan(1.0) ! 파이
    real,    parameter :: 
    integer, parameter :: nbin   = 200           ! bin의 개수       (#)


end module constants


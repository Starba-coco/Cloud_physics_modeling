module constants
    implicit none

    real,    parameter :: rmin   = 1.00E-06      ! Drop의 최소 크기 (m)
    real,    parameter :: rmax   = 1.00E-02      ! Drop의 최대 크기 (m)
    real,    parameter :: N0     = 1.00E8        ! 수농도           (#/m^3)
    real,    parameter :: sigma  = 1.2           ! 표준편차
    real,    parameter :: rm     = 5.00E-04
    real(8), parameter :: pi     = 4.0*atan(1.0) ! 파이
    integer, parameter :: nbin   = 200           ! bin의 개수       (#)

end module constants
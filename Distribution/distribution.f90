program bin_model
    use constants
    implicit none
    integer              :: i, unit_num
    real(8)              :: delta_logr, log_rmin, log_rmax
    real(8)              :: pdf_value, dr, total_drops
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:)
    character(len=100)   :: filename

    unit_num = 20
    filename = 'result_distribution.txt'

    allocate(r(nbin+1))
    allocate(r_center(nbin))
    allocate(n_bin(nbin))
    allocate(log_r(nbin+1))

    log_rmin   = log(rmin)
    log_rmax   = log(rmax)
    delta_logr = (log_rmax - log_rmin)/nbin

    ! 각 Bin의 Size를 log scale로 계산
    do i = 1, nbin+1
        log_r(i) = log_rmin + (i-1)*delta_logr
        r(i)     = exp(log_r(i))
    end do
    
    ! 빈 중심 계산 (기하평균)
    do i = 1, nbin
        r_center(i) = sqrt(r(i)*r(i+1))
    end do

    ! 각 bin마다 포함되는 drop의 개수 계산
    total_drops = 0.0
    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr       = r(i+1) - r(i)
        n_bin(i) = N0 * pdf_value * dr
        total_drops = total_drops + n_bin(i)
    end do

    open(unit=unit_num, file=filename, status='replace')
    write(unit_num, '(A)') 'Bin Center (m), Drop Count (#/m³)'

    ! 각 bin 결과 파일에 쓰기
    do i = 1, nbin
        write(unit_num, '(F15.8, 2X,",", E15.8)') r_center(i), n_bin(i)
    end do

    print *, total_drops

    close(unit_num)

end program bin_model

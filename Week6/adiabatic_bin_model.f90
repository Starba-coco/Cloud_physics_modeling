program adiabatic_bin_model
    use constants
    implicit none

    integer              :: i, i_vant
    real(8)              :: dt, time, z, i_time(2), w_val(2)
    real(8)              :: Ms, rho_s, S, w, T, p, rho, qv
    real(8), allocatable :: r(:), r_center(:), n_bin(:), log_r(:), &
                            r_drop(:), r_center_drop(:), n_bin_drop(:), log_r_drop(:), &
                            m(:), dm(:), m_new(:), r_new(:), n_bin_new(:), Vt(:)
    real(8), allocatable :: k(:,:), ec(:,:)
    real(8)              :: pdf_value, dr, total_aerosols
    character(len=20)    :: aerosol_type
    real(8)              :: log_rmin, log_rmax, delta_logr
    real(8)              :: log_rmin_drop, log_rmax_drop, delta_logr_drop
    real(8)              :: rmin, rmax, rmin_drop, rmax_drop

    namelist /input_params/ aerosol_type, w_val, T, z, p, i_time, &
                            rmin, rmax, qv, rmin_drop, rmax_drop

    dt           = 0.1d0
    time         = 0.0d0
    aerosol_type = 'NaCl'

    open(unit=15, file='input.nml', status='old', action='read')
    read(15, nml=input_params)
    close(15)

    call cal_rhoa(p, T, rho)

    ! Allocate arrays
    allocate(m(nbin_drop), Vt(nbin_drop), m_new(nbin_drop), dm(nbin_drop))
    allocate(r(nbin+1), r_center(nbin), n_bin(nbin), n_bin_new(nbin_drop), log_r(nbin+1))
    allocate(r_drop(nbin_drop+1), r_new(nbin_drop), r_center_drop(nbin_drop), n_bin_drop(nbin_drop), log_r_drop(nbin_drop+1))
    allocate(k(nbin_drop, nbin_drop))
    allocate(ec(nbin_drop, nbin_drop))

    log_rmin   = log(rmin)
    log_rmax   = log(rmax)
    delta_logr = (log_rmax - log_rmin) / nbin

    log_rmin_drop   = log(rmin_drop)
    log_rmax_drop   = log(rmax_drop)
    delta_logr_drop = (log_rmax_drop - log_rmin_drop) / nbin_drop

    call set_aerosol(aerosol_type, Ms, rho_s, i_vant)
    call cal_bin(r, r_center, log_r, log_rmin, delta_logr, &
                 r_drop, r_center_drop, log_r_drop, log_rmin_drop, delta_logr_drop)

    call terminal_velocity(r_center_drop, rho, T, P, Vt)

    do i = 1, nbin
        call lognormal(r_center(i), pdf_value)
        dr       = r(i + 1) - r(i)
        n_bin(i) = pdf_value * dr
    end do

    ! #/m3 -> #/kg 변환
    n_bin = n_bin / rho
    total_aerosols = sum(n_bin)

    open(10, file='output.txt', status='unknown', action='write')
    write(10, '(A7,3X,A6,3X,A10,3X,A10,3X,A8,3X,A12,3X,A12,3X,A12)') &
         'Time(s)', 'w(m/s)', 'T(K)', 'p(Pa)', 'z(m)', 'RH(%)', 'Nd(#/kg)', 'qv(g/kg)'

    ! ---- 효율 계산을 위한 호출 (예시) ----
    ! 여기서 r_center_drop와 동일 차원의 r를 사용해야 함.
    ! 만약 effic가 r, r_center_drop 둘 다 같은 의미의 반경 배열을 기대한다면
    ! r_center_drop를 r로 사용하거나, 별도 배열 r_eff를 만들어 복제 가능.
    !
    ! 가정: effic에서는 두 배열이 같은 의미이므로 둘 다 r_center_drop 사용
    !
    call effic(r_center_drop, r_center_drop, ec)

    do
        if (time <= i_time(1)) then
            w = w_val(1)
        else
            w = w_val(2)
        end if
        if (time > i_time(2)) exit

        call adiabatic_process(z, T, p, rho, w, dt)
        call cal_rhoa(p, T, rho)
        call update_saturation(p, T, qv, S)

        call activation(p, T, S, nbin, nbin_drop, r, n_bin, n_bin_drop, Ms, rho_s, i_vant)
        call condensation(T, p, S, dt, m, m_new, rho, r_new, n_bin_drop, r_center_drop, qv)
        call redistribution(m, m_new, n_bin_drop, nbin_drop, r_new, r_drop, n_bin_new)

        if (time == 0.0d0 .or. time == 600.0d0 .or. time == 1200.0d0 .or. time == 1800.0d0 .or. &
            time == 2400.0d0 .or. time == 3000.0d0 .or. time == 3600.0d0 ) then
            call write_drop_distribution(nbin_drop, r_center_drop, n_bin_drop, delta_logr_drop, time)
        end if 

        write(10, '(F7.2,3X,F6.2,3X,F10.2,3X,F10.2,3X,F8.2,3X,E12.4,3X,E12.4,3X,F12.6)') &
              time, w, T, p, z, (S+1.0d0)*100.0d0, sum(n_bin_drop), qv

        if (time >= 205.0d0 .and. time <= 220.0d0) then
            call write_aerosol_distribution(nbin, r_center, n_bin, delta_logr, time)
        end if

        call collision_kernel(r_center_drop, Vt, ec, k)

        ! call collision(dt, rho, r_center_drop, n_bin_drop, k)

        time = time + dt    
    end do

    close(10)

    deallocate(r, r_center, n_bin, log_r, r_drop, r_center_drop, n_bin_drop, log_r_drop)
    deallocate(m, m_new, dm, r_new, n_bin_new, Vt, k, ec)
end program adiabatic_bin_model

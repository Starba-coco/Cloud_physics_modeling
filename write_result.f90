subroutine write_results(unit_num, time, T)
    implicit none
    integer, intent(in) :: unit_num
    real,    intent(in) :: time, T

    write(unit_num, '(F10.1, ",", F10.6)') time, T
end subroutine write_results

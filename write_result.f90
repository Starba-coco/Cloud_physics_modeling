subroutine write_results(unit_num, time, T)
    implicit none
    integer, intent(in) :: unit_num
    real, intent(in)    :: time, T

    ! 파일에 출력 (CSV 형식)
    write(unit_num, '(F10.2, ",", F10.2)') time, T
end subroutine write_results

subroutine open_file(unit_num, filename)
    implicit none
    integer, intent(in) :: unit_num
    character(len=*), intent(in) :: filename

    ! 파일 열기
    open(unit=unit_num, file=filename, status='unknown')

    ! 헤더 쓰기
    write(unit_num, '(A)') 'Time(s), Temperature(K)'
end subroutine open_file

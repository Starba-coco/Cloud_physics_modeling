subroutine lognormal(r, pdf_value)
    use constants
    implicit none
    real(8), intent(in)  :: r
    real(8), intent(out) :: pdf_value
    real(8)              :: term1, term2

    ! 로그 정규 분포 계산
    term1     = - (log(r / rm))**2 / (2.0d0 * sigma**2)
    term2     = (r * sigma * sqrt(2.0d0 * pi))
    pdf_value = (1.0d0 / term2) * exp(term1)

end subroutine lognormal
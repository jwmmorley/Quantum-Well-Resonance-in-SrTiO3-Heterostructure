module spectral
    implicit none
    
    public spectral_function

contains
    function spectral_function(values, weights, point, broadening) result(a)
        real*8, intent(in) :: weights(:), values(:), point, broadening
        integer            :: i           
        real*8             :: a
        a = 0d0
        do i = 1, size(values)
            a = a - dimag(weights(i) &
                / (point - values(i) + dcmplx(0d0, broadening)))
        end do
        
    end function spectral_function

end module spectral

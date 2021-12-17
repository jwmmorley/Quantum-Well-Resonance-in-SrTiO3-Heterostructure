module spectral
    implicit none
    
    public spectral_function

contains
    function spectral_function(values, weights, point, broadening) result(a)
        real*8, intent(in) :: weights, values, point, broadening
        integer            :: i           
        real*8             :: a
        a = 0d0
        a = a - dimag(weights &
            / (point - values + dcmplx(0d0, broadening)))
        
    end function spectral_function

end module spectral

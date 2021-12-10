module bulk
    implicit none
    
    private :: bulk_identity_matrix
    public  :: bulk_build, bulk_add_potential

contains
    function bulk_build(h, num_layers) result(hb)
    integer,    intent(in) :: num_layers
    complex*16, intent(in) :: h(:, :, :)
    complex*16             :: hb(num_layers * size(H(1, :, 1)), &
                                 num_layers * size(H(1, :, 1)))
    integer                :: num_bands, max_hopping
    integer                :: i, j
    num_bands = size(h(1, 1, :))
    max_hopping = (size(h(:, 1, 1)) - 1) / 2
    hb = dcmplx(0d0, 0d0)
    do j = 1, num_layers
        do i = 1, num_layers
            if (abs(i - j) <= max_hopping) then
                hb(num_bands * (i - 1) + 1:num_bands * i, num_bands * (j - 1) + 1:num_bands * j) &
                    = h(i - j + max_hopping + 1, :, :)
            end if
        end do
    end do
    
    end function bulk_build
    
    subroutine bulk_add_potential(hb, potential)
    complex*16, intent(inout) :: hb(:, :)
    real*8,     intent(in)    :: potential(:)
    integer                   :: num_layers, num_bands
    integer                   :: i, j
    num_layers = size(potential)
    num_bands = size(hb(1, :)) / num_layers
    
    do j = 1, num_layers
        do i = 1, num_layers
            if (i == j) then
                hb(num_bands * (i - 1) + 1:num_bands * i, num_bands * (j - 1) + 1:num_bands * j) &
                    = hb(num_bands * (i - 1) + 1:num_bands * i, num_bands * (j - 1) + 1:num_bands * j) &
                    + bulk_identity_matrix(num_bands, potential(i))
            end if
        end do
    end do
    
    end subroutine bulk_add_potential
    
    function bulk_identity_matrix(length, scaler) result(matrix)
        integer, intent(in) :: length
        real*8,  intent(in) :: scaler
        real*8              :: matrix(length, length)
        integer             :: i
        matrix = 0d0
        do i = 1, length
            matrix(i, i) = scaler
        end do
    
    end function bulk_identity_matrix

end module bulk

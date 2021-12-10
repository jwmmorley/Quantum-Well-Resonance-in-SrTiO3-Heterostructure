module transform
    implicit none
    
    public:: transform_r_to_kz

contains
        function transform_r_to_kz(t, w, kx, ky) result(h)
            integer,    intent(in) :: w(:, :, :)
            real*8,     intent(in) :: kx, ky
            complex*16, intent(in) :: t(:, :, :, :, :)
            complex*16             :: h(size(t(1, 1, :, 1, 1)), size(t(1, 1, 1, 1, :)), size(t(1, 1, 1, 1, :)))
            real*8                 :: phase
            integer                :: num_bands, max_hopping, x, y, z, i, j
            num_bands   = size(t(1, 1, 1, 1, :))
            max_hopping = (size(t(1, 1, :, 1, 1)) - 1) / 2
            h           = dcmplx(0d0, 0d0)
            do x = 1, 2 * max_hopping + 1
                do y = 1, 2 * max_hopping + 1
                    do z = 1, 2 * max_hopping + 1
                        do i = 1, num_bands
                            do j = 1, num_bands
                                phase = -1d0 * (dble(x - max_hopping - 1) * kx + dble(y - max_hopping - 1) * ky)
                                h(z, i, j) = h(z, i, j) &
                                    + t(x, y, z, i, j) * dcmplx(cos(phase), sin(phase)) / dble(w(x, y, z))
                            end do
                        end do
                    end do
                end do
            end do
                        
        end function transform_r_to_kz

end module transform

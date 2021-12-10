module import
    implicit none

    private :: import_meta_data, import_weight_data, import_tunneling_data
    public  :: import_data

contains
    subroutine import_meta_data(file_number, num_r, num_bands)
        integer, intent(in)  :: file_number
        integer, intent(out) :: num_bands, num_r
        read(file_number, *) ! Skip Header
        read(file_number, *) num_bands
        read(file_number, *) num_r

    end subroutine import_meta_data

    subroutine import_weight_data(file_number, w, num_r)
        integer, intent(in)               :: file_number, num_r
        integer                           :: i, iso
        integer, intent(out), allocatable :: w(:)
        allocate(w(num_r))
        do i = 1, num_r
            read(file_number, fmt = "(I5.1)", iostat = iso, advance = "no") w(i)
            if (iso /= 0) then
                read(file_number, fmt = "(I5.1)", iostat = iso, advance = "no") w(i)
            end if
        end do
        read(file_number, *) ! Go to new line
        
    end subroutine import_weight_data


    subroutine import_tunneling_data(file_number, rx, ry, rz, bi, bj, t, max_hopping, num_r, num_bands)
        integer,    intent(in)               :: file_number, num_r, num_bands
        integer                              :: i
        complex*16, intent(out), allocatable :: t(:)
        integer,    intent(out), allocatable :: rx(:), ry(:), rz(:), bi(:), bj(:)
        integer,    intent(out)              :: max_hopping
        max_hopping = 0
        allocate(rx(num_r * num_bands**2))
        allocate(ry(num_r * num_bands**2))
        allocate(rz(num_r * num_bands**2))
        allocate(bi(num_r * num_bands**2))
        allocate(bj(num_r * num_bands**2))
        allocate(t(num_r * num_bands**2))
        do i = 1, num_r * num_bands**2
            read(file_number, &
                fmt = "(I5.1, I5.1, I5.1, I5.1, I5.1, F12.6, F12.6)") &
                   rx(i), ry(i), rz(i), bi(i), bj(i), t(i)
            if (abs(rx(i)) > max_hopping) then
                max_hopping = abs(rx(i))
            end if
            if (abs(ry(i)) > max_hopping) then
                max_hopping = abs(ry(i))
            end if
            if (abs(rz(i)) > max_hopping) then
                max_hopping = abs(rz(i))
            end if
        end do         

    end subroutine import_tunneling_data

    subroutine import_data(file_name, hr, hrw, max_hopping, num_bands)
        character(256), intent(in)               :: file_name
        complex*16,     allocatable              :: t(:)
        integer,        allocatable              :: rx(:), ry(:), rz(:), bi(:), bj(:), w(:)
        integer                                  :: file_number, num_r, i, n
        integer,        intent(out)              :: max_hopping, num_bands
        integer,        intent(out), allocatable :: hrw(:, :, :)
        complex*16,     intent(out), allocatable :: hr(:, :, :, :, :)
        open(newunit = file_number, file = file_name)
        call import_meta_data(file_number, num_r, num_bands)
        call import_weight_data(file_number, w, num_r)
        call import_tunneling_data(file_number, rx, ry, rz, bi, bj, t, max_hopping, num_r, num_bands)
        close(file_number)
        allocate(hr(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1, num_bands, num_bands))
        allocate(hrw(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1))
        hr  = dcmplx(0d0, 0d0)
        hrw = -1
        n   = 0
        do i = 1, num_r * num_bands**2
            hr(rx(i) + max_hopping + 1, ry(i) + max_hopping + 1, rz(i) + max_hopping + 1, &
                bi(i), bj(i)) = t(i)
            if (hrw(rx(i) + max_hopping + 1, ry(i) + max_hopping + 1, rz(i) + max_hopping + 1) == -1) then
                n = n + 1
                hrw(rx(i) + max_hopping + 1, ry(i) + max_hopping + 1, rz(i) + max_hopping + 1) = w(n)
            end if
        end do

    end subroutine import_data

end module import

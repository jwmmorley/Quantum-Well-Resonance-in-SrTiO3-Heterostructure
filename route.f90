module route
    use cpu
    use spglib_f08
    implicit none
        
        private :: route_count_unique_integer, route_loc_integer
        public  :: route_range_double, route_build

contains
    function route_range_double(start_val, end_val, num_points) result(array)
        integer, intent(in) :: num_points
        real*8,  intent(in) :: start_val, end_val
        real*8              :: array(num_points)
        integer             :: i
        array = start_val
        do i = 2, num_points
            array(i) = array(i) + (i - 1) * (end_val - start_val) / (num_points - 1)
        end do
    
    end function route_range_double
    
    function route_count_unique_integer(array) result(num)
        integer, intent(in) :: array(:)
        integer             :: num
        integer             :: unique(size(array))
        integer             :: i, j
        num = 0
        unique = 0
        if (any(array == 0)) then
            num = num + 1
        end if
        do i = 1, size(array)
            if (array(i) /= 0 .and. .not. any(unique == array(i))) then
                num = num + 1
                do j = 1, size(unique)
                    if (unique(j) == 0) then
                        unique(j) = array(i)
                        exit
                    end if
                end do
            end if
        end do
    
    end function route_count_unique_integer
    
    function route_loc_integer(array, val) result(lo)
        integer, intent(in) :: array(:)
        integer, intent(in) :: val
        integer             :: lo
        integer             :: i
        do i = 1, size(array)
            if (array(i) == val) then
                lo = i
                exit
            end if
        end do
    
    end function route_loc_integer
    
    subroutine route_build(Kxs, Kys, weight, k_path_map, num_length, length_scale, lattice, positions, atom_types)
        integer, intent(in)               :: num_length, atom_types(:)
        real*8,  intent(in)               :: lattice(:, :), positions(:, :), length_scale
        real*8,  intent(out), allocatable :: kxs(:), kys(:)
        integer, intent(out), allocatable :: weight(:)
        integer, intent(out), allocatable :: k_path_map(:)
        integer                           :: map(floor(dble(num_length) / length_scale)**3), &
                                             grid_point(3, floor(dble(num_length) / length_scale)**3), & 
                                             num_ir_grid, i, j, n, used_map_length
        integer, allocatable              :: used_map_to_map(:), used_grid_point(:, :), used_map(:)
        num_ir_grid = spg_get_ir_reciprocal_mesh(grid_point, map, &
                      (/ floor(dble(num_length) / length_scale), floor(dble(num_length) / length_scale), &
                      floor(dble(num_length) / length_scale) /), (/ 0, 0, 0 /), 1, &
                      transpose(lattice), transpose(positions), atom_types, size(atom_types), 1D-5)
        used_map_length = num_length
        allocate(used_map(used_map_length**2))
        used_map = -1
        do i = 1, used_map_length
            do j = 1, used_map_length
                if (i <= used_map_length / 2 .and. j >= used_map_length / 2 + mod(used_map_length, 2)) then
                    used_map((i - 1) * used_map_length + j) = map(&
                        ((i + floor(dble(num_length) / length_scale) - used_map_length / 2) - 1) &
                            * floor(dble(num_length / length_scale)) &
                        + (j - used_map_length / 2 - mod(used_map_length, 2))) + 1
                end if
                if (i >= used_map_length / 2 + mod(used_map_length, 2) .and. j <= used_map_length / 2) then
                    used_map((i - 1) * used_map_length + j) = map(&
                        ((i - used_map_length / 2) - mod(used_map_length, 2)) &
                            * floor(dble(num_length / length_scale)) &
                        + (j + floor(dble(num_length) / length_scale) - used_map_length / 2)) + 1
                end if
                if (i <= used_map_length / 2 .and. j <= used_map_length / 2) then
                    used_map((i - 1) * used_map_length + j) = map(&
                        ((i + floor(dble(num_length) / length_scale) - used_map_length / 2) - 1) &
                            * floor(dble(num_length / length_scale)) &
                        + (j + floor(dble(num_length) / length_scale) - used_map_length / 2)) + 1
                end if
                if (i > used_map_length / 2 .and. j > used_map_length / 2) then
                    used_map((i - 1) * used_map_length + j) = map(&
                        ((i - used_map_length / 2) - 1) * floor(dble(num_length / length_scale)) + (j - used_map_length / 2)) + 1
                end if
            end do
        end do
        allocate(used_map_to_map(route_count_unique_integer(used_map)))
        allocate(weight(size(used_map_to_map)))
        allocate(used_grid_point(3, size(used_map_to_map)))
        allocate(kxs(size(used_map_to_map)))
        allocate(kys(size(used_map_to_map)))
        used_map_to_map = -1
        weight = 0
        n = 0
        do i = 1, size(used_map)
            if (used_map(i) == -1) then
                cycle
            end if
            if (any(used_map_to_map == used_map(i))) then
                used_map(i) = route_loc_integer(used_map_to_map, used_map(i))
            else ! Unique
                n = n + 1
                used_map_to_map(n)    = used_map(i)
                used_grid_point(:, n) = grid_point(:, used_map(i))
                used_map(i) = n
            end if
        end do
        do i = 1, size(used_map)
            if (used_map(i) == -1) then
                cycle
            end if
            weight(used_map(i)) = weight(used_map(i)) + 1
        end do
        do i = 1, size(used_grid_point(1, :))
            kxs(i) = dble(used_grid_point(1, i)) / dble(floor(dble(num_length) / length_scale))
            kys(i) = dble(used_grid_point(2, i)) / dble(floor(dble(num_length) / length_scale))
        end do
        allocate(k_path_map(2 * floor((used_map_length - 1) / 2d0) + 1))
        do i = 1, floor((used_map_length - 1) / 2d0)
            k_path_map(i) = used_map((used_map_length - i) * used_map_length + used_map_length - i + 1)
        end do
        do i = floor((used_map_length - 1) / 2d0) + 2, size(k_path_map)
            k_path_map(i) = used_map(&
                (i + mod(used_map_length, floor((used_map_length - 1) / 2d0)) - 2) * used_map_length &
                + used_map_length - floor((used_map_length - 1) / 2d0))
        end do
        k_path_map(floor((used_map_length - 1) / 2d0) + 1) = used_map(&
            (used_map_length - floor((used_map_length - 1) / 2d0) - 1) * used_map_length &
            + used_map_length - floor((used_map_length - 1) / 2d0))
        
    end subroutine route_build

end module route

module potential
    implicit none
    
    public :: potential_add_well
    
contains
    subroutine potential_add_well(array, start_ly, end_ly , initial_depth, final_depth)
    real*8,  intent(inout) :: array(:)
    integer, intent(in)    :: start_ly, end_ly
    real*8,  intent(in)    :: initial_depth, final_depth
    real*8                 :: potential(end_ly - start_ly + 1)
    integer                :: i
    potential = initial_depth
    do i = 0, end_ly - start_ly 
        potential(i + 1) = potential(i + 1) + (final_depth - initial_depth) * i / (end_ly - start_ly)
    end do
    do i = 1, size(potential)
        if (start_ly >= 1 .and. end_ly <= size(array)) then 
            array(start_ly + i -1 ) = array(start_ly + i -1) + potential(i)
        else 
            print*, ' The layer no. inputted is out of bound for potential. '
        end if
    end do
    
    end subroutine potential_add_well    
    
end module potential

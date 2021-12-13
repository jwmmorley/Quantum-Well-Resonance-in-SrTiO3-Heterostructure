module potential
    implicit none
    
    public :: potential_add_well
    
contains
    subroutine potential_add_well(array, start_ly, end_ly, initial_depth, final_depth)
    real*8,  intent(inout) :: array(:)
    integer, intent(in)    :: start_ly, end_ly
    real*8,  intent(in)    :: initial_depth, final_depth
    real*8                 :: potential(abs(end_ly - start_ly) + 1)
    integer                :: i, start, finish
    start = min(start_ly, end_ly)
    finish = max(start_ly, end_ly)
    potential = initial_depth
    if (finish /= start) then
        do i = 1, finish - start + 1 
            potential(i) = potential(i) + (final_depth - initial_depth) * (i - 1) / (finish - start)
        end do
    end if
    do i = 1, size(potential)
        if (1 <= start + i - 1 .and. start + i - 1 <= size(array)) then
            array(start + i - 1) = array(start + i -1) + potential(i)
        else 
            print*, " Warning! The layer no. inputted is out of bound for potential. The well has been cropped. "
        end if
    end do
    
    end subroutine potential_add_well    
    
end module potential

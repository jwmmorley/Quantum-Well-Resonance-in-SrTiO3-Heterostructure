module timer
    implicit none
    
    public :: timer_start, timer_split
    
    integer, private, save :: timer_start_tic       = -1
    integer, private, save :: timer_last_tic        = -1
    real*8,  private, save :: timer_tics_per_second

contains
    subroutine timer_start()
        call system_clock(timer_start_tic, timer_tics_per_second)
        
    end subroutine timer_start
    
    subroutine timer_split()
        integer :: timer_tic
        if (timer_start_tic == -1) then
            call timer_start()
        end if
        if (timer_last_tic == -1) then
            timer_last_tic = timer_start_tic
        end if
        call system_clock(timer_tic)
        write (*, fmt = "(A22, f8.1, A9, f8.1)") &
            "Ellapsed Time; Total: ", dble(timer_tic - timer_start_tic) / timer_tics_per_second, &
            ", Split: ", dble(timer_tic - timer_last_tic) / timer_tics_per_second
        timer_last_tic = timer_tic
        
    end subroutine timer_split
    
end module timer

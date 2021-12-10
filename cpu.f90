module cpu
    use mpi
    implicit none
    
    private :: cpu_sum_double, &
               cpu_broadcast_integer, cpu_broadcast_double, cpu_broadcast_double_complex
    public  :: cpu_start, cpu_stop, &
               cpu_get_id, cpu_get_num, cpu_get_master_id, cpu_is_master, &
               cpu_broadcast, cpu_sum
    
    logical, private, save      :: cpu_started   = .false.
    integer, private, parameter :: cpu_master_id = 0
    integer, private, save      :: cpu_id, cpu_num
    integer, private            :: cpu_error
    
    interface cpu_sum
        module procedure cpu_sum_double
    end interface cpu_sum
    interface cpu_broadcast
        module procedure cpu_broadcast_integer, cpu_broadcast_double, cpu_broadcast_double_complex
    end interface cpu_broadcast
    
contains
    subroutine cpu_start()
        call MPI_INIT(cpu_error)
        call MPI_COMM_RANK(MPI_COMM_WORLD, cpu_id, cpu_error)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, cpu_num, cpu_error)
        cpu_started = .true.
    
    end subroutine cpu_start
    
    subroutine cpu_stop()
        if (.not. cpu_started) then
            call cpu_start()
        end if
        call MPI_FINALIZE(cpu_error)
        
    end subroutine cpu_stop
    
    subroutine cpu_split_work(start, stop, ammount)
        integer, intent(in)  :: ammount
        integer, intent(out) :: start, stop
        if (.not. cpu_started) then
            call cpu_start()
        end if
        if (cpu_num - mod(cpu_id, cpu_num) < cpu_id) then
            start = ammount - (cpu_num - cpu_id) * (floor(dble(ammount) / dble(cpu_num)) + 1) + 1
            stop = start + floor(dble(ammount) / dble(cpu_num))
        else
            start = cpu_id * floor(dble(ammount) / dble(cpu_num)) + 1
            stop = start + floor(dble(ammount) / dble(cpu_num)) - 1
        end if
        
    end subroutine cpu_split_work
        
    subroutine cpu_sum_double(v_in, v_out)
        real*8, intent(in)  :: v_in(..)
        real*8, intent(out) :: v_out(..)
        call MPI_ALLREDUCE(v_in, v_out, size(v_in), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, cpu_error)
    
    end subroutine cpu_sum_double
        
    subroutine cpu_broadcast_integer(v, length, target_id)
        integer, intent(in)    :: length
        integer, intent(inout) :: v(..)
        integer, optional      :: target_id
        integer                :: id = cpu_master_id
        if (.not. cpu_started) then
            call cpu_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_INT, id, MPI_COMM_WORLD, cpu_error)
    
    end subroutine cpu_broadcast_integer
    
    subroutine cpu_broadcast_double(v, length, target_id)
        integer, intent(in)    :: length
        real*16, intent(inout) :: v(..)
        integer, optional      :: target_id
        integer                :: id = cpu_master_id 
        if (.not. cpu_started) then
            call cpu_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_DOUBLE_PRECISION, id, MPI_COMM_WORLD, cpu_error)
    
    end subroutine cpu_broadcast_double
    
    subroutine cpu_broadcast_double_complex(v, length, target_id)
        integer,    intent(in)    :: length
        complex*16, intent(inout) :: v(..)
        integer,    optional      :: target_id
        integer                   :: id = cpu_master_id 
        if (.not. cpu_started) then
            call cpu_start()
        end if
        if (present(target_id)) then
            id = target_id
        end if
        call MPI_BCAST(v, length, MPI_DOUBLE_COMPLEX, id, MPI_COMM_WORLD, cpu_error)
    
    end subroutine cpu_broadcast_double_complex
    
    function cpu_get_id()
        integer :: cpu_get_id 
        if (.not. cpu_started) then
            call cpu_start()
        end if
        cpu_get_id = cpu_id
    
    end function cpu_get_id
    
    function cpu_get_num()
        integer :: cpu_get_num
        if (.not. cpu_started) then
            call cpu_start()
        end if
        cpu_get_num = cpu_num
    
    end function cpu_get_num
    
    function cpu_get_master_id()
        integer :: cpu_get_master_id
        if (.not. cpu_started) then
            call cpu_start()
        end if
        cpu_get_master_id = cpu_master_id
        
    end function cpu_get_master_id
    
    function cpu_is_master()
        logical :: cpu_is_master
        if (.not. cpu_started) then
            call cpu_start()
        end if
        cpu_is_master = (cpu_get_id() == cpu_get_master_id())
        
    end function cpu_is_master

end module cpu

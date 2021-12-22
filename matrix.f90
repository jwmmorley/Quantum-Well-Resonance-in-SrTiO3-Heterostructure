module matrix
    implicit none
    
    private :: matrix_get_eigen_range, matrix_get_eigen_bound
    public  :: matrix_get_eigen
    
    integer, private, save :: range_lwork  = 1
    integer, private, save :: range_lrwork = 1
    integer, private, save :: range_liwork = 1
    integer, private, save :: bound_lwork  = 1
    integer, private, save :: bound_lrwork = 1
    integer, private, save :: bound_liwork = 1
    
    interface matrix_get_eigen
        module procedure matrix_get_eigen_range, matrix_get_eigen_bound
    end interface matrix_get_eigen
            

contains
    recursive subroutine matrix_get_eigen_range(matrix, val, vec, num_found, min_val, max_val)
        complex*16, intent(in)  :: matrix(:, :)
        real*8,     intent(in)  :: min_val, max_val
        complex*16, intent(out) :: vec(size(matrix(:, 1)), size(matrix(1, :)))
        real*8,     intent(out) :: val(size(matrix(:, 1)))
        integer,    intent(out) :: num_found
        complex*16              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                 :: info
        integer                 :: isuppz(2 * size(matrix(:, 1)))
        complex*16              :: work(range_lwork)
        real*8                  :: rwork(range_lrwork)
        integer                 :: iwork(range_liwork)
        copy_matrix = matrix
        if (range_lwork < 2 * size(matrix(:, 1)) .or. &
            range_lrwork < 24 * size(matrix(:, 1)) .or. &
            range_liwork < 10 * size(matrix(:, 1))) then
            range_lwork  = -1
            range_lrwork = -1
            range_liwork = -1
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), min_val, max_val, 0, 0, 0d0, &
                num_found, val, vec, size(matrix(:, 1)), isuppz, work, range_lwork, rwork, range_lrwork, iwork, range_liwork, info)
            range_lwork  = int(work(1))
            range_lrwork = int(rwork(1))
            range_liwork = int(iwork(1))
            call matrix_get_eigen(matrix, val, vec, num_found, min_val, max_val)
        else
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), min_val, max_val, 0, 0, 0d0, &
                num_found, val, vec, size(matrix(:, 1)), isuppz, work, range_lwork, rwork, range_lrwork, iwork, range_liwork, info)
        end if
        
    end subroutine matrix_get_eigen_range
    
    recursive subroutine matrix_get_eigen_bound(matrix, val, vec, num_found, min_val, max_val)
        complex*16, intent(in)  :: matrix(:, :)
        integer,    intent(in)  :: min_val, max_val
        complex*16, intent(out) :: vec(size(matrix(:, 1)), size(matrix(1, :)))
        real*8,     intent(out) :: val(size(matrix(:, 1)))
        integer,    intent(out) :: num_found
        complex*16              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                 :: info
        integer                 :: isuppz(2 * size(matrix(:, 1)))
        complex*16              :: work(bound_lwork)
        real*8                  :: rwork(bound_lrwork)
        integer                 :: iwork(bound_liwork)
        copy_matrix = matrix
        if (bound_lwork < 2 * size(matrix(:, 1)) .or. &
            bound_lrwork < 24 * size(matrix(:, 1)) .or. &
            bound_liwork < 10 * size(matrix(:, 1))) then
            bound_lwork  = -1
            bound_lrwork = -1
            bound_liwork = -1
            call zheevr('V', 'I', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), 0d0, 0d0, min_val, max_val, 0d0, &
                num_found, val, vec, size(matrix(:, 1)), isuppz, work, bound_lwork, rwork, bound_lrwork, iwork, bound_liwork, info)
            bound_lwork  = int(work(1))
            bound_lrwork = int(rwork(1))
            bound_liwork = int(iwork(1))
            call matrix_get_eigen(matrix, val, vec, num_found, min_val, max_val)
        else
            call zheevr('V', 'I', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), 0d0, 0d0, min_val, max_val, 0d0, &
                num_found, val, vec, size(matrix(:, 1)), isuppz, work, bound_lwork, rwork, bound_lrwork, iwork, bound_liwork, info)
        end if
        
    end subroutine matrix_get_eigen_bound

end module matrix

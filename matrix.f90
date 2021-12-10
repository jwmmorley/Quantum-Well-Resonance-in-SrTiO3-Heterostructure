module matrix
    implicit none
    
    public  :: matrix_get_eigen
    
    integer, private, save :: lwork  = 1
    integer, private, save :: lrwork = 1
    integer, private, save :: liwork = 1
            

contains
    recursive subroutine matrix_get_eigen(matrix, val, vec, num_found, min_val, max_val)
        complex*16, intent(in)  :: matrix(:, :)
        real*8,     intent(in)  :: min_val, max_val
        complex*16, intent(out) :: vec(size(matrix(:, 1)), size(matrix(1, :)))
        real*8,     intent(out) :: val(size(matrix(:, 1)))
        integer,    intent(out) :: num_found
        complex*16              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                 :: info
        integer,    allocatable :: isuppz(:)
        complex*16              :: work(lwork)
        real*8                  :: rwork(lrwork)
        integer                 :: iwork(liwork)
        copy_matrix = matrix
        if (lwork < 2 * size(matrix(:, 1)) .or. lrwork < 24 * size(matrix(:, 1)) .or. liwork < 10 * size(matrix(:, 1))) then
            lwork  = -1
            lrwork = -1
            liwork = -1
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), min_val, max_val, 1, 1, 0d0, &
                num_found, val, vec, size(matrix(:, 1)), isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
            lwork  = int(work(1))
            lrwork = int(rwork(1))
            liwork = int(iwork(1))
            call matrix_get_eigen(matrix, val, vec, num_found, min_val, max_val)
        else
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), min_val, max_val, 1, 1, 0d0, &
                num_found, val, vec, size(matrix(:, 1)), isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
        end if
        
    end subroutine matrix_get_eigen

end module matrix

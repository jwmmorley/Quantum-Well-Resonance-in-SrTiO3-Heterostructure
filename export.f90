module export
    implicit none
    
    private :: export_to_string_integer, export_open_append, &
               export_data0_double   , export_data1_double   , export_data2_double, &
               export_data0_character, export_data1_character
    public  :: export_meta_data, export_data, &
               export_create_dir, export_to_string
    
    interface export_data
        module procedure export_data0_character, export_data1_character, &
                         export_data0_double   , export_data1_double   , export_data2_double
    end interface export_data
    interface export_to_string
        module procedure export_to_string_integer
    end interface export_to_string


contains
    subroutine export_meta_data(path, num_k_length, num_layers, num_k_path, num_energy, &
        broadening, crystal_length, max_crystal_length, &
        min_potential, max_potential, min_density, max_density, energy_min, energy_max, min_heatmap, max_heatmap)
        character(*), intent(in) :: path
        integer,      intent(in) :: num_k_length, num_layers, num_k_path, num_energy
        real*8,       intent(in) :: broadening, crystal_length, max_crystal_length
        real*8,       intent(in) :: min_potential, max_potential, min_density, max_density, &
                                    energy_min, energy_max, min_heatmap, max_heatmap
        integer                  :: file_number
        open(newunit = file_number, file = trim(path), status = "new", action = "write")
        write(file_number, fmt = "(4A16)")      "# Num K Len" , "Num Layers"     , "Num Path K"      , "Num Path Energy"
        write(file_number, fmt = "(4I16.1)")    num_k_length  , num_layers       , num_k_path        , num_energy
        write(file_number, fmt = "(3A16)")      "# Broadening", "Crystal Len (a)", "K Len [2pi*a]"
        write(file_number, fmt = "(3F16.8)")    broadening    , crystal_length   , max_crystal_length
        write(file_number, fmt = "(2A16)")      "# Min"       , "Max"
        write(file_number, fmt = "(2F16.8, A)") min_potential , max_potential    , " # Potential"
        write(file_number, fmt = "(2F16.8, A)") min_density   , max_density      , " # Density"
        write(file_number, fmt = "(2F16.8, A)") energy_min    , energy_max       , " # Heatmap y"
        write(file_number, fmt = "(2F16.8, A)") min_heatmap   , max_heatmap      , " # Heatmap c"
        close(file_number)
    
    end subroutine export_meta_data
    
    subroutine export_data0_character(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        format = "(A)"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        write(file_number, fmt = format) data
        close(file_number)
    
    end subroutine export_data0_character
    
    subroutine export_data1_character(path, data, fmt)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: data(:)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i
        format = "(A)"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data)
            write(file_number, fmt = format) data(i)
        end do
        close(file_number)
    
    end subroutine export_data1_character
    
    subroutine export_data0_double(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        format = "(F16.8)"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        write(file_number, fmt = format) data
        close(file_number)
    
    end subroutine export_data0_double
    
    subroutine export_data1_double(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data(:)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i
        format = "(F16.8)"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data)
            write(file_number, fmt = format) data(i)
        end do
        close(file_number)
    
    end subroutine export_data1_double
    
    subroutine export_data2_double(path, data, fmt)
        character(*), intent(in)  :: path
        real*8,       intent(in)  :: data(:, :)
        character(*), optional    :: fmt
        character(:), allocatable :: format
        integer                   :: file_number
        integer                   :: i, j
        format = "(F16.8)"
        if (present(fmt)) then
            format = fmt
        end if
        call export_open_append(path, file_number)
        do i = 1, size(data(:, 1))
            do j = 1, size(data(1, :))
                write(file_number, fmt = format, advance = "no") data(i, j)
            end do
            write(file_number, *)
        end do
        close(file_number)
    
    end subroutine export_data2_double
    
    subroutine export_open_append(path, file_number)
        character(*), intent(in)  :: path
        integer,      intent(out) :: file_number
        logical                   :: exists
        inquire(file = trim(path), exist = exists)
        if (exists) then
            open(newunit = file_number, file = trim(path), status = "old", position = "append", action = "write")
        else
            open(newunit = file_number, file = trim(path), status = "new", action = "write")
        end if
    
    end subroutine export_open_append
    
    function export_create_dir(path, prefix) result(dir)
        character(*), intent(in)  :: path
        character(*), intent(in)  :: prefix
        integer                   :: file_label
        logical                   :: directory_exists
        character(:), allocatable :: dir
        inquire(file = path//prefix//"/.", exist = directory_exists)
        if (.not. directory_exists) then
            dir = path//prefix//"/"
            call execute_command_line("mkdir -p """//dir//"""")
        else
            file_label       = 0
            do while (directory_exists)
                file_label = file_label + 1
                inquire(file = path//prefix//export_to_string(file_label)//"/.", exist = directory_exists)
            end do 
            dir = path//prefix//export_to_string(file_label)//"/"
            call execute_command_line("mkdir -p """//dir//"""")
        end if
        
    
    end function export_create_dir
    
    function export_to_string_integer(value) result(string)
        integer,      intent(in)  :: value
        character(256)            :: buffer
        character(:), allocatable :: string
        write(buffer, "(I6.1)") value
        buffer = adjustl(buffer)
        string = trim(buffer)
    
    end function export_to_string_integer



end module export

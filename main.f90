program main
    use timer
    use cpu
    use import
    use route
    use potential
    use transform
    use bulk
    use matrix
    use spectral
    use density
    use export
    implicit none
    
    integer,        parameter   :: num_variation      = 40
    integer,        parameter   :: num_layers         = 125
    
    integer                     :: well_1_start       = 30
    integer                     :: well_1_stop        = 50
    real*8                      :: well_1_start_depth = -0.18d0
    real*8                      :: well_1_stop_depth  = -0.18d0
    
    integer                     :: well_2_start       = 50
    integer                     :: well_2_stop        = 55
    real*8                      :: well_2_start_depth = -0.03d0 
    real*8                      :: well_2_stop_depth  = -0.03d0
    
    integer,        parameter   :: num_k_length       = 50
    integer,        parameter   :: num_energy         = 50
    real*8,         parameter   :: energy_min         = -0.375d0
    real*8,         parameter   :: energy_max         = 0d0
    real*8,         parameter   :: length_scale       = 0.25d0
    real*8,         parameter   :: broadening         = 0.005d0
    real*8,         parameter   :: zheevr_tolerance   = 0.1d0
    
    character(256), parameter   :: input_file         = "SrTiO3_hr.dat"
    character(256), parameter   :: out_dir            = "./out/"
    character(256), parameter   :: all_denity_file    = "all_density.dat"
    character(256), parameter   :: meta_file          = "meta.dat"
    character(256), parameter   :: heatmap_file       = "heatmap.dat"
    character(256), parameter   :: potential_file     = "potential.dat"
    character(256), parameter   :: density_file       = "density.dat"
    real*8,         parameter   :: crystal_length     = 3.905d-10
    real*8,         parameter   :: pi                 = 3.141592653589793d0
    
    complex*16,     allocatable :: hr(:, :, :, :, :)
    integer,        allocatable :: hrw(:, :, :)
    integer                     :: max_hopping, num_bands
    
    real*8,         allocatable :: kx(:), ky(:)
    integer,        allocatable :: kw(:)
    integer,        allocatable :: kp(:)
    real*8                      :: lattice(3, 3)   = reshape((/ 1d0, 0d0, 0d0, &
                                                                0d0, 1d0, 0d0, &
                                                                0d0, 0d0, 1d0 /), &
                                                                shape(lattice), order = (/ 2, 1 /))
    real*8                      :: positions(5, 3) = reshape((/ 0d0  , 0d0  , 0d0  , &
                                                                0.5d0, 0.5d0, 0.5d0, &
                                                                0.5d0, 0.5d0, 0d0  , &
                                                                0.5d0, 0d0  , 0.5d0, &
                                                                0d0  , 0.5d0, 0.5d0 /), &
                                                               shape(positions), order = (/ 2, 1 /))
    integer                     :: atom_types(size(positions(:, 1))) = (/ 1, 2, 3 ,3, 3/)
    
    real*8                      :: pot(num_layers)
    
    integer                     :: k_start, k_stop
    
    complex*16,     allocatable :: hb(:, :)
    
    real*8,         allocatable :: spec_slice(:, :, :)
    real*8,         allocatable :: spec(:, :, :)
    
    real*8,         allocatable :: energy(:), weight(:, :)
    complex*16,     allocatable :: eigen_vec(:, :)
    integer                     :: num_found
    
    real*8                      :: cbm
    
    real*8                      :: energy_range(num_energy)
    
    real*8,         allocatable :: heatmap(:, :)
    real*8                      :: den(num_layers)
    
    character(256)              :: path, variation_path
    
    real*8                      :: output_density
    
    integer                     :: i, j, k, l, m

    call cpu_start()
    if (cpu_is_master()) then
        call timer_start()
        print *, "Innit..."
    end if
    
    ! Extract Data
    if (cpu_is_master()) then
        call import_data(input_file, hr, hrw, max_hopping, num_bands)
    end if
    call cpu_broadcast(max_hopping, 1)
    call cpu_broadcast(num_bands, 1)
    if (.not. cpu_is_master()) then
        allocate(hr(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1, num_bands, num_bands))
        allocate(hrw(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1))
    end if
    call cpu_broadcast(hr, size(hr))
    call cpu_broadcast(hrw, size(hrw))
    
    ! Build Route
    call route_build(kx, ky, kw, kp, num_k_length, length_scale, &
        lattice, positions, atom_types)
    kx = kx * 2 * pi
    ky = ky * 2 * pi
    
    ! Allocate Dynamic Arrays
    allocate(energy(num_layers * num_bands))
    allocate(eigen_vec(num_layers * num_bands, num_layers * num_bands))
    allocate(weight(num_layers, num_layers * num_bands))
    allocate(spec_slice(size(kw), num_layers, num_energy))
    allocate(spec(size(kw), num_layers, num_energy))
    allocate(heatmap(size(kp), num_energy))
    
    ! Determine Conductiong Band Minimum
    hb = bulk_build(transform_r_to_kz(hr, hrw, 0d0, 0d0), num_layers)
    call matrix_get_eigen(hb, energy, eigen_vec, num_found, -100d0, 100d0)
    cbm = minval(energy(:num_found))
    
    ! Set Up Energy Range
    energy_range = route_range_double(energy_min, energy_max, num_energy)
    
    ! Set up Output Folder
    if (cpu_is_master()) then
        path = export_create_dir(out_dir)
        call export_data(trim(path)//trim(all_denity_file), &
            "#          Total  Quantum Well 1  Quantum Well 2")
    end if
    
    ! Analyse
    if (cpu_is_master()) then
        call timer_split()
        print *, "Begin Analysis..."
    end if
    do i = 1, num_variation
        pot        = 0d0
        spec_slice = 0d0
        spec       = 0d0
        
        ! Manipulate Wells
        well_2_start = well_2_start + 1
        well_2_stop  = well_2_stop + 1
        
        ! Set Up Potential
        call potential_add_well(pot, well_1_start, well_1_stop, well_1_start_depth, well_1_stop_depth)
        call potential_add_well(pot, well_2_start, well_2_stop, well_2_start_depth, well_2_stop_depth)
        
        ! Determine Spectral Function
        call cpu_split_work(k_start, k_stop, size(kw))
        do k = k_start, k_stop
            hb = bulk_build(transform_r_to_kz(hr, hrw, kx(k), ky(k)), num_layers)
            call bulk_add_potential(hb, pot)
            call matrix_get_eigen(hb, energy, eigen_vec, num_found, &
                energy_min + cbm - zheevr_tolerance, energy_max + cbm + zheevr_tolerance)
            energy = energy - cbm
            weight = 0d0
            do l = 1, num_layers
                do j = 1, num_found
                    do m = 1, num_bands
                        weight(l, j) = weight(l, j) + abs(eigen_vec((l - 1) * num_bands + m, j))**2
                    end do
                end do
            end do
            do l = 1, num_layers
                do j = 1, num_energy
                    spec_slice(k, l, j) = spectral_function(energy(:num_found), dble(weight(l, :num_found)), &
                        energy_range(j), broadening)
                end do
            end do
        end do
        call cpu_sum(spec_slice, spec)
        
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            ! Determine Heatmap
            do k = 1, size(kp)
                do j = 1, num_energy
                    heatmap(k, j) = sum(spec(kp(k), :, j))
                end do
            end do
            where (heatmap == 0d0)
                heatmap = minval(heatmap, mask = heatmap /= 0d0)
            end where
            heatmap = log(heatmap)
            
            ! Determine Density
            do l = 1, num_layers
                den(l) = density_function(spec(:, l, :), kw, crystal_length, energy_min, energy_max, num_energy)
            end do
            den = den * 1d-18
            
            ! Export Data
            variation_path = export_create_dir(path)
            call export_meta_data(trim(variation_path)//trim(meta_file), &
                num_k_length, num_layers, broadening, crystal_length, length_scale, &
                minval(pot), maxval(pot), minval(den), maxval(den), energy_min, energy_max, minval(heatmap), maxval(heatmap))
            call export_data(trim(variation_path)//trim(heatmap_file), transpose(heatmap))
            call export_data(trim(variation_path)//trim(potential_file), pot)
            call export_data(trim(variation_path)//trim(density_file), den)
            call export_data(trim(path)//trim(all_denity_file), reshape( &
                (/ sum(den), sum(den(well_1_start:well_1_stop)), sum(den(well_2_start:well_2_stop)) /), (/ 1, 3 /)))
            
            ! Produce Graphs
            call execute_command_line("gnuplot -c heatmap.p "//trim(variation_path))
            call execute_command_line("gnuplot -c density.p "//trim(variation_path))
            if (any(pot /= 0d0)) then
                call execute_command_line("gnuplot -c potential.p "//trim(variation_path))
            end if
            if (i > 2) then
                call execute_command_line("gnuplot -c total_density.p "//trim(path)// &
                    " Seperation[a] "//export_to_string(well_2_start - well_1_stop - i + 1))
                call execute_command_line("gnuplot -c reservoir_density.p "//trim(path)// &
                    " Seperation[a] "//export_to_string(well_2_start - well_1_stop - i + 1))
                call execute_command_line("gnuplot -c transport_density.p "//trim(path)// &
                    " Seperation[a] "//export_to_string(well_2_start - well_1_stop - i + 1))
                call execute_command_line("gnuplot -c reservoir_plus_transport_density.p "//trim(path)// &
                    " Seperation[a] "//export_to_string(well_2_start - well_1_stop - i + 1))
            end if
            
            write(*, fmt = "(A10, I4.1, A9)") "Variation ", i, " complete"
            call timer_split()
        end if
    end do

    call cpu_stop()
    

end program main
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
    use graph
    implicit none
    
    integer,        parameter   :: num_variation      = 3
    integer,        parameter   :: num_layers         = 40
    
    integer                     :: well_1_start       = 1!30
    integer                     :: well_1_stop        = 1!66
    real*8                      :: well_1_start_depth = 0d0!-0.22d0
    real*8                      :: well_1_stop_depth  = 0d0!-0.22d0
    
    integer                     :: well_2_start       = 1!18!66
    integer                     :: well_2_stop        = 10!21!66 + 4
    real*8                      :: well_2_start_depth = -0.28d0!-0.03d0 
    real*8                      :: well_2_stop_depth  = 0d0!-0.2d0!0d0!-0.03d0
    
    integer,        parameter   :: num_k_length       = 50
    integer,        parameter   :: num_energy         = 50
    integer,        Parameter   :: max_num_base_bands = 3
    real*8,         parameter   :: energy_min         = -0.3d0
    real*8,         parameter   :: energy_max         = 0d0
    real*8,         parameter   :: length_scale       = 0.5d0
    real*8,         parameter   :: broadening         = 0.0025d0
    real*8,         parameter   :: zheevr_tolerance   = 1d0
    
    integer                     :: x_offset           = 0
    character(*),   parameter   :: x_label            = "Reservoir and Transport Seperation"
    
    character(*),   parameter   :: input_file         = "SrTiO3_hr.dat"
    character(*),   parameter   :: out_dir            = "./out/"
    character(*),   parameter   :: all_denity_file    = "all_density"
    character(*),   parameter   :: meta_file          = "meta"
    character(*),   parameter   :: heatmap_file       = "heatmap"
    character(*),   parameter   :: potential_file     = "potential"
    character(*),   parameter   :: density_file       = "density"
    real*8,         parameter   :: crystal_length     = 3.905d-10
    real*8,         parameter   :: pi                 = 3.141592653589793d0
    
    complex*16,     allocatable :: hr(:, :, :, :, :)
    integer,        allocatable :: hrw(:, :, :)
    integer                     :: max_hopping, num_bands, num_base_bands
    
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
    
    real*8,         allocatable :: spec_slice(:, :, :, :, :)
    real*8,         allocatable :: spec(:, :, :, :, :)
    
    real*8,         allocatable :: energy(:), weight(:, :, :)
    complex*16,     allocatable :: eigen_vec(:, :)
    integer                     :: num_found
    
    real*8                      :: cbm
    
    real*8                      :: energy_range(num_energy)
    
    real*8,         allocatable :: heatmap(:, :, :, :)
    real*8,         allocatable :: den(:, :, :)
    
    character(256)              :: path, variation_path
    
    real*8                      :: output_density
    
    integer                     :: i, j, k, l, m, e

    call cpu_start()
    if (cpu_is_master()) then
        call timer_start()
        print *, "Initilise..."
    end if
    
    ! Extract Data
    if (cpu_is_master()) then
        call import_data(input_file, hr, hrw, max_hopping, num_base_bands)
    end if
    call cpu_broadcast(max_hopping, 1)
    call cpu_broadcast(num_base_bands, 1)
    num_bands = num_base_bands * num_layers
    if (.not. cpu_is_master()) then
        allocate(hr(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1, num_base_bands, num_base_bands))
        allocate(hrw(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1))
        hr = dcmplx(0d0, 0d0)
        hrw = 0
    end if
    call cpu_broadcast(hr, size(hr))
    call cpu_broadcast(hrw, size(hrw))
    
    ! Build Route
    call route_build(kx, ky, kw, kp, num_k_length, length_scale, &
        lattice, positions, atom_types)
    kx = kx * 2 * pi
    ky = ky * 2 * pi
    
    ! Allocate Dynamic Arrays
    allocate(energy(num_bands))
    allocate(eigen_vec(num_bands, num_bands))
    allocate(weight(num_layers, num_base_bands, max_num_base_bands))
    allocate(spec_slice(size(kw), num_layers, num_base_bands, max_num_base_bands, num_energy))
    allocate(spec(size(kw), num_layers, num_base_bands, max_num_base_bands, num_energy))
    allocate(heatmap(size(kp), num_base_bands, max_num_base_bands, num_energy))
    allocate(den(num_layers, num_base_bands, max_num_base_bands))
    
    ! Determine Conductiong Band Minimum
    hb = bulk_build(transform_r_to_kz(hr, hrw, 0d0, 0d0), num_layers)
    energy = 0d0
    eigen_vec = dcmplx(0d0, 0d0)
    call matrix_get_eigen(hb, energy, eigen_vec, num_found, -100d0, 100d0)
    cbm = minval(energy(:num_found))
    
    ! Set Up Energy Range
    energy_range = route_range_double(energy_min, energy_max, num_energy)
    
    ! Set up Output Folder
    if (cpu_is_master()) then
        path = export_create_dir(out_dir)
        call export_data(trim(path)//trim(all_denity_file)//".dat", &
            "#          Total  Quantum Well 1  Quantum Well 2             Sum"// &
            "       Avg Total      Avg Well 1      Avg Well 2    Avg Well Sum")
    end if
    
    ! Analyse
    if (cpu_is_master()) then
        call timer_split()
        print *, "Begin Analysis..."
    end if
    do i = 1, num_variation
        ! Manipulate Wells
        !well_2_start = well_2_start + 1
        !well_2_stop  = well_2_start + 2
        
        ! Set Up Potential
        pot        = 0d0
        call potential_add_well(pot, well_1_start, well_1_stop, well_1_start_depth, well_1_stop_depth)
        call potential_add_well(pot, well_2_start, well_2_stop, well_2_start_depth, well_2_stop_depth)
        
        ! Determine Spectral Function
        spec_slice = 0d0
        call cpu_split_work(k_start, k_stop, size(kw))
        do k = k_start, k_stop
            hb = bulk_build(transform_r_to_kz(hr, hrw, kx(k), ky(k)), num_layers)
            call bulk_add_potential(hb, pot)
            energy = 0d0
            eigen_vec = dcmplx(0d0, 0d0)
            num_found = 0
            call matrix_get_eigen(hb, energy, eigen_vec, num_found, &
                energy_min + cbm - zheevr_tolerance, energy_max + cbm + zheevr_tolerance)
            energy(:num_found) = energy(:num_found) - cbm
            weight = 0d0
            do l = 1, num_bands
                do j = 1, num_base_bands
                    do m = 1, num_base_bands
                        if (m <= max_num_base_bands) then
                            weight(l, j, m) = weight(l, j, m) + abs(eigen_vec((l - 1) * num_base_bands + j, m))**2
                        end if
                    end do
                end do
            end do
            do l = 1, num_layers
                do j = 1, num_base_bands
                    do m = 1, num_found
                        if (m <= max_num_base_bands) then
                            do e = 1, num_energy
                                spec_slice(k, l, j, m, e) = spectral_function(energy(m), dble(weight(l, j, m)), &
                                    energy_range(e), broadening)
                            end do
                        end if
                    end do
                end do
            end do
        end do
        spec       = 0d0
        call cpu_sum(spec_slice, spec)
        
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            ! Determine Density
            den = 0d0
            do l = 1, num_layers
                do j = 1, num_base_bands
                    do m = 1, max_num_base_bands
                        den(l, j, m) = density_function(spec(:, l, j, m, :), &
                            kw, crystal_length, energy_min, energy_max, num_energy)
                    end do
                end do
            end do
            den = den * 1d-18
            
            ! Determine Heatmap
            heatmap = 0d0
            do k = 1, size(kp)
                do j = 1, num_base_bands
                    do m = 1, max_num_base_bands
                        do e = 1, num_energy
                            heatmap(k, j, m, e) = sum(spec(kp(k), :, j, m, e))
                        end do
                    end do
                end do
            end do
            do j = 1, num_base_bands
                do m = 1, max_num_base_bands
                    where (heatmap(:, j, m, :) == 0d0)
                        heatmap(:, j, m, :) = minval(heatmap(:, j, m, :), mask = heatmap(:, j, m, :) /= 0d0)
                    end where
                end do
            end do
            
            ! Export Meta Data
            variation_path = export_create_dir(path)
            call export_meta_data(trim(variation_path)//trim(meta_file)//".dat", &
                num_k_length, num_layers, broadening, crystal_length, length_scale, &
                minval(pot), maxval(pot), minval(den), maxval(den), &
                energy_min, energy_max, minval(log(heatmap)), maxval(log(heatmap)))
                
            ! Export Density Data
            call export_data(trim(variation_path)//trim(potential_file)//".dat", pot)
            call export_data(trim(variation_path)//trim(density_file)//".dat", sum(sum(den, 2), 2))
            do j = 1, num_base_bands
                do m = 1, max_num_base_bands
                    call export_data(trim(variation_path)//export_to_string(j)//trim(density_file)//export_to_string(m)//".dat", &
                        den(:, j, m))
                end do
            end do
            call export_data(trim(path)//trim(all_denity_file)//".dat", reshape((/ &
                sum(den), &
                sum(den(well_1_start:well_1_stop, :, :)), &
                sum(den(well_2_start:well_2_stop, :, :)), &
                sum(den(well_1_start:well_1_stop, :, :)) + &
                    sum(den(well_2_start:well_2_stop, :, :)), &
                sum(den) / size(den), &
                sum(den(well_1_start:well_1_stop, :, :)) / (abs(well_1_stop - well_1_start) + 1), &
                sum(den(well_2_start:well_2_stop, :, :)) / (abs(well_2_stop - well_2_start) + 1), &
                sum(den(well_1_start:well_1_stop, :, :)) + &
                    sum(den(well_2_start:well_2_stop, :, :)) / (abs(well_1_stop + well_2_stop - well_1_start - well_2_start) + 2) &
                /), (/ 1, 8 /)))
                
            ! Export Heatmap Data
            call export_data(trim(variation_path)//trim(heatmap_file)//".dat", &
                    transpose(log(sum(sum(heatmap(:, :, :, :), 2), 2))))
            do j = 1, num_base_bands
                do m = 1, max_num_base_bands
                    call export_data(trim(variation_path)//export_to_string(j)//trim(heatmap_file)//export_to_string(m)//".dat", &
                        transpose(log(heatmap(:, j, m, :))))
                end do
            end do
            
            ! Produce Graphs
            call graph_heatmap_plot(heatmap_file, heatmap_file, trim(variation_path))
            do j = 1, num_base_bands
                do m = 1, max_num_base_bands
                    call graph_heatmap_plot(export_to_string(j)//heatmap_file//export_to_string(m), &
                        export_to_string(j)//heatmap_file//export_to_string(m), trim(variation_path))
                end do
            end do
            call graph_basic_plot(density_file, density_file, &
                1, "Layer", "Carrier Density [nm^{-2}]", 1, variation_path)
            do j = 1, num_base_bands
                do m = 1, max_num_base_bands
                    call graph_basic_plot(export_to_string(j)//density_file//export_to_string(m), &
                        export_to_string(j)//density_file//export_to_string(m), &
                        1, "Layer", "Carrier Density [nm^{-2}]", 1, variation_path)
                end do
            end do
            call graph_basic_plot(potential_file, "potential", &
                1, "Layer", "Potential [eV]", 1, variation_path)
            if (i > 2) then
                call graph_basic_plot(all_denity_file, "Total_Density", &
                    1, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "reservoir_density", &
                    2, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "transport_density", &
                    3, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "reservoir_plus_transport_density", &
                    4, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "Average_Total_Density", &
                    5, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "Average_reservoir_density", &
                    6, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "Average_transport_density", &
                    7, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
                call graph_basic_plot(all_denity_file, "Average_reservoir_plus_transport_density", &
                    8, x_label, "Carrier Density [nm^{-2}]", x_offset, path)
            end if
            
            write(*, fmt = "(A10, I4.1, A9)") "Variation ", i, " complete"
            call timer_split()
        end if
    end do

    call cpu_stop()

end program main

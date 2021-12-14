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
    
    integer,        parameter   :: num_variation      = 25
    integer,        parameter   :: num_layers         = 200
    
    integer                     :: well_1_start       = 30
    integer                     :: well_1_stop        = 66
    real*8                      :: well_1_start_depth = -0.22d0
    real*8                      :: well_1_stop_depth  = -0.22d0
    
    integer                     :: well_2_start       = 66
    integer                     :: well_2_stop        = 66 + 4
    real*8                      :: well_2_start_depth = -0.03d0 
    real*8                      :: well_2_stop_depth  = -0.03d0
    
    integer,        parameter   :: num_k_length       = 80
    integer,        parameter   :: num_energy         = 80
    real*8,         parameter   :: energy_min         = -0.3d0
    real*8,         parameter   :: energy_max         = 0d0
    real*8,         parameter   :: length_scale       = 0.5d0
    real*8,         parameter   :: broadening         = 0.0025d0
    real*8,         parameter   :: zheevr_tolerance   = 0.1d0
    
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
    
    real*8,         allocatable :: spec_slice(:, :, :, :)
    real*8,         allocatable :: spec(:, :, :, :)
    
    real*8,         allocatable :: energy(:), weight(:, :, :)
    complex*16,     allocatable :: eigen_vec(:, :)
    integer                     :: num_found
    
    real*8                      :: cbm
    
    real*8                      :: energy_range(num_energy)
    
    real*8,         allocatable :: heatmap(:, :, :)
    real*8,         allocatable :: den(:, :)
    
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
    allocate(weight(num_layers, num_bands, num_layers * num_bands))
    allocate(spec_slice(size(kw), num_layers, num_bands, num_energy))
    allocate(spec(size(kw), num_layers, num_bands, num_energy))
    allocate(heatmap(size(kp), num_bands, num_energy))
    allocate(den(num_layers, num_bands))
    
    ! Determine Conductiong Band Minimum
    hb = bulk_build(transform_r_to_kz(hr, hrw, 0d0, 0d0), num_layers)
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
        pot        = 0d0
        spec_slice = 0d0
        spec       = 0d0
        
        ! Manipulate Wells
        well_2_start = well_2_start + 1
        well_2_stop  = well_2_start + 2
        
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
                do j = 1, num_bands
                    do m = 1, num_found
                        weight(l, j, m) = weight(l, j, m) + abs(eigen_vec((l - 1) * num_bands + j, m))**2
                    end do
                end do
            end do
            do l = 1, num_layers
                do j = 1, num_bands
                    do e = 1, num_energy
                        spec_slice(k, l, j, e) = spectral_function(energy(:num_found), dble(weight(l, j, :num_found)), &
                            energy_range(e), broadening)
                    end do
                end do
            end do
        end do
        call cpu_sum(spec_slice, spec)
        
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            ! Determine Density
            do l = 1, num_layers
                do j = 1, num_bands
                    den(l, j) = density_function(spec(:, l, j, :), kw, crystal_length, energy_min, energy_max, num_energy)
                end do
            end do
            den = den * 1d-18
            
            ! Determine Heatmap
            do k = 1, size(kp)
                do j = 1, num_bands
                    do e = 1, num_energy
                        heatmap(k, j, e) = sum(spec(kp(k), :, j, e))
                    end do
                end do
            end do
            do j = 1, num_bands
                where (heatmap(:, j, :) == 0d0)
                    heatmap(:, j, :) = minval(heatmap(:, j, :), mask = heatmap(:, j, :) /= 0d0)
                end where
            end do
            
            ! Export Meta Data
            variation_path = export_create_dir(path)
            call export_meta_data(trim(variation_path)//trim(meta_file)//".dat", &
                num_k_length, num_layers, broadening, crystal_length, length_scale, &
                minval(pot), maxval(pot), minval(den), maxval(den), &
                energy_min, energy_max, minval(log(heatmap)), maxval(log(heatmap)))
                
            ! Export Density Data
            call export_data(trim(variation_path)//trim(potential_file)//".dat", pot)
            call export_data(trim(variation_path)//trim(density_file)//".dat", sum(den, 2))
            do j = 1, num_bands
                call export_data(trim(variation_path)//export_to_string(j)//trim(density_file)//".dat", &
                    den(:, j))
            end do
            call export_data(trim(path)//trim(all_denity_file)//".dat", reshape((/ &
                sum(den), &
                sum(den(well_1_start:well_1_stop, :)), &
                sum(den(well_2_start:well_2_stop, :)), &
                sum(den(well_1_start:well_1_stop, :)) + &
                    sum(den(well_2_start:well_2_stop, :)), &
                sum(den) / size(den), &
                sum(den(well_1_start:well_1_stop, :)) / (abs(well_1_stop - well_1_start) + 1), &
                sum(den(well_2_start:well_2_stop, :)) / (abs(well_2_stop - well_2_start) + 1), &
                sum(den(well_1_start:well_1_stop, :)) + &
                    sum(den(well_2_start:well_2_stop, :)) / (abs(well_1_stop + well_2_stop - well_1_start - well_2_start) + 2) &
                /), (/ 1, 8 /)))
                
            ! Export Heatmap Data
            call export_data(trim(variation_path)//trim(heatmap_file)//".dat", &
                    transpose(log(sum(heatmap(:, :, :), 2))))
            do j = 1, num_bands
                call export_data(trim(variation_path)//export_to_string(j)//trim(heatmap_file)//".dat", &
                    transpose(log(heatmap(:, j, :))))
            end do
            
            ! Produce Graphs
            call execute_command_line("gnuplot -e """ &
                    //"data_file='"//heatmap_file//"';" &
                    //"output_file='"//heatmap_file//"';" &
                    //"path='"//trim(variation_path)//"'" &
                    //""" basic_heatmap.p")
            do j = 1, num_bands
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//export_to_string(j)//heatmap_file//"';" &
                    //"output_file='"//export_to_string(j)//heatmap_file//"';" &
                    //"path='"//trim(variation_path)//"'" &
                    //""" basic_heatmap.p")
            end do
            call execute_command_line("gnuplot -e """ &
                //"data_file='"//density_file//"';" &
                //"output_file='density';" &
                //"data_column=1;" &
                //"x_label='Layer';" &
                //"y_label='Carrier Density [nm^{-2}]';" &
                //"x_offset=1;" &
                //"path='"//trim(variation_path)//"'" &
                //""" basic_plot.p")
            do j = 1, num_bands
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//export_to_string(j)//density_file//"';" &
                    //"output_file='"//export_to_string(j)//density_file//"';" &
                    //"data_column=1;" &
                    //"x_label='Layer';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset=1;" &
                    //"path='"//trim(variation_path)//"'" &
                    //""" basic_plot.p")
            end do
            call execute_command_line("gnuplot -e """ &
                //"data_file='"//potential_file//"';" &
                //"output_file='potential';" &
                //"data_column=1;" &
                //"x_label='Layer';" &
                //"y_label='Potential [eV]';" &
                //"x_offset=1;" &
                //"path='"//trim(variation_path)//"'" &
                //""" basic_plot.p")
            if (i > 2) then
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='Total_Density';" &
                    //"data_column=1;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='reservoir_density';" &
                    //"data_column=2;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='transport_density';" &
                    //"data_column=3;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='reservoir_plus_transport_density';" &
                    //"data_column=4;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='Average_Total_Density';" &
                    //"data_column=5;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='Average_reservoir_density';" &
                    //"data_column=6;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='Average_transport_density';" &
                    //"data_column=7;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
                call execute_command_line("gnuplot -e """ &
                    //"data_file='"//all_denity_file//"';" &
                    //"output_file='Average_reservoir_plus_transport_density';" &
                    //"data_column=8;" &
                    //"x_label='"//x_label//"';" &
                    //"y_label='Carrier Density [nm^{-2}]';" &
                    //"x_offset="//export_to_string(x_offset)//";" &
                    //"path='"//trim(path)//"'" &
                    //""" basic_plot.p")
            end if
            
            write(*, fmt = "(A10, I4.1, A9)") "Variation ", i, " complete"
            call timer_split()
        end if
    end do

    call cpu_stop()
    

end program main

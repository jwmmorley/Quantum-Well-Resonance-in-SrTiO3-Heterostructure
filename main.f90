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
    
    integer,        parameter   :: num_variation      = 36
    integer,        parameter   :: num_layers         = 60
    
    integer                     :: well_1_start       = 1!30
    integer                     :: well_1_stop        = 7!66
    real*8                      :: well_1_start_depth = -0.1d0!-0.22d0
    real*8                      :: well_1_stop_depth  = 0d0!-0.22d0
    
    integer                     :: well_2_start       = 10!18!66
    integer                     :: well_2_stop        = 12!21!66 + 4
    real*8                      :: well_2_start_depth = -0.28d0!-0.28d0!-0.03d0 
    real*8                      :: well_2_stop_depth  = -0.28d0!-0.2d0!0d0!-0.03d0
    
    integer,        parameter   :: num_k_length       = 256
    integer,        parameter   :: num_energy         = 256
    integer,        parameter   :: output_bands(3)    = (/ 1, 2, 3 /)
    integer,        parameter   :: output_states(*)   = (/ 1, 2, 3 /)
    real*8                      :: energy_min         = -1.3d0
    real*8                      :: energy_max         = 0d0
    real*8,         parameter   :: length_scale       = 1d0
    real*8,         parameter   :: broadening         = 0.0025d0
    
    real*8,         parameter   :: kcbmx              = 0d0
    real*8,         parameter   :: kcbmy              = 0d0
    
    integer                     :: x_offset           = 0
    character(*),   parameter   :: x_label            = "Surface Layer Doping [-0.1eV]"
    
    logical,        parameter   :: auto_energy_min    = .true.
    logical,        parameter   :: auto_energy_max    = .false.
    
    character(*),   parameter   :: input_file         = "SrTiO3_hr.dat"
    character(*),   parameter   :: out_dir            = "./out/"
    real*8,         parameter   :: crystal_length     = 3.905d-10
    real*8,         parameter   :: pi                 = 3.141592653589793d0
    
    real*8,         parameter   :: lattice(3, 3)      = reshape((/ 1d0, 0d0, 0d0, &
                                                                0d0, 1d0, 0d0, &
                                                                0d0, 0d0, 1d0 /), &
                                                                shape(lattice))
    real*8,         parameter   :: positions(3, 5)    = reshape((/ 0d0  , 0d0  , 0d0  , &
                                                                0.5d0, 0.5d0, 0.5d0, &
                                                                0.5d0, 0.5d0, 0d0  , &
                                                                0.5d0, 0d0  , 0.5d0, &
                                                                0d0  , 0.5d0, 0.5d0 /), &
                                                               shape(positions))
    integer,        parameter   :: atom_types(*)      = (/ 1, 2, 3 ,3, 3/)
    
    complex*16,     allocatable :: hr(:, :, :, :, :)
    integer,        allocatable :: hrw(:, :, :)
    
    integer                     :: max_hopping, num_states, max_num_states(num_variation), num_bands
    
    real*8,         allocatable :: kx(:), ky(:)
    integer,        allocatable :: kw(:)
    integer,        allocatable :: kp(:), kl(:, :)
    integer                     :: k_start, k_stop
    
    real*8                      :: pot(num_variation, num_layers)
    
    complex*16,     allocatable :: hb(:, :)
    real*8,         allocatable :: energy(:)
    complex*16,     allocatable :: weight(:, :)
    integer                     :: num_found
    real*8                      :: cbm
    
    real*8,         allocatable :: spec(:, :, :)
    real*8                      :: energy_range(num_variation, num_energy)
    
    real*8,         allocatable :: heatmap_bands_cpu(:, :, :), heatmap_states_cpu(:, :, :), heatmap_total_cpu(:, :), &
                                   heatmap_bands(:, :, :), heatmap_states(:, :, :), heatmap_total(:, :)
    real*8,         allocatable :: den_bands_cpu(:, :), den_states_cpu(:, :), den_total_cpu(:), &
                                   den_bands(:, :), den_states(:, :), den_total(:)
    
    character(256)              :: path, variation_dir, band_dir, energy_state_dir
    
    integer                     :: i, j, k, l, m, e, n

    call cpu_start()
    if (cpu_is_master()) then
        call timer_start()
        print *, "Initilise..."
    end if
    
    ! Set Up potentials
    pot = 0d0
    do i = 1, num_variation
        ! Manipulate Wells
        ! -------------------------------------------------------------------------------------------------------------------------
        well_1_start_depth = (i - 1) * -0.1d0
        ! -------------------------------------------------------------------------------------------------------------------------
        
        ! Set Up Potential
        call potential_add_well(pot(i, :), well_1_start, well_1_stop, well_1_start_depth, well_1_stop_depth)
        call potential_add_well(pot(i, :), well_2_start, well_2_stop, well_2_start_depth, well_2_stop_depth)
    end do
    
    ! Extract Data
    if (cpu_is_master()) then
        call import_data(input_file, hr, hrw, max_hopping, num_bands)
    end if
    call cpu_broadcast(max_hopping, 1)
    call cpu_broadcast(num_bands, 1)
    num_states = num_bands * num_layers
    if (.not. cpu_is_master()) then
        allocate(hr(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1, num_bands, num_bands))
        allocate(hrw(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1))
    end if
    call cpu_broadcast(hr, size(hr))
    call cpu_broadcast(hrw, size(hrw))
    
    ! Build Route
    if (cpu_is_master()) then
        call route_build(kx, ky, kw, kp, num_k_length, length_scale, &
            lattice, positions, atom_types)
        i = size(kw)
        j = size(kp)
    end if
    call cpu_broadcast(i, 1)
    call cpu_broadcast(j, 1)
    if (.not. cpu_is_master()) then
        allocate(kx(i))
        allocate(ky(i))
        allocate(kw(i))
        allocate(kp(j))
    end if
    call cpu_broadcast(kx, size(kx))
    call cpu_broadcast(ky, size(ky))
    call cpu_broadcast(kw, size(kw))
    call cpu_broadcast(kp, size(kp))
    
    ! Determine Conductiong Band Minimum
    allocate(energy(num_states))
    allocate(weight(num_states, num_states))
    hb = bulk_build(transform_r_to_kz(hr, hrw, kcbmx, kcbmy), num_layers)
    energy = 0d0
    weight = dcmplx(0d0, 0d0)
    call matrix_get_eigen(hb, energy, weight, num_found, 1, 1)
    cbm = energy(1)
    
    ! Determine Maximum number of bands
    do i = 1, num_variation
        hb = bulk_build(transform_r_to_kz(hr, hrw, kcbmx, kcbmy), num_layers)
        call bulk_add_potential(hb, pot(i, :))
        if (auto_energy_min) then
            energy = 0d0
            weight = dcmplx(0d0, 0d0)
            call matrix_get_eigen(hb, energy, weight, num_found, 1, 1)
            energy_min = (energy(1) - cbm) - sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
        end if
        if (auto_energy_max) then
            energy = 0d0
            weight = dcmplx(0d0, 0d0)
            call matrix_get_eigen(hb, energy, weight, num_found, (num_states - 1), num_states)
            energy_max = (energy(1) - cbm) + sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
        end if
        if (energy_min == energy_max) then
            if (auto_energy_min .and. .not. auto_energy_max) then
                energy_min = energy_min - sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
            else if (.not. auto_energy_min .and. auto_energy_max) then
                energy_max = energy_max + sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
            else
                energy_min = energy_min - sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
                energy_max = energy_max + sqrt((broadening / 10d0) - broadening**2) - sqrt(broadening - broadening**2)
            end if
        end if
        energy    = 0d0
        weight    = dcmplx(0d0, 0d0)
        num_found = 0
        call matrix_get_eigen(hb, energy, weight, num_found, energy_min + cbm, energy_max + cbm)
        max_num_states(i) = num_found
        if (max_num_states(i) == 0) then
            max_num_states(i) = 1
        end if
        ! Set Up Energy Range
        energy_range(i, :) = route_range_double(energy_min, energy_max, num_energy)
    end do
    
    ! Allocate Dynamic Arrays
    allocate(spec(num_bands, maxval(max_num_states), num_energy))
    allocate(den_bands_cpu(num_layers, size(output_bands)))
    allocate(den_states_cpu(num_layers, size(output_states)))
    allocate(den_total_cpu(num_layers))
    allocate(heatmap_bands_cpu(size(kp), size(output_bands), num_energy))
    allocate(heatmap_states_cpu(size(kp), size(output_states), num_energy))
    allocate(heatmap_total_cpu(size(kp), num_energy))
    if (cpu_is_master()) then
        allocate(den_bands(num_layers, size(output_bands)))
        allocate(den_states(num_layers, size(output_states)))
        allocate(den_total(num_layers))
        allocate(heatmap_bands(size(kp), size(output_bands), num_energy))
        allocate(heatmap_states(size(kp), size(output_states), num_energy))
        allocate(heatmap_total(size(kp), num_energy))
    end if
    
    ! Set up Output Folder
    if (cpu_is_master()) then
        path = export_create_dir(out_dir, "Run ")
        call export_data(trim(path)//"density.dat", &
            "#          Total  Quantum Well 1  Quantum Well 2             Sum"// &
            "       Avg Total      Avg Well 1      Avg Well 2    Avg Well Sum")
    end if
    
    ! Determine CPU Work
    call cpu_split_work(k_start, k_stop, size(kw))
    
    ! Map Total K to K Path
    allocate(kl(size(kw), sort_mode(kp)))
    kl = -1
    do k = k_start, k_stop
        i = 1
        do j = 1, size(kp)
            if (kp(j) == k) then
                kl(k, i) = j
                i = i + 1
            end if
        end do
    end do
    
    ! Analyse
    if (cpu_is_master()) then
        call timer_split()
        print *, "Begin Analysis..."
    end if
    do i = 1, num_variation
        den_bands_cpu      = 0d0
        den_states_cpu     = 0d0
        den_total_cpu      = 0d0
        heatmap_bands_cpu  = 0d0
        heatmap_states_cpu = 0d0
        heatmap_total_cpu  = 0d0
        if (cpu_is_master()) then
            den_bands      = 0d0
            den_states     = 0d0
            den_total      = 0d0
            heatmap_bands  = 0d0
            heatmap_states = 0d0
            heatmap_total  = 0d0
        end if
        
        ! Spectral into Density and Heatmap
        do k = k_start, k_stop
            hb = bulk_build(transform_r_to_kz(hr, hrw, kx(k), ky(k)), num_layers)
            call bulk_add_potential(hb, pot(i, :))
            energy = 0d0
            weight = dcmplx(0d0, 0d0)
            num_found = 0
            call matrix_get_eigen(hb, energy, weight, num_found, 1, max_num_states(i))
            energy(:num_found) = energy(:num_found) - cbm
            if (kl(k, 1) == -1) then
                do l = 1, num_layers
                    do j = 1, num_bands
                        do m = 1, num_found
                            do e = 1, num_energy
                                spec(j, m, e) = spectral_function( &
                                    energy(m), dble(abs(weight((l - 1) * num_bands + j, m))**2), energy_range(i, e), broadening)
                            end do
                        end do
                    end do
                    den_total_cpu(l) = den_total_cpu(l) + sum(spec(:, :num_found, :))
                    do j = 1, size(output_bands)
                        den_bands_cpu(l, j) = den_bands_cpu(l, j) + sum(spec(output_bands(j), :num_found, :))
                    end do
                    do m = 1, size(output_states)
                        den_states_cpu(l, m) = den_states_cpu(l, m) + sum(spec(:, output_states(m), :))
                    end do
                end do
            else
                do l = 1, num_layers
                    do j = 1, num_bands
                        do m = 1, num_found
                            do e = 1, num_energy
                                spec(j, m, e) = spectral_function( &
                                    energy(m), dble(abs(weight((l - 1) * num_bands + j, m))**2), energy_range(i, e), broadening)
                            end do
                        end do
                    end do
                    den_total_cpu(l) = den_total_cpu(l) + sum(spec(:, :num_found, :))
                    do n = 1, size(kl(1, :))
                        if (kl(k, n) == -1) then
                            exit
                        end if
                        heatmap_total_cpu(kl(k, n), :) = heatmap_total_cpu(kl(k, n), :) + sum(sum(spec(:, :num_found, :), 1), 1)
                    end do
                    do j = 1, size(output_bands)
                        den_bands_cpu(l, j) = den_bands_cpu(l, j) + sum(spec(output_bands(j), :num_found, :))
                        do n = 1, size(kl(1, :))
                            if (kl(k, n) == -1) then
                                exit
                            end if
                            heatmap_bands_cpu(kl(k, n), j, :) = heatmap_bands_cpu(kl(k, n), j, :) &
                                + sum(spec(output_bands(j), :num_found, :), 1)
                        end do
                    end do
                    do m = 1, size(output_states)
                        den_states_cpu(l, m) = den_states_cpu(l, m) + sum(spec(:, output_states(m), :))
                        do n = 1, size(kl(1, :))
                            if (kl(k, n) == -1) then
                                exit
                            end if
                            heatmap_states_cpu(kl(k, n), m, :) = heatmap_states_cpu(kl(k, n), m, :) &
                                + sum(spec(:, output_states(m), :), 1)
                        end do
                    end do
                end do
            end if
        end do
        den_bands_cpu  = den_bands_cpu * (1d0 / pi) &
                          * ((energy_range(i, 2) - energy_range(i, 1)) / dble(num_energy)) &
                          * (1d0 / crystal_length**2) &
                          * (1d0 / dble(sum(kw))) &
                          * (1d-18)
        den_states_cpu = den_states_cpu * (1d0 / pi) &
                          * ((energy_range(i, 2) - energy_range(i, 1)) / dble(num_energy)) &
                          * (1d0 / crystal_length**2) &
                          * (1d0 / dble(sum(kw))) &
                          * (1d-18)
        den_total_cpu  = den_total_cpu * (1d0 / pi) &
                          * ((energy_range(i, 2) - energy_range(i, 1)) / dble(num_energy)) &
                          * (1d0 / crystal_length**2) &
                          * (1d0 / dble(sum(kw))) &
                          * (1d-18)
        call cpu_sum(den_bands_cpu, den_bands)
        call cpu_sum(den_states_cpu, den_states)
        call cpu_sum(den_total_cpu, den_total)
        call cpu_sum(heatmap_bands_cpu, heatmap_bands)
        call cpu_sum(heatmap_states_cpu, heatmap_states)
        call cpu_sum(heatmap_total_cpu, heatmap_total)
        
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            ! Export Meta Data
            variation_dir = export_create_dir(trim(path), "Variation "//export_to_string(i))
            call export_meta_data(trim(variation_dir)//"meta.dat", &
                num_k_length, num_layers, size(kp), num_energy, &
                broadening, crystal_length, length_scale, &
                minval(pot(i, :)), maxval(pot(i, :)), minval(den_total), maxval(den_total), &
                minval(energy_range(i, :)), maxval(energy_range(i, :)), minval(log(heatmap_total)), maxval(log(heatmap_total)))
                
            ! Export and Plot Variation Data
            call export_data(trim(variation_dir)//"potential.dat",      pot(i, :))
            call export_data(trim(variation_dir)//"density.dat",        den_total)
            call export_data(trim(variation_dir)//"band structure.dat", transpose(log(heatmap_total)))
            call graph_basic_plot("potential", "potential", &
                1, "Layer", "Potential [eV]", 1, trim(variation_dir))
            call graph_basic_plot("density", "density", &
                1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(variation_dir))
            call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(variation_dir))
            
            ! Export and Plot Variation Band Data
            do j = 1, size(output_bands)
                band_dir = export_create_dir(trim(variation_dir)//"Bands/", "Band "//export_to_string(output_bands(j)))
                call export_data(trim(band_dir)//"density"//".dat", &
                    den_bands(:, j))
                call export_data(trim(band_dir)//"band structure"//".dat", &
                    transpose(log(heatmap_bands(:, j, :))))
                call graph_basic_plot("density", "density", &
                    1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(band_dir))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(band_dir))
            end do
            
            ! Export and Plot Variation Energy State Data
            do j = 1, size(output_states)
                energy_state_dir = export_create_dir(trim(variation_dir)//"Energy States/", &
                    "Energy State "//export_to_string(output_states(j)))
                call export_data(trim(energy_state_dir)//"density"//".dat", &
                    den_states(:, j))
                call export_data(trim(energy_state_dir)//"band structure"//".dat", &
                    transpose(log(heatmap_states(:, j, :))))
                call graph_basic_plot("density", "density", &
                    1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(energy_state_dir))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(energy_state_dir))
            end do
            
            ! Export and Plot Data
            call export_data(trim(path)//"density.dat", reshape((/ &
                sum(den_total), &
                sum(den_total(well_1_start:well_1_stop)), &
                sum(den_total(well_2_start:well_2_stop)), &
                sum(den_total(well_1_start:well_1_stop)) + &
                    sum(den_total(well_2_start:well_2_stop)), &
                sum(den_total) / size(den_total), &
                sum(den_total(well_1_start:well_1_stop)) / (abs(well_1_stop - well_1_start) + 1), &
                sum(den_total(well_2_start:well_2_stop)) / (abs(well_2_stop - well_2_start) + 1), &
                sum(den_total(well_1_start:well_1_stop)) + &
                    sum(den_total(well_2_start:well_2_stop)) / (abs(well_1_stop + well_2_stop - well_1_start - well_2_start) + 2) &
                /), (/ 1, 8 /)))
            if (i > 2) then
                call graph_basic_plot("density", "Total_Density", &
                    1, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "reservoir_density", &
                    2, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "transport_density", &
                    3, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "reservoir_plus_transport_density", &
                    4, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_Total_Density", &
                    5, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_reservoir_density", &
                    6, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_transport_density", &
                    7, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_reservoir_plus_transport_density", &
                    8, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
            end if
            
            write(*, fmt = "(A10, I4.1, A9)") "Variation ", i, " complete"
            call timer_split()
        end if
    end do
    
    ! Deallocate Arrays
    deallocate(hr)
    deallocate(hrw)
    deallocate(hb)
    deallocate(kx)
    deallocate(ky)
    deallocate(kw)
    deallocate(kp)
    deallocate(spec)
    deallocate(energy)
    deallocate(weight)
    deallocate(den_bands_cpu)
    deallocate(den_states_cpu)
    deallocate(den_total_cpu)
    deallocate(heatmap_bands_cpu)
    deallocate(heatmap_states_cpu)
    deallocate(heatmap_total_cpu)
    if (cpu_is_master()) then
        deallocate(den_bands)
        deallocate(den_states)
        deallocate(den_total)
        deallocate(heatmap_bands)
        deallocate(heatmap_states)
        deallocate(heatmap_total)
    end if
    
    call cpu_stop()
    
contains
    subroutine sort_bubble(array)
        integer, intent(inout) :: array(:)
        integer                :: temp
        integer                :: i, j
        logical                :: swapped
        do j = size(array) - 1, 1, -1
            swapped = .FALSE.
            do i = 1, j
                if (array(i) > array(i + 1)) then
                    temp       = array(i)
                    array(i)   = array(i + 1)
                    array(i+1) = temp
                    swapped    = .true.
                end if
            end do
            if (.not. swapped) then
                return
            end if
        end do
      
    end subroutine sort_bubble
    
    function sort_mode(array) result(val)
        integer, intent(in) :: array(:)
        integer             :: val
        integer             :: temp_array(size(array))
        integer             :: max_count, count, i
        temp_array = array
        call sort_bubble(temp_array) 
        val        = temp_array(1)
        max_count  = 1
        count      = 1
        do i = 2, size(temp_array)
            if (temp_array(i) == temp_array(i - 1)) then
                count = count + 1
            else
                if (count > max_count) then
                    max_count = count
                    val       = temp_array(i - 1)
                end if
                count = 1
            end if
        end do
        if (count > max_count) then
            val = temp_array(size(temp_array))
        end if
    
    end function sort_mode

end program main

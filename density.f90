module density
    implicit none
    
    public :: density_function

contains
    function density_function(spec, kw, crystal_length, energy_min, energy_max, num_energy) result(density)
        real*8,  intent(in) :: spec(:, :), crystal_length
        real*8,  intent(in) :: energy_min, energy_max
        integer, intent(in) :: kw(:)
        integer, intent(in) :: num_energy
        real*8              :: density
        integer             :: k
        density = 0d0
        do k  = 1, size(kw)
            density = density + dble(kw(k)) * sum(spec(k, :))
        end do
        density = density &
            * (1d0 / 3.14d0) &
            * ((energy_max - energy_min) / dble(num_energy)) &
            * (1d0 / crystal_length**2) &
            * (1d0 / dble(sum(kw)))
    
    end function density_function

end module density

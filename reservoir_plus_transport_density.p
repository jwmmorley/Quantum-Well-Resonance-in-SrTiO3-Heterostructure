path = ARG1
x_label = ARG2
x_offset = ARG3
cd path
set print "-"
data_file = "all_density.dat"

set terminal png size 1680,1050
set output "reservoir_plus_transport_density.png"
set key off
set tics scale 0

set title "Reservoir Plus Transport Density"

set xlabel x_label

set ylabel "n [nm-2]"
set ytics

plot data_file using (column(0) + x_offset):(column(2) + column(3)) with lines

path = ARG1
cd path
set print "-"
data_file = "density.dat"
meta_file = "meta.dat"

getValue(row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
num_col = system("awk 'NR==1{print NF}' ".data_file)
num_rows = system("awk 'END{ print NR }' ".data_file)

set terminal png size 1680,1050
set output "density.png"
set key off
set tics scale 0

set title "Density"

set xlabel "Layer"

set ylabel "n [nm-2]"
set yrange [getValue(7, 1, meta_file): getValue(7, 2, meta_file)]
set ytics

plot data_file using (column(0) + 1):1 with lines

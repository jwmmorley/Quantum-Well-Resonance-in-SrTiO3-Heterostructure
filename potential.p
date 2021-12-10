path = ARG1
cd path
set print "-"
data_file = "potential.dat"
meta_file = "meta.dat"

getValue(row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
num_col = system("awk 'NR==1{print NF}' ".data_file)
num_rows = system("awk 'END{ print NR }' ".data_file)

set terminal png size 1680,1050
set output "potential.png"
set key off

set xlabel "Layer"

set ylabel "E [ev]"
set yrange [getValue(6, 1, meta_file):getValue(6, 2, meta_file)]
set ytics

plot data_file using (column(0) + 1):1 with lines

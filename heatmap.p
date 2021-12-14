path = ARG1
cd path
set print "-"
data_file = "heatmap.dat"
meta_file = "meta.dat"

getValue(row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
num_bounds = system("awk 'NR==3{print NF}' ".meta_file)
num_col = system("awk 'NR==1{print NF}' ".data_file)
num_rows = system("awk 'END{ print NR }' ".data_file)

set terminal png size 1680,1050
set output "heatmap.png"
set palette rgbformulae 21, 22, 23
set key off
set tics front

set title "Band Structure"

set cbrange [0:getValue(9, 2, meta_file)] 

set xlabel "K [pi/a]"
set xtics (sprintf("%3.2f", -1.00 * getValue(4, 3, meta_file)) 0, \
    sprintf("%3.2f", -0.50 * getValue(4, 3, meta_file)) (num_col / 4), \
    "0.00" (num_col / 2), \
    sprintf("%3.2f", 0.50 * getValue(4, 3, meta_file)) ((3 * num_col) / 4), \
    sprintf("%3.2f", 1.00 * getValue(4, 3, meta_file)) (num_col - 1))
        
num_y_tics = 8
y_min = getValue(8, 1, meta_file)
y_max = getValue(8, 2, meta_file)
y_difference = y_max - y_min
y_zero = (num_rows - 1) * (-1 * y_min) / y_difference    
   
set ylabel "E - E_cbm [ev]"
set ytics ("0" y_zero)
do for [y = 1 : num_y_tics] {
    set ytics add (sprintf("%3.2f", y * y_difference / num_y_tics) (y_zero + y * (num_rows - 1) / num_y_tics))
    set ytics add (sprintf("%3.2f", -1 * y * y_difference / num_y_tics) (y_zero - y * (num_rows - 1) / num_y_tics))
}
plot data_file matrix with image

set print "-"

getValue(row, col, filename) = system("awk "."'{if (NR == '".row."') print $'".col."'}'"." '".filename."'")

num_row = getValue(2, 3, meta_file.".dat")
num_col = getValue(2, 4, meta_file.".dat")

set terminal png size 1680,1050;
set output path.output_file.".png"
set key off
set tics front
set palette rgbformulae 21, 22, 23

set title "Band Structure"

set xlabel "K [pi/a]"
set xtics (sprintf("%3.2f", -1.00 * getValue(4, 3, meta_file.".dat")) 0, \
    sprintf("%3.2f", -0.50 * getValue(4, 3, meta_file.".dat")) (num_col / 4), \
    "0.00" (num_col / 2), \
    sprintf("%3.2f", 0.50 * getValue(4, 3, meta_file.".dat")) ((3 * num_col) / 4), \
    sprintf("%3.2f", 1.00 * getValue(4, 3, meta_file.".dat")) (num_col - 1))
        
num_y_tics = 8
y_min = getValue(8, 1, meta_file.".dat")
y_max = getValue(8, 2, meta_file.".dat")
y_difference = y_max - y_min
y_zero = (num_row - 1) * (-1 * y_min) / y_difference    
   
set ylabel "E - E_cbm [ev]"
set ytics ("0" y_zero)
do for [y = 1 : num_y_tics] {
    set ytics add (sprintf("%3.2f", y * y_difference / num_y_tics) (y_zero + y * (num_row - 1) / num_y_tics))
    set ytics add (sprintf("%3.2f", -1 * y * y_difference / num_y_tics) (y_zero - y * (num_row - 1) / num_y_tics))
}
plot path.data_file.".dat" matrix with image

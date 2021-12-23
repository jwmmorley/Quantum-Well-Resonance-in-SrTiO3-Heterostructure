set print "-";

getValue(row, col, filename) = system("awk "."'{if (NR == '".row."') print $'".col."'}'"." '".filename."'");

num_row = getValue(2, 3, meta_file.".dat");
num_col = getValue(2, 4, meta_file.".dat");

set terminal png size 1680,1050;
set encoding utf8;
set output path.output_file.".png";
set key off;
set tics front;
set palette rgbformulae 21, 22, 23;

#set title "Band Structure";

unset cbtics
set colorbox user
set style line 1001 linecolor rgb "white"
set colorbox border 1001 
set colorbox origin 0.1, 0.2
set colorbox size 0.03, 0.33
set label 1 "Max" at screen 0.115, 0.54 centre textcolor rgb "white";
set label 2 "Min" at screen 0.115, 0.185 centre textcolor rgb "white";
set cbrange [getValue(9, 1, meta_file.".dat") to getValue(9, 2, meta_file.".dat")]

set label 3 "K_{M Γ X} [π/a]" at screen (2 - sqrt(2)), 0.02 centre;
set xtics (sprintf("%3.2f", sqrt(2) * getValue(4, 3, meta_file.".dat")) 0, \
    sprintf("%3.2f", 1.00 * getValue(4, 3, meta_file.".dat")) (num_col / 2) - ((sqrt(2) / 2) * (num_col / 2)), \
    sprintf("%3.2f", 0.50 * getValue(4, 3, meta_file.".dat")) (num_col / 2) - ((sqrt(2) / 4) * (num_col / 2)), \
    "0.00" (num_col / 2), \
    sprintf("%3.2f", 0.50 * getValue(4, 3, meta_file.".dat")) ((3 * num_col) / 4), \
    sprintf("%3.2f", 1.00 * getValue(4, 3, meta_file.".dat")) (num_col - 2));
        
num_y_tics = 4;
y_min = getValue(8, 1, meta_file.".dat") + 0; # + 0 Converts from string to number
y_max = getValue(8, 2, meta_file.".dat") + 0;
y_difference = y_max - y_min;
set format y "%4.3f"
if (y_max < 0) { \
    set ytics (sprintf("%4.3f", y_max) (num_row - 1));
    y_zero = num_row;
    y_zero_value = y_max;
};
if (0 < y_min){ \
    set ytics (sprintf("%4.3f", y_min) 0);
    y_zero = 0;
    y_zero_value = y_min;
};
if (y_min < 0 && 0 < y_max) { \
    set ytics ("0" y_zero);
    y_zero = (num_col - 1) * (-1 * y_min) / y_difference; 
    y_zero_value = 0;
};
do for [y = 1 : num_y_tics] {
    set ytics add (sprintf("%4.3f", y_zero_value + y * (y_difference + y_zero_value) / num_y_tics) \
        (y_zero + y * (num_col - 1) / num_y_tics));
    set ytics add (sprintf("%4.3f", y_zero_value + -1 * y * (y_difference + y_zero_value) / num_y_tics) \
        (y_zero - y * (num_col - 1) / num_y_tics));
}; 
   
set ylabel "E - E_{cbm} [ev]";


set multiplot layout 1, 2;
set size (2 - sqrt(2)), (1 - 0.03);
set origin 0, 0.03
set rmargin 0
plot [0:(num_row / 2)] path.data_file.".dat" matrix with image;

set size (sqrt(2) - 1), (1 - 0.03);
set lmargin 0
unset rmargin
set origin (2 - sqrt(2)), 0.03;
unset ylabel
set ytics font ",0"
unset colorbox
plot [((num_row / 2) + 1):num_row] path.data_file.".dat" matrix with image;
unset multiplot

cd path;
set print "-";

set terminal svg size 1680,1050 name output_file font "Helvetica,16" background rgb 'white';
set output output_file.".svg";
set key off;
set tics scale 0;

set xlabel x_label;

set ylabel "Carrier Density [nm^{-2}]";

plot data_file.".dat" using (column(0) + x_offset):data_column with lines linewidth 2 linecolor rgbcolor 'blue';

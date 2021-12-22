module graph
    use export
    implicit none

    public :: graph_basic_plot, graph_heatmap_plot

contains
    subroutine graph_basic_plot(data_file, output_file, data_column, x_label, y_label, x_offset, path)
        character(*), intent(in) :: data_file, output_file, x_label, y_label, path
        integer,      intent(in) :: data_column, x_offset
        call execute_command_line("gnuplot -e """ &
            //"data_file='"//data_file//"';" &
            //"output_file='"//output_file//"';" &
            //"data_column="//export_to_string(data_column)//";" &
            //"x_label='"//x_label//"';" &
            //"y_label='"//y_label//"';" &
            //"x_offset="//export_to_string(x_offset)//";" &
            //"path='"//trim(path)//"'" &
            //""" basic_plot.p")
    end subroutine graph_basic_plot
    
    subroutine graph_heatmap_plot(data_file, output_file, meta_file, path)
        character(*), intent(in) :: data_file, output_file, meta_file, path
        call execute_command_line("gnuplot -e """ &
            //"data_file='"//data_file//"';" &
            //"output_file='"//output_file//"';" &
            //"meta_file='"//meta_file//"';" &
            //"path='"//path//"'" &
            //""" basic_heatmap.p")
    end subroutine graph_heatmap_plot

end module graph

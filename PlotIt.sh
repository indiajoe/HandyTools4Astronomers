#!/bin/bash
#This is just quickly plot repetedly X vs Y of two column you copy paste into input in terminal.
#Usage: ./PlotIt.sh ["w lp"]
#....................... indiajoe@gmail.com
cat > plotit.gnuplotTEMP << EOF
plot 'plotit.dataTEMP' $1 title ''
pause -1
EOF
while :
do
    echo "Enter the X Y in column to plot"
    cat > plotit.dataTEMP 
    gnuplot plotit.gnuplotTEMP
done
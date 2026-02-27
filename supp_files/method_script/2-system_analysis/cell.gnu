set terminal pngcairo size 1600,1200 font "sans,35"
set encoding iso_8859_1
set output 'cell.png'
set border 31 lw 4
set key top center horizontal font " ,25"
set xlabel "Frame" 
set ylabel "XY or Z cell dimensions ({\305})"
#plot "rmsd.dat" u 1:2 w l lw 2 lc rgb "black" title "RMSD"
plot "cell.dat" u 1:2 every 5 w l lw 2 lc rgb "black" notitle, "" u 1:4 every 5 w l lw 2 lc rgb "black" notitle

quit


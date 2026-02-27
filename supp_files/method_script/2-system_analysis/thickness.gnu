set terminal pngcairo size 1600,1200 font "sans,35"
set encoding iso_8859_1
set output 'thickness.png'
set border 31 lw 4
set key top center horizontal font " ,25"
set xlabel "Frame" 
set ylabel "Membrane thickness ({\305})"
#plot "rmsd.dat" u 1:2 w l lw 2 lc rgb "black" title "RMSD"
plot "thickness.dat" u 1:2 every 5 w l lw 2 lc rgb "black" notitle

quit


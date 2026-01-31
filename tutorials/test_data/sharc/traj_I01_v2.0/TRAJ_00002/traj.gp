unset key
unset colorbox
set title "Energies in diagonal basis"
set xlabel "Time t in fs"
set ylabel "Energy in eV"
set cbrange [0:17]
set palette defined (0.0 "gray90", 1e-5 "gray60", 1e-4 "gray30", 1e-3 "orange", 1e-2 "red", 1e-1 "magenta", 1e-0 "blue", 10 "blue", 11 "green", 12 "red", 13 "turquoise", 14 "orange", 15 "cyan", 16 "brown", 17 "skyblue")

plot "output_data/expec.out" u 1:($4)               title "Total Energy" lw   0.50 lc rgbcolor "#000000" w l, \
""               u 1:($5):(abs($8)+10) title "State 1"     lw   4.50 pal w l, \
""               u 1:($5):(abs($11))    title "State 1"     lw   3.50 pal w l, \
""               u 1:($6):(abs($9)+10) title "State 2"     lw   4.50 pal w l, \
""               u 1:($6):(abs($12))    title "State 2"     lw   3.50 pal w l, \
""               u 1:($7):(abs($10)+10) title "State 3"     lw   4.50 pal w l, \
""               u 1:($7):(abs($13))    title "State 3"     lw   3.50 pal w l, \
""               u 1:($3)               title "Trajectory"   lw   1.00 lc rgbcolor "#000000" pt 6 w p

pause -1

set key
set yrange [0:1]
set ylabel "Wavefunction Amplitude"
set title "MCH Quantum Amplitudes"
plot "output_data/coeff_MCH.out" 	u 1:2  	title "Sum of Amplitudes" 	lw   3.00 	lc rgbcolor "#000000" 	w l, \
"" 			u 1:($3**2+$4**2)  	title "S 0, 0" 	lw   1.00 	lc rgbcolor "#FF0000" 	w l, \
"" 			u 1:($5**2+$6**2)  	title "S 1, 0" 	lw   1.00 	lc rgbcolor "#00FF66" 	w l, \
"" 			u 1:($7**2+$8**2)  	title "S 2, 0" 	lw   1.00 	lc rgbcolor "#3200FF" 	w l

pause -1

set yrange [0:1]
set ylabel "Wavefunction Amplitude"
set title "Diag Quantum Amplitudes"
plot "output_data/coeff_diag.out" 	u 1:2  	title "Sum of Amplitudes" 	lw   3.00 	lc rgbcolor "#000000" 	w l, \
"" 			u 1:($3**2+$4**2)  	title "State 1" 	lw   2.00 	lc rgbcolor "#FF0000" 	w l, \
"" 			u 1:($5**2+$6**2)  	title "State 2" 	lw   2.00 	lc rgbcolor "#00FF66" 	w l, \
"" 			u 1:($7**2+$8**2)  	title "State 3" 	lw   2.00 	lc rgbcolor "#3200FF" 	w l

pause -1

set xrange[0:*]
set ylabel "Hopping Probability"
set title "Diag Hopping Probablities and Random Number"
set style fill solid 0.25 border
plot "output_data/prob.out"	u 1:5  	title "State 3" 	lw   1.00 	lc rgbcolor "#3200FF" 	w boxes, \
""		u 1:4  	title "State 2" 	lw   1.00 	lc rgbcolor "#00FF66" 	w boxes, \
""		u 1:3  	title "State 1" 	lw   1.00 	lc rgbcolor "#FF0000" 	w boxes, \
""		u 1:($1>0 ? $2 : 1/0)	title "Random number" 			lc rgbcolor "black"	w lp

pause -1



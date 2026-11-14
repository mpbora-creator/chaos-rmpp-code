root=system("pwd")
hist=root."/data/control_m05/time_history.tsv"
out=root."/build/no_debris_M05_equilibrium_diagnostics_tikz.tex"

set terminal tikz standalone color size 7.2in,5.0in font "\\normalsize" \
    preamble "\\usepackage{amsmath}\\usepackage[T1]{fontenc}"
set output out
set datafile separator "\t"
set border lw 0.8
set grid back lc rgb "#dddddd" lw 0.5
set tics in nomirror scale 0.55
set xtics 200
set xrange [0:1000]
set key opaque samplen 1.8 spacing 1.0

set multiplot

# Geometry matches the reservoir-control figure exactly.
left_x0=0.1
left_x1=0.47
right_x0=0.58
right_x1=0.95
bottom_y0=0.14
bottom_y1=0.49
top_y0=0.55
top_y1=0.90

# Upper left: volume-averaged electron and ion densities.
set lmargin at screen left_x0
set rmargin at screen left_x1
set bmargin at screen top_y0
set tmargin at screen top_y1
set format x ""
unset xlabel
set yrange [0.9997:1.0003]
set ytics 0.0001
set format y "%g"
set ylabel '$\langle n_s\rangle/n_0$' offset 0.5,0
set key at graph 0.98,0.97 right top
set label 90 '(a)' at graph 0.02,0.96 left front textcolor rgb "black"
plot hist every 5 using 2:16 with lines lw 0.8 lc rgb "#0072B2" title '$\langle n_e\rangle$', \
     hist every 5 using 2:17 with lines lw 0.8 lc rgb "#D55E00" title '$\langle n_i\rangle$', \
     1 with lines dt 2 lw 0.6 lc rgb "#444444" notitle

# Upper right: local ion flow, with its late-time mean.
set lmargin at screen right_x0
set rmargin at screen right_x1
set bmargin at screen top_y0
set tmargin at screen top_y1
set format x ""
unset xlabel
set yrange [0.0108:0.0126]
set ytics 0.0003
set format y "%.4f"
set ylabel '$\langle v_i\rangle/v_{te0}$' offset 0.5,0
set key at graph 0.98,0.97 right top
set label 90 '(b)' at graph 0.02,0.96 left front textcolor rgb "black"
plot hist every 5 using 2:23 with lines lw 0.8 lc rgb "#D55E00" title 'ion mean', \
     0.011699474126171385 with lines dt 2 lw 0.6 lc rgb "#444444" title 'late-time mean'

# Lower left: electrostatic field energy.
set lmargin at screen left_x0
set rmargin at screen left_x1
set bmargin at screen bottom_y0
set tmargin at screen bottom_y1
set format x "%g"
set xlabel '$t\omega_{pe}$' offset 0,0.4
unset logscale y
set yrange [0:0.10]
set ytics 0.02
set format y "%.2f"
set ylabel '$W_E$' offset 0.5,0
unset y2tics
unset y2label
unset key
set label 90 '(c)' at graph 0.02,0.96 left front textcolor rgb "black"
plot hist every 5 using 2:9 with lines lw 0.8 lc rgb "#009E73" notitle

# Lower right: total-energy variation relative to the initial value.
set lmargin at screen right_x0
set rmargin at screen right_x1
set bmargin at screen bottom_y0
set tmargin at screen bottom_y1
unset logscale y
unset y2tics
unset y2label
set format x "%g"
set xlabel '$t\omega_{pe}$' offset 0,0.4
set yrange [-0.4:0.6]
set ytics (-0.4,-0.2,0.0,0.2,0.4,0.6)
set format y "%.1f"
set ylabel '$100[\mathcal{E}(t)/\mathcal{E}(0)-1]$' offset 0.5,0
unset key
set label 90 '(d)' at graph 0.02,0.96 left front textcolor rgb "black"
plot hist every 5 using 2:(100.0*($10/220.31618311932633-1.0)) with lines lw 0.8 lc rgb "#6A3D9A" notitle, \
     0 with lines dt 2 lw 0.6 lc rgb "#444444" notitle

unset multiplot
unset output

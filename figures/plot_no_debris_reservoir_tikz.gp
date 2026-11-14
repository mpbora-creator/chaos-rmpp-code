root=system("pwd")
hist=root."/data/control_m05/time_history.tsv"
reg=root."/data/control_m05/time_history_density_regulator.tsv"
out=root."/build/no_debris_M05_reservoir_balance_tikz.tex"

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

# Common, exactly aligned plot-box geometry.
left_x0=0.09
left_x1=0.47
right_x0=0.57
right_x1=0.95
bottom_y0=0.14
bottom_y1=0.49
top_y0=0.55
top_y1=0.90

# Upper left: retain x ticks but suppress their labels.
set lmargin at screen left_x0
set rmargin at screen left_x1
set bmargin at screen top_y0
set tmargin at screen top_y1
set format x ""
unset xlabel
set yrange [-1.5:1.5]
set ytics 0.5
set ylabel '$\Delta N_s/N_{\mathrm{tot}}\;(10^{-4})$' offset 0.5,0
set key at graph 0.98,0.97 right top
set label 90 '(a)' at graph 0.02,0.96 left front textcolor rgb "black"
plot hist using 2:(1e4*($3-429779)/859558.0) ev 10 with lines lw 0.8 lc rgb "#0072B2" title '$(N_e-N_{e0})/N_{\mathrm{tot}}$', \
     hist using 2:(1e4*($4-429779)/859558.0) ev 10 with lines lw 0.8 lc rgb "#D55E00" title '$(N_i-N_{i0})/N_{\mathrm{tot}}$'

# Upper right: identical plot-box size and vertical placement.
set lmargin at screen right_x0
set rmargin at screen right_x1
set bmargin at screen top_y0
set tmargin at screen top_y1
set format x ""
unset xlabel
set yrange [-1.5:1.5]
set ytics 0.5
set ylabel '$(N_e-N_i)/N_{\mathrm{tot}}\;(10^{-4})$' offset 0.5,0
unset key
set label 90 '(b)' at graph 0.02,0.96 left front textcolor rgb "black"
plot hist using 2:(1e4*($3-$4)/859558.0) ev 10 with lines lw 0.8 lc rgb "#CC3311" notitle, \
     0 with lines dt 2 lw 0.6 lc rgb "#444444" notitle

# Lower left: moved upward and aligned exactly below the upper-left box.
set lmargin at screen left_x0
set rmargin at screen left_x1
set bmargin at screen bottom_y0
set tmargin at screen bottom_y1
set format x "%g"
set xlabel '$t\omega_{pe}$' offset 0,0.4
set yrange [0:6]
set ytics 1
set ylabel '$2N_{\mathrm{pair}}/N_{\mathrm{tot}}\;(10^{-5})$' offset 0.5,0
set key at graph 0.98,0.97 right top
set label 90 '(c)' at graph 0.02,0.96 left front textcolor rgb "black"
plot reg using 2:(1e5*2.0*$5/859558.0) ev 10 with lines lw 0.8 lc rgb "#6A3D9A" title 'pair deficit', \
     reg using 2:(1e5*2.0*$8/859558.0) ev 10 with lines lw 0.8 lc rgb "#777777" title 'fractional remainder'

# Lower right: identical plot-box size and placement to lower left.
set lmargin at screen right_x0
set rmargin at screen right_x1
set bmargin at screen bottom_y0
set tmargin at screen bottom_y1
set format x "%g"
set xlabel '$t\omega_{pe}$' offset 0,0.4
set yrange [0:10]
set ytics 2
set ylabel '$N_{\mathrm{pair}}^{\mathrm{cum}}$' offset 0.5,0
unset key
set label 90 '(d)' at graph 0.02,0.96 left front textcolor rgb "black"
plot reg using 2:7 with lines lw 0.8 lc rgb "#009E73" notitle

unset multiplot
unset output

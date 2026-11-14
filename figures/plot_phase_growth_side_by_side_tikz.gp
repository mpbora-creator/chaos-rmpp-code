root = system("pwd")
growth = root."/data/derived/iicsi_growth_history_M05.tsv"
out = root."/build/phase_growth_side_by_side_M05_tikz.tex"
raster = "ion_phase_scatter_naked_M05_t1400.png"

set terminal tikz standalone color size 8.4in,3.4in font "\\normalsize" \
    preamble "\\usepackage{amsmath}\\usepackage[T1]{fontenc}\\usepackage{graphicx}"
set output out
set datafile separator "\t"
set border lw 0.8
set grid back lc rgb "#dddddd" lw 0.5
set tics in nomirror scale 0.55
unset key

set multiplot

# Left panel: empty TikZ axes and color bar.  The transparent PNG is inserted
# as a single LaTeX graphic at the exact data-rectangle dimensions.
set lmargin at screen 0.06
set rmargin at screen 0.38
set bmargin at screen 0.18
set tmargin at screen 0.92
set xrange [0:350]
set yrange [-0.08:0.08]
set xtics 50
set ytics 0.02
set format x "%g"
set format y "%.2f"
set xlabel '$x/\lambda_{De}$' offset 0,0.4
set ylabel '$v_i/v_{te0}$' offset 0.5,0
unset grid
set cbrange [-0.08:0.08]
set cbtics 0.02
set format cb "%.2f"
set cblabel '$v_i/v_{te0}$' offset 0.8,0
set palette defined (-0.08 "#440154", -0.04 "#3B528B", 0.0 "#21918C", 0.04 "#5EC962", 0.08 "#FDE725")
set colorbox user origin 0.40,0.18 size 0.025,0.74
set label 10 sprintf('\raisebox{-1.258in}{\includegraphics[width=2.688in,height=2.516in]{%s}}',raster) \
    at screen 0.22,0.55 center front
set arrow 10 from 105,-0.08 to 105,0.08 nohead dt 2 lw 0.8 lc rgb "#CC3311" front
set label 90 '(a)' at graph 0.02,0.96 left front textcolor rgb "black"
# The invalid palette point creates the color scale without emitting points.
plot '+' using 1:(1/0):(0) with points palette notitle

unset label 10
unset label 90
unset arrow 10
unset colorbox

# Right panel: effective IICSI growth rate over the extended record.
set lmargin at screen 0.60
set rmargin at screen 0.965
set bmargin at screen 0.18
set tmargin at screen 0.92
set grid back lc rgb "#dddddd" lw 0.5
set xrange [0:5000]
set yrange [-0.02:0.02]
set xtics 1000
set ytics 0.004
set format x "%g"
set format y "%.3f"
set xlabel '$t\omega_{pe}$' offset 0,0.4
set ylabel '$\gamma_{\mathrm{eff}}/\omega_{pe}$' offset 0.5,0
set arrow 1 from 1220,graph 0 to 1220,graph 1 nohead dt 3 lw 0.7 lc rgb "#777777" front
set arrow 2 from 1460,graph 0 to 1460,graph 1 nohead dt 3 lw 0.7 lc rgb "#777777" front
set label 90 '(b)' at graph 0.02,0.96 left front textcolor rgb "black"
plot growth using 1:3 with lines lw 1.0 lc rgb "#D55E00" notitle, \
     0 with lines dt 2 lw 0.7 lc rgb "#444444" notitle

unset multiplot
unset output

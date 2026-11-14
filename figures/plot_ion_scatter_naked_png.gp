root = system("pwd")
ions = root."/data/derived/ion_particles_M05_t1400.tsv"
out = root."/build/ion_phase_scatter_naked_M05_t1400.png"

# Aspect ratio exactly matches the left TikZ data rectangle.
set terminal pngcairo transparent enhanced size 2535,2375
set output out
set datafile separator "\t"
unset border
unset tics
unset xlabel
unset ylabel
unset key
unset colorbox
unset object
set lmargin at screen 0
set rmargin at screen 1
set bmargin at screen 0
set tmargin at screen 1
set xrange [0:350]
set yrange [-0.08:0.08]
set cbrange [-0.08:0.08]
set palette defined (-0.08 "#440154", -0.04 "#3B528B", 0.0 "#21918C", 0.04 "#5EC962", 0.08 "#FDE725")
plot ions using 1:2:2 with dots palette notitle

unset output

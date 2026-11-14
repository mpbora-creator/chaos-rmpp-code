#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p build

gnuplot figures/plot_no_debris_equilibrium_tikz.gp
gnuplot figures/plot_ion_scatter_naked_png.gp
gnuplot figures/plot_phase_growth_side_by_side_tikz.gp
gnuplot figures/plot_no_debris_reservoir_tikz.gp

for tex in \
  no_debris_M05_equilibrium_diagnostics_tikz.tex \
  phase_growth_side_by_side_M05_tikz.tex \
  no_debris_M05_reservoir_balance_tikz.tex
do
  (cd build && pdflatex -interaction=nonstopmode -halt-on-error "$tex" >/dev/null)
done

printf 'Built paper figures in %s/build\n' "$repo_root"

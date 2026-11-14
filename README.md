# Code and figure data for *A generic route to chaos in self-excited plasma oscillations*

This repository is the archival software release accompanying the paper
*A generic route to chaos in self-excited plasma oscillations*. It contains
the 1D–1V electrostatic particle-in-cell (PIC) solver, the exact parameter
files used for the reported `M=0.5` debris and debris-free control cases,
the analysis used to obtain the time-dependent IICSI growth rate, and the
source data and Gnuplot/TikZ scripts for the principal PIC figures.

## Repository contents

- `src/`: Julia PIC implementation, open reservoirs, debris charging,
  photoemission, secondary-electron emission, and diagnostics.
- `config/paper_m05_debris.jl`: archived debris-run parameters.
- `config/paper_m05_no_debris.jl`: matched debris-free control parameters.
- `scripts/run_pic.jl`: simulation launcher.
- `scripts/extract_iicsi_growth_history.jl`: reconstruction of the coherent
  ion-density-mode envelope and effective growth rate.
- `figures/`: portable Gnuplot/TikZ figure sources.
- `data/control_m05/`: debris-free time histories used in Figures 3 and 5.
- `data/derived/`: compact derived data used in Figure 4.
- `data/theory/`: compact theoretical trajectory data.

The large restart checkpoints and raw Vlasov archives are not committed to
GitHub because individual files exceed GitHub's normal file-size limit. The
compact data included here reproduce the deposited figures. Raw archives may
be deposited separately as a Zenodo dataset if full re-analysis from particle
snapshots is required.

## Requirements

- Julia 1.9 or later (the release was tested with Julia 1.12)
- Gnuplot with the `tikz` and `pngcairo` terminals
- A LaTeX installation containing TikZ, `amsmath`, and `graphicx`
- Poppler only if PNG previews are desired

Install the Julia environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick verification

```bash
julia --project=. scripts/quick_check.jl
julia --project=. scripts/quick_check_open.jl
```

## Run the paper configurations

Debris case:

```bash
julia --project=. scripts/run_pic.jl --parameters-file config/paper_m05_debris.jl
```

Matched debris-free control:

```bash
julia --project=. scripts/run_pic.jl --parameters-file config/paper_m05_no_debris.jl
```

These are computationally substantial runs. Output is written below
`diagnostics/` and is excluded from version control.

## Reproduce the deposited PIC figures

From the repository root:

```bash
./scripts/reproduce_figures.sh
```

The command creates `build/` and compiles:

- `no_debris_M05_equilibrium_diagnostics_tikz.pdf` (Figure 3),
- `phase_growth_side_by_side_M05_tikz.pdf` (Figure 4), and
- `no_debris_M05_reservoir_balance_tikz.pdf` (Figure 5).

Figure 4 uses a high-resolution transparent PNG containing only the ion
macroparticles, embedded in a TikZ-generated frame. This avoids TeX-memory
failure while retaining LaTeX fonts, axes, tick marks, labels, and colour bar.

## Recompute the IICSI growth history

The included `data/derived/iicsi_growth_history_M05.tsv` is ready for plotting.
To recompute it from a raw serialized Vlasov archive:

```bash
CHAOS_VLASOV_FILE=/path/to/vlasov_snapshots.jls \
  julia --project=. scripts/extract_iicsi_growth_history.jl
```

The effective rate is

```math
\gamma_{\mathrm{eff}}(t)=\frac{d}{dt}\ln A(t),
```

where `A(t)` is the analytic envelope of the leading coherent ion-density
mode in the band `0.015 <= omega/omega_pe <= 0.05`.

## Reproducibility notes

- The reported random seed is `42`.
- The open-grid run uses `dt*omega_pe = 0.05`, `L/lambda_De = 350`, and
  domain-scaled `N_x = 2634`.
- Restart files preserve particles, fields, RNG streams, reservoir state,
  debris state, and diagnostic accumulators.
- `Manifest.toml` records the complete Julia dependency resolution.

## Citation

Use the metadata in `CITATION.cff`. After the first GitHub release is archived
by Zenodo, replace the DOI placeholder in the manuscript's Code Availability
statement with the DOI displayed by Zenodo.

## License

The software is released under the MIT License. The compact numerical data in
`data/` may be reused under the same terms with citation of the associated
paper and software DOI.

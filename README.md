# ENSO-conditioned evolution of global mean surface temperature

Everything needed — and nothing more — to reproduce the figures and tables of:

> Tippett, M. K. and Becker, E. J.: *ENSO-conditioned evolution of global mean
> surface temperature.* Preprint: <https://eartharxiv.org/repository/view/13005/>

**What `make all` regenerates:** all 12 figures and 6 tables used by the
manuscript, from the frozen data in `data/` alone.

## Headline forecast: February 2027 GMST (unverified — the date hasn't happened yet)

The paper's central forward-looking result, shown together in
`figures/forecast.pdf` (manuscript Fig. 3e):

- **NMME multi-model forecast** (June 2026 initialization): global mean
  surface temperature peaks in **February 2027 at about 1.84°C** above the
  1850–1900 baseline.
- **This paper's ENSO-conditioned model**, using only the observed recent
  12-month GMST level (+1.30°C, June 2025–May 2026) and the NMME-forecast
  December 2026 Niño-3.4 value (+3.3°C) as inputs: **February 2027 peak of
  about 1.68°C**, 95% prediction interval +1.42°C to +1.94°C.

Both numbers are genuine predictions of a date that has not yet occurred —
rerunning this capsule once observations catch up to February 2027 will show
how they verified.

## Browse the code and figures (no clone required)

Executed notebooks with every figure inline (retina resolution), viewable
straight on GitHub:

- [`notebooks/global_temperature_enso-prediction.ipynb`](notebooks/global_temperature_enso-prediction.ipynb)
- [`notebooks/nmme_comparison.ipynb`](notebooks/nmme_comparison.ipynb)

These are a viewing convenience, **not** something `make all` regenerates:
the authoritative source is `scripts/*.py`, and `make all` is what regenerates
the manuscript outputs. Refresh them with `make notebooks` (see
`scripts/build_notebooks.py`).

## Quick start

```bash
conda env create -f environment.yml    # or: mamba env create -f environment.yml
conda activate enso-gmt-capsule

make all          # 12 figures -> figures/ , 6 tables -> tables/
make manuscript   # compile manuscript/manuscript.pdf from those outputs
make verify       # rebuild and check against the committed reference outputs
```

Make targets: `all` (figures + tables), `figures`, `tables`, `manuscript`,
`notebooks` (committed inline-figure notebooks under `notebooks/`), `verify`,
`clean`. Set `PYTHON=...` to point at a specific interpreter.

## Environment

`environment.yml` is a version-pinned environment covering the analysis stack
(numpy, pandas, xarray, scipy, statsmodels, netCDF4, matplotlib), the
jupytext/nbconvert notebook build, and `poppler` (for the `pdftoppm` figure
comparison in `make verify`). No cartopy or gridded fields are required.

## Repository layout

```
data/         frozen input data (Zenodo v2) + SHA256SUMS
scripts/      analysis in jupytext percent format + verify.py
figures/      the manuscript's 12 figures    (committed reference outputs)
tables/       the manuscript's 6 LaTeX tables (committed reference outputs)
manuscript/   self-contained LaTeX source + cited-only bibliography
notebooks/    executed notebooks, figures inline (built by `make notebooks`,
              generated from scripts/, never hand-edited; committed for browsing)
environment.yml  Makefile
```

## Reproduction map

Two scripts produce every output; each writes both figures and tables.

### `scripts/global_temperature_enso-prediction.py` — 9 figures, 4 tables

| Output | Kind |
|---|---|
| `eofs_pcs.pdf` | figure |
| `reconstruction_2pc.pdf` | figure |
| `pc_correlations.pdf` | figure |
| `fit_all.pdf` | figure |
| `variance_explained_dc.pdf` | figure |
| `simplified_pred.pdf` | figure |
| `corr_dec_n34_by_month.pdf` | figure |
| `fit_idea2.pdf` | figure |
| `forecast.pdf` | figure |
| `regression_table_full.tex` | table |
| `regression_table_r1.tex` | table |
| `regression_table_dc2.tex` | table |
| `regression_table_dc2_roni.tex` | table |

### `scripts/nmme_comparison.py` — 3 figures, 2 tables

| Output | Kind |
|---|---|
| `nmme_skill_comparison.pdf` | figure |
| `nmme_trajectory_examples.pdf` | figure |
| `nmme_traj_dc_pca.pdf` | figure |
| `nmme_trends_table.tex` | table |
| `nmme_significance_table.tex` | table |

The scripts are jupytext percent-format `.py` files (cells + markdown
narration). `make notebooks` renders and executes them into
`notebooks/*.ipynb` with figures inline; the `.ipynb` files are build products
and are never hand-edited.

## Verification

`make verify` rebuilds every figure and table into a scratch tree (`.verify/`)
and compares against the committed reference copies in `figures/` and `tables/`:

- **tables** must be character-identical (`diff`-clean);
- **figures** are rendered with `pdftoppm` and compared pixel-by-pixel.

On the machine and BLAS/LAPACK build that generated the reference outputs the
match is exact (18/18 outputs). On other hardware, **visual identity** is the
standard (sub-percent pixel shifts in the sub-leading modes are possible — see
below); `verify` reports each output as `EXACT`, `WITHIN-TOL`, or `FAIL`.

The figure PDFs are written with `savefig(..., metadata={'CreationDate': None})`
so their bytes carry no build timestamp: `make all` is byte-reproducible and
leaves a clean `git status` on any machine that reproduces the pixels exactly.

### Known limitation — `regression_table_full.tex` (float32 SVD, sub-leading modes)

The archive stores the observed Niño-3.4 and GMST series as **float32** (the
analysis's native in-memory dtype — the two are bit-identical). The PCA
therefore runs in float32. The 4-PC benchmark table `regression_table_full.tex`
uses the low-variance **GMST-PC2** mode, whose singular-value gap is small
(s₂→s₃ ≈ 0.12; variance 3.1 % → 1.3 %). A float32 singular vector on such a gap
is sensitive at the ~10⁻³ level to the exact LAPACK/BLAS implementation, so two
cells of this table can differ by one unit in the last printed digit across
environments (e.g. a coefficient of `-0.11899` rounding to `-0.119` vs `-0.118`,
and its adjusted R²). This is a numerical-library artifact, **not** a data or
porting difference:

- the archived float32 inputs are exact representations of the underlying data
  (max |diff| = 0);
- the other three regression tables and all figures that depend only on the
  dominant leading modes reproduce exactly;
- only PC2-derived quantities move (this table, plus sub-2 % pixel shifts in the
  PC2 curves and correlation-value annotations of `eofs_pcs.pdf`,
  `reconstruction_2pc.pdf`, `simplified_pred.pdf`, and `pc_correlations.pdf`).

The analysis deliberately keeps the float32 SVD (its native in-memory dtype);
casting to float64 would improve conditioning but change the results.
`make verify` treats this table under the visual-identity standard (reported as
`WARN`, not `FAIL`, if it moves).

## Data

`data/` holds byte-identical copies of the Zenodo archive
(v2, DOI [10.5281/zenodo.21432181](https://doi.org/10.5281/zenodo.21432181);
concept DOI [10.5281/zenodo.20514350](https://doi.org/10.5281/zenodo.20514350)):

- `obs_n34_GMSTa_185001-202605.nc` — observed Niño-3.4 SST (from Jan 1979) and
  global mean surface temperature anomaly (from Jan 1850), monthly through May 2026
- `nmme_n34a_GMSTa_199106-202606.nc` — NMME forecasts, June starts 1991–2026:
  per-model ensemble means and members, observed Niño-3.4 on the forecast grid,
  per-model metadata
- `Rnino34.ascii.txt` — NOAA CPC relative Niño-3.4 (RONI) index
  ([source](https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt),
  retrieved 2026-04-03)
- `SHA256SUMS` — checksums; run `sha256sum -c data/SHA256SUMS` to check integrity

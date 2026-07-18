# Reproduction capsule — ENSO-conditioned evolution of global mean surface temperature

**Status: under construction (Phase 1 of 6).** See `CLAUDE.md` for the build
plan and `docs/phase0_dependency_trace.md` for the code/data dependency map.

This repository will contain exactly what is needed — and no more — to
reproduce the figures and tables of:

> Tippett, M. K.: *ENSO-conditioned evolution of global mean surface
> temperature.* Preprint: <https://eartharxiv.org/repository/view/13005/>

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
- `make_zenodo_archive.py` — provenance: how the netCDFs were derived from
  ERSSTv5, NOAA GlobalTemp v6.1.0, and NMME (not runnable without those sources)
- `SHA256SUMS` — checksums of the above

## Planned layout

```
data/        frozen input data (Zenodo v2)          scripts/    analysis (percent-format)
figures/     the manuscript's 12 figures            tables/     the 6 LaTeX tables
manuscript/  LaTeX source + trimmed bibliography    notebooks/  executed walkthroughs (built)
docs/        dependency trace, build notes
```

## Reproduction and verification

Outputs are verified against the source-project reference copies: tables must
`diff` character-identical, figures must render (`pdftoppm`) and pixel-compare
clean. These gates are exact on the machine and environment that generated the
reference artifacts; on other hardware or BLAS/LAPACK builds, **visual identity**
is the standard.

### Known limitation — `regression_table_full.tex` (float32 SVD, sub-leading modes)

The archive stores the observed Niño-3.4 and GMST series as **float32** (matching
the source project's in-memory dtype — the two are bit-identical). The PCA
therefore runs in float32. The 4-PC benchmark table `regression_table_full.tex`
uses the low-variance **GMST-PC2** mode, whose singular-value gap is small
(s₂→s₃ ≈ 0.12; variance 3.1 % → 1.3 %). A float32 singular vector on such a gap
is sensitive at the ~10⁻³ level to the exact LAPACK/BLAS implementation, so two
cells of this table can differ by one unit in the last printed digit across
environments (e.g. a coefficient of `-0.11899` rounding to `-0.119` vs `-0.118`,
and its adjusted R²). This is a numerical-library artifact, **not** a data or
porting difference:

- the archived inputs are bit-identical (max |diff| = 0) to the source raw data;
- the other three regression tables and all figures that depend only on the
  dominant leading modes reproduce exactly;
- only PC2-derived quantities move (this table, plus sub-2 % pixel shifts in the
  PC2 curves and correlation-value annotations of `eofs_pcs.pdf`,
  `reconstruction_2pc.pdf`, `simplified_pred.pdf`, and `pc_correlations.pdf`).

The port deliberately keeps the float32 SVD to match the source code exactly;
casting to float64 would improve conditioning but diverge from the reference.
This table is accepted under the visual-identity standard.
```

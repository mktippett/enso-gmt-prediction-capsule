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

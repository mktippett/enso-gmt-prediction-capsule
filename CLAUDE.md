# CLAUDE.md — enso-gmt-prediction-capsule

Reproduction capsule for the manuscript "ENSO-conditioned evolution of global
mean surface temperature" (EarthArXiv: https://eartharxiv.org/repository/view/13005/).
Contract: `make all` regenerates every figure and table used by the manuscript
from the data in `data/` alone.

## Ground rules

1. **The source project `/Users/tippett/claude/enso-gmt-prediction/` is READ-ONLY.**
   It is the verification oracle. Never write, delete, or move anything there.
   All surgery happens on copies inside this repo; all outputs go to this repo's
   `figures/` and `tables/`.
2. **Verification gates** (definition of done for the script ports):
   - the 6 tables must `diff` character-identical against the oracle's `tex/` copies;
   - the 12 figures must image-diff clean (`pdftoppm` render + pixel compare)
     against the oracle's `plots/` copies.
   Gates are exact on this machine; on other hardware/BLAS, visual identity is
   the standard (SVD last-digit differences possible).
3. Scripts are written in jupytext percent format (`# %%` cells + markdown
   narration); `notebooks/*.ipynb` are build products (jupytext → execute),
   never hand-edited.
4. Data files in `data/` are frozen: byte-identical to Zenodo v2
   (DOI 10.5281/zenodo.21432181; concept DOI 10.5281/zenodo.20514350;
   verified via `data/SHA256SUMS` and Zenodo md5 on 2026-07-18).
   `data/make_zenodo_archive.py` is provenance documentation (how the netCDFs
   were derived from ERSSTv5 / NOAA GlobalTemp / NMME pickles); it is not
   runnable without the raw private data.

## Plan and current state

Phases (full context: `docs/phase0_dependency_trace.md`):

- [x] Phase 0 — dependency trace (see docs/)
- [x] Phase 1 — scaffold, data copies + checksums (this commit)
- [ ] Phase 2 — port the two scripts (copy → retarget to `data/` → trim to the
      mapped outputs → percent-format narration). Trim strictly by the trace's
      KEEP/CUT lists; mind the "trim traps" section.
- [ ] Phase 3 — run + verify against the oracle (gates above)
- [ ] Phase 4 — manuscript copy: relativize `\includegraphics`/`\input` paths
      (Edit tool, never sed on .tex), extract cited-only bib entries from the
      `.aux`, compile, PDF-diff against an oracle compile
- [ ] Phase 5 — environment lock (from `pangeo-2025`, pruned), Makefile
      (`all`/`figures`/`tables`/`manuscript`/`notebooks`/`verify`), README with
      the reproduction map, notebook build
- [ ] Phase 6 — clean-machine test (fresh clone + locked env), GitHub publish,
      Zenodo-GitHub integration for a capsule DOI, optional Binder badge

## Source locations (oracle, read-only)

- Scripts: `../enso-gmt-prediction/scripts/global_temperature_enso-prediction.py`
  (1957 lines), `nmme_comparison.py` (1289 lines). Neither imports `config.py`.
- Reference outputs: `../enso-gmt-prediction/plots/*.pdf`, `../enso-gmt-prediction/tex/*.tex`
- Manuscript: `~/tex/papers/enso-global-mean-temperature/prediction/revision/manuscript-revision1-v2.tex`
  (12 figures + 6 tables, absolute paths); bibliography `~/tex/bib/all.bib`
  (cites the archive as `ENSO-conditioned-data`, concept DOI).

## Environment

- Python: `mamba run -n pangeo-2025 python <script>.py` until Phase 5 produces
  the capsule's own lock file.
- The ported scripts need numpy, pandas, xarray, netCDF4, scipy, statsmodels,
  matplotlib. NOT needed (dropped with cut blocks): cartopy, gridded ERSSTv5,
  NOAA gridded file, PMM.txt.

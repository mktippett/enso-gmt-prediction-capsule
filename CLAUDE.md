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
- [x] Phase 2 — ported both scripts to `scripts/*.py` (jupytext percent):
      retargeted to `data/`, trimmed by the trace's KEEP/CUT lists, trim traps
      kept. `nmme_comparison.py` aliases the one archive NMME file to both
      original pickles (`init_time→S`, `member→M`, `lead→L`, `obs_n34→obs`).
- [x] Phase 3 — ran + verified against the oracle. 5/6 tables diff
      character-identical; all 12 figures within the visual-identity standard
      (most pixel-identical). **Known exception:** `regression_table_full.tex`
      differs by one last-digit unit in two GMST-PC2 cells — the archive stores
      the obs series as float32 (bit-identical to the oracle raw data), so the
      SVD runs in float32 and the sub-leading GMST-PC2 mode (small s₂→s₃ gap) is
      LAPACK/BLAS-sensitive. Accepted under the visual-identity standard;
      documented in README. Port keeps float32 to match the oracle exactly.
- [x] Phase 4 — manuscript copy in `manuscript/`: `manuscript.tex` with all 18
      figure/table paths relativized to `../figures/` + `../tables/` (19 lines
      changed, nothing else); local `agujournal2019.cls`. Cited-only bib (38
      entries) extracted **verbatim** from `all.bib` via `manuscript/extract_bib.py`
      (provenance only — needs external `all.bib`, not in `make all`); `references.bib`
      is a frozen committed source. NB `bibexport` was rejected: it silently drops
      non-standard `day` fields and needs `@string` resolution — verbatim copy of
      each entry + all 96 `@string` macros yields a `.bbl` byte-identical to a
      full-`all.bib` compile. Compiles self-contained; PDF-diff vs a fresh oracle
      compile is pixel-identical on all text/layout/reference pages. Remaining
      figure/table page diffs (7,9,12,15,23,25) are the same float32 SVD exception
      as Phase 3 — surfacing as the **N34-PC2** mode in the figures (all GMST-PC
      panels clean) and GMST-PC2 in `regression_table_full`. Build products
      (`manuscript.{aux,bbl,blg,log,pdf}`) gitignored. `references.bib` reflects
      the current `all.bib` (concept DOI `...20514350`, matching this file).
- [x] Phase 5 — build system + docs. `environment.yml`: version-pinned, pruned
      from `pangeo-2025` (analysis stack + jupytext/nbconvert + `poppler` for the
      verify render; env name `enso-gmt-capsule`). `Makefile` targets
      `all`/`figures`/`tables`/`manuscript`/`notebooks`/`verify`/`clean`
      (`PYTHON ?=` override; per-script stamps in `.build/` so `all` runs each
      script once; `manuscript` = `latexmk`; `notebooks` = `$(PYTHON) -m jupytext
      --execute`). `scripts/verify.py` rebuilds into `.verify/` and compares to
      the committed `figures/`+`tables/` — tables char-identical, figures
      `pdftoppm`@150dpi + per-pixel diff (matplotlib reader, no ImageMagick dep);
      reports EXACT/WITHIN-TOL/WARN/FAIL, `regression_table_full` float32 case is
      WARN-not-FAIL. README rewritten with the reproduction map. **On-machine
      result: 18/18 EXACT, `make all` leaves figures/tables git-clean.**
      Three script changes (both scripts, behaviour-preserving in script mode):
      (1) `import os`; (2) output dirs read `CAPSULE_FIGURES_DIR`/
      `CAPSULE_TABLES_DIR` (default unchanged) so verify builds to scratch;
      (3) `_ROOT` now falls back from `__file__` to a CWD-ancestor search, since
      `__file__` is undefined under `jupytext --execute`. Bug found + fixed:
      `references.bib`'s header comment contained the literal `@string`, which
      BibTeX parses (it has no `%` comments) → non-fatal error that still broke
      `latexmk`; reworded to "string-macro" in both `references.bib` and
      `extract_bib.py`. `.bbl` content is unchanged (comment is not an entry).
- [x] Phase 6 — publish (**done**; the two items below are permanently optional,
      not blockers). Clean-machine test done via **fresh clone** (temp dir,
      `pangeo-2025` interpreter; a fresh locked-env create was deferred by user
      choice): `make all` + `make verify` = 18/18 EXACT. The clone test exposed
      that "make all leaves git-clean" (Phase 5) held only because the gitignored
      `.build/*.stamp` short-circuits a rebuild on the origin machine — a genuine
      fresh clone regenerated all 12 figure PDFs, and matplotlib stamps each with
      the build-time `/CreationDate`, so they came back byte-modified (render
      pixel-identical; tables text-only, unaffected). **Fix (user-approved):**
      added `metadata={'CreationDate': None}` to all 12 `savefig` calls — a
      deliberate metadata-only deviation from the read-only oracle (which emits
      the timestamp). Verified: figures now byte-identical across runs (12/12),
      no `CreationDate` field, and pixel-identical to the pre-change oracle-matching
      references (18/18 EXACT). `make all` on a clean clone is now truly byte-clean.
      **Published:** public GitHub repo `mktippett/enso-gmt-prediction-capsule`
      (topics set); `CITATION.cff` + `.zenodo.json` added (both paper
      authors as creators; record linked to the preprint and dataset concept DOI
      `10.5281/zenodo.20514350`).

- [x] Browse notebooks (post-Phase-6, 2026-07-18) — committed inline-figure
      notebooks under `docs/browse/` so the code + all 12 figures are viewable on
      GitHub without cloning (Option C: a browsing convenience, deliberately NOT
      in the reproduction contract or `make verify`). `scripts/build_browse.py`
      reads each `scripts/*.py` **unmodified**, prepends a markdown banner + a
      shim cell (`%matplotlib inline` + a no-op override of `matplotlib.use` to
      neutralize the scripts' `Agg` switch), executes via nbclient, then
      `normalize()`s away the three nondeterminism sources (random cell ids →
      `cell-{i}`; drops `metadata.execution` timing; coalesces same-stream
      outputs) so a rebuild is **byte-identical** (verified both notebooks).
      Safe to commit / refresh without churn because the scripts are
      deterministic: main script has no RNG; `nmme_comparison`'s only stochastic
      step (`block_bootstrap_trend_decade`) hard-codes `seed=42`. `make browse`
      target (not a dep of `all`); README "Browse the code and figures" section.
      Distinct from `make notebooks` (gitignored `notebooks/*.ipynb`, Agg, no
      inline output — the plain executed walkthroughs).
      **Permanently optional (deferred by user 2026-07-18 — capsule is complete
      without these; pick up only if desired):**
        1. **License** — no `LICENSE`/`.zenodo.json` license field (repo is
           all-rights-reserved by default). Add if/when a choice is made;
           trivially reversible even after a DOI exists.
        2. **Zenodo capsule DOI** — GitHub↔Zenodo webhook not enabled. If pursued:
           log into zenodo.org via "Log in with GitHub", toggle the repo ON at
           zenodo.org/account/settings/github/, THEN cut a `v1.0.0` GitHub release
           (webhook only archives releases created after it is enabled) → Zenodo
           mints concept + version DOIs → add the concept-DOI badge to README.
        3. **Binder** — badge + "Run in the browser" section were removed
           2026-07-18: the mybinder.org build fails resolving `environment.yml`
           (libmamba "Cannot find a valid extracted directory cache for
           `zipp-3.23.0`" / "Package cache error" during `mamba env update` —
           a package-cache/link failure in the Binder image build, not a bad
           pin per se). The committed `docs/browse/` notebooks now cover the
           "see it without installing" use case. If Binder is revisited, the
           lead is the mamba cache-link step in the Binder Dockerfile; try
           loosening the `zipp` pin or a `pip`-based `requirements.txt`/
           `postBuild` instead of the full conda solve.

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

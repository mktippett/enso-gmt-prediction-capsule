# Phase 0 — Dependency trace (2026-07-18)

Line-by-line trace of the two analysis scripts, mapping the manuscript's
12 figures + 6 tables to producing code blocks, with KEEP/CUT closure and
data-input audit. Line numbers refer to the source-project scripts as of
2026-07-18 (oracle copies; sha256 recorded at Phase 2 fork time).

Archive-adequacy gaps C1–C3 found by this trace were CLOSED by Zenodo v2
(DOI 10.5281/zenodo.21432181, published 2026-07-18): GMSTa extended to Jan 1850
(C1: preindustrial shift derivable, +0.7840 °C), n34/GMSTa from Jan 1979
(C2: climatology months), member-level ssta/trefa + obs_n34 + Mmax/Lmax
(C3: DC-PCA figure and trends table). `Rnino34.ascii.txt` included in deposit.

## Mapping: output → script → blocks → data

### global_temperature_enso-prediction.py (9 figures + 4 tables)

| Output | savefig/write at | KEEP blocks (line ranges) | Data |
|---|---|---|---|
| eofs_pcs.pdf | 292 | setup 1–75; data load 92–130; JM matrices 139–176; PCA 188–225; fig 263–293 | n34, GMSTa |
| reconstruction_2pc.pdf | 337 | + N34 1991–2020 shift 160–171; fig 296–338 | n34, GMSTa |
| pc_correlations.pdf | 383 | corr diagnostics 343–384 | derived |
| fit_all.pdf | 918 | 4-PC regression 386–420; `time_fit` at 425; 2-PC model 453–470; fig 897–919 | derived |
| variance_explained_dc.pdf | 600 | Dec-proxy 499–510; DC decomposition 512–559; fig 561–601 | derived |
| simplified_pred.pdf | 683 | 604–608 + 626–628 (see trim traps); fig 646–684 | derived |
| corr_dec_n34_by_month.pdf | 730 | RONI load 686–708; fig 710–731 | + Rnino34.ascii.txt |
| fit_idea2.pdf | 893 | Simplified regression 854–877; fig 879–895 | derived |
| forecast.pdf | 1386 | forecast+PI 1228–1288; NMME Jun-2026 traj 1293–1313; shared_ylim 1315–1340; ENSO-response betas 1190–1209; gmt_fit_2c 1060–1066; fig 1342–1387 | + NMME GMSTa (ens. means, incl. 2026 start) |
| regression_table_full.tex | 1005 | make_reg_table 985–998; OLS 1000–1008 | derived |
| regression_table_r1.tex | 1014 | OLS 1010–1017 | derived |
| regression_table_dc2.tex | 1042 | Simplified OLS 1019–1031; 1038–1043 | derived |
| regression_table_dc2_roni.tex | 1057 | RONI OLS 1045–1058 | + Rnino34.ascii.txt |

Retargeting notes:
- Lines 92–122 (gridded ERSST/NOAAGMT area means + download fallback) are
  replaced by opening `data/obs_n34_GMSTa_185001-202605.nc`:
  `n34 = ds.n34.sel(time=PERIOD)` (drop NaN pre-1979 rows implicitly via sel),
  `gmt_all = ds.GMSTa` (full 1850–2026 record → `GMT_PREIND_SHIFT` computed
  exactly as in the original, line 119).
- Lines 1294–1297 (NMME GMT pickle load) are replaced by the archive's `GMSTa`
  (model, init_time, lead) with `rename({'init_time':'S','lead':'L'})`; the
  downstream `.sel(S=June).mean('M')` step becomes unnecessary (already
  ensemble means) — verify the June-selection and per-model anomaly logic
  (1298–1305) reproduces identically.
- `dec_n34_forecast = 3.3` (line 1245) stays hard-coded.

### nmme_comparison.py (3 figures + 2 tables)

| Output | savefig/write at | KEEP blocks | Data |
|---|---|---|---|
| nmme_skill_comparison.pdf | 525 | helpers 39–78; scaffold 85–172; NMME load 179–184; Dec N34 187–204; GMT anom 206–220; alignment 222–273; predictions 276–303; skill 306–354; fig 478–527 | n34, GMSTa, NMME n34a+GMSTa |
| nmme_trajectory_examples.pdf | 779 | Wilcoxon 356–371; fig 708–781 (uses GMT_PREIND_SHIFT) | + preind shift (from archive GMSTa) |
| nmme_traj_dc_pca.pdf | 967 | `_gmt_p`/`_n34_p`/`_June_p`/`_June_p_n34` at 791–794; fig 858–969 | **member-level trefa + ssta; obs_n34; Mmax/Lmax** |
| nmme_trends_table.tex | 1241 | bootstrap fn 54–78; table 1160–1243 (seed=42, N_BOOT=2000 — deterministic) | member-level + obs_n34 |
| nmme_significance_table.tex | 1285 | sign tests 373–417; table 1246–1286 | derived |

Retargeting notes:
- Scaffold data load 91–115 → archive obs file (same as main script; this
  script's n34 anomaly uses the 1991–2020 climatology, lines 122–124 — fully
  covered).
- Pickle loads 179–184 → archive NMME file; rename
  `{'init_time':'S','member':'M','lead':'L'}` on load so downstream selection
  code (`M=slice(0, Mmax)`, `L=6.5`) is unchanged. `ds.obs` → archive `obs_n34`.
  Where the original computes `.mean('M')` from member data (196, 215), the
  archive's `n34a`/`GMSTa` can be used directly OR recomputed from members —
  recomputation reproduces exactly (verified: |mean(trefa)−GMSTa|max = 0).
- `Mmax`/`Lmax` come from the archive variables instead of pickle attributes.

## Cut lists

### global_temperature_enso-prediction.py — CUT
227–258 (variance_explained), 428–442 (fit_full_model fig), 472–488
(fit_reduction1), 610–624 + 629–644 (eofs_pcs_dc fig), 733–852 (composites,
double-dip diagnostics, n34_pc2_pc1_scatter), 921–977 (DC corr prints + 2
heatmaps), 1033–1036 (regression_table_dc — NOT in manuscript), 1069–1182
(rmse_by_month, residuals_by_year), 1211–1226 (enso_response fig only),
1290–1291 (fcast_months), 1389–1478 (hindcast_panels, hindcast_neutral),
1480–1760 (SST correlation maps + PMM — removes cartopy + gridded data deps),
1762–1880 (E/C indices + index_correlations.tex — not in manuscript),
1882–1957 (fit_forecast_combined).
Also cuttable: `ols_pc2` (1024) — feeds only the cut table_dc and a print.

### nmme_comparison.py — CUT
419–443 (direction checks), 446–475 (nmme_n34_scatter), 530–618 (fig2b
trend-baseline), 621–705 (fig2c LOO), 795–855 (nmme_traj_pca fig body —
KEEP 791–794), `_cs_*` accumulators inside the Fig-5 loop (894–895, 945–959),
972–1097 (coupling-skill scatters), 1100–1157 (nmme_traj_pca_n34).

## TRIM TRAPS — intermediates inside cut blocks that MUST be kept

1. `time_fit` (main:425, inside cut fit_full_model) → fit_all, fit_idea2.
2. `dec_n34_z`, `beta_dec_n34` (main:606–608, inside cut eofs_pcs_dc) →
   corr_dec_n34_by_month.
3. `pc_level`, `r_level_pc1`, `r_dcpc1_pc2` (main:626–628, same cut block) →
   plotted in simplified_pred.
4. ENSO-response computations main:1190–1209 (inside cut enso_response) →
   left panel of forecast.pdf.
5. `gmt_fit_2c` (main:1062–1064), `lev_insample` (main:1318), and the whole
   shared_ylim block main:1315–1340: derived from the six CUT hindcast years
   but sets forecast.pdf's y-limits. Cutting changes the published axis.
6. nmme: `_gmt_p`, `_n34_p`, `_June_p`, `_June_p_n34` (791–794, inside cut
   Fig 4) and `_June_p_n34obs` (885, inside KEEP Fig 5) → DC-PCA figure and
   trends table.

## Reproducibility notes

- All four regression tables, variance_explained_dc, and corr_dec_n34_by_month
  are provably invariant to the climatology base period (column-centering
  cancels it); the plotted GMST trajectories are NOT — hence C2 (Jan–May 1979)
  in the archive.
- The trends-table bootstrap is seeded (`seed=42`); Wilcoxon/sign tests are
  deterministic. SVD is deterministic on one machine; last-digit differences
  possible across BLAS implementations (see CLAUDE.md gate definitions).
- matplotlib figures save via Agg backend; PDF byte-diffs will differ by
  timestamps — compare rendered pixels (`pdftoppm`), not bytes.

## Unresolved / decisions for later phases

- Figure-display accommodation for notebook builds: original scripts save PDFs
  and `plt.show()` under Agg; inline rendering in executed notebooks may need
  figures left open (Phase 5).
- `ols_pc2` cut is optional (trivial either way).

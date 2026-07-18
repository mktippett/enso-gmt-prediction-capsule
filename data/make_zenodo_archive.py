#!/usr/bin/env python
# coding: utf-8

# make_zenodo_archive.py
#
# Produces two netCDF files for Zenodo archive (v2 of record 10.5281/zenodo.20514351):
#
#   observations/obs_n34_GMSTa_185001-202605.nc
#       n34   — ERSSTv5 area-mean SST (absolute, degC), 5S-5N, 170W-120W;
#               provided from Jan 1979 (NaN before)
#       GMSTa — NOAA GlobalTemp v6.1 global area-mean temperature anomaly
#                (1991-2020 reference, degC); full record from Jan 1850 so the
#                1850-1900 preindustrial baseline shift is derivable
#       time  — monthly, January 1850 - May 2026
#
#   observations/nmme_n34a_GMSTa_199106-202606.nc
#       n34a    — NMME ssta ensemble mean per model (degC)
#       GMSTa   — NMME trefa ensemble mean per model (degC)
#       ssta    — NMME Nino-3.4 SST anomaly per member (degC)          [v2]
#       trefa   — NMME GMST anomaly per member (degC)                   [v2]
#       obs_n34 — observed Nino-3.4 on the NMME grid (init_time, lead)  [v2]
#       Mmax, Lmax — per-model ensemble size / forecast length          [v2]
#       dims  — (model, init_time, [member,] lead)
#               init_time: June 1 dates, 1991-2026
#               lead: 0.5-11.5 (months after init_time; 6.5 = December)
#
# The Zenodo deposit also includes observations/Rnino34.ascii.txt (not produced
# here): NOAA CPC monthly relative Nino-3.4 (RONI) index, copied unmodified from
# https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt (retrieved
# 2026-04-03). The manifest check at the end of this script reports its span.
#
# v2 changes (2026-07-18): obs record extended (GMSTa to Jan 1850, n34 to Jan
# 1979), member-level NMME fields + obs_n34 + Mmax/Lmax added, and
# Rnino34.ascii.txt included in the deposit, so that every figure and table of
# the manuscript is reproducible from the archive alone.

import pickle
from pathlib import Path

import numpy as np
import xarray as xr

# ── Paths ──────────────────────────────────────────────────────────────────
ERSST_PATH    = '/Users/tippett/notebooks/data/ERSSTv5.sst.mnmean.nc'
NOAAGMT_PATH  = ('/Users/tippett/notebooks/data/'
                 'NOAAGlobalTemp_v6.1.0_gridded_s185001_e202605_c20260608T115341.nc')
NMME_N34_PATH = '/Users/tippett/notebooks/NMME/n34_nmme_saved_on_disk'
NMME_GMT_PATH = '/Users/tippett/notebooks/NMME/GMT_nmme_saved_on_disk-2026-06'
OUT_DIR       = Path(__file__).parent.parent / 'observations'

N34_SLICE        = slice('1979-01-01', '2026-05-01')   # analysis PERIOD start (Jan 1979)
GMST_SLICE       = slice('1850-01-01', '2026-05-01')   # full record: preindustrial baseline
NMME_START_YEARS = np.arange(1991, 2027)   # 36 June starts (1991–2026)

print(f"NOAA GlobalTemp: {Path(NOAAGMT_PATH).name}")

# =============================================================================
# OBSERVATIONS
# =============================================================================

# n34: ERSSTv5 absolute SST, area-weighted mean 5S-5N, 170W-120W
ds_sst = xr.open_dataset(ERSST_PATH).rename({'lat': 'latitude', 'lon': 'longitude'})
sst    = ds_sst.sst.sortby('latitude')
w_n34  = np.cos(np.deg2rad(sst.latitude))
n34    = (sst
          .sel(longitude=slice(190, 240), latitude=slice(-5, 5))
          .weighted(w_n34)
          .mean(['longitude', 'latitude'])
          .sel(time=N34_SLICE)
          .load())

# GMSTa: NOAA GlobalTemp anom, cosine-weighted global mean (1991-2020 base)
ds_gmt = (xr.open_dataset(NOAAGMT_PATH)
            .rename({'lat': 'latitude', 'lon': 'longitude'})
            .squeeze('z', drop=True))
w_gmt  = np.cos(np.deg2rad(ds_gmt.latitude))
GMSTa  = (ds_gmt.anom
          .weighted(w_gmt)
          .mean(['longitude', 'latitude'])
          .sel(time=GMST_SLICE)
          .load())

# Put n34 on the full GMSTa time axis (NaN before Jan 1979)
n34_full = n34.reindex(time=GMSTa.time)
assert int(n34_full.notnull().sum()) == len(n34.time), "n34 reindex lost values"

n34_full.attrs = {
    'long_name': 'Nino-3.4 sea surface temperature',
    'units': 'degC',
    'source': 'NOAA ERSSTv5',
    'region': '5S-5N, 170W-120W (lon 190-240E), cosine-latitude weighted area mean',
    'note': 'Absolute SST, not anomaly. Provided from Jan 1979; NaN before.',
}
GMSTa.attrs = {
    'long_name': 'Global mean surface temperature anomaly',
    'units': 'degC',
    'reference_period': '1991-2020',
    'source': f'NOAA GlobalTemp v6.1.0 ({Path(NOAAGMT_PATH).name})',
    'region': 'global cosine-latitude weighted area mean of anom variable',
    'note': ('Full record from Jan 1850: the preindustrial baseline shift is '
             'GMT_PREIND_SHIFT = -mean(GMSTa, 1850-1900); add to GMSTa to '
             'express degC above 1850-1900'),
}

ds_obs = xr.Dataset({'n34': n34_full.rename('n34'), 'GMSTa': GMSTa.rename('GMSTa')})
ds_obs.attrs = {
    'title': 'Observed Nino-3.4 SST and global mean surface temperature anomaly',
    'period': 'January 1850 - May 2026 (n34 from January 1979)',
    'created': '2026-07-18',
    'version': 'v2: GMSTa extended to Jan 1850 (preindustrial baseline), n34 to Jan 1979',
}

obs_path = OUT_DIR / 'obs_n34_GMSTa_185001-202605.nc'
ds_obs.to_netcdf(obs_path)
print(f"Saved: {obs_path}")
print(f"  n34:   {dict(n34.sizes)} non-NaN  {float(n34.min()):.2f} to {float(n34.max()):.2f} degC")
print(f"  GMSTa: {dict(GMSTa.sizes)}  {float(GMSTa.min()):.3f} to {float(GMSTa.max()):.3f} degC")
print(f"  GMT_PREIND_SHIFT check: {-float(GMSTa.sel(time=slice('1850', '1900')).mean()):+.4f} degC")

# =============================================================================
# NMME
# =============================================================================

print("\nLoading NMME pickles ...")
with open(NMME_N34_PATH, 'rb') as f:
    ds_n34_nmme = pickle.load(f)
with open(NMME_GMT_PATH, 'rb') as f:
    ds_gmt_nmme = pickle.load(f)
print("Loaded.")

# n34a: ssta ensemble mean per model, June starts 1991-2025
June_n34  = ds_n34_nmme.ssta.S.dt.month == 6
n34a_raw  = ds_n34_nmme.ssta.sel(S=June_n34).mean('M')   # (model, S_june, L)
n34a_sel  = n34a_raw.sel(S=n34a_raw.S.dt.year.isin(NMME_START_YEARS))

# GMSTa: trefa ensemble mean per model, June starts 1991-2025
June_gmt  = ds_gmt_nmme.trefa.S.dt.month == 6
GMSTa_raw = ds_gmt_nmme.trefa.sel(S=June_gmt).mean('M')  # (model, S_june, L)
GMSTa_sel = GMSTa_raw.sel(S=GMSTa_raw.S.dt.year.isin(NMME_START_YEARS))

assert list(n34a_sel.S.dt.year.values)  == list(NMME_START_YEARS), "n34a start years mismatch"
assert list(GMSTa_sel.S.dt.year.values) == list(NMME_START_YEARS), "GMSTa start years mismatch"
assert list(ds_n34_nmme.model.values) == list(ds_gmt_nmme.model.values), "model lists differ"

# ── v2: member-level fields, observed N34 on the NMME grid, model metadata ──
# Member-level ssta/trefa are required for the trajectory PCA figure
# (nmme_traj_dc_pca.pdf) and the trends table (block bootstrap over members).
ssta_mem  = ds_n34_nmme.ssta.sel(S=June_n34)
ssta_mem  = ssta_mem.sel(S=ssta_mem.S.dt.year.isin(NMME_START_YEARS))
trefa_mem = ds_gmt_nmme.trefa.sel(S=June_gmt)
trefa_mem = trefa_mem.sel(S=trefa_mem.S.dt.year.isin(NMME_START_YEARS))

# Observed N34 (S, L) — verification reference in the trends table and DC-PCA figure
obs_n34_sel = ds_n34_nmme.obs.sel(S=June_n34)
obs_n34_sel = obs_n34_sel.sel(S=obs_n34_sel.S.dt.year.isin(NMME_START_YEARS))

# Per-model ensemble size and forecast length (used as M/L selection bounds)
Mmax = ds_gmt_nmme.Mmax
Lmax = ds_gmt_nmme.Lmax

# Rename S → init_time, M → member, L → lead; keep coordinate values as-is.
# Drop any 'units'/'calendar' attrs inherited from the pickle — xarray will
# re-encode init_time as CF datetime on write.
n34a_out    = n34a_sel.rename({'S': 'init_time', 'L': 'lead'})
GMSTa_out   = GMSTa_sel.rename({'S': 'init_time', 'L': 'lead'})
ssta_out    = ssta_mem.rename({'S': 'init_time', 'M': 'member', 'L': 'lead'})
trefa_out   = trefa_mem.rename({'S': 'init_time', 'M': 'member', 'L': 'lead'})
obs_n34_out = obs_n34_sel.rename({'S': 'init_time', 'L': 'lead'})
for da in [n34a_out, GMSTa_out, ssta_out, trefa_out, obs_n34_out]:
    for attr in ('units', 'calendar'):
        da['init_time'].attrs.pop(attr, None)

n34a_out.attrs = {
    'long_name': 'NMME forecast Nino-3.4 SST anomaly, ensemble mean per model',
    'units': 'degC',
    'source_variable': 'ssta from n34_nmme_saved_on_disk',
    'source': 'IRI Data Library NMME, https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/',
    'anomaly_reference': ('pre-computed; two-constant climatology for '
                          'COLA-RSMAS-CCSM4, COLA-RSMAS-CESM1, NCEP-CFSv2; '
                          '1991-2020 for CanESM5, GEM5.2-NEMO, GFDL-SPEAR, NASA-GEOSS2S'),
    'processing': 'ensemble mean over members (M) per model; no further anomaly removal',
}
GMSTa_out.attrs = {
    'long_name': 'NMME forecast global mean surface temperature anomaly, ensemble mean per model',
    'units': 'degC',
    'source_variable': 'trefa from GMT_nmme_saved_on_disk-2026-06',
    'source': 'IRI Data Library NMME, https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/',
    'processing': 'ensemble mean over members (M) per model; no further anomaly removal',
}

ssta_out.attrs = {
    'long_name': 'NMME forecast Nino-3.4 SST anomaly, per member',
    'units': 'degC',
    'source_variable': 'ssta from n34_nmme_saved_on_disk',
    'anomaly_reference': n34a_out.attrs['anomaly_reference'],
    'note': 'NaN-padded beyond each model Mmax members / Lmax leads',
}
trefa_out.attrs = {
    'long_name': 'NMME forecast global mean surface temperature anomaly, per member',
    'units': 'degC',
    'source_variable': 'trefa from GMT_nmme_saved_on_disk-2026-06',
    'note': 'NaN-padded beyond each model Mmax members / Lmax leads',
}
obs_n34_out.attrs = {
    'long_name': 'Observed Nino-3.4 SST anomaly on the NMME (init_time, lead) grid',
    'units': 'degC',
    'source_variable': 'obs from n34_nmme_saved_on_disk',
    'note': 'Verification reference; NaN where the target month is not yet observed',
}

ds_nmme = xr.Dataset({'n34a': n34a_out, 'GMSTa': GMSTa_out,
                      'ssta': ssta_out, 'trefa': trefa_out,
                      'obs_n34': obs_n34_out,
                      'Mmax': Mmax, 'Lmax': Lmax})
ds_nmme['Mmax'].attrs = {'long_name': 'last valid member coordinate per model'}
ds_nmme['Lmax'].attrs = {'long_name': 'number of valid leads per model'}
ds_nmme['init_time'].attrs = {'long_name': 'June 1 initialization date'}
ds_nmme['lead'].attrs = {
    'long_name': 'forecast lead (months after init_time)',
    'note': 'lead=0.5 is June, lead=6.5 is December, lead=11.5 is May of init_year+1',
}
ds_nmme.attrs = {
    'title': 'NMME forecast Nino-3.4 and GMST anomaly, June starts 1991-2026',
    'models': ', '.join(str(m) for m in ds_n34_nmme.model.values),
    'period': 'June initializations 1991-2026 (36 start years), 12 leads',
    'created': '2026-07-18',
    'version': 'v2: member-level ssta/trefa, obs_n34, and Mmax/Lmax added',
}

nmme_path = OUT_DIR / 'nmme_n34a_GMSTa_199106-202606.nc'
ds_nmme.to_netcdf(nmme_path,
    encoding={'init_time': {'units': 'days since 1970-01-01', 'calendar': 'standard'}})
print(f"\nSaved: {nmme_path}")
print(f"  n34a:    {dict(n34a_out.sizes)}")
print(f"  GMSTa:   {dict(GMSTa_out.sizes)}")
print(f"  ssta:    {dict(ssta_out.sizes)}")
print(f"  trefa:   {dict(trefa_out.sizes)}")
print(f"  obs_n34: {dict(obs_n34_out.sizes)}")
print(f"  Mmax: {dict(zip(ds_nmme.model.values, ds_nmme.Mmax.values))}")
print(f"  Lmax: {dict(zip(ds_nmme.model.values, ds_nmme.Lmax.values))}")
print(f"  models: {list(ds_n34_nmme.model.values)}")
print(f"  init_time: {n34a_out.init_time.values[0]} to {n34a_out.init_time.values[-1]}")
print(f"  lead: {n34a_out.lead.values[0]} to {n34a_out.lead.values[-1]}")

# =============================================================================
# DEPOSIT MANIFEST — third file uploaded as-is (not produced by this script)
# =============================================================================
RONI_PATH   = OUT_DIR / 'Rnino34.ascii.txt'
RONI_SOURCE = 'https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt'
roni = np.loadtxt(RONI_PATH, skiprows=1, usecols=(0, 1, 2))
print(f"\nDeposit alongside (unmodified copy): {RONI_PATH.name}")
print(f"  source: {RONI_SOURCE}")
print(f"  {int(roni[0, 0])}-{int(roni[0, 1]):02d} to "
      f"{int(roni[-1, 0])}-{int(roni[-1, 1]):02d}, {len(roni)} monthly rows")

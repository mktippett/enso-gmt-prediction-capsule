#!/usr/bin/env python
"""Verify freshly-built outputs against the committed reference copies.

The Makefile's ``verify`` target rebuilds figures/ and tables/ into a scratch
tree (``.verify/``) by pointing the analysis scripts at it through
``CAPSULE_FIGURES_DIR`` / ``CAPSULE_TABLES_DIR``, then runs this script to
compare that tree against the committed ``figures/`` and ``tables/``:

  * tables  — character-identical text diff (deterministic; must match exactly);
  * figures — ``pdftoppm`` render at a fixed DPI + per-pixel comparison.

Cross-machine BLAS/LAPACK differences perturb the float32 SVD's sub-leading
modes (see README, "Known limitation"). PC2-derived outputs may therefore
differ from the reference at the sub-percent pixel level, and one cell of
``regression_table_full.tex`` may move by a last digit. Such cases are reported
as WITHIN-TOL / WARN rather than failures — the capsule's cross-machine
standard is visual identity, not bytewise identity.

Usage:
    python scripts/verify.py <built_root>

<built_root> must contain figures/ and tables/ subdirectories. Exit status is
0 if every output is EXACT or WITHIN-TOL, 1 if any output fails.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np

_ROOT = Path(__file__).parent.parent
REF_FIGURES = _ROOT / 'figures'
REF_TABLES = _ROOT / 'tables'

# Fraction of pixels allowed to differ before a figure is flagged. Cross-BLAS
# PC2 perturbations move well under 2 % of pixels (README); a clean match is
# essentially 0.
PIXEL_FRAC_TOL = 0.02
DPI = 150

# Known cross-machine numerical exception (float32 SVD, GMST-PC2). A last-digit
# move in these outputs is a LAPACK/BLAS artifact, not a porting difference.
KNOWN_FLOAT32_TABLES = {'regression_table_full.tex'}

GREEN, YELLOW, RED, RESET = '\033[32m', '\033[33m', '\033[31m', '\033[0m'


def _render(pdf: Path, out_prefix: Path) -> Path:
    """Render a single-page PDF to PNG at DPI; return the PNG path."""
    subprocess.run(
        ['pdftoppm', '-png', '-r', str(DPI), '-singlefile',
         str(pdf), str(out_prefix)],
        check=True, capture_output=True,
    )
    # pdftoppm -singlefile appends ".png" to the full prefix.
    return out_prefix.with_name(out_prefix.name + '.png')


def compare_figure(name: str, ref: Path, built: Path, tmp: Path):
    """Return (status, detail) for one figure PDF."""
    if not built.exists():
        return 'FAIL', 'not built'
    ref_png = _render(ref, tmp / f'{name}.ref')
    built_png = _render(built, tmp / f'{name}.new')
    a = mpimg.imread(ref_png)
    b = mpimg.imread(built_png)
    if a.shape != b.shape:
        return 'FAIL', f'shape {a.shape} vs {b.shape}'
    diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
    # A pixel "differs" if any channel moves by more than 1/255.
    per_pixel = diff.max(axis=-1) if diff.ndim == 3 else diff
    frac = float((per_pixel > 1.0 / 255).mean())
    if frac == 0.0:
        return 'EXACT', 'pixel-identical'
    if frac < PIXEL_FRAC_TOL:
        return 'WITHIN-TOL', f'{frac:.4%} pixels differ (max {diff.max():.3f})'
    return 'FAIL', f'{frac:.4%} pixels differ (> {PIXEL_FRAC_TOL:.0%})'


def compare_table(name: str, ref: Path, built: Path):
    """Return (status, detail) for one LaTeX table."""
    if not built.exists():
        return 'FAIL', 'not built'
    if ref.read_bytes() == built.read_bytes():
        return 'EXACT', 'character-identical'
    if name in KNOWN_FLOAT32_TABLES:
        return 'WARN', 'differs (known float32 SVD / GMST-PC2 exception)'
    return 'FAIL', 'differs'


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    built_root = Path(sys.argv[1])
    built_figures = built_root / 'figures'
    built_tables = built_root / 'tables'

    results = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        print('── Tables (character-identical) ' + '─' * 33)
        for ref in sorted(REF_TABLES.glob('*.tex')):
            status, detail = compare_table(ref.name, ref, built_tables / ref.name)
            results.append((status, ref.name))
            print(f'  {_fmt(status):<22} {ref.name:<34} {detail}')
        print('── Figures (pdftoppm ' + f'{DPI} dpi + pixel diff) ' + '─' * 22)
        for ref in sorted(REF_FIGURES.glob('*.pdf')):
            name = ref.stem
            status, detail = compare_figure(name, ref, built_figures / ref.name, tmp)
            results.append((status, ref.name))
            print(f'  {_fmt(status):<22} {ref.name:<34} {detail}')

    fails = [n for s, n in results if s == 'FAIL']
    exact = sum(1 for s, _ in results if s == 'EXACT')
    tol = sum(1 for s, _ in results if s in ('WITHIN-TOL', 'WARN'))
    print('─' * 64)
    print(f'  {exact} exact, {tol} within tolerance, {len(fails)} failed '
          f'(of {len(results)} outputs)')
    if fails:
        print(f'{RED}VERIFY FAILED:{RESET} ' + ', '.join(fails))
        return 1
    print(f'{GREEN}VERIFY PASSED{RESET} '
          '(all outputs exact or within the visual-identity tolerance)')
    return 0


def _fmt(status: str) -> str:
    color = {'EXACT': GREEN, 'WITHIN-TOL': YELLOW, 'WARN': YELLOW,
             'FAIL': RED}.get(status, '')
    return f'{color}{status}{RESET}'


if __name__ == '__main__':
    raise SystemExit(main())

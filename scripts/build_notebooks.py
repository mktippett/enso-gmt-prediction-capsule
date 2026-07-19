#!/usr/bin/env python
"""Build executed notebooks with figures rendered inline (committed to notebooks/).

GitHub renders ``.ipynb`` with embedded outputs, so a reader can see the code
and every figure without cloning or running anything. These notebooks are NOT
what ``make all`` / ``make verify`` regenerate: the authoritative source is
``scripts/*.py``.

Usage:
    python scripts/build_notebooks.py <script.py> <out.ipynb>

The committed scripts force ``matplotlib.use('Agg')`` so that ``make all``
writes figures to disk without an interactive display. For notebook viewing we
want the opposite: figures inline at retina resolution. Rather than fork the
scripts, we prepend a shim cell that (1) switches on the inline backend at
retina resolution and (2) neutralizes the script's own ``matplotlib.use('Agg')``
call. The scripts themselves are read unmodified.
"""
import sys

import jupytext
import nbformat
from nbclient import NotebookClient

BANNER = (
    "## Executed notebook\n"
    "\n"
    "This notebook is generated from `{script}` with every figure rendered "
    "inline, viewable on GitHub without cloning or running anything.\n"
    "\n"
    "It is **not** what `make all` regenerates. The authoritative "
    "source is the script under `scripts/`; regenerate the manuscript figures "
    "and tables with `make all`. Outputs here may carry machine-specific "
    "last-digit differences (the float32-SVD caveat documented in the README)."
)

SHIM = (
    "# notebook-only shim (not in the committed scripts): render figures inline\n"
    "%matplotlib inline\n"
    "%config InlineBackend.figure_format = 'retina'\n"
    "import matplotlib\n"
    "matplotlib.use = lambda *a, **k: None  # neutralize the script's Agg switch\n"
)


def normalize(nb) -> None:
    """Strip execution-timing nondeterminism so a rebuild is byte-stable.

    nbclient stamps each cell with random ids, wall-clock ``execution``
    metadata, and timing-dependent stream chunking. None of it affects what a
    reader sees; all of it churns the committed file on every refresh. Pin the
    ids to the cell index, drop the timing metadata, and coalesce consecutive
    same-stream outputs into one.
    """
    for i, cell in enumerate(nb.cells):
        cell["id"] = f"cell-{i}"
        cell.get("metadata", {}).pop("execution", None)
        merged = []
        for out in cell.get("outputs", []):
            if (out.get("output_type") == "stream" and merged
                    and merged[-1].get("output_type") == "stream"
                    and merged[-1].get("name") == out.get("name")):
                merged[-1]["text"] += out["text"]
            else:
                merged.append(out)
        if "outputs" in cell:
            cell["outputs"] = merged


def main(script: str, out: str) -> None:
    nb = jupytext.read(script)  # py:percent -> notebook object (unexecuted)
    nb.cells.insert(0, nbformat.v4.new_code_cell(SHIM))
    nb.cells.insert(0, nbformat.v4.new_markdown_cell(BANNER.format(script=script)))

    NotebookClient(nb, timeout=1800, kernel_name="python3").execute()
    normalize(nb)
    nbformat.write(nb, out)
    print(f"wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python scripts/build_notebooks.py <script.py> <out.ipynb>")
    main(sys.argv[1], sys.argv[2])

# Reproduction capsule — "ENSO-conditioned evolution of global mean surface
# temperature". `make all` regenerates every figure and table used by the
# manuscript from the frozen data in data/ alone.
#
#   make all         figures + tables
#   make figures     the 12 figures  -> figures/
#   make tables      the 6 LaTeX tables -> tables/
#   make manuscript  compile manuscript/manuscript.pdf (needs figures + tables)
#   make notebooks   committed inline-figure notebooks -> notebooks/*.ipynb
#   make verify      rebuild into .verify/ and compare against committed outputs
#   make clean       remove build intermediates (keeps all committed outputs)
#
# Run inside the locked environment (see environment.yml):
#   conda env create -f environment.yml && conda activate enso-gmt-capsule

PYTHON ?= python

MAIN_SCRIPT := scripts/global_temperature_enso-prediction.py
NMME_SCRIPT := scripts/nmme_comparison.py
DATA        := $(wildcard data/*.nc) data/Rnino34.ascii.txt

STAMP_DIR   := .build
MAIN_STAMP  := $(STAMP_DIR)/main.stamp
NMME_STAMP  := $(STAMP_DIR)/nmme.stamp

.PHONY: all figures tables manuscript notebooks verify clean help

all: figures tables

help:
	@sed -n 's/^# *//p' $(MAKEFILE_LIST) | sed -n '1,20p'

# Each script reads data/ and writes both figures/ and tables/; a stamp per
# script keeps `make all` from running either one twice.
$(MAIN_STAMP): $(MAIN_SCRIPT) $(DATA)
	$(PYTHON) $(MAIN_SCRIPT)
	@mkdir -p $(STAMP_DIR) && touch $@

$(NMME_STAMP): $(NMME_SCRIPT) $(DATA)
	$(PYTHON) $(NMME_SCRIPT)
	@mkdir -p $(STAMP_DIR) && touch $@

figures tables: $(MAIN_STAMP) $(NMME_STAMP)

manuscript: figures tables
	cd manuscript && latexmk -pdf -interaction=nonstopmode manuscript.tex

# Committed notebooks with figures rendered inline (retina), so the code +
# outputs are viewable on GitHub without cloning. Deliberately NOT a dependency
# of `all` or `verify` — a viewing convenience, not something `make all`
# regenerates (see scripts/build_notebooks.py).
notebooks: notebooks/global_temperature_enso-prediction.ipynb notebooks/nmme_comparison.ipynb

notebooks/%.ipynb: scripts/%.py scripts/build_notebooks.py
	@mkdir -p notebooks
	$(PYTHON) scripts/build_notebooks.py $< $@

# Build fresh outputs into a scratch tree and diff against the committed
# figures/ and tables/ (tables: character-identical; figures: pixel compare).
verify:
	@rm -rf .verify && mkdir -p .verify/figures .verify/tables
	CAPSULE_FIGURES_DIR=$(abspath .verify/figures) CAPSULE_TABLES_DIR=$(abspath .verify/tables) $(PYTHON) $(MAIN_SCRIPT)
	CAPSULE_FIGURES_DIR=$(abspath .verify/figures) CAPSULE_TABLES_DIR=$(abspath .verify/tables) $(PYTHON) $(NMME_SCRIPT)
	$(PYTHON) scripts/verify.py .verify

clean:
	rm -rf $(STAMP_DIR) .verify
	-cd manuscript && latexmk -C manuscript.tex >/dev/null 2>&1

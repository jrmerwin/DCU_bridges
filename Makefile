PYTHON ?= python
.PHONY: reproduce test figures paper check all
reproduce:
	$(PYTHON) scripts/reproduce.py --suite compact
test:
	$(PYTHON) -m unittest discover -s tests -v
figures:
	$(PYTHON) scripts/make_figures.py
paper:
	$(PYTHON) scripts/build_paper.py
check:
	$(PYTHON) scripts/verify_manifest.py
all: reproduce test figures paper

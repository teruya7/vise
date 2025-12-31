# Minimal makefile for Sphinx documentation

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help clean html dirhtml singlehtml livehtml

help:
	@echo "Sphinx documentation build targets:"
	@echo "  html       Build standalone HTML files"
	@echo "  dirhtml    Build HTML files named index.html in directories"
	@echo "  singlehtml Build a single large HTML file"
	@echo "  livehtml   Build HTML and auto-rebuild on changes"
	@echo "  clean      Remove build directory"

clean:
	rm -rf $(BUILDDIR)

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(SOURCEDIR) $(BUILDDIR)/dirhtml $(SPHINXOPTS)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(SOURCEDIR) $(BUILDDIR)/singlehtml $(SPHINXOPTS)
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

livehtml:
	sphinx-autobuild $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS)

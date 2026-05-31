PKG := $(shell Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Package"])')
VERSION := $(shell Rscript -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])')
TARBALL := $(PKG)_$(VERSION).tar.gz

.PHONY: help document test build check install install.tarball install.fresh quick clean clean-check distclean

help:
	@echo "Targets:"
	@echo "  make document        # regenerate roxygen documentation"
	@echo "  make test            # run package-aware testthat tests"
	@echo "  make build           # build source tarball without building vignettes"
	@echo "  make check           # run R CMD check on built tarball without building vignettes"
	@echo "  make install         # install current package source"
	@echo "  make install.tarball # build tarball without building vignettes and install it"
	@echo "  make install.fresh   # document, test, build, check, install tarball"
	@echo "  make quick           # document, test, install from source"
	@echo "  make clean           # remove build/check artefacts"
	@echo "  make clean-check     # remove R CMD check directories"
	@echo "  make distclean       # clean + remove generated docs"

document:
	Rscript -e 'if (!requireNamespace("roxygen2", quietly = TRUE)) stop("Install the roxygen2 package first: install.packages(\"roxygen2\")"); roxygen2::roxygenise(".", roclets = c("rd","namespace"))'

test:
	Rscript -e 'if (!requireNamespace("testthat", quietly = TRUE)) stop("Install the testthat package first: install.packages(\"testthat\")"); testthat::test_local(".", reporter = "summary")'

build:
	R CMD build --no-build-vignettes .

check: build
	R CMD check --no-manual --no-build-vignettes "$(TARBALL)"

install:
	R CMD INSTALL .

install.tarball: build
	R CMD INSTALL "$(TARBALL)"

install.fresh: document test check
	R CMD INSTALL "$(TARBALL)"

quick: document test install

clean-check:
	rm -rf "$(PKG).Rcheck"
	rm -rf ..Rcheck

clean: clean-check
	rm -f "$(TARBALL)"
	rm -f *~ src/*~ R/*~ man/*~ tests/testthat/*~

distclean: clean
	rm -f NAMESPACE
	rm -f man/*.Rd

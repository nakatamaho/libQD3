ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src include tests fortran

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = qd.pc

dist_doc_DATA = README README.md
EXTRA_DIST = docs/qd.pdf docs/td.pdf docs/edd.pdf
dist_noinst_DATA = docs/Makefile \
  docs/qd.tex docs/td.tex docs/edd.tex docs/qd.bib \
  docs/nine-two-sum.eps docs/nine-two-sum.fig \
  docs/normal_sum_prod.eps docs/normal_sum_prod.fig \
  docs/qd_add.eps docs/qd_add.fig \
  docs/qd_add_proof.eps docs/qd_add_proof.fig \
  docs/qd_add_qd_d.eps docs/qd_add_qd_d.fig \
  docs/qd_mul_accum.eps docs/qd_mul_accum.fig \
  docs/qd_mul_qd_d.eps docs/qd_mul_qd_d.fig \
  docs/quick-two-sum.eps docs/quick-two-sum.fig \
  docs/six-three-sum.eps docs/six-three-sum.fig \
  docs/three-sum-2.eps docs/three-sum-2.fig \
  docs/three-sum-3.eps docs/three-sum-3.fig \
  docs/three-sum.eps docs/three-sum.fig \
  docs/two-prod.eps docs/two-prod.fig \
  docs/two-sum.eps docs/two-sum.fig

doc: docs/qd.pdf docs/td.pdf docs/edd.pdf

docs/qd.pdf:
	$(MAKE) -C docs qd.pdf

docs/td.pdf:
	$(MAKE) -C docs td.pdf

docs/edd.pdf:
	$(MAKE) -C docs edd.pdf

cpp-demo:
	$(MAKE) -C tests demo

if HAVE_FORTRAN

fortran-demo:
	$(MAKE) -C fortran demo

demo: cpp-demo fortran-demo

else

fortran-demo:
	@echo "You need a Fortran 95 compiler to build fortran demo programs." && false

demo: cpp-demo

endif

time:
	$(MAKE) -C tests time

bin_SCRIPTS=qd-config

if ENABLE_COVERAGE
coverage:
	$(LCOV) --zerocounters --directory "$(abs_top_builddir)"
	$(MAKE) check
	$(LCOV) --capture --directory "$(abs_top_builddir)" --output-file "$(abs_top_builddir)/coverage.info"
	$(LCOV) --extract "$(abs_top_builddir)/coverage.info" "$(abs_top_srcdir)/src/*" "$(abs_top_srcdir)/include/qd/*_inline.h" --output-file "$(abs_top_builddir)/coverage.filtered.info" --ignore-errors unused
	$(GENHTML) "$(abs_top_builddir)/coverage.filtered.info" --output-directory "$(abs_top_builddir)/coverage-html"
endif

.PHONY: cpp-demo fortran-demo demo time doc coverage

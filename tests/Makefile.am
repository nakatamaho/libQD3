LDADD = $(top_builddir)/src/libqd.la
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/include -I$(top_srcdir)/include

TESTS = qd_test pslq_test c_test huge
check_PROGRAMS = qd_test pslq_test c_test huge
EXTRA_PROGRAMS = qd_timer quadt_test

dist_noinst_DATA = coeff.dat

if HAVE_FORTRAN
TESTS += f_test
check_PROGRAMS += f_test
AM_FCFLAGS = @ax_f90_modflag@$(top_builddir)/fortran
AM_FFLAGS = $(AM_FCFLAGS)

f_test_SOURCES = f_test.f
f_test_LINK=$(CXXLINK)
f_test_LDADD = $(top_builddir)/fortran/libqdmod.la \
               $(top_builddir)/fortran/libqd_f_main.la \
               $(LDADD) $(top_builddir)/src/libqd.la $(FCLIBS)
endif

if HAVE_MPFR
TESTS += oracle_test_arith oracle_test_unary oracle_test_binary \
         oracle_test_special oracle_test_identities oracle_test_io \
         oracle_test_capi oracle_test_rounding oracle_test_bf_commutativity
check_PROGRAMS += oracle_test_arith oracle_test_unary oracle_test_binary \
                  oracle_test_special oracle_test_identities oracle_test_io \
                  oracle_test_capi oracle_test_rounding oracle_test_bf_commutativity

ORACLE_COMMON_SOURCES = oracle/mpfr_oracle.cpp oracle/mpfr_oracle.h \
                        oracle/qd_rng.cpp oracle/qd_rng.h \
                        oracle/tap.h oracle/fn_registry.h

oracle_test_arith_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_arith_SOURCES = oracle/test_arith.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_arith_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_unary_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle -I$(top_srcdir)/src
oracle_test_unary_SOURCES = oracle/test_unary.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_unary_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_binary_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_binary_SOURCES = oracle/test_binary.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_binary_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_special_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_special_SOURCES = oracle/test_special.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_special_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_identities_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_identities_SOURCES = oracle/test_identities.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_identities_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_io_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_io_SOURCES = oracle/test_io.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_io_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_capi_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_capi_SOURCES = oracle/test_capi.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_capi_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_rounding_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_rounding_SOURCES = oracle/test_rounding_corners.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_rounding_LDADD = $(LDADD) @MPFR_LIBS@

oracle_test_bf_commutativity_CPPFLAGS = $(AM_CPPFLAGS) @MPFR_CFLAGS@ -I$(srcdir)/oracle
oracle_test_bf_commutativity_SOURCES = oracle/test_bf_commutativity.cpp $(ORACLE_COMMON_SOURCES)
oracle_test_bf_commutativity_LDADD = $(LDADD) @MPFR_LIBS@
endif

CLEANFILES=qd_timer quadt_test huge

qd_test_SOURCES = qd_test.cpp
pslq_test_SOURCES = pslq.h pslq_test.cpp tictoc.cpp tictoc.h
qd_timer_SOURCES = qd_timer.cpp tictoc.cpp tictoc.h
quadt_test_SOURCES = quadt_test.cpp tictoc.cpp quadt.h tictoc.h
huge_SOURCES = huge.cpp

c_test_SOURCES = c_test.c
c_test_LINK = $(CXXLINK)

time: qd_timer$(EXEEXT)
	./qd_timer$(EXEEXT)

demo: $(check_PROGRAMS) $(EXTRA_PROGRAMS)

.PHONY: time demo

AM_FCFLAGS = $(FDFLAGS)
INCLUDES = $(FFTW_INC)
LDADD = $(FFTW_LIB)

lib_LTLIBRARIES = libp3dfft.la
libp3dfft_la_LDFLAGS = -avoid-version
libp3dfft_la_SOURCES = fft_spec.F90 module.F90 fft_init.F90 fft_exec.F90 wrap.F90

module.o: setup.F90 init_plan.F90 ftran.F90 btran.F90 reorder.F90 fcomm1.F90 \
fcomm2_trans.F90 bcomm1_trans.F90 fcomm2.F90 bcomm1.F90 bcomm2.F90

all-local:
	-@[ -e "p3dfft.mod" ] && mv -f p3dfft.mod ../include
	-@[ -e ".libs/libp3dfft.a" ] && ln -s .libs/libp3dfft.a libp3dfft.a

clean-local:
	-test -z "*.mod" || rm -f *.mod

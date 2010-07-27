AM_FCFLAGS = $(DFLAGS) # $(FCFLAGS)
INCLUDES = $(FFTW_INC)
LDADD = $(FFTW_LIB)

lib_LIBRARIES = libp3dfft.a

libp3dfft_a_SOURCES = fft_spec.F90 module.F90 fft_init.F90 fft_exec.F90 wrap.F90

module.o: setup.F90 init_plan.F90 ftran.F90 btran.F90 reorder.F90 fcomm1.F90 \
fcomm2_trans.F90 bcomm1_trans.F90 fcomm2.F90 bcomm1.F90 bcomm2.F90 ghost_cell.F90

all-local:
	-@[ -e "p3dfft.mod" ] && mv -f p3dfft.mod ../include

clean-local:
	-test -z "*.mod" || rm -f *.mod

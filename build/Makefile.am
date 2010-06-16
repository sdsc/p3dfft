AM_FCFLAGS = $(DFLAGS) # $(FCFLAGS)
INCLUDES = $(FFTW_INC)
LDADD = $(FFTW_LIB)

lib_LIBRARIES = libp3dfft.a

libp3dfft_a_SOURCES = fft_spec.F90 module.F90 fft_init.F90 fft_exec.F90 wrap.f90

all-local:
	-mv -f p3dfft.mod ../include

clean-local:
	-test -z "*.mod" || rm -f *.mod

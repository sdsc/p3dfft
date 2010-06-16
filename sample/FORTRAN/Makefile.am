AM_FCFLAGS = $(DFLAGS) $(INCLUDES) # $(FCFLAGS)
INCLUDES = -I$(top_builddir)/include/
LDADD = $(top_builddir)/build/libp3dfft.a $(FFTW_LIB) $(FFTWF) 

fsampledir = $(datadir)/p3dfft_samples/
fsample_PROGRAMS =  test_sine_f.x test_sine_inplace_f.x test_rand_f.x test_spec_f.x test_inverse_f.x

test_sine_f_x_SOURCES = driver_sine.f90

test_inverse_f_x_SOURCES = driver_inverse.f90

test_sine_inplace_f_x_SOURCES =  driver_sine_inplace.f90

test_rand_f_x_SOURCES = driver_rand.f90

test_spec_f_x_SOURCES = driver_spec.f90

clean-local:
	-test -z "*.x" || rm -f *.x

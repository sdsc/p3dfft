AM_FCFLAGS = $(FDFLAGS) $(INCLUDES) 
INCLUDES = -I$(top_builddir)/include/
LDADD = $(top_builddir)/build/libp3dfft.a $(FFTW_LIB) $(FFTWF) $(ESSL_LIB) 

fsampledir = $(datadir)/p3dfft-samples/
fsample_PROGRAMS =  test_sine_many_f.x test_sine_f.x test_sine_inplace_f.x test_rand_f.x test_spec_f.x test_inverse_f.x test_cheby_f.x test_noop_f.x test_sine_inplace_many_f.x test_rand_many_f.x

test_sine_many_f_x_SOURCES = driver_sine_many.F90
test_sine_f_x_SOURCES = driver_sine.F90

test_inverse_f_x_SOURCES = driver_inverse.F90

test_sine_inplace_f_x_SOURCES =  driver_sine_inplace.F90
test_sine_inplace_many_f_x_SOURCES = driver_sine_inplace_many.F90

test_rand_f_x_SOURCES = driver_rand.F90
test_rand_many_f_x_SOURCES = driver_rand_many.F90

test_spec_f_x_SOURCES = driver_spec.F90

test_cheby_f_x_SOURCES = driver_cheby.F90

test_noop_f_x_SOURCES = driver_noop.F90

clean-local:
	-test -z "*.x" || rm -f *.x

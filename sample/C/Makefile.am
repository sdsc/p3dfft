AM_CPPFLAGS = $(DFLAGS) 
INCLUDES = -I$(top_builddir)/include/ 
LDADD = $(top_builddir)/build/libp3dfft.a $(FFTW_LIB) $(FFTWF)

csampledir = $(datadir)/p3dfft_samples/
csample_PROGRAMS =  test_sine_c.x test_sine_inplace_c.x test_rand_c.x test_spec_c.x test_inverse_c.x

test_sine_c_x_SOURCES = driver_sine.c

test_inverse_c_x_SOURCES = driver_inverse.c

test_sine_inplace_c_x_SOURCES =  driver_sine_inplace.c

test_rand_c_x_SOURCES = driver_rand.c

test_spec_c_x_SOURCES = driver_spec.c

clean-local:
	-test -z "*.x" || rm -f *.x

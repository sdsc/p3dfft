#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(P3DFFT, 2.4, dmitry@sdsc.edu)
AM_INIT_AUTOMAKE(no-dependencies)
AC_CONFIG_SRCDIR([build/module.F90])
AC_CONFIG_HEADER([config.h])

# check compilers

# initialize array index
N="0"

# save arguments passed in as a string
AC_DEFINE_UNQUOTED([CONFIGURE_ARGS],["$ac_configure_args"],[arguments passed to configure script])

# check whether to compile using IBM
AC_MSG_CHECKING([whether to use IBM compiler])
AC_ARG_ENABLE(ibm, [AC_HELP_STRING([--enable-ibm], [compile P3DFFT using IBM compiler])], enableibm=$enableval, enableibm=no)
AC_MSG_RESULT([$enableibm])
if test "$enableibm" = "yes"; then
        AC_DEFINE(IBM, 1, [Define if you want to compile P3DFFT using IBM compiler])
	AC_DEFINE(CC, "mpfort", [overriding mpicc to link C programs (only for IBM)])
	CC="mpfort"
	eval "ARRAY${N}='-DIBM'"
        N=`expr $N + 1`
fi

# check whether to compile using Intel
AC_MSG_CHECKING([whether to use Intel compiler])
AC_ARG_ENABLE(intel, [AC_HELP_STRING([--enable-intel], [compile P3DFFT using Intel compiler])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(INTEL, 1, [Define if you want to compile P3DFFT using Intel compiler])
	eval "ARRAY${N}='-DINTEL'"
	N=`expr $N + 1`
fi

# check whether to compile using PGI
AC_MSG_CHECKING([whether to use PGI compiler])
AC_ARG_ENABLE(pgi, [AC_HELP_STRING([--enable-pgi], [compile P3DFFT using PGI compiler])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
	AC_SUBST(LDFLAGS,"$LDFLAGS -pgf90libs")
        AC_DEFINE(PGI, 1, [Define if you want to compile P3DFFT using PGI compiler])
	eval "ARRAY${N}='-DPGI'"
        N=`expr $N + 1`
fi

# check whether to compile using GCC
AC_MSG_CHECKING([whether to use GCC compiler])
AC_ARG_ENABLE(gcc, [AC_HELP_STRING([--enable-gcc], [compile P3DFFT using GCC compiler])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(GCC, 1, [Define if you want to compile P3DFFT using GCC compiler])
	eval "ARRAY${N}='-DGCC'"
        N=`expr $N + 1`
fi

# check whether to enable single precision
AC_MSG_CHECKING([whether to enable single precision])
AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single], [compile P3DFFT in single precision (default is double precision)])], spval=$enableval, spval=no)
AC_MSG_RESULT([$spval])
if test "$spval" = "yes"; then
	AC_DEFINE(SINGLE_PREC, 1, [Define if you want to compile P3DFFT in single precision])
	eval "ARRAY${N}='-DSINGLE_PREC'"
        N=`expr $N + 1`
fi

# check whether to enable 1D decomposition
AC_MSG_CHECKING([whether to enable 1D decomposition])
AC_ARG_ENABLE(oned, [AC_HELP_STRING([--enable-oned], [for 1D decomposition (default is 2D but it can be made 1D through setting up grid 1 x n when running the code) ])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(ONED, 1, [Define if you want 1D decomposition])
	eval "ARRAY${N}='-DONED'"
	N=`expr $N + 1`
fi

# check whether to enable estimation
AC_MSG_CHECKING([whether to enable estimation])
AC_ARG_ENABLE(estimate, [AC_HELP_STRING([--enable-estimate], [for initialization thoroughness in setting up FFTW calls (this option only goes with --enable-fftw) .])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(ESTIMATE, 1, [Define whether you want to enable estimation])
	eval "ARRAY${N}='-DESTIMATE'"
        N=`expr $N + 1`
fi

# check whether to enable the measure algorithm
AC_MSG_CHECKING([whether to enable the measure algorithm])
AC_ARG_ENABLE(measure, [AC_HELP_STRING([--enable-measure], [for search-once-for-the-fast-algorithm (takes more time on p3dfft_setup()).])], ok=$enableval, ok=yes)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(MEASURE, 1, [Define if you want to enable the measure algorithm])
	eval "ARRAY${N}='-DMEASURE'"
        N=`expr $N + 1`
fi

# check whether to enable the patient algorithm
AC_MSG_CHECKING([whehter to enable the patient algorithm])
AC_ARG_ENABLE(patient, [AC_HELP_STRING([--enable-patient], [for search-once-for-the-fastest-algorithm (takes much more time on p3dfft_setup()).])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(PATIENT, 1, [Define if you want to enable the patient algorithm])
	eval "ARRAY${N}='-DPATIENT'"
        N=`expr $N + 1`
fi

# check whether to enable C convention for processor dimensions
AC_MSG_CHECKING([whether to enable C convention for processor dimensions])
AC_ARG_ENABLE(dimsc, [AC_HELP_STRING([--enable-dimsc], [to assign processor rows and columns according to C convention, rather than "Fortran", which is default.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(DIMS_C, 1, [Define if you want to enable C convention for processor dimensions])
	eval "ARRAY${N}='-DDIMS_C'"
        N=`expr $N + 1`
fi

# check whether to use MPI_Alltoall instead of MPI_Alltotallv
AC_MSG_CHECKING([whether to use MPI_Alltoall instead of MPI_Alltotallv])
AC_ARG_ENABLE(useeven, [AC_HELP_STRING([--enable-useeven], [for using MPI_Alltoall instead of MPI_Alltotallv. This method pads the send buffers with zeros to make them equal size. This options is not needed on most architectures.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(USE_EVEN, 1, [Define if you want to MPI_Alltoall instead of MPI_Alltotallv])
	eval "ARRAY${N}='-DUSE_EVEN'"
        N=`expr $N + 1`
fi

# check whether to enable stride-1 data structures
AC_MSG_CHECKING([whether to enable stride-1 data structures])
AC_ARG_ENABLE(stride1, [AC_HELP_STRING([--enable-stride1], [to enable stride-1 data structures on output (this may in some cases give some advantage in performance). You can define loop blocking factors NBL_X and NBL_Y to experiment, otherwise they are set to default values.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(STRIDE1, 1, [Define if you want to enable stride-1 data structures])
	eval "ARRAY${N}='-DSTRIDE1'"
        N=`expr $N + 1`
fi

# check whether to override default value of the NBL_X
AC_MSG_CHECKING([whether to override default value of NBL_X])
AC_ARG_ENABLE(nblx, [AC_HELP_STRING([--enable-nblx], [to define loop blocking factor NBL_X])], nblval=$enableval, nblval="")
if test "$nblval" != ""; then
        AC_MSG_RESULT([yes])
        AC_DEFINE(NBL_X, $nblval, [Define if you want to override the default value of NBL_X]) 
	eval "ARRAY${N}='-DNBL_X=$nblval'"
        N=`expr $N + 1`
else
        AC_MSG_RESULT([no])
fi

# check whether to override default value of NBL_Y
AC_MSG_CHECKING([whether to override the default value of NBL_Y])
AC_ARG_ENABLE(nbly, [AC_HELP_STRING([--enable-nbly], [to define loop blocking factor NBL_Y])], nblval=$enableval, nblval="")
if test "$nblval" != ""; then
        AC_MSG_RESULT([yes])
        AC_DEFINE(NBL_Y, $nblval, [Define if you want to override the default value of NBL_Y])
	eval "ARRAY${N}='-DNBL_VAL=$nblval'"
        N=`expr $N + 1`
else
        AC_MSG_RESULT([no])
fi

# check whether to use OpenMP
AC_MSG_CHECKING([whether to use OpenMP library])
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp], [for using the OpenMP library (disabled by default)])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_DEFINE(OPENMP, 1, [Define if you want to use the OpenMP library])
	eval "ARRAY${N}='-DOPENMP'"
        N=`expr $N + 1`
fi

# check whether to use FFTW
AC_MSG_CHECKING([whether to use FFTW library])
AC_ARG_ENABLE(fftw, [AC_HELP_STRING([--enable-fftw], [for using the FFTW library])], fftwval=$enableval, fftwval=no)
AC_MSG_RESULT([$fftwval])
if test "$fftwval" = "yes"; then
        AC_DEFINE(FFTW, 1, [Define if you want to use the FFTW library])
	eval "ARRAY${N}='-DFFTW'"
        N=`expr $N + 1`
fi

# check whether to use ESSL
AC_MSG_CHECKING([whether to use ESSL library])
AC_ARG_ENABLE(essl, [AC_HELP_STRING([--enable-essl], [for using ESSL (on IBM) in 1D FFTs.])], esslval=$enableval, esslval=no)
AC_MSG_RESULT([$esslval])
if test "$esslval" = "yes"; then
        AC_DEFINE(ESSL, 1, [Define if you want to use the ESSL library instead of FFTW])
	AC_SUBST(LDFLAGS,"$LDFLAGS -lessl")
	eval "ARRAY${N}='-DESSL'"
        N=`expr $N + 1`
fi

# check to see if FFTW or ESSL has been selected
if test "$fftwval" == "no" && test "$esslval" == "no"; then
        AC_MSG_ERROR([You must choose to use either ESSL or FFTW!])
fi

# check to see if FFTW and ESSL have both been selected
if test "$fftwval" == "yes" && test "$esslval" == "yes"; then
        AC_MSG_ERROR([You must only choose either ESSL or FFTW!])
fi

# check for fftw3
AC_MSG_CHECKING([FFTW directory])
AC_ARG_WITH(fftw, [AC_HELP_STRING([--with-fftw], [path to location of FFTW library])], withfftw=$withval, withfftw="")
if test "$withfftw" != "yes" && test "$withfftw" != ""; then
	AC_MSG_RESULT([$withfftw])
	if test -e [$withfftw/lib/libfftw3.a ]; then
		AC_SUBST(FFTW_INC, "-I$withfftw/include")
		AC_SUBST(FFTW_LIB, "$withfftw/lib/libfftw3.a")
	else
		AC_MSG_ERROR([libfftw3.a was not found in given location!])
	fi

	if test "$spval" == "yes"; then
		AC_SUBST(FFTWF, "$withfftw/lib/libfftw3f.a")
	fi
else
        AC_MSG_RESULT([none])
	if test "$fftwval" == "yes"; then
		AC_MSG_ERROR([location to FFTW was not given!])
	fi
fi

# construct compiler flags
N=`expr $N - 1`
TMP0=
eval TMP0="\$ARRAY0"
for i in $(seq 1 $N)
do
	eval TMP0="\$TMP0' '\$ARRAY${i}"
done

# set C compiler flags
AC_SUBST(CDFLAGS, $TMP0)

# create a new set of flags if compiler is IBM
if test "$enableibm" != "yes"; then
	AC_SUBST(FDFLAGS, $TMP0)
else
	TMP1=
	eval TMP1="-WF,-P"
	for i in $(seq 0 $N)
	do
		eval TMP1="\$TMP1','\$ARRAY${i}"
	done
	AC_SUBST(FDFLAGS, $TMP1)
fi

# Custom strings
AC_SUBST(SUCCESS, "Build successful!")

# Checks for programs.
AC_PROG_FC(mpif90)
AC_CHECK_PROGS(CC, mpicc, $CC)
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile
                 build/Makefile
	 	 sample/Makefile
	 	 sample/C/Makefile
	 	 sample/FORTRAN/Makefile])

AC_SUBST(TESTLIBS, $TESTLIBS)
AC_OUTPUT()
cp config.h include/config.h
AC_MSG_NOTICE([configure script complete!])


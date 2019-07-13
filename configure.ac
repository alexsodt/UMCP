#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([src/hd.C])

AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT(..)

AC_CHECK_FILE("./qhull-2011.1/lib/libqhullstatic.a",	
	[AC_MSG_NOTICE("Found qhull library. Assuming this is built for the proper architecture")],
	[AC_MSG_ERROR("Could not find file qhull-2011.1/lib/libqhullstatic.a. Make this file first")])

AC_ARG_ENABLE(parallel, [AC_HELP_STRING([--parallel],[compile using MPI])], parallel=yes, parallel=no)
AC_ARG_ENABLE(cuda, [AC_HELP_STRING([--cuda],[Use CUDA GPU for special tasks])], cuda=yes, cuda=no)

if test -z $CXXFLAGS; then
    CXXFLAGS=''
fi
if test -z $CFLAGS; then
    CFLAGS=''
fi

c_compilers="clang gcc"
cxx_compilers="clang g++"

if test "$parallel" = "yes"; then
	AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
	AC_DEFINE(PARALLEL,1,[Define to compile in parallel.])
	c_compilers="mpicc $c_compilers"
	cxx_compilers="mpic++ mpicxx $cxx_compilers"

if test "a$cuda" = "ayes"; then
	AC_CHECK_PROG(NVCC, nvcc, nvcc)
	AC_DEFINE(USE_CUDA,1,[Define to use CUDA GPU for special tasks.])
	AM_CUDA="-lcudart"
fi

AC_PROG_CC([$c_compilers])
AC_PROG_CXX([$cxx_compilers])
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h unistd.h])
AC_CHECK_HEADERS([gsl/gsl_nan.h],
	[],
	[AC_MSG_ERROR("Couldn't find GSL header file. Is its path in your CPPFLAGS?")])
AC_CHECK_HEADERS([fftw3.h],
	[],
	[AC_MSG_ERROR("Couldn't find fftw3.h header file. Is its path in your CPPFLAGS?")])
	


AM_CONDITIONAL(MPI, test "$parallel" = "yes" )
AM_CONDITIONAL(PARALLEL, test "$parallel" = "yes" )
AM_CONDITIONAL(CUDA, test "$cuda" = "yes" )

linearAlgebraString="no"
usingUnderscore="yes"

save_LDFLAGS=$LDFLAGS
AC_MSG_CHECKING([for Accelerate framework.])
LDFLAGS="$LDFLAGS -framework Accelerate"
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [dsyev_()])],[AC_MSG_RESULT(yes) 
		usingMacAcclerate="yes"
		linearAlgebraString="-framework Accelerate"], [AC_MSG_RESULT(no)]) 
LDFLAGS=$save_LDFLAGS

if test "$linearAlgebraString" = "no"; then 
AC_MSG_CHECKING([for vecLib framework.])
LDFLAGS="$LDFLAGS -framework vecLib"
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [dsyev()])],[AC_MSG_RESULT(yes) 
		usingMacvecLib="yes"
		linearAlgebraString="-framework vecLib"], [AC_MSG_RESULT(no)]) 
LDFLAGS=$save_LDFLAGS
fi

if test "$linearAlgebraString" = "no"; then 
# test for ATLAS linear algebra, optional.
	AC_CHECK_LIB(lapack,dsyev,[
		usingLapack="yes"
		usingUnderscore="no"
		linearAlgebraString="-llapack -lblas"])
	AC_CHECK_LIB(lapack,dsyev_,[
		usingLapack="yes"
		usingUnderscore="yes"
		linearAlgebraString="-llapack -lblas"])

	atlasTest="dgemm"

	if test "$usingUnderscore" = "yes"; then
		atlasTest=$atlasTest"_"
	fi

	if test "$usingLapack" = "yes"; then
		LDFLAGS="$LDFLAGS -llapack -lblas"
		AC_CHECK_LIB(atlas,$atlasTest,[
			usingATLAS="yes"
			linearAlgebraString="-latlas -llapack -lblas"])
		LDFLAGS=$save_LDFLAGS
	fi
fi

if test "$linearAlgebraString" = "no"; then
	AC_MSG_ERROR("Couldn't find a LAPACK/BLAS package to which to link.")
fi

if test "$usingUnderscore" = "yes"; then
	AC_MSG_NOTICE([Using a trailing underscore for lapack and blas.])
	AC_DEFINE(LAPACK_UNDERSCORE,1,[Using underscore for lapack/blas.])
else
	AC_MSG_NOTICE("Not using a trailing underscore for lapack and blas.")
fi

AM_LDADD="$AM_LDADD $linearAlgebraString"

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lqhullcpp':
AC_CHECK_LIB([qhullcpp], [main])
# FIXME: Replace `main' with a function in `-lqhullstatic':
AC_CHECK_LIB([qhullstatic], [main])
# FIXME: Replace `main' with a function in `-lqhullstatic_p':
AC_CHECK_LIB([qhullstatic_p], [main])

AC_CHECK_LIB(gsl,main,
	[],
	[AC_MSG_ERROR("Couldn't find or use the GSL libraries.")])
AC_CHECK_LIB(fftw3,main,
	[],
	[AC_MSG_ERROR("Couldn't find or use the FFTW v3 library.")])


AC_CHECK_HEADERS([gsl/gsl_nan.h], [my_found_gsl_headers=yes; break;])

AS_IF([test "x$my_found_gsl_headers" != "xyes"],
	[AC_MSG_ERROR([Unable to find GNU scientific library headers.])])

AC_CHECK_HEADERS([fftw3.h], [my_found_fftw_headers=yes; break;])

AS_IF([test "x$my_found_fftw_headers" != "xyes"],
	[AC_MSG_ERROR([Unable to find FFTW (v3) headers.])])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

AM_QHULL="../qhull-2011.1/lib/libqhullstatic.a"
AM_INCLUDE="-I../qhull-2011.1/src/libqhull"

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor gettimeofday memset pow sqrt strcasecmp strchr strncasecmp strpbrk strrchr strstr strtol])

AC_SUBST(AM_LDADD)
AC_SUBST(AM_INCLUDE)
AC_SUBST(AM_QHULL)
AC_SUBST(AM_CUDA)

AC_CONFIG_FILES([Makefile optimized/Makefile debug/Makefile examples/test/Makefile])
AC_OUTPUT

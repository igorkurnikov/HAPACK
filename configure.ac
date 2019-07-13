# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([HARLEM],[0.1])
#AC_CONFIG_SRCDIR([config/config.mk])

AC_CONFIG_MACRO_DIR([config])
m4_include([config/wxwin.m4])
#AC_CONFIG_AUX_DIR([config])
AC_CANONICAL_TARGET

dnl ========================================================
dnl =
dnl = Don't change the following two lines.  Doing so breaks:
dnl =
dnl = CFLAGS="-foo" ./configure
dnl = ok
dnl ========================================================
CFLAGS="${CFLAGS=}"
CPPFLAGS="${CPPFLAGS=}"
CXXFLAGS="${CXXFLAGS=}"
LDFLAGS="${LDFLAGS=}"

#AM_INIT_AUTOMAKE([1.10.1 foreign subdir-objects])
AM_INIT_AUTOMAKE([1.10.1 foreign ])
LT_PREREQ([2.2.6])
LT_INIT

# config header with general definitions and parameters 
AC_CONFIG_HEADER(config.h)
# Define HAPACK variable to show files are compiled as part of the package,
AC_DEFINE([HAPACK],,[Define HAPACK variable to show compilation in the package])

### PMEMD Library
AC_ARG_WITH([libpmemd],
	[AS_HELP_STRING([--with-libpmemd],[link to modified pmemd library])],
	AS_IF([test "$with_libpmemd" == "yes"],WITH_LIB_PMEMD=1),
	WITH_LIB_PMEMD=1)

if test -n "$WITH_LIB_PMEMD"; then
	SANDER_DEFS="-DWITH_LIB_PMEMD" 
fi

AM_CONDITIONAL([WITH_LIB_PMEMD], test -n "$WITH_LIB_PMEMD")

MODS_DEFS="-DHAOGL -DBLAS -DUSE_IPACK -DPYTHON_COMMAND -DINT_JUMNA -DINT_DELPHI -DWX ${SANDER_DEFS}"

datarootdir='${prefix}/share/${PACKAGE_TARNAME}'
docdir='${datarootdir}/doc'
bindir='${prefix}/bin'
datarootdir='${prefix}'

AC_SUBST(pkgscriptdir,'${datarootdir}/scripts')

AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[enable debug compilation]),enable_debug=yes,enable_debug=no)

AM_CONDITIONAL(WITH_DEBUG, test ${enable_debug} == "yes")

AC_SUBST(CXXLD,'$(CXX)')	
AC_SUBST(HARLEM_MOD_DEFS,["-DHARLEM_MOD -Dlinux -DGNU -DGCC_NEW  $MODS_DEFS "]) 	

AC_PROG_CC([gcc icc])
AC_PROG_CXX([g++ icpc])
AC_PROG_F77([g77 gfortran ifort])
AC_PROG_FC([g90 gfortran ifort])

AM_CONDITIONAL(WITH_GFORTRAN, test ${F77} == "gfortran")
#AM_CONDITIONAL(WITH_GFORTRAN, test ${F77} == "gfortran-4.5")

AC_PATH_PROGS(AS, $AS as, $CC)
AC_CHECK_PROGS(AR, ar, :)
AC_CHECK_PROGS(LD, ld, :)
AC_CHECK_PROGS(STRIP, strip, :)
AC_CHECK_PROGS(WINDRES, windres, :)	
AC_F77_LIBRARY_LDFLAGS

##############Library Check
#zlib check
AC_CHECK_HEADERS(zlib.h,,
	AC_MSG_ERROR([Cannot find zlib.h]))
AC_CHECK_LIB(z, zlibVersion,,
	AC_MSG_ERROR([Cannot find zlib]))
AC_SUBST(ZLIB_LIBS,['-lz -lm'])
AC_SUBST(ZLIB_CPPFLAGS,"")

###BOOST
AC_ARG_WITH([boostlibs],
	  AS_HELP_STRING([--with-boostlibs=DIR],[build with boost libraries located at DIR]),
	  [
	    AS_IF([test "$with_boostlibs" != yes],
	      [
		AS_IF([test "$BOOST_LIBS" != ""],
			AC_SUBST(BOOST_LIBS,$BOOST_LIBS),
			AC_SUBST(BOOST_LIBS,"-L[$with_boostlibs] -lboost_filesystem -lboost_system  "))
	      ],
	      [
		AC_MSG_NOTICE(["Path to boost libs is not specified"])
	      ]
	    )
	  ]
)

AC_ARG_WITH([boostinc],
	  AS_HELP_STRING([--with-boostinc=DIR],[build with boost includes located in DIR]),
	  [
	    AS_IF([test "$with_boostinc" != yes],
	      [
		AS_IF([test "$BOOST_INCLUDE" != ""],
			AC_SUBST(BOOST_INCLUDE,$BOOST_INCLUDE),
			AC_SUBST(BOOST_INCLUDE,"-I[$with_boostinc]"))
	      ],
	      [
		AC_MSG_NOTICE(["Path to boost include files is not specified"])
	      ]
	    )
	  ]
)
##AC_SUBST(BOOST_INCLUDE,"-I../../local/include")
##AC_SUBST(BOOST_LIBS,"-L ../../local/lib -lboost_filesystem -lboost_system ")

###SWIG
#AC_PROG_SWIG
#SWIG_ENABLE_CXX
#SWIG_MULTI_MODULE_SUPPORT
AC_SUBST(SWIG,"swig")

###PYTHON
AC_SUBST(PYTHON_PREFIX,`python3-config --prefix`)
AC_SUBST(PYTHON_CPPFLAGS,`python3-config --cflags`)
AC_SUBST(PYTHON_LDFLAGS, `python3-config --ldflags`)
#AC_SUBST(PYTHON_LDFLAGS, "-L${PYTHON_PREFIX}/lib ${PYTHON_LDFLAGS}")

###MPI
mpi_include_cxx=`mpicxx -showme:compile`
#mpi_include_fort=`mpifort -showme:compile`
mpi_include_fort=`mpif90 -showme:compile`
mpi_libs_cxx=`mpicxx -showme:link`
#mpi_libs_fort=`mpifort -showme:link`
mpi_libs_fort=`mpif90 -showme:link`
AC_SUBST(MPI_INCLUDE,"$mpi_include_cxx $mpi_include_fort")
AC_SUBST(MPI_LIBS,"$mpi_libs_cxx $mpi_libs_fort")

##############wxWidgets

AM_OPTIONS_WXCONFIG

AM_PATH_WXCONFIG([3.0.1],[wxWin=1],[wxWin=0],,)

if test $wxWin != 1; then 
	  AC_MSG_ERROR([
		wxWidgets must be installed on your system. 
		Please check that wx-config is in path, the directory
		where wxWidgets libraries are installed 
		returned by
		'wx-config --libs' or 'wx-config --static --libs' command
		is in LD_LIBRARY_PATH or equivalent variable and
		wxWidgets version is 3.0.1 or above])
fi

##############mkl
AC_ARG_WITH([mkl],
	AS_HELP_STRING([--with-mkl=DIR],[build with intel math kernel library, if DIR is specified will use this directory as Root Directory of mkl package]),
	[
		AS_IF([test "$with_mkl" != yes],
			[
				AC_MSG_NOTICE([Path to intel math kernel library: $with_mkl])
				MKL_FLAG=
				MKL_INCLUDE=
				MKL_LIBS="[$with_mkl]/lib/32"
#				BLAS_LIBS="${MKL_LIBS} -lmkl -lmkl_vml_p4"
				BLAS_LIBS="-i-static -L${MKL_LIBS} ${MKL_LIBS}/libmkl_intel.a ${MKL_LIBS}/libmkl_core.a "
#				LAPACK_LIBS="-lmkl_lapack -lmkl_ia32 -lguide"
				LAPACK_LIBS=" ${MKL_LIBS}/libmkl_lapack.a ${MKL_LIBS}/libmkl_ia32.a ${MKL_LIBS}/libguide.a "
				
				AC_SUBST(MKL_FLAG)
				AC_SUBST(MKL_INCLUDE)
				AC_SUBST(MKL_LIBS)
				AC_SUBST(BLAS_LIBS)
				AC_SUBST(LAPACK_LIBS)
			],
			[
				AC_MSG_ERROR([Path to mkl is not specified])
			]
		)
	]
	)
	
AC_ARG_WITH([mkl64],
	AS_HELP_STRING([--with-mkl64=DIR],[build with 64-bit intel math kernel library, if DIR is specified will use this directory as Root Directory of mkl package]),
	[
		AS_IF([test "$with_mkl" != yes],
			[
				AC_MSG_NOTICE([Path to intel math kernel library: $with_mkl64])
				MKL_FLAG= 
				MKL_INCLUDE=
#				MKL_LIBS="-L[$with_mkl64]/lib/intel64"
#				MKL_LIBS="-L[$with_mkl64]/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
#				MKL_LIBS="-L[$with_mkl64]/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
				MKL_LIBS=" "
				BLAS_LIBS=" "
#				LAPACK_LIBS=" -mkl "
#				LAPACK_LIBS="[$with_mkl64]/lib/intel64/libmkl_lapack95_ilp64.a [$with_mkl64]/lib/intel64/libmkl_blas95_ilp64.a"
				LAPACK_LIBS="-L[$with_mkl64]/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
				
				AC_SUBST(MKL_FLAG)
				AC_SUBST(MKL_INCLUDE)
				AC_SUBST(MKL_LIBS)
				AC_SUBST(BLAS_LIBS)
				AC_SUBST(LAPACK_LIBS)
			],
			[
				AC_MSG_ERROR([Path to mkl64 is not specified])
			]
		)
	]
	)
	

##############CUDA
AC_ARG_WITH([cuda],
  AS_HELP_STRING([--with-cuda],[build with CUDA support]),
  [
		AC_MSG_NOTICE([building with CUDA])
		AC_SUBST(NVCC,"nvcc")
		AS_IF([test "$NVCCFLAGS" != ""],
			AC_SUBST(NVCCFLAGS),
			AC_SUBST(NVCCFLAGS,"-O3 --ptxas-options=-v -arch=sm_13"))
		AC_SUBST(NVCCLIBS,"-lcudart -lcuda -L/usr/local/cuda/lib")
		AC_SUBST(NVCCHOSTFLAGS,"-I/usr/local/cuda/include")
		#AC_SUBST(NVCCFLAGS,"NVCCFLAGS = -g --ptxas-options=-v --device-emulation")
  ],
  [
    AC_MSG_NOTICE([building without CUDA])
  ]
  )
AM_CONDITIONAL(WITH_CUDA, test "$with_cuda" = yes)
##############AC_CONFIG_FILES
AC_CONFIG_FILES([Makefile])
#Config PNPS
AC_CONFIG_FILES([PNPS/Makefile])
#AS_IF([test "$enable_mpi" == "yes"],
#    AC_CONFIG_FILES([PNPS/mpi/Makefile]))
#AC_CONFIG_FILES([PNPS/pnps/Makefile])
#AC_CONFIG_FILES([PNPS/pnps_d/Makefile])
#AC_CONFIG_FILES([PNPS/wxUtil/Makefile])
AC_CONFIG_FILES([PNPS/libpnp/Makefile])
#AC_CONFIG_FILES([PNPGUI/Makefile])
#HARLEM
#HARLEM Dependencies
AC_CONFIG_FILES([MY_VFLIB/Makefile])
AC_CONFIG_FILES([DELPHI_LIB/Makefile])
AC_CONFIG_FILES([WOLFGANG/Makefile WOLFGANG/DO_LIB/Makefile WOLFGANG/IPACK/Makefile ])
AC_CONFIG_FILES([GEOBALL/Makefile])
AC_CONFIG_FILES([JUMNA_LIB/Makefile])
AC_CONFIG_FILES([HARLEM/Makefile])
AC_CONFIG_FILES([HARLEM_NOGUI/Makefile])
AC_CONFIG_FILES([HARLEMLL/Makefile])
#AS_IF([test -n "$WITH_LIB_PMEMD"],
AC_CONFIG_FILES([PMEMD_LIB/Makefile])
#Documentation
AC_CONFIG_FILES([doc/Makefile doc/API/Makefile])
#diffrent stuff
#AC_CONFIG_FILES([misc/Makefile misc/polgrid/Makefile])

AC_OUTPUT

##############Print OUTPUT variables: 
AC_MSG_RESULT([Configuration Summary:])

AC_MSG_RESULT(AS_HELP_STRING([CC],$CC))
AC_MSG_RESULT(AS_HELP_STRING([CXX],$CXX))
AC_MSG_RESULT(AS_HELP_STRING([CXXLD],$CXXLD))
AC_MSG_RESULT(AS_HELP_STRING([CPP],$CPP))
AC_MSG_RESULT(AS_HELP_STRING([CCFLAGS],$CCFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([CXXFLAGS],$CXXFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([LDFLAGS],$LDFLAGS))

AC_MSG_RESULT(AS_HELP_STRING([F77],$F77))
AC_MSG_RESULT(AS_HELP_STRING([FFLAGS],$FFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([FLIBS],$FLIBS))
AC_MSG_RESULT(AS_HELP_STRING([FC],$FC))
AC_MSG_RESULT(AS_HELP_STRING([FCFLAGS],$FCFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([FCLIBS],$FCLIBS))

#AC_MSG_RESULT([HARLEM flags:])
AC_MSG_RESULT(AS_HELP_STRING([HARLEM_MOD_DEFS],$HARLEM_MOD_DEFS))

#AC_MSG_RESULT([MPI flags:])
#AC_MSG_RESULT(AS_HELP_STRING([MPI_FLAG],$MPI_FLAG))
AC_MSG_RESULT(AS_HELP_STRING([MPI_INCLUDE],$MPI_INCLUDE))
AC_MSG_RESULT(AS_HELP_STRING([MPI_LIBS],$MPI_LIBS))

AC_MSG_RESULT(AS_HELP_STRING([BOOST_INCLUDE],$BOOST_INCLUDE))
AC_MSG_RESULT(AS_HELP_STRING([BOOST_LIBS],$BOOST_LIBS))

#AC_MSG_RESULT([BLAS and LAPACK flags:])
AC_MSG_RESULT(AS_HELP_STRING([BLAS_LIBS],$BLAS_LIBS))
AC_MSG_RESULT(AS_HELP_STRING([LAPACK_LIBS],$LAPACK_LIBS))
#AC_MSG_RESULT(AS_HELP_STRING([MKL_FLAG],$MKL_FLAG))
#AC_MSG_RESULT(AS_HELP_STRING([MKL_INCLUDE],$MKL_INCLUDE))
#AC_MSG_RESULT(AS_HELP_STRING([MKL_LIBS],$MKL_LIBS))

#AC_MSG_RESULT([PYTHON flags:])
AC_MSG_RESULT(AS_HELP_STRING([PYTHON_CPPFLAGS],$PYTHON_CPPFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([PYTHON_LDFLAGS],$PYTHON_LDFLAGS))
#AC_MSG_RESULT(AS_HELP_STRING([PYTHON_SITE_PKG],$PYTHON_SITE_PKG))
#AC_MSG_RESULT(AS_HELP_STRING([PYTHON_EXTRA_LIBS],$PYTHON_EXTRA_LIBS))
#AC_MSG_RESULT(AS_HELP_STRING([PYTHON_EXTRA_LDFLAGS],$PYTHON_EXTRA_LDFLAGS))

#AC_MSG_RESULT([wxWidgets flags:])
AC_MSG_RESULT(AS_HELP_STRING([WX_CPPFLAGS],$WX_CPPFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([WX_CFLAGS],$WX_CFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([WX_CXXFLAGS],$WX_CXXFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([WX_CFLAGS_ONLY],$WX_CFLAGS_ONLY))
AC_MSG_RESULT(AS_HELP_STRING([WX_CXXFLAGS_ONLY],$WX_CXXFLAGS_ONLY))
AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS],$WX_LIBS))
AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS_GL],$WX_LIBS_GL))
AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS_BASE],$WX_LIBS_BASE))
AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS_STATIC],$WX_LIBS_STATIC))
AC_MSG_RESULT(AS_HELP_STRING([WX_VERSION],$WX_VERSION))
AC_MSG_RESULT(AS_HELP_STRING([WX_RESCOMP],$WX_RESCOMP))
AC_MSG_RESULT(AS_HELP_STRING([WX_RC],$WX_RC))

AC_MSG_RESULT(AS_HELP_STRING([WITH_LIB_PMEMD],$WITH_LIB_PMEMD))

AC_MSG_RESULT([CUDA flags:])
AC_MSG_RESULT(AS_HELP_STRING([WITH_CUDA],$WITH_CUDA))
AC_MSG_RESULT(AS_HELP_STRING([NVCC],$NVCC))
AC_MSG_RESULT(AS_HELP_STRING([NVCCFLAGS],$NVCCFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([NVCCHOSTFLAGS],$NVCCHOSTFLAGS))

AC_MSG_RESULT([Installation paths:])
AC_MSG_RESULT(AS_HELP_STRING([prefix],$prefix))
AC_MSG_RESULT(AS_HELP_STRING([bindir],$bindir))
AC_MSG_RESULT(AS_HELP_STRING([datarootdir],$datarootdir))
AC_MSG_RESULT(AS_HELP_STRING([datadir],$datadir))
AC_MSG_RESULT(AS_HELP_STRING([pkgscriptdir],$pkgscriptdir))
AC_MSG_RESULT(AS_HELP_STRING([pkgbasisdir],$pkgbasisdir))
AC_MSG_RESULT(AS_HELP_STRING([pkgresiduesdbdir],$pkgresiduesdbdir))

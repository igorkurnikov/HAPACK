# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([HARLEM],[0.1])
#AC_CONFIG_SRCDIR([config/config.mk])

AC_CONFIG_MACRO_DIR([config])
m4_include([config/wxwin.m4])
#AC_CONFIG_AUX_DIR([config])
AC_CANONICAL_TARGET

dnl ========================================================
dnl = This is for:
dnl = CFLAGS="-foo" ./configure
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
	WITH_LIB_PMEMD=0)

#if test -n "$WITH_LIB_PMEMD"; then
#	SANDER_DEFS="-DWITH_LIB_PMEMD" 
#fi

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
AC_PROG_F77([gfortran g77 ifort])
AC_PROG_FC([gfortran g90 ifort])

AM_CONDITIONAL(WITH_GFORTRAN, test ${F77} == "gfortran")
#AM_CONDITIONAL(WITH_GFORTRAN, test ${F77} == "gfortran-4.5")

AC_PATH_PROGS(AS, $AS as, $CC)
AC_CHECK_PROGS(AR, ar, :)
AC_CHECK_PROGS(LD, ld, :)
AC_CHECK_PROGS(STRIP, strip, :)
AC_CHECK_PROGS(WINDRES, windres, :)	
AC_F77_LIBRARY_LDFLAGS

#### ZLIB
AC_SUBST(ZLIB_LIBS,['-lz -lm'])

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
AC_SUBST(BOOST_LIBS," -lboost_filesystem -lboost_system ")

#PLPLOT
#AC_SUBST(PLPLOT_LIBS,`pkg-config plplot --libs`)
#AC_SUBST(PLPLOTCPP_LIBS,`pkg-config plplot-c++ --libs`)
#AC_SUBST(PLPLOTWX_LIBS,`pkg-config plplot-wxwidgets --libs`)
#AS_IF( test "x$PLPLOT_LIBS" = "x",
#AC_SUBST(PLPLOT_LIBS,`pkg-config plplotd --libs`)
#AC_SUBST(PLPLOTCPP_LIBS,`pkg-config plplotd-c++ --libs`)
#AC_SUBST(PLPLOTWX_LIBS,`pkg-config plplotd-wxwidgets --libs`) 
#)
###SWIG
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

# BLAS and LAPACK libraries
#LAPACK_LIBS="-lblas -llapack"
#AC_SUBST(LAPACK_LIBS,`pkg-config lapack-netlib --libs --static`) 
#AC_SUBST(LAPACK_LIBS,`pkg-config lapack --libs --static`) 
AC_SUBST(LAPACK_LIBS,`pkg-config openblas --libs --static`) 
#AC_SUBST(LAPACK_LIBS,"/usr/lib/x86_64-linux-gnu/liblapack.a  /usr/lib/x86_64-linux-gnu/libblas.a")
##############mkl
AC_ARG_WITH([mkl64],
	AS_HELP_STRING([--with-mkl64=DIR],[build with 64-bit intel math kernel library, if DIR is specified will use this directory as Root Directory of mkl package]),
	[
		AS_IF([test "$with_mkl" != yes],
			[
				AC_MSG_NOTICE([Path to intel math kernel library: $with_mkl64])
#				LAPACK_LIBS="-L[$with_mkl64]/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_sequential -lmkl_core -lmkl_def -lmkl_avx -lmkl_avx2  -lpthread -lm -ldl"
				LAPACK_LIBS="-L[$with_mkl64]/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_sequential -lgomp -lmkl_core   -lpthread -lm -ldl"
#				LAPACK_LIBS=`pkg-config mkl-static-lp64-seq --libs`
#				LAPACK_LIBS=`pkg-config mkl-dynamic-lp64-seq --libs`
				
			],
			[
				AC_MSG_ERROR([Path to mkl64 is not specified])
			]
		)
	]
	)
AC_SUBST(LAPACK_LIBS)

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
#HARLEM Dependencies
AC_CONFIG_FILES([MY_VFLIB/Makefile])
AC_CONFIG_FILES([DELPHI_LIB/Makefile])
AC_CONFIG_FILES([WOLFGANG/Makefile WOLFGANG/DO_LIB/Makefile WOLFGANG/IPACK/Makefile ])
AC_CONFIG_FILES([GEOBALL/Makefile])
AC_CONFIG_FILES([JUMNA_LIB/Makefile])
AC_CONFIG_FILES([HARLEMLL/Makefile])
#AS_IF([test -n "$WITH_LIB_PMEMD"],
#AC_CONFIG_FILES([PMEMD_LIB/Makefile])
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
AC_MSG_RESULT(AS_HELP_STRING([MPI_INCLUDE],$MPI_INCLUDE))
AC_MSG_RESULT(AS_HELP_STRING([MPI_LIBS],$MPI_LIBS))

AC_MSG_RESULT(AS_HELP_STRING([BOOST_INCLUDE],$BOOST_INCLUDE))
AC_MSG_RESULT(AS_HELP_STRING([BOOST_LIBS],$BOOST_LIBS))

AC_MSG_RESULT(AS_HELP_STRING([PLPLOT_LIBS],$PLPLOT_LIBS))
AC_MSG_RESULT(AS_HELP_STRING([PLPLOTCPP_LIBS],$PLPLOTCPP_LIBS))
AC_MSG_RESULT(AS_HELP_STRING([PLPLOTWX_LIBS],$PLPLOTWX_LIBS))

#AC_MSG_RESULT([BLAS and LAPACK flags:])
AC_MSG_RESULT(AS_HELP_STRING([LAPACK_LIBS],$LAPACK_LIBS))

#AC_MSG_RESULT([PYTHON flags:])
AC_MSG_RESULT(AS_HELP_STRING([PYTHON_CPPFLAGS],$PYTHON_CPPFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([PYTHON_LDFLAGS],$PYTHON_LDFLAGS))

#AC_MSG_RESULT([wxWidgets flags:])
AC_MSG_RESULT(AS_HELP_STRING([WX_CPPFLAGS],$WX_CPPFLAGS))
AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS],$WX_LIBS))
AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS_GL],$WX_LIBS_GL))
#AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS_BASE],$WX_LIBS_BASE))
#AC_MSG_RESULT(AS_HELP_STRING([WX_LIBS_STATIC],$WX_LIBS_STATIC))
AC_MSG_RESULT(AS_HELP_STRING([WX_VERSION],$WX_VERSION))
# AC_MSG_RESULT(AS_HELP_STRING([WX_RESCOMP],$WX_RESCOMP))
# AC_MSG_RESULT(AS_HELP_STRING([WX_RC],$WX_RC))

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

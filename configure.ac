#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(PVMOS, 0.84, b.pieters@fz-juelich.de)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/mesh2d.c])
AC_CONFIG_HEADERS([src/config.h])
AC_PREFIX_DEFAULT([/usr])
m4_include([ax_blas.m4])
# Checks for programs.

if test -z $CFLAGS; then
    CFLAGS='-Ofast -Wall -fPIC -flto'
fi

AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MKDIR_P

AC_ARG_ENABLE([openblas-multithreaded],[AS_HELP_STRING([--enable-openblas-multithreaded],[Enable multi-threading in openblas (as of OpenBLAS 2.15 this is a bad idea)])],
[case "${enableval}" in
  yes) AC_DEFINE(OPENBLAS, 0,  [1: disable multi-threading for openblas, 0: do not]);;
  no)  AC_DEFINE(OPENBLAS, 1,  [1: disable multi-threading for openblas, 0: do not]) ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openblas-singlethreaded]) ;;
esac],[AC_DEFINE(OPENBLAS, 0,  [1: disable multi-threading for openblas, 0: do not])])

AX_BLAS()
if [test ! x$BLAS_LIBS = x-lopenblas]; then
	AC_DEFINE(OPENBLAS, 0)
fi

AC_ARG_ENABLE([generate-octave-pkg],[AS_HELP_STRING([--enable-generate-octave-pkg],[Create the pvmos-mesh package for GNU octave])],
[case "${enableval}" in
  yes) pvmosmesh=yes ;;
  no)  pvmosmesh=no;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-generate-octave-pkg]) ;;
esac],[pvmosmesh=yes])
AM_CONDITIONAL([PVMOSMESH], [test x$pvmosmesh = xyes])


AC_ARG_ENABLE([dmalloc],[AS_HELP_STRING([--enable-dmalloc],[use the dmalloc memory allocator (recommended when compiling with mingw)])],
[case "${enableval}" in
  yes) dmalloc=malloc-2.8.6.o ;;
  no)  dmalloc="";;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-dmalloc]) ;;
esac],[dmalloc=""])

AM_CONDITIONAL([DMALLOC], [test x$dmalloc = xmalloc-2.8.6.o])
# Checks for libraries.
AC_CHECK_LIB(m, sin, [], [
        echo "Error! math library not functional."
        exit -1
        ])
# Check for CHOLMOD
AC_CHECK_LIB(cholmod,cholmod_start,[],[echo "Error! missing cholmod."
        exit -1])


AC_CHECK_LIB(matheval, evaluator_create,[],[])

AC_ARG_ENABLE([inkscape-extension],[AS_HELP_STRING([--enable-inkscape-extension],[Attempt to install the inkscape extension to extract polygons])],
[case "${enableval}" in
  yes) inkscape_extension=yes ;;
  no)  inkscape_extenion=no;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-inkscape-extension]) ;;
esac],[inkscape_extension=yes])
INKSCAPE_EXT_PATH="~/.config/inkscape/extensions"

if [test x$inkscape_extension = xyes ]; then
	AC_CHECK_PROG(INKSCAPE, inkscape, yes, no)
	if test x$INKSCAPE = xno; then
		echo "Inkscape not available, cannot install extension"
		inkscape_extension=no
	fi
	REAL_INKSCAPE_EXT_PATH=$(inkscape -x)
	INKSCAPE_EXT_PATH=`echo $REAL_INKSCAPE_EXT_PATH |sed -e 's|^.*\(/inkscape.*\)|\1|g'`
	INKSCAPE_PREFIX=`echo $REAL_INKSCAPE_EXT_PATH |sed -e 's|^\(.*\)/inkscape.*|\1|g'`
fi
AM_CONDITIONAL([INKSCAPE_EXT], [test x$inkscape_extension = xyes ])
AC_SUBST([inkscape_ext_path], [$INKSCAPE_EXT_PATH])

	


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE


# default is "no" except on OS X
# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow putenv sqrt strchr])
LIBS="$BLAS_LIBS $LIBS $FLIBS"

AC_CHECK_PROG(tarp, tar, yes, no)
if test x$tarp = xno; then
    echo tar not found
    exit -1
fi

AC_CHECK_PROG(build_doc, pdflatex, yes, no)
AM_CONDITIONAL([BUILDDOC], [test x$build_doc = xyes])
if test x$build_doc = xno; then
    echo Cannot build the pdf.
fi
AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile examples/Makefile extractpoly/Makefile])
AC_OUTPUT()

DIR=$bindir
while [test x`expr index "$DIR" \$` = x1]
do
	ECHOD="echo $DIR"
	DIR=`eval $ECHOD`
done
BINDIR=$DIR

DIR=$datadir
while [test x`expr index "$DIR" \$` = x1]
do
	ECHOD="echo $DIR"
	DIR=`eval $ECHOD`
done
DATADIR=$DIR

DIR=$docdir
while [test x`expr index "$DIR" \$` = x1]
do
	ECHOD="echo $DIR"
	DIR=`eval $ECHOD`
done
DOCDIR=$DIR


echo "--------------------------------------------------------------------------"
echo "Configuration Summary  For $PACKAGE $VERSION"
echo "--------------------------------------------------------------------------"
echo "bindir			$BINDIR"
echo "datadir			$DATADIR"
echo "docdir			$DOCDIR"
echo "--------------------------------------------------------------------------"
if [test x$pvmosmesh = xyes]; then
	echo "After installing $PACKAGE you may install the pvmos-mesh package"
	echo "in GNU Octave by executing (in Octave):"
	echo "\"pkg install $DATADIR/pvmos-mesh-$VERSION.tar.gz\""
	echo
fi
if [test  x$inkscape_extension = xyes ]; then
	echo "The extract poly extension will be installed to"
	echo "$DATADIR$INKSCAPE_EXT_PATH"
	echo
	if [test ! x$INKSCAPE_PREFIX = x$DATADIR ]; then
	echo "--> Warning: this is not the path returned by \"inkscape -x\""
	echo "    To fix this pass \"--datadir=$INKSCAPE_PREFIX\" to configure"
	fi
else
	echo "The extract poly extension will not be installed"
fi
echo "To install $PACKAGE" proceed with
echo "\"make\""
echo 
echo "Followed by (as root)"
echo "\"make install\""
echo 
echo "--------------------------------------------------------------------------"

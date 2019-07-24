#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT([UMCP-docs], [0.1], [alexander.sodt@nih.gov])
AC_CONFIG_SRCDIR([doc.tex])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([config.h])

# Find latex program, probably pdflatex
AC_PROG_LATEX

if [test "x$latex" == "xpdflatex" ]; then
	AC_DEFINE(HAVE_PDFLATEX,1,"Using pdflatex")	
fi

# To the best of our ability, test for all required packages.

ac_latex_err=0

AC_TEST_LATEX(footmisc,HAVE_LATEX_FOOTMISC)
AC_TEST_LATEX(amsmath,HAVE_LATEX_AMSMATH)
AC_TEST_LATEX(graphicx,HAVE_LATEX_GRAPHICX)
AC_TEST_LATEX(xcolor,HAVE_LATEX_XCOLOR)
AC_TEST_LATEX(comment,HAVE_LATEX_COMMENT)
AC_TEST_LATEX(etoolbox,HAVE_LATEX_ETOOLBOX)
AC_TEST_LATEX(tcolorbox,HAVE_LATEX_TCOLORBOX)
AC_TEST_LATEX(helvet,HAVE_LATEX_HELVET)
AC_TEST_LATEX(setspace,HAVE_LATEX_SETSPACE)
AC_TEST_LATEX(wrapfig,HAVE_LATEX_WRAPFIG)
AC_TEST_LATEX(newfloat,HAVE_LATEX_NEWFLOAT)
AC_TEST_LATEX(mdframed,HAVE_LATEX_MDFRAMED)
AC_TEST_LATEX(cleveref,HAVE_LATEX_CLEVEREF)
AC_TEST_LATEX(calc,HAVE_LATEX_CALC)
AC_TEST_LATEX(regexpatch,HAVE_LATEX_REGEXPATCH)
AC_TEST_LATEX(filecontents,HAVE_LATEX_FILECONTENTS)
AC_TEST_LATEX(ifxetex,HAVE_LATEX_IFXETEX )
AC_TEST_LATEX(ifluate,HAVE_LATEX_IFLUATEX)

if [test "x$ac_latex_err" == "x1"]; then
	$as_echo "The following files were missing: $ac_latex_missing"
	$as_echo "Failed to build Makefile. Install missing latex packages."
else
	$as_echo "Found and built with all (known) required latex packages."
fi
#

AC_CONFIG_FILES(Makefile pdfs/Makefile)
AC_OUTPUT

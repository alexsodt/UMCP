# from Matthieu Boretti (boretti@@bss-network.com) GPL/AC exception
AC_DEFUN([AC_PROG_LATEX],[
AC_CHECK_PROGS(latex,[pdflatex],no)
export latex;
if test $latex = "no" ;
then
	AC_MSG_ERROR([Unable to find a LaTeX application]);
fi
AC_SUBST(latex)
])

# written by AJS
# first argument is name of packate to try
AC_DEFUN([AC_TEST_LATEX],[
cat > conftest.tex << _ACEOF
\documentclass{article}
\usepackage{$1}
\begin{document}
Hello, world!
\end{document}
_ACEOF
$latex -interaction=nonstopmode conftest >& conftest_lat.out
# did we write conftest.pdf
if test -e conftest.pdf; then 
	$as_echo "Found required latex package $1"
	AC_DEFINE($2,1,[Defining latex package available])
else
	$as_echo "Could not find required latex package $1"
	ac_latex_err="1"
	ac_latex_missing="$ac_latex_missing $1"
fi
ac_clean_files="$ac_clean_files conftest.tex conftest.aux conftest_lat.out"
])

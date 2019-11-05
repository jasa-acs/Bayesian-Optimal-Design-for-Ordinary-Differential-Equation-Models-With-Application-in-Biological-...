#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP compartmentalnaelcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compartmentalnsel4cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP compartmentalsigcpp(SEXP, SEXP, SEXP);
extern SEXP compsolver2cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fitzhughnaelcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fitzhughnsel4cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP fitzhughsigcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP fitzsolver2cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP jakstatnaelcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP jakstatnsel4cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP jakstatsigcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP jakstatsolver2cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP modcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP naelcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nsel4cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP placsolver1cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"compartmentalnaelcpp",  (DL_FUNC) &compartmentalnaelcpp,   5},
    {"compartmentalnsel4cpp", (DL_FUNC) &compartmentalnsel4cpp,  4},
    {"compartmentalsigcpp",   (DL_FUNC) &compartmentalsigcpp,    3},
    {"compsolver2cpp",        (DL_FUNC) &compsolver2cpp,         7},
    {"fitzhughnaelcpp",       (DL_FUNC) &fitzhughnaelcpp,        5},
    {"fitzhughnsel4cpp",      (DL_FUNC) &fitzhughnsel4cpp,       4},
    {"fitzhughsigcpp",        (DL_FUNC) &fitzhughsigcpp,         4},
    {"fitzsolver2cpp",        (DL_FUNC) &fitzsolver2cpp,         7},
    {"jakstatnaelcpp",        (DL_FUNC) &jakstatnaelcpp,         5},
    {"jakstatnsel4cpp",       (DL_FUNC) &jakstatnsel4cpp,        4},
    {"jakstatsigcpp",         (DL_FUNC) &jakstatsigcpp,          4},
    {"jakstatsolver2cpp",     (DL_FUNC) &jakstatsolver2cpp,     12},
    {"modcpp",                (DL_FUNC) &modcpp,                 4},
    {"naelcpp",               (DL_FUNC) &naelcpp,                5},
    {"nsel4cpp",              (DL_FUNC) &nsel4cpp,               4},
    {"placsolver1cpp",        (DL_FUNC) &placsolver1cpp,         9},
    {NULL, NULL, 0}
};

void R_init_aceodes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

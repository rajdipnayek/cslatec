/* defehl.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* DECK DEFEHL */
/* Subroutine */ int defehl_(S_fp f, integer *neq, real *t, real *y, real *
	h__, real *yp, real *f1, real *f2, real *f3, real *f4, real *f5, real 
	*ys, real *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer k;
    static real ch;

/* ***BEGIN PROLOGUE  DEFEHL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DEFEHL-S, DFEHL-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     Fehlberg Fourth-Fifth order Runge-Kutta Method */
/* ********************************************************************** */

/*    DEFEHL integrates a system of NEQ first order */
/*    ordinary differential equations of the form */
/*               dU/DX = F(X,U) */
/*    over one step when the vector Y(*) of initial values for U(*) and */
/*    the vector YP(*) of initial derivatives, satisfying  YP = F(T,Y), */
/*    are given at the starting point X=T. */

/*    DEFEHL advances the solution over the fixed step H and returns */
/*    the fifth order (sixth order accurate locally) solution */
/*    approximation at T+H in the array YS(*). */
/*    F1,---,F5 are arrays of dimension NEQ which are needed */
/*    for internal storage. */
/*    The formulas have been grouped to control loss of significance. */
/*    DEFEHL should be called with an H not smaller than 13 units of */
/*    roundoff in T so that the various independent arguments can be */
/*    distinguished. */

/*    This subroutine has been written with all variables and statement */
/*    numbers entirely compatible with DERKFS. For greater efficiency, */
/*    the call to DEFEHL can be replaced by the module beginning with */
/*    line 222 and extending to the last line just before the return */
/*    statement. */

/* ********************************************************************** */

/* ***SEE ALSO  DERKF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891009  Removed unreferenced statement label.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DEFEHL */



/* ***FIRST EXECUTABLE STATEMENT  DEFEHL */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --ys;
    --f5;
    --f4;
    --f3;
    --f2;
    --f1;
    --yp;
    --y;

    /* Function Body */
    ch = *h__ / 4.f;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L230: */
	ys[k] = y[k] + ch * yp[k];
    }
    r__1 = *t + ch;
    (*f)(&r__1, &ys[1], &f1[1], &rpar[1], &ipar[1]);

    ch = *h__ * 3.f / 32.f;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L240: */
	ys[k] = y[k] + ch * (yp[k] + f1[k] * 3.f);
    }
    r__1 = *t + *h__ * 3.f / 8.f;
    (*f)(&r__1, &ys[1], &f2[1], &rpar[1], &ipar[1]);

    ch = *h__ / 2197.f;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L250: */
	ys[k] = y[k] + ch * (yp[k] * 1932.f + (f2[k] * 7296.f - f1[k] * 
		7200.f));
    }
    r__1 = *t + *h__ * 12.f / 13.f;
    (*f)(&r__1, &ys[1], &f3[1], &rpar[1], &ipar[1]);

    ch = *h__ / 4104.f;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L260: */
	ys[k] = y[k] + ch * (yp[k] * 8341.f - f3[k] * 845.f + (f2[k] * 
		29440.f - f1[k] * 32832.f));
    }
    r__1 = *t + *h__;
    (*f)(&r__1, &ys[1], &f4[1], &rpar[1], &ipar[1]);

    ch = *h__ / 20520.f;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L270: */
	ys[k] = y[k] + ch * (yp[k] * -6080.f + (f3[k] * 9295.f - f4[k] * 
		5643.f) + (f1[k] * 41040.f - f2[k] * 28352.f));
    }
    r__1 = *t + *h__ / 2.f;
    (*f)(&r__1, &ys[1], &f5[1], &rpar[1], &ipar[1]);

/*     COMPUTE APPROXIMATE SOLUTION AT T+H */

    ch = *h__ / 7618050.f;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L290: */
	ys[k] = y[k] + ch * (yp[k] * 902880.f + (f3[k] * 3855735.f - f4[k] * 
		1371249.f) + (f2[k] * 3953664.f + f5[k] * 277020.f));
    }

    return 0;
} /* defehl_ */


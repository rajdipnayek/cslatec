/* dfehl.f -- translated by f2c (version 12.02.01).
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

/* DECK DFEHL */
/* Subroutine */ int dfehl_(S_fp df, integer *neq, doublereal *t, doublereal *
	y, doublereal *h__, doublereal *yp, doublereal *f1, doublereal *f2, 
	doublereal *f3, doublereal *f4, doublereal *f5, doublereal *ys, 
	doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k;
    static doublereal ch;

/* ***BEGIN PROLOGUE  DFEHL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (DEFEHL-S, DFEHL-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     Fehlberg Fourth-Fifth Order Runge-Kutta Method */
/* ********************************************************************** */

/*    DFEHL integrates a system of NEQ first order */
/*    ordinary differential equations of the form */
/*               DU/DX = DF(X,U) */
/*    over one step when the vector Y(*) of initial values for U(*) and */
/*    the vector YP(*) of initial derivatives, satisfying  YP = DF(T,Y), */
/*    are given at the starting point X=T. */

/*    DFEHL advances the solution over the fixed step H and returns */
/*    the fifth order (sixth order accurate locally) solution */
/*    approximation at T+H in the array YS(*). */
/*    F1,---,F5 are arrays of dimension NEQ which are needed */
/*    for internal storage. */
/*    The formulas have been grouped to control loss of significance. */
/*    DFEHL should be called with an H not smaller than 13 units of */
/*    roundoff in T so that the various independent arguments can be */
/*    distinguished. */

/*    This subroutine has been written with all variables and statement */
/*    numbers entirely compatible with DRKFS. For greater efficiency, */
/*    the call to DFEHL can be replaced by the module beginning with */
/*    line 222 and extending to the last line just before the return */
/*    statement. */

/* ********************************************************************** */

/* ***SEE ALSO  DDERKF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DFEHL */


/* ***FIRST EXECUTABLE STATEMENT  DFEHL */
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
    ch = *h__ / 4.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	ys[k] = y[k] + ch * yp[k];
/* L10: */
    }
    d__1 = *t + ch;
    (*df)(&d__1, &ys[1], &f1[1], &rpar[1], &ipar[1]);

    ch = *h__ * 3. / 32.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	ys[k] = y[k] + ch * (yp[k] + f1[k] * 3.);
/* L20: */
    }
    d__1 = *t + *h__ * 3. / 8.;
    (*df)(&d__1, &ys[1], &f2[1], &rpar[1], &ipar[1]);

    ch = *h__ / 2197.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	ys[k] = y[k] + ch * (yp[k] * 1932. + (f2[k] * 7296. - f1[k] * 7200.));
/* L30: */
    }
    d__1 = *t + *h__ * 12. / 13.;
    (*df)(&d__1, &ys[1], &f3[1], &rpar[1], &ipar[1]);

    ch = *h__ / 4104.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	ys[k] = y[k] + ch * (yp[k] * 8341. - f3[k] * 845. + (f2[k] * 29440. - 
		f1[k] * 32832.));
/* L40: */
    }
    d__1 = *t + *h__;
    (*df)(&d__1, &ys[1], &f4[1], &rpar[1], &ipar[1]);

    ch = *h__ / 20520.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	ys[k] = y[k] + ch * (yp[k] * -6080. + (f3[k] * 9295. - f4[k] * 5643.) 
		+ (f1[k] * 41040. - f2[k] * 28352.));
/* L50: */
    }
    d__1 = *t + *h__ / 2.;
    (*df)(&d__1, &ys[1], &f5[1], &rpar[1], &ipar[1]);

/*     COMPUTE APPROXIMATE SOLUTION AT T+H */

    ch = *h__ / 7618050.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	ys[k] = y[k] + ch * (yp[k] * 902880. + (f3[k] * 3855735. - f4[k] * 
		1371249.) + (f2[k] * 3953664. + f5[k] * 277020.));
/* L60: */
    }

    return 0;
} /* dfehl_ */


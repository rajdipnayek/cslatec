/* ddatrp.f -- translated by f2c (version 12.02.01).
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

/* DECK DDATRP */
/* Subroutine */ int ddatrp_(doublereal *x, doublereal *xout, doublereal *
	yout, doublereal *ypout, integer *neq, integer *kold, doublereal *phi,
	 doublereal *psi)
{
    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__, d__;
    static integer i__, j;
    static doublereal temp1, gamma;
    static integer koldp1;

/* ***BEGIN PROLOGUE  DDATRP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Interpolation routine for DDASSL. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      DOUBLE PRECISION (SDATRP-S, DDATRP-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------------- */
/*     THE METHODS IN SUBROUTINE DDASTP USE POLYNOMIALS */
/*     TO APPROXIMATE THE SOLUTION. DDATRP APPROXIMATES THE */
/*     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING */
/*     ONE OF THESE POLYNOMIALS, AND ITS DERIVATIVE,THERE. */
/*     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM */
/*     DDASTP, SO DDATRP CANNOT BE USED ALONE. */

/*     THE PARAMETERS ARE: */
/*     X     THE CURRENT TIME IN THE INTEGRATION. */
/*     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED */
/*     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT */
/*           (THIS IS OUTPUT) */
/*     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT */
/*           (THIS IS OUTPUT) */
/*     NEQ   NUMBER OF EQUATIONS */
/*     KOLD  ORDER USED ON LAST SUCCESSFUL STEP */
/*     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y */
/*     PSI   ARRAY OF PAST STEPSIZE HISTORY */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/* ***END PROLOGUE  DDATRP */



/* ***FIRST EXECUTABLE STATEMENT  DDATRP */
    /* Parameter adjustments */
    --yout;
    --ypout;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --psi;

    /* Function Body */
    koldp1 = *kold + 1;
    temp1 = *xout - *x;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yout[i__] = phi[i__ + phi_dim1];
/* L10: */
	ypout[i__] = 0.;
    }
    c__ = 1.;
    d__ = 0.;
    gamma = temp1 / psi[1];
    i__1 = koldp1;
    for (j = 2; j <= i__1; ++j) {
	d__ = d__ * gamma + c__ / psi[j - 1];
	c__ *= gamma;
	gamma = (temp1 + psi[j - 1]) / psi[j];
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    yout[i__] += c__ * phi[i__ + j * phi_dim1];
/* L20: */
	    ypout[i__] += d__ * phi[i__ + j * phi_dim1];
	}
/* L30: */
    }
    return 0;

/* ------END OF SUBROUTINE DDATRP------ */
} /* ddatrp_ */


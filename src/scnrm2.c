/* scnrm2.f -- translated by f2c (version 12.02.01).
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

/* DECK SCNRM2 */
doublereal scnrm2_(integer *n, complex *cx, integer *incx)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;
    static real cutlo = 4.441e-16f;
    static real cuthi = 1.304e19f;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_90[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real ret_val, r__1;

    /* Local variables */
    static integer i__, nn;
    static real sum;
    static logical imag;
    static real absx, xmax;
    static integer next;
    static logical scale;
    static real hitest;

    /* Assigned format variables */
    static char *next_fmt;

/* ***BEGIN PROLOGUE  SCNRM2 */
/* ***PURPOSE  Compute the unitary norm of a complex vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A3B */
/* ***TYPE      COMPLEX (SNRM2-S, DNRM2-D, SCNRM2-C) */
/* ***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2, */
/*             LINEAR ALGEBRA, UNITARY, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       CX  complex vector with N elements */
/*     INCX  storage spacing between elements of CX */

/*     --Output-- */
/*   SCNRM2  single precision result (zero if N .LE. 0) */

/*     Unitary norm of the complex N-vector stored in CX with storage */
/*     increment INCX. */
/*     If N .LE. 0, return with result = 0. */
/*     If N .GE. 1, then INCX must be .GE. 1 */

/*     Four phase method using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         CUTLO = maximum of  SQRT(U/EPS)  over all known machines. */
/*         CUTHI = minimum of  SQRT(V)      over all known machines. */
/*     where */
/*         EPS = smallest no. such that EPS + 1. .GT. 1. */
/*         U   = smallest positive no.   (underflow limit) */
/*         V   = largest  no.            (overflow  limit) */

/*     Brief outline of algorithm. */

/*     Phase 1 scans zero components. */
/*     Move to phase 2 when a component is nonzero and .LE. CUTLO */
/*     Move to phase 3 when a component is .GT. CUTLO */
/*     Move to phase 4 when a component is .GE. CUTHI/M */
/*     where M = N for X() real and M = 2*N for complex. */

/*     Values for CUTLO and CUTHI. */
/*     From the environmental parameters listed in the IMSL converter */
/*     document the limiting values are as follows: */
/*     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are */
/*                   Univac and DEC at 2**(-103) */
/*                   Thus CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC. */
/*                   Thus CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC. */
/*                   Thus CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/ */
/*     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/ */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SCNRM2 */
    /* Parameter adjustments */
    --cx;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  SCNRM2 */
    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;

/*                                                 BEGIN MAIN LOOP */

    i__1 = nn;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	i__3 = i__;
	absx = (r__1 = cx[i__3].r, dabs(r__1));
	imag = FALSE_;
	switch (next) {
	    case 0: goto L30;
	    case 1: goto L50;
	    case 2: goto L70;
	    case 3: goto L110;
	    case 4: goto L90;
	}
L30:
	if (absx > cutlo) {
	    goto L85;
	}
	next = 1;
	next_fmt = fmt_50;
	scale = FALSE_;

/*                        PHASE 1.  SUM IS ZERO */

L50:
	if (absx == zero) {
	    goto L200;
	}
	if (absx > cutlo) {
	    goto L85;
	}

/*                                PREPARE FOR PHASE 2. */

	next = 2;
	next_fmt = fmt_70;
	goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
	next = 3;
	next_fmt = fmt_110;
	sum = sum / absx / absx;
L105:
	scale = TRUE_;
	xmax = absx;
	goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
	if (absx > cutlo) {
	    goto L75;
	}

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

L110:
	if (absx <= xmax) {
	    goto L115;
	}
/* Computing 2nd power */
	r__1 = xmax / absx;
	sum = one + sum * (r__1 * r__1);
	xmax = absx;
	goto L200;

L115:
/* Computing 2nd power */
	r__1 = absx / xmax;
	sum += r__1 * r__1;
	goto L200;

/*                  PREPARE FOR PHASE 3. */

L75:
	sum = sum * xmax * xmax;

L85:
	next = 4;
	next_fmt = fmt_90;
	scale = FALSE_;

/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

	hitest = cuthi / *n;

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

L90:
	if (absx >= hitest) {
	    goto L100;
	}
/* Computing 2nd power */
	r__1 = absx;
	sum += r__1 * r__1;
L200:

/*                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS. */

	if (imag) {
	    goto L210;
	}
	absx = (r__1 = r_imag(&cx[i__]), dabs(r__1));
	imag = TRUE_;
	switch (next) {
	    case 0: goto L30;
	    case 1: goto L50;
	    case 2: goto L70;
	    case 3: goto L110;
	    case 4: goto L90;
	}

L210:
	;
    }

/*              END OF MAIN LOOP. */
/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = sqrt(sum);
    if (scale) {
	ret_val *= xmax;
    }
L300:
    return ret_val;
} /* scnrm2_ */


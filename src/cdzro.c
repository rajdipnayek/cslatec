/* cdzro.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__0 = 0;
static real c_b5 = 1.f;

/* DECK CDZRO */
/* Subroutine */ int cdzro_(real *ae, E_fp f, real *h__, integer *n, integer *
	nq, integer *iroot, real *re, real *t, complex *yh, real *uround, 
	real *b, real *c__, real *fb, real *fc, complex *y)
{
    /* System generated locals */
    integer yh_dim1, yh_offset;
    real r__1;

    /* Local variables */
    static real a, p, q, fa;
    static integer ic;
    static real er, rw, cmb, tol, acmb, acbs;
    extern /* Subroutine */ int cdntp_(real *, integer *, integer *, integer *
	    , real *, real *, complex *, complex *);
    static integer kount;

/* ***BEGIN PROLOGUE  CDZRO */
/* ***SUBSIDIARY */
/* ***PURPOSE  CDZRO searches for a zero of a function F(N, T, Y, IROOT) */
/*            between the given values B and C until the width of the */
/*            interval (B, C) has collapsed to within a tolerance */
/*            specified by the stopping criterion, */
/*              ABS(B - C) .LE. 2.*(RW*ABS(B) + AE). */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDZRO-S, DDZRO-D, CDZRO-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*     This is a special purpose version of ZEROIN, modified for use with */
/*     the CDRIV package. */

/*     Sandia Mathematical Program Library */
/*     Mathematical Computing Services Division 5422 */
/*     Sandia Laboratories */
/*     P. O. Box 5800 */
/*     Albuquerque, New Mexico  87115 */
/*     Control Data 6600 Version 4.5, 1 November 1971 */

/*     PARAMETERS */
/*        F     - Name of the external function, which returns a */
/*                real result.  This name must be in an */
/*                EXTERNAL statement in the calling program. */
/*        B     - One end of the interval (B, C).  The value returned for */
/*                B usually is the better approximation to a zero of F. */
/*        C     - The other end of the interval (B, C). */
/*        RE    - Relative error used for RW in the stopping criterion. */
/*                If the requested RE is less than machine precision, */
/*                then RW is set to approximately machine precision. */
/*        AE    - Absolute error used in the stopping criterion.  If the */
/*                given interval (B, C) contains the origin, then a */
/*                nonzero value should be chosen for AE. */

/* ***REFERENCES  L. F. Shampine and H. A. Watts, ZEROIN, a root-solving */
/*                 routine, SC-TM-70-631, Sept 1970. */
/*               T. J. Dekker, Finding a zero by means of successive */
/*                 linear interpolation, Constructive Aspects of the */
/*                 Fundamental Theorem of Algebra, edited by B. Dejon */
/*                 and P. Henrici, 1969. */
/* ***ROUTINES CALLED  CDNTP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  CDZRO */
/* ***FIRST EXECUTABLE STATEMENT  CDZRO */
    /* Parameter adjustments */
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --y;

    /* Function Body */
    er = *uround * 4.f;
    rw = dmax(*re,er);
    ic = 0;
    acbs = (r__1 = *b - *c__, dabs(r__1));
    a = *c__;
    fa = *fc;
    kount = 0;
/*                                                    Perform interchange */
L10:
    if (dabs(*fc) < dabs(*fb)) {
	a = *b;
	fa = *fb;
	*b = *c__;
	*fb = *fc;
	*c__ = a;
	*fc = fa;
    }
    cmb = (*c__ - *b) * .5f;
    acmb = dabs(cmb);
    tol = rw * dabs(*b) + *ae;
/*                                                Test stopping criterion */
    if (acmb <= tol) {
	return 0;
    }
    if (kount > 50) {
	return 0;
    }
/*                                    Calculate new iterate implicitly as */
/*                                    B + P/Q, where we arrange P .GE. 0. */
/*                         The implicit form is used to prevent overflow. */
    p = (*b - a) * *fb;
    q = fa - *fb;
    if (p < 0.f) {
	p = -p;
	q = -q;
    }
/*                          Update A and check for satisfactory reduction */
/*                          in the size of our bounding interval. */
    a = *b;
    fa = *fb;
    ++ic;
    if (ic >= 4) {
	if (acmb * 8.f >= acbs) {
/*                                                                 Bisect */
	    *b = (*c__ + *b) * .5f;
	    goto L20;
	}
	ic = 0;
    }
    acbs = acmb;
/*                                            Test for too small a change */
    if (p <= dabs(q) * tol) {
/*                                                 Increment by tolerance */
	*b += r_sign(&tol, &cmb);
/*                                               Root ought to be between */
/*                                               B and (C + B)/2. */
    } else if (p < cmb * q) {
/*                                                            Interpolate */
	*b += p / q;
    } else {
/*                                                                 Bisect */
	*b = (*c__ + *b) * .5f;
    }
/*                                             Have completed computation */
/*                                             for new iterate B. */
L20:
    cdntp_(h__, &c__0, n, nq, t, b, &yh[yh_offset], &y[1]);
    *fb = (*f)(n, b, &y[1], iroot);
    if (*n == 0) {
	return 0;
    }
    if (*fb == 0.f) {
	return 0;
    }
    ++kount;

/*             Decide whether next step is interpolation or extrapolation */

    if (r_sign(&c_b5, fb) == r_sign(&c_b5, fc)) {
	*c__ = a;
	*fc = fa;
    }
    goto L10;
} /* cdzro_ */


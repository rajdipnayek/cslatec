/* cpqr79.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__3 = 3;
static complex c_b12 = {1.f,0.f};

/* DECK CPQR79 */
/* Subroutine */ int cpqr79_(integer *ndeg, complex *coeff, complex *root, 
	integer *ierr, real *work)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1, q__2;

    /* Local variables */
    static complex c__;
    static integer k, kj, km1, kad, khi, khr, kwi, kwr;
    static complex scale;
    extern /* Subroutine */ int comqr_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, integer *), xermsg_(
	    char *, char *, char *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);

/* ***BEGIN PROLOGUE  CPQR79 */
/* ***PURPOSE  Find the zeros of a polynomial with complex coefficients. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F1A1B */
/* ***TYPE      COMPLEX (RPQR79-S, CPQR79-C) */
/* ***KEYWORDS  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS */
/* ***AUTHOR  Vandevender, W. H., (SNLA) */
/* ***DESCRIPTION */

/*   Abstract */
/*       This routine computes all zeros of a polynomial of degree NDEG */
/*       with complex coefficients by computing the eigenvalues of the */
/*       companion matrix. */

/*   Description of Parameters */
/*       The user must dimension all arrays appearing in the call list */
/*            COEFF(NDEG+1), ROOT(NDEG), WORK(2*NDEG*(NDEG+1)) */

/*    --Input-- */
/*      NDEG    degree of polynomial */

/*      COEFF   COMPLEX coefficients in descending order.  i.e., */
/*              P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1) */

/*      WORK    REAL work array of dimension at least 2*NDEG*(NDEG+1) */

/*   --Output-- */
/*      ROOT    COMPLEX vector of roots */

/*      IERR    Output Error Code */
/*           - Normal Code */
/*          0  means the roots were computed. */
/*           - Abnormal Codes */
/*          1  more than 30 QR iterations on some eigenvalue of the */
/*             companion matrix */
/*          2  COEFF(1)=0.0 */
/*          3  NDEG is invalid (less than or equal to 0) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  COMQR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   911010  Code reworked and simplified.  (RWC and WRB) */
/* ***END PROLOGUE  CPQR79 */
/* ***FIRST EXECUTABLE STATEMENT  CPQR79 */
    /* Parameter adjustments */
    --work;
    --root;
    --coeff;

    /* Function Body */
    *ierr = 0;
    if (c_abs(&coeff[1]) == 0.f) {
	*ierr = 2;
	xermsg_("SLATEC", "CPQR79", "LEADING COEFFICIENT IS ZERO.", &c__2, &
		c__1, (ftnlen)6, (ftnlen)6, (ftnlen)28);
	return 0;
    }

    if (*ndeg <= 0) {
	*ierr = 3;
	xermsg_("SLATEC", "CPQR79", "DEGREE INVALID.", &c__3, &c__1, (ftnlen)
		6, (ftnlen)6, (ftnlen)15);
	return 0;
    }

    if (*ndeg == 1) {
	q__2.r = -coeff[2].r, q__2.i = -coeff[2].i;
	c_div(&q__1, &q__2, &coeff[1]);
	root[1].r = q__1.r, root[1].i = q__1.i;
	return 0;
    }

    c_div(&q__1, &c_b12, &coeff[1]);
    scale.r = q__1.r, scale.i = q__1.i;
    khr = 1;
    khi = khr + *ndeg * *ndeg;
    kwr = khi + khi - khr;
    kwi = kwr + *ndeg;

    i__1 = kwr;
    for (k = 1; k <= i__1; ++k) {
	work[k] = 0.f;
/* L10: */
    }

    i__1 = *ndeg;
    for (k = 1; k <= i__1; ++k) {
	kad = (k - 1) * *ndeg + 1;
	i__2 = k + 1;
	q__1.r = scale.r * coeff[i__2].r - scale.i * coeff[i__2].i, q__1.i = 
		scale.r * coeff[i__2].i + scale.i * coeff[i__2].r;
	c__.r = q__1.r, c__.i = q__1.i;
	work[kad] = -c__.r;
	kj = khi + kad - 1;
	work[kj] = -r_imag(&c__);
	if (k != *ndeg) {
	    work[kad + k] = 1.f;
	}
/* L20: */
    }

    comqr_(ndeg, ndeg, &c__1, ndeg, &work[khr], &work[khi], &work[kwr], &work[
	    kwi], ierr);

    if (*ierr != 0) {
	*ierr = 1;
	xermsg_("SLATEC", "CPQR79", "NO CONVERGENCE IN 30 QR ITERATIONS.", &
		c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)35);
	return 0;
    }

    i__1 = *ndeg;
    for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	i__2 = k;
	i__3 = kwr + km1;
	i__4 = kwi + km1;
	q__1.r = work[i__3], q__1.i = work[i__4];
	root[i__2].r = q__1.r, root[i__2].i = q__1.i;
/* L30: */
    }
    return 0;
} /* cpqr79_ */


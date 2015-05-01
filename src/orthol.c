/* orthol.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;
static real c_b9 = 10.f;
static integer c__0 = 0;

/* DECK ORTHOL */
/* Subroutine */ int orthol_(real *a, integer *m, integer *n, integer *nrda, 
	integer *iflag, integer *irank, integer *iscale, real *diag, integer *
	kpivot, real *scales, real *cols, real *cs)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer j, k, l;
    static real as, sc;
    static integer mk, kp;
    static real acc, akk, sad, sig, dum, css, uro;
    static integer jcol;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real diagk, asave, sigma, anorm, sruro;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int cscale_(real *, integer *, integer *, integer 
	    *, real *, real *, real *, real *, real *, real *, integer *, 
	    integer *), xermsg_(char *, char *, char *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  ORTHOL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (ORTHOL-S) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   Reduction of the matrix A to upper triangular form by a sequence of */
/*   orthogonal HOUSEHOLDER transformations pre-multiplying A */

/*   Modeled after the ALGOL codes in the articles in the REFERENCES */
/*   section. */

/* ********************************************************************** */
/*   INPUT */
/* ********************************************************************** */

/*     A -- Contains the matrix to be decomposed, must be dimensioned */
/*           NRDA by N */
/*     M -- Number of rows in the matrix, M greater or equal to N */
/*     N -- Number of columns in the matrix, N greater or equal to 1 */
/*     IFLAG -- Indicates the uncertainty in the matrix data */
/*             = 0 when the data is to be treated as exact */
/*             =-K when the data is assumed to be accurate to about */
/*                 K digits */
/*     ISCALE -- Scaling indicator */
/*               =-1 if the matrix A is to be pre-scaled by */
/*               columns when appropriate. */
/*               Otherwise no scaling will be attempted */
/*     NRDA -- Row dimension of A, NRDA greater or equal to M */
/*     DIAG,KPIVOT,COLS -- Arrays of length at least n used internally */
/*         ,CS,SCALES */

/* ********************************************************************** */
/*   OUTPUT */
/* ********************************************************************** */

/*     IFLAG - Status indicator */
/*            =1 for successful decomposition */
/*            =2 if improper input is detected */
/*            =3 if rank of the matrix is less than N */
/*     A -- Contains the reduced matrix in the strictly upper triangular */
/*          part and transformation information in the lower part */
/*     IRANK -- Contains the numerically determined matrix rank */
/*     DIAG -- Contains the diagonal elements of the reduced */
/*             triangular matrix */
/*     KPIVOT -- Contains the pivotal information, the column */
/*               interchanges performed on the original matrix are */
/*               recorded here. */
/*     SCALES -- Contains the column scaling parameters */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***REFERENCES  G. Golub, Numerical methods for solving linear least */
/*                 squares problems, Numerische Mathematik 7, (1965), */
/*                 pp. 206-216. */
/*               P. Businger and G. Golub, Linear least squares */
/*                 solutions by Householder transformations, Numerische */
/*                 Mathematik  7, (1965), pp. 269-276. */
/* ***ROUTINES CALLED  CSCALE, R1MACH, SDOT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900402  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  ORTHOL */

/* ********************************************************************** */

/*     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED */
/*     BY THE FUNCTION R1MACH. */

/* ***FIRST EXECUTABLE STATEMENT  ORTHOL */
    /* Parameter adjustments */
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --diag;
    --kpivot;
    --scales;
    --cols;
    --cs;

    /* Function Body */
    uro = r1mach_(&c__3);

/* ********************************************************************** */

    if (*m >= *n && *n >= 1 && *nrda >= *m) {
	goto L1;
    }
    *iflag = 2;
    xermsg_("SLATEC", "ORTHOL", "INVALID INPUT PARAMETERS.", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;

L1:
    acc = uro * 10.f;
    if (*iflag < 0) {
/* Computing MAX */
	r__1 = acc, r__2 = pow_ri(&c_b9, iflag);
	acc = dmax(r__1,r__2);
    }
    sruro = sqrt(uro);
    *iflag = 1;
    *irank = *n;

/*     COMPUTE NORM**2 OF JTH COLUMN AND A MATRIX NORM */

    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	kpivot[j] = j;
	cols[j] = sdot_(m, &a[j * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	cs[j] = cols[j];
	anorm += cols[j];
/* L2: */
    }

/*     PERFORM COLUMN SCALING ON A WHEN SPECIFIED */

    cscale_(&a[a_offset], nrda, m, n, &cols[1], &cs[1], &dum, &dum, &anorm, &
	    scales[1], iscale, &c__0);

    anorm = sqrt(anorm);


/*     CONSTRUCTION OF UPPER TRIANGULAR MATRIX AND RECORDING OF */
/*     ORTHOGONAL TRANSFORMATIONS */


    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	mk = *m - k + 1;
	if (k == *n) {
	    goto L25;
	}
	kp = k + 1;

/*        SEARCHING FOR PIVOTAL COLUMN */

	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
	    if (cols[j] >= sruro * cs[j]) {
		goto L5;
	    }
	    cols[j] = sdot_(&mk, &a[k + j * a_dim1], &c__1, &a[k + j * a_dim1]
		    , &c__1);
	    cs[j] = cols[j];
L5:
	    if (j == k) {
		goto L7;
	    }
	    if (sigma >= cols[j] * .99f) {
		goto L10;
	    }
L7:
	    sigma = cols[j];
	    jcol = j;
L10:
	    ;
	}
	if (jcol == k) {
	    goto L25;
	}

/*        PERFORM COLUMN INTERCHANGE */

	l = kpivot[k];
	kpivot[k] = kpivot[jcol];
	kpivot[jcol] = l;
	cols[jcol] = cols[k];
	cols[k] = sigma;
	css = cs[k];
	cs[k] = cs[jcol];
	cs[jcol] = css;
	sc = scales[k];
	scales[k] = scales[jcol];
	scales[jcol] = sc;
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    asave = a[l + k * a_dim1];
	    a[l + k * a_dim1] = a[l + jcol * a_dim1];
/* L20: */
	    a[l + jcol * a_dim1] = asave;
	}

/*        CHECK RANK OF THE MATRIX */

L25:
	sig = sdot_(&mk, &a[k + k * a_dim1], &c__1, &a[k + k * a_dim1], &c__1)
		;
	diagk = sqrt(sig);
	if (diagk > acc * anorm) {
	    goto L30;
	}

/*        RANK DEFICIENT PROBLEM */
	*iflag = 3;
	*irank = k - 1;
	xermsg_("SLATEC", "ORTHOL", "RANK OF MATRIX IS LESS THAN THE NUMBER "
		"OF COLUMNS.", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)50);
	return 0;

/*        CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A */

L30:
	akk = a[k + k * a_dim1];
	if (akk > 0.f) {
	    diagk = -diagk;
	}
	diag[k] = diagk;
	a[k + k * a_dim1] = akk - diagk;
	if (k == *n) {
	    goto L50;
	}
	sad = diagk * akk - sig;
	i__2 = *n;
	for (j = kp; j <= i__2; ++j) {
	    as = sdot_(&mk, &a[k + k * a_dim1], &c__1, &a[k + j * a_dim1], &
		    c__1) / sad;
	    i__3 = *m;
	    for (l = k; l <= i__3; ++l) {
/* L35: */
		a[l + j * a_dim1] += as * a[l + k * a_dim1];
	    }
/* L40: */
/* Computing 2nd power */
	    r__1 = a[k + j * a_dim1];
	    cols[j] -= r__1 * r__1;
	}
L50:
	;
    }


    return 0;
} /* orthol_ */


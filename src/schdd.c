/* schdd.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;

/* DECK SCHDD */
/* Subroutine */ int schdd_(real *r__, integer *ldr, integer *p, real *x, 
	real *z__, integer *ldz, integer *nz, real *y, real *rho, real *c__, 
	real *s, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real a, b;
    static integer i__, j;
    static real t;
    static integer ii;
    static real xx, zeta;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real norm;
    extern doublereal snrm2_(integer *, real *, integer *);
    static real alpha, scale, azeta;

/* ***BEGIN PROLOGUE  SCHDD */
/* ***PURPOSE  Downdate an augmented Cholesky decomposition or the */
/*            triangular factor of an augmented QR decomposition. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D7B */
/* ***TYPE      SINGLE PRECISION (SCHDD-S, DCHDD-D, CCHDD-C) */
/* ***KEYWORDS  CHOLESKY DECOMPOSITION, DOWNDATE, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     SCHDD downdates an augmented Cholesky decomposition or the */
/*     triangular factor of an augmented QR decomposition. */
/*     Specifically, given an upper triangular matrix R of order P, a */
/*     row vector X, a column vector Z, and a scalar Y, SCHDD */
/*     determines an orthogonal matrix U and a scalar ZETA such that */

/*                        (R   Z )     (RR  ZZ) */
/*                    U * (      )  =  (      ) , */
/*                        (0 ZETA)     ( X   Y) */

/*     where RR is upper triangular.  If R and Z have been obtained */
/*     from the factorization of a least squares problem, then */
/*     RR and ZZ are the factors corresponding to the problem */
/*     with the observation (X,Y) removed.  In this case, if RHO */
/*     is the norm of the residual vector, then the norm of */
/*     the residual vector of the downdated problem is */
/*     SQRT(RHO**2 - ZETA**2). SCHDD will simultaneously downdate */
/*     several triplets (Z,Y,RHO) along with R. */
/*     For a less terse description of what SCHDD does and how */
/*     it may be applied, see the LINPACK guide. */

/*     The matrix U is determined as the product U(1)*...*U(P) */
/*     where U(I) is a rotation in the (P+1,I)-plane of the */
/*     form */

/*                       ( C(I)     -S(I)     ) */
/*                       (                    ) . */
/*                       ( S(I)       C(I)    ) */

/*     The rotations are chosen so that C(I) is real. */

/*     The user is warned that a given downdating problem may */
/*     be impossible to accomplish or may produce */
/*     inaccurate results.  For example, this can happen */
/*     if X is near a vector whose removal will reduce the */
/*     rank of R.  Beware. */

/*     On Entry */

/*         R      REAL(LDR,P), where LDR .GE. P. */
/*                R contains the upper triangular matrix */
/*                that is to be downdated.  The part of  R */
/*                below the diagonal is not referenced. */

/*         LDR    INTEGER. */
/*                LDR is the leading dimension of the array R. */

/*         P      INTEGER. */
/*                P is the order of the matrix R. */

/*         X      REAL(P). */
/*                X contains the row vector that is to */
/*                be removed from R.  X is not altered by SCHDD. */

/*         Z      REAL(LDZ,NZ), where LDZ .GE. P. */
/*                Z is an array of NZ P-vectors which */
/*                are to be downdated along with R. */

/*         LDZ    INTEGER. */
/*                LDZ is the leading dimension of the array Z. */

/*         NZ     INTEGER. */
/*                NZ is the number of vectors to be downdated */
/*                NZ may be zero, in which case Z, Y, and RHO */
/*                are not referenced. */

/*         Y      REAL(NZ). */
/*                Y contains the scalars for the downdating */
/*                of the vectors Z.  Y is not altered by SCHDD. */

/*         RHO    REAL(NZ). */
/*                RHO contains the norms of the residual */
/*                vectors that are to be downdated. */

/*     On Return */

/*         R */
/*         Z      contain the downdated quantities. */
/*         RHO */

/*         C      REAL(P). */
/*                C contains the cosines of the transforming */
/*                rotations. */

/*         S      REAL(P). */
/*                S contains the sines of the transforming */
/*                rotations. */

/*         INFO   INTEGER. */
/*                INFO is set as follows. */

/*                   INFO = 0  if the entire downdating */
/*                             was successful. */

/*                   INFO =-1  if R could not be downdated. */
/*                             In this case, all quantities */
/*                             are left unaltered. */

/*                   INFO = 1  if some RHO could not be */
/*                             downdated.  The offending RHOs are */
/*                             set to -1. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SDOT, SNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SCHDD */


/*     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT */
/*     IN THE ARRAY S. */

/* ***FIRST EXECUTABLE STATEMENT  SCHDD */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    --rho;
    --c__;
    --s;

    /* Function Body */
    *info = 0;
    s[1] = x[1] / r__[r_dim1 + 1];
    if (*p < 2) {
	goto L20;
    }
    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	s[j] = x[j] - sdot_(&i__2, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	s[j] /= r__[j + j * r_dim1];
/* L10: */
    }
L20:
    norm = snrm2_(p, &s[1], &c__1);
    if (norm < 1.f) {
	goto L30;
    }
    *info = -1;
    goto L120;
L30:
/* Computing 2nd power */
    r__1 = norm;
    alpha = sqrt(1.f - r__1 * r__1);

/*        DETERMINE THE TRANSFORMATIONS. */

    i__1 = *p;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *p - ii + 1;
	scale = alpha + (r__1 = s[i__], dabs(r__1));
	a = alpha / scale;
	b = s[i__] / scale;
/* Computing 2nd power */
	r__1 = a;
/* Computing 2nd power */
	r__2 = b;
	norm = sqrt(r__1 * r__1 + r__2 * r__2);
	c__[i__] = a / norm;
	s[i__] = b / norm;
	alpha = scale * norm;
/* L40: */
    }

/*        APPLY THE TRANSFORMATIONS TO R. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx = 0.f;
	i__2 = j;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = j - ii + 1;
	    t = c__[i__] * xx + s[i__] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = c__[i__] * r__[i__ + j * r_dim1] - s[i__] 
		    * xx;
	    xx = t;
/* L50: */
	}
/* L60: */
    }

/*        IF REQUIRED, DOWNDATE Z AND RHO. */

    if (*nz < 1) {
	goto L110;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	zeta = y[j];
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__ + j * z_dim1] = (z__[i__ + j * z_dim1] - s[i__] * zeta) / 
		    c__[i__];
	    zeta = c__[i__] * zeta - s[i__] * z__[i__ + j * z_dim1];
/* L70: */
	}
	azeta = dabs(zeta);
	if (azeta <= rho[j]) {
	    goto L80;
	}
	*info = 1;
	rho[j] = -1.f;
	goto L90;
L80:
/* Computing 2nd power */
	r__1 = azeta / rho[j];
	rho[j] *= sqrt(1.f - r__1 * r__1);
L90:
/* L100: */
	;
    }
L110:
L120:
    return 0;
} /* schdd_ */


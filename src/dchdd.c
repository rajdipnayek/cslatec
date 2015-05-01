/* dchdd.f -- translated by f2c (version 12.02.01).
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

/* DECK DCHDD */
/* Subroutine */ int dchdd_(doublereal *r__, integer *ldr, integer *p, 
	doublereal *x, doublereal *z__, integer *ldz, integer *nz, doublereal 
	*y, doublereal *rho, doublereal *c__, doublereal *s, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, b;
    static integer i__, j;
    static doublereal t;
    static integer ii;
    static doublereal xx;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal zeta, norm;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal alpha, scale, azeta;

/* ***BEGIN PROLOGUE  DCHDD */
/* ***PURPOSE  Downdate an augmented Cholesky decomposition or the */
/*            triangular factor of an augmented QR decomposition. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D7B */
/* ***TYPE      DOUBLE PRECISION (SCHDD-S, DCHDD-D, CCHDD-C) */
/* ***KEYWORDS  CHOLESKY DECOMPOSITION, DOWNDATE, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     DCHDD downdates an augmented Cholesky decomposition or the */
/*     triangular factor of an augmented QR decomposition. */
/*     Specifically, given an upper triangular matrix R of order P,  a */
/*     row vector X, a column vector Z, and a scalar Y, DCHDD */
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
/*     SQRT(RHO**2 - ZETA**2).  DCHDD will simultaneously downdate */
/*     several triplets (Z,Y,RHO) along with R. */
/*     For a less terse description of what DCHDD does and how */
/*     it may be applied, see the LINPACK guide. */

/*     The matrix U is determined as the product U(1)*...*U(P) */
/*     where U(I) is a rotation in the (P+1,I)-plane of the */
/*     form */

/*                       ( C(I)     -S(I)     ) */
/*                       (                    ) . */
/*                       ( S(I)       C(I)    ) */

/*     The rotations are chosen so that C(I) is double precision. */

/*     The user is warned that a given downdating problem may */
/*     be impossible to accomplish or may produce */
/*     inaccurate results.  For example, this can happen */
/*     if X is near a vector whose removal will reduce the */
/*     rank of R.  Beware. */

/*     On Entry */

/*         R      DOUBLE PRECISION(LDR,P), where LDR .GE. P. */
/*                R contains the upper triangular matrix */
/*                that is to be downdated.  The part of  R */
/*                below the diagonal is not referenced. */

/*         LDR    INTEGER. */
/*                LDR is the leading dimension of the array R. */

/*         P      INTEGER. */
/*                P is the order of the matrix R. */

/*         X      DOUBLE PRECISION(P). */
/*                X contains the row vector that is to */
/*                be removed from R.  X is not altered by DCHDD. */

/*         Z      DOUBLE PRECISION(LDZ,N)Z), where LDZ .GE. P. */
/*                Z is an array of NZ P-vectors which */
/*                are to be downdated along with R. */

/*         LDZ    INTEGER. */
/*                LDZ is the leading dimension of the array Z. */

/*         NZ     INTEGER. */
/*                NZ is the number of vectors to be downdated */
/*                NZ may be zero, in which case Z, Y, and RHO */
/*                are not referenced. */

/*         Y      DOUBLE PRECISION(NZ). */
/*                Y contains the scalars for the downdating */
/*                of the vectors Z.  Y is not altered by DCHDD. */

/*         RHO    DOUBLE PRECISION(NZ). */
/*                RHO contains the norms of the residual */
/*                vectors that are to be downdated. */

/*     On Return */

/*         R */
/*         Z      contain the downdated quantities. */
/*         RHO */

/*         C      DOUBLE PRECISION(P). */
/*                C contains the cosines of the transforming */
/*                rotations. */

/*         S      DOUBLE PRECISION(P). */
/*                S contains the sines of the transforming */
/*                rotations. */

/*         INFO   INTEGER. */
/*                INFO is set as follows. */

/*                   INFO = 0  if the entire downdating */
/*                             was successful. */

/*                   INFO =-1  if R could not be downdated. */
/*                             in this case, all quantities */
/*                             are left unaltered. */

/*                   INFO = 1  if some RHO could not be */
/*                             downdated.  The offending RHO's are */
/*                             set to -1. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DDOT, DNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DCHDD */


/*     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT */
/*     IN THE ARRAY S. */

/* ***FIRST EXECUTABLE STATEMENT  DCHDD */
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
	s[j] = x[j] - ddot_(&i__2, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	s[j] /= r__[j + j * r_dim1];
/* L10: */
    }
L20:
    norm = dnrm2_(p, &s[1], &c__1);
    if (norm < 1.) {
	goto L30;
    }
    *info = -1;
    goto L120;
L30:
/* Computing 2nd power */
    d__1 = norm;
    alpha = sqrt(1. - d__1 * d__1);

/*        DETERMINE THE TRANSFORMATIONS. */

    i__1 = *p;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *p - ii + 1;
	scale = alpha + (d__1 = s[i__], abs(d__1));
	a = alpha / scale;
	b = s[i__] / scale;
/* Computing 2nd power */
	d__1 = a;
/* Computing 2nd power */
	d__2 = b;
	norm = sqrt(d__1 * d__1 + d__2 * d__2);
	c__[i__] = a / norm;
	s[i__] = b / norm;
	alpha = scale * norm;
/* L40: */
    }

/*        APPLY THE TRANSFORMATIONS TO R. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx = 0.;
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
	azeta = abs(zeta);
	if (azeta <= rho[j]) {
	    goto L80;
	}
	*info = 1;
	rho[j] = -1.;
	goto L90;
L80:
/* Computing 2nd power */
	d__1 = azeta / rho[j];
	rho[j] *= sqrt(1. - d__1 * d__1);
L90:
/* L100: */
	;
    }
L110:
L120:
    return 0;
} /* dchdd_ */


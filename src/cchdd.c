/* cchdd.f -- translated by f2c (version 12.02.01).
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

/* DECK CCHDD */
/* Subroutine */ int cchdd_(complex *r__, integer *ldr, integer *p, complex *
	x, complex *z__, integer *ldz, integer *nz, complex *y, real *rho, 
	real *c__, complex *s, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static real a;
    static complex b;
    static integer i__, j;
    static complex t;
    static integer ii;
    static complex xx, zeta;
    static real norm, alpha, scale;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    static real azeta;
    extern doublereal scnrm2_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CCHDD */
/* ***PURPOSE  Downdate an augmented Cholesky decomposition or the */
/*            triangular factor of an augmented QR decomposition. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D7B */
/* ***TYPE      COMPLEX (SCHDD-S, DCHDD-D, CCHDD-C) */
/* ***KEYWORDS  CHOLESKY DECOMPOSITION, DOWNDATE, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     CCHDD downdates an augmented Cholesky decomposition or the */
/*     triangular factor of an augmented QR decomposition. */
/*     Specifically, given an upper triangular matrix R of order P,  a */
/*     row vector X, a column vector Z, and a scalar Y, CCHDD */
/*     determines a unitary matrix U and a scalar ZETA such that */

/*                        (R   Z )     (RR  ZZ) */
/*                    U * (      )  =  (      ) , */
/*                        (0 ZETA)     ( X   Y) */

/*     where RR is upper triangular.  If R and Z have been obtained */
/*     from the factorization of a least squares problem, then */
/*     RR and ZZ are the factors corresponding to the problem */
/*     with the observation (X,Y) removed.  In this case, if RHO */
/*     is the norm of the residual vector, then the norm of */
/*     the residual vector of the downdated problem is */
/*     SQRT(RHO**2 - ZETA**2).  CCHDD will simultaneously downdate */
/*     several triplets (Z,Y,RHO) along with R. */
/*     For a less terse description of what CCHDD does and how */
/*     it may be applied, see the LINPACK Guide. */

/*     The matrix U is determined as the product U(1)*...*U(P) */
/*     where U(I) is a rotation in the (P+1,I)-plane of the */
/*     form */

/*                       ( C(I)  -CONJG(S(I)) ) */
/*                       (                    ) . */
/*                       ( S(I)       C(I)    ) */

/*     the rotations are chosen so that C(I) is real. */

/*     The user is warned that a given downdating problem may */
/*     be impossible to accomplish or may produce */
/*     inaccurate results.  For example, this can happen */
/*     if X is near a vector whose removal will reduce the */
/*     rank of R.  Beware. */

/*     On Entry */

/*         R      COMPLEX(LDR,P), where LDR .GE. P. */
/*                R contains the upper triangular matrix */
/*                that is to be downdated.  The part of R */
/*                below the diagonal is not referenced. */

/*         LDR    INTEGER. */
/*                LDR is the leading dimension of the array R. */

/*         p      INTEGER. */
/*                P is the order of the matrix R. */

/*         X      COMPLEX(P). */
/*                X contains the row vector that is to */
/*                be removed from R.  X is not altered by CCHDD. */

/*         Z      COMPLEX(LDZ,NZ), where LDZ .GE. P. */
/*                Z is an array of NZ P-vectors which */
/*                are to be downdated along with R. */

/*         LDZ    INTEGER. */
/*                LDZ is the leading dimension of the array Z. */

/*         NZ     INTEGER. */
/*                NZ is the number of vectors to be downdated */
/*                NZ may be zero, in which case Z, Y, and RHO */
/*                are not referenced. */

/*         Y      COMPLEX(NZ). */
/*                Y contains the scalars for the downdating */
/*                of the vectors Z.  Y is not altered by CCHDD. */

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

/*         S      COMPLEX(P). */
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
/* ***ROUTINES CALLED  CDOTC, SCNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CCHDD */


/*     SOLVE THE SYSTEM CTRANS(R)*A = X, PLACING THE RESULT */
/*     IN THE ARRAY S. */

/* ***FIRST EXECUTABLE STATEMENT  CCHDD */
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
    r_cnjg(&q__2, &x[1]);
    r_cnjg(&q__3, &r__[r_dim1 + 1]);
    c_div(&q__1, &q__2, &q__3);
    s[1].r = q__1.r, s[1].i = q__1.i;
    if (*p < 2) {
	goto L20;
    }
    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	r_cnjg(&q__2, &x[j]);
	i__3 = j - 1;
	cdotc_(&q__3, &i__3, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = j;
	r_cnjg(&q__2, &r__[j + j * r_dim1]);
	c_div(&q__1, &s[j], &q__2);
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
/* L10: */
    }
L20:
    norm = scnrm2_(p, &s[1], &c__1);
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
	scale = alpha + c_abs(&s[i__]);
	a = alpha / scale;
	i__2 = i__;
	q__1.r = s[i__2].r / scale, q__1.i = s[i__2].i / scale;
	b.r = q__1.r, b.i = q__1.i;
/* Computing 2nd power */
	r__1 = a;
/* Computing 2nd power */
	r__2 = b.r;
/* Computing 2nd power */
	r__3 = r_imag(&b);
	norm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	c__[i__] = a / norm;
	i__2 = i__;
	r_cnjg(&q__2, &b);
	q__1.r = q__2.r / norm, q__1.i = q__2.i / norm;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	alpha = scale * norm;
/* L40: */
    }

/*        APPLY THE TRANSFORMATIONS TO R. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx.r = 0.f, xx.i = 0.f;
	i__2 = j;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = j - ii + 1;
	    i__3 = i__;
	    q__2.r = c__[i__3] * xx.r, q__2.i = c__[i__3] * xx.i;
	    i__4 = i__;
	    i__5 = i__ + j * r_dim1;
	    q__3.r = s[i__4].r * r__[i__5].r - s[i__4].i * r__[i__5].i, 
		    q__3.i = s[i__4].r * r__[i__5].i + s[i__4].i * r__[i__5]
		    .r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + j * r_dim1;
	    i__4 = i__;
	    i__5 = i__ + j * r_dim1;
	    q__2.r = c__[i__4] * r__[i__5].r, q__2.i = c__[i__4] * r__[i__5]
		    .i;
	    r_cnjg(&q__4, &s[i__]);
	    q__3.r = q__4.r * xx.r - q__4.i * xx.i, q__3.i = q__4.r * xx.i + 
		    q__4.i * xx.r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    r__[i__3].r = q__1.r, r__[i__3].i = q__1.i;
	    xx.r = t.r, xx.i = t.i;
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
	i__2 = j;
	zeta.r = y[i__2].r, zeta.i = y[i__2].i;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * z_dim1;
	    i__4 = i__ + j * z_dim1;
	    r_cnjg(&q__4, &s[i__]);
	    q__3.r = q__4.r * zeta.r - q__4.i * zeta.i, q__3.i = q__4.r * 
		    zeta.i + q__4.i * zeta.r;
	    q__2.r = z__[i__4].r - q__3.r, q__2.i = z__[i__4].i - q__3.i;
	    i__5 = i__;
	    q__1.r = q__2.r / c__[i__5], q__1.i = q__2.i / c__[i__5];
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = i__;
	    q__2.r = c__[i__3] * zeta.r, q__2.i = c__[i__3] * zeta.i;
	    i__4 = i__;
	    i__5 = i__ + j * z_dim1;
	    q__3.r = s[i__4].r * z__[i__5].r - s[i__4].i * z__[i__5].i, 
		    q__3.i = s[i__4].r * z__[i__5].i + s[i__4].i * z__[i__5]
		    .r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    zeta.r = q__1.r, zeta.i = q__1.i;
/* L70: */
	}
	azeta = c_abs(&zeta);
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
} /* cchdd_ */


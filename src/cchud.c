/* cchud.f -- translated by f2c (version 12.02.01).
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

/* DECK CCHUD */
/* Subroutine */ int cchud_(complex *r__, integer *ldr, integer *p, complex *
	x, complex *z__, integer *ldz, integer *nz, complex *y, real *rho, 
	real *c__, complex *s)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, j;
    static complex t, xj;
    static integer jm1;
    static complex zeta;
    static real scale, azeta;
    extern /* Subroutine */ int crotg_(complex *, complex *, real *, complex *
	    );

/* ***BEGIN PROLOGUE  CCHUD */
/* ***PURPOSE  Update an augmented Cholesky decomposition of the */
/*            triangular part of an augmented QR decomposition. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D7B */
/* ***TYPE      COMPLEX (SCHUD-S, DCHUD-D, CCHUD-C) */
/* ***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             UPDATE */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     CCHUD updates an augmented Cholesky decomposition of the */
/*     triangular part of an augmented QR decomposition.  Specifically, */
/*     given an upper triangular matrix R of order P, a row vector */
/*     X, a column vector Z, and a scalar Y, CCHUD determines a */
/*     unitary matrix U and a scalar ZETA such that */


/*                              (R  Z)     (RR   ZZ ) */
/*                         U  * (    )  =  (        ) , */
/*                              (X  Y)     ( 0  ZETA) */

/*     where RR is upper triangular.  If R and Z have been */
/*     obtained from the factorization of a least squares */
/*     problem, then RR and ZZ are the factors corresponding to */
/*     the problem with the observation (X,Y) appended.  In this */
/*     case, if RHO is the norm of the residual vector, then the */
/*     norm of the residual vector of the updated problem is */
/*     SQRT(RHO**2 + ZETA**2).  CCHUD will simultaneously update */
/*     several triplets (Z,Y,RHO). */

/*     For a less terse description of what CCHUD does and how */
/*     it may be applied see the LINPACK Guide. */

/*     The matrix U is determined as the product U(P)*...*U(1), */
/*     where U(I) is a rotation in the (I,P+1) plane of the */
/*     form */

/*                       (     (CI)      S(I) ) */
/*                       (                    ) . */
/*                       ( -CONJG(S(I))  (CI) ) */

/*     The rotations are chosen so that C(I) is real. */

/*     On Entry */

/*         R      COMPLEX(LDR,P), where LDR .GE. P. */
/*                R contains the upper triangular matrix */
/*                that is to be updated.  The part of R */
/*                below the diagonal is not referenced. */

/*         LDR    INTEGER. */
/*                LDR is the leading dimension of the array R. */

/*         P      INTEGER. */
/*                P is the order of the matrix R. */

/*         X      COMPLEX(P). */
/*                X contains the row to be added to R.  X is */
/*                not altered by CCHUD. */

/*         Z      COMPLEX(LDZ,NZ), where LDZ .GE. P. */
/*                Z is an array containing NZ P-vectors to */
/*                be updated with R. */

/*         LDZ    INTEGER. */
/*                LDZ is the leading dimension of the array Z. */

/*         NZ     INTEGER. */
/*                NZ is the number of vectors to be updated */
/*                NZ may be zero, in which case Z, Y, and RHO */
/*                are not referenced. */

/*         Y      COMPLEX(NZ). */
/*                Y contains the scalars for updating the vectors */
/*                Z.  Y is not altered by CCHUD. */

/*         RHO    REAL(NZ). */
/*                RHO contains the norms of the residual */
/*                vectors that are to be updated.  If RHO(J) */
/*                is negative, it is left unaltered. */

/*     On Return */

/*         RC */
/*         RHO    contain the updated quantities. */
/*         Z */

/*         C      REAL(P). */
/*                C contains the cosines of the transforming */
/*                rotations. */

/*         S      COMPLEX(P). */
/*                S contains the sines of the transforming */
/*                rotations. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CROTG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CCHUD */


/*     UPDATE R. */

/* ***FIRST EXECUTABLE STATEMENT  CCHUD */
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
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	xj.r = x[i__2].r, xj.i = x[i__2].i;

/*        APPLY THE PREVIOUS ROTATIONS. */

	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + j * r_dim1;
	    q__2.r = c__[i__3] * r__[i__4].r, q__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = i__;
	    q__3.r = s[i__5].r * xj.r - s[i__5].i * xj.i, q__3.i = s[i__5].r *
		     xj.i + s[i__5].i * xj.r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__;
	    q__2.r = c__[i__3] * xj.r, q__2.i = c__[i__3] * xj.i;
	    r_cnjg(&q__4, &s[i__]);
	    i__4 = i__ + j * r_dim1;
	    q__3.r = q__4.r * r__[i__4].r - q__4.i * r__[i__4].i, q__3.i = 
		    q__4.r * r__[i__4].i + q__4.i * r__[i__4].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    xj.r = q__1.r, xj.i = q__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
/* L10: */
	}
L20:

/*        COMPUTE THE NEXT ROTATION. */

	crotg_(&r__[j + j * r_dim1], &xj, &c__[j], &s[j]);
/* L30: */
    }

/*     IF REQUIRED, UPDATE Z AND RHO. */

    if (*nz < 1) {
	goto L70;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	zeta.r = y[i__2].r, zeta.i = y[i__2].i;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + j * z_dim1;
	    q__2.r = c__[i__3] * z__[i__4].r, q__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = i__;
	    q__3.r = s[i__5].r * zeta.r - s[i__5].i * zeta.i, q__3.i = s[i__5]
		    .r * zeta.i + s[i__5].i * zeta.r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__;
	    q__2.r = c__[i__3] * zeta.r, q__2.i = c__[i__3] * zeta.i;
	    r_cnjg(&q__4, &s[i__]);
	    i__4 = i__ + j * z_dim1;
	    q__3.r = q__4.r * z__[i__4].r - q__4.i * z__[i__4].i, q__3.i = 
		    q__4.r * z__[i__4].i + q__4.i * z__[i__4].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    zeta.r = q__1.r, zeta.i = q__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
/* L40: */
	}
	azeta = c_abs(&zeta);
	if (azeta == 0.f || rho[j] < 0.f) {
	    goto L50;
	}
	scale = azeta + rho[j];
/* Computing 2nd power */
	r__1 = azeta / scale;
/* Computing 2nd power */
	r__2 = rho[j] / scale;
	rho[j] = scale * sqrt(r__1 * r__1 + r__2 * r__2);
L50:
/* L60: */
	;
    }
L70:
    return 0;
} /* cchud_ */


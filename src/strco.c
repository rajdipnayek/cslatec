/* strco.f -- translated by f2c (version 12.02.01).
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

/* DECK STRCO */
/* Subroutine */ int strco_(real *t, integer *ldt, integer *n, real *rcond, 
	real *z__, integer *job)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer j, k, l;
    static real s, w;
    static integer i1, j1, j2;
    static real ek;
    static integer kk;
    static real sm, wk, wkm;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    extern doublereal sasum_(integer *, real *, integer *);
    static logical lower;
    static real tnorm, ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  STRCO */
/* ***PURPOSE  Estimate the condition number of a triangular matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2A3 */
/* ***TYPE      SINGLE PRECISION (STRCO-S, DTRCO-D, CTRCO-C) */
/* ***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK, */
/*             TRIANGULAR MATRIX */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     STRCO estimates the condition of a real triangular matrix. */

/*     On Entry */

/*        T       REAL(LDT,N) */
/*                T contains the triangular matrix.  The zero */
/*                elements of the matrix are not referenced, and */
/*                the corresponding elements of the array can be */
/*                used to store other information. */

/*        LDT     INTEGER */
/*                LDT is the leading dimension of the array T. */

/*        N       INTEGER */
/*                N is the order of the system. */

/*        JOB     INTEGER */
/*                = 0         T  is lower triangular. */
/*                = nonzero   T  is upper triangular. */

/*     On Return */

/*        RCOND   REAL */
/*                an estimate of the reciprocal condition of  T . */
/*                For the system  T*X = B , relative perturbations */
/*                in  T  and  B  of size  EPSILON  may cause */
/*                relative perturbations in  X  of size  EPSILON/RCOND . */
/*                If  RCOND  is so small that the logical expression */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                is true, then  T  may be singular to working */
/*                precision.  In particular,  RCOND  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        Z       REAL(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  T  is close to a singular matrix, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SASUM, SAXPY, SSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  STRCO */

/* ***FIRST EXECUTABLE STATEMENT  STRCO */
    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --z__;

    /* Function Body */
    lower = *job == 0;

/*     COMPUTE 1-NORM OF T */

    tnorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	if (lower) {
	    l = *n + 1 - j;
	}
	i1 = 1;
	if (lower) {
	    i1 = j;
	}
/* Computing MAX */
	r__1 = tnorm, r__2 = sasum_(&l, &t[i1 + j * t_dim1], &c__1);
	tnorm = dmax(r__1,r__2);
/* L10: */
    }

/*     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E . */
/*     TRANS(T)  IS THE TRANSPOSE OF T . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF Y . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE TRANS(T)*Y = E */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L20: */
    }
    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = kk;
	if (lower) {
	    k = *n + 1 - kk;
	}
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= (r__2 = t[k + k * t_dim1], 
		dabs(r__2))) {
	    goto L30;
	}
	s = (r__1 = t[k + k * t_dim1], dabs(r__1)) / (r__2 = ek - z__[k], 
		dabs(r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	if (t[k + k * t_dim1] == 0.f) {
	    goto L40;
	}
	wk /= t[k + k * t_dim1];
	wkm /= t[k + k * t_dim1];
	goto L50;
L40:
	wk = 1.f;
	wkm = 1.f;
L50:
	if (kk == *n) {
	    goto L90;
	}
	j1 = k + 1;
	if (lower) {
	    j1 = 1;
	}
	j2 = *n;
	if (lower) {
	    j2 = k - 1;
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    sm += (r__1 = z__[j] + wkm * t[k + j * t_dim1], dabs(r__1));
	    z__[j] += wk * t[k + j * t_dim1];
	    s += (r__1 = z__[j], dabs(r__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	w = wkm - wk;
	wk = wkm;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    z__[j] += w * t[k + j * t_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     SOLVE T*Z = Y */

    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *n + 1 - kk;
	if (lower) {
	    k = kk;
	}
	if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = t[k + k * t_dim1], dabs(
		r__2))) {
	    goto L110;
	}
	s = (r__1 = t[k + k * t_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(
		r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L110:
	if (t[k + k * t_dim1] != 0.f) {
	    z__[k] /= t[k + k * t_dim1];
	}
	if (t[k + k * t_dim1] == 0.f) {
	    z__[k] = 1.f;
	}
	i1 = 1;
	if (lower) {
	    i1 = k + 1;
	}
	if (kk >= *n) {
	    goto L120;
	}
	w = -z__[k];
	i__2 = *n - kk;
	saxpy_(&i__2, &w, &t[i1 + k * t_dim1], &c__1, &z__[i1], &c__1);
L120:
/* L130: */
	;
    }
/*     MAKE ZNORM = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (tnorm != 0.f) {
	*rcond = ynorm / tnorm;
    }
    if (tnorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* strco_ */


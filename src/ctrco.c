/* ctrco.f -- translated by f2c (version 12.02.01).
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

/* DECK CTRCO */
/* Subroutine */ int ctrco_(complex *t, integer *ldt, integer *n, real *rcond,
	 complex *z__, integer *job)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer j, k, l;
    static real s;
    static complex w;
    static integer i1, j1, j2;
    static complex ek;
    static integer kk;
    static real sm;
    static complex wk, wkm;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static logical lower;
    static real tnorm, ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CTRCO */
/* ***PURPOSE  Estimate the condition number of a triangular matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C3 */
/* ***TYPE      COMPLEX (STRCO-S, DTRCO-D, CTRCO-C) */
/* ***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK, */
/*             TRIANGULAR MATRIX */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CTRCO estimates the condition of a complex triangular matrix. */

/*     On Entry */

/*        T       COMPLEX(LDT,N) */
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

/*        Z       COMPLEX(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  T  is close to a singular matrix, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CSSCAL, SCASUM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CTRCO */


/* ***FIRST EXECUTABLE STATEMENT  CTRCO */
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
	r__1 = tnorm, r__2 = scasum_(&l, &t[i1 + j * t_dim1], &c__1);
	tnorm = dmax(r__1,r__2);
/* L10: */
    }

/*     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  CTRANS(T)*Y = E . */
/*     CTRANS(T)  IS THE CONJUGATE TRANSPOSE OF T . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF Y . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE CTRANS(T)*Y = E */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
/* L20: */
    }
    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = kk;
	if (lower) {
	    k = *n + 1 - kk;
	}
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) != 0.f) {
	    i__3 = k;
	    q__2.r = -z__[i__3].r, q__2.i = -z__[i__3].i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4)
		    );
	    r__8 = (r__5 = q__1.r, dabs(r__5)) + (r__6 = r_imag(&q__1), dabs(
		    r__6));
	    q__4.r = q__1.r / r__8, q__4.i = q__1.i / r__8;
	    q__3.r = r__7 * q__4.r, q__3.i = r__7 * q__4.i;
	    ek.r = q__3.r, ek.i = q__3.i;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = k + k * t_dim1;
	if ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(r__2)) 
		<= (r__3 = t[i__3].r, dabs(r__3)) + (r__4 = r_imag(&t[k + k * 
		t_dim1]), dabs(r__4))) {
	    goto L30;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = k + k * t_dim1;
	s = ((r__1 = t[i__3].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * 
		t_dim1]), dabs(r__2))) / ((r__3 = q__1.r, dabs(r__3)) + (r__4 
		= r_imag(&q__1), dabs(r__4)));
	csscal_(n, &s, &z__[1], &c__1);
	q__2.r = s, q__2.i = 0.f;
	q__1.r = q__2.r * ek.r - q__2.i * ek.i, q__1.i = q__2.r * ek.i + 
		q__2.i * ek.r;
	ek.r = q__1.r, ek.i = q__1.i;
L30:
	i__2 = k;
	q__1.r = ek.r - z__[i__2].r, q__1.i = ek.i - z__[i__2].i;
	wk.r = q__1.r, wk.i = q__1.i;
	q__2.r = -ek.r, q__2.i = -ek.i;
	i__2 = k;
	q__1.r = q__2.r - z__[i__2].r, q__1.i = q__2.i - z__[i__2].i;
	wkm.r = q__1.r, wkm.i = q__1.i;
	s = (r__1 = wk.r, dabs(r__1)) + (r__2 = r_imag(&wk), dabs(r__2));
	sm = (r__1 = wkm.r, dabs(r__1)) + (r__2 = r_imag(&wkm), dabs(r__2));
	i__2 = k + k * t_dim1;
	if ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * t_dim1]
		), dabs(r__2)) == 0.f) {
	    goto L40;
	}
	r_cnjg(&q__2, &t[k + k * t_dim1]);
	c_div(&q__1, &wk, &q__2);
	wk.r = q__1.r, wk.i = q__1.i;
	r_cnjg(&q__2, &t[k + k * t_dim1]);
	c_div(&q__1, &wkm, &q__2);
	wkm.r = q__1.r, wkm.i = q__1.i;
	goto L50;
L40:
	wk.r = 1.f, wk.i = 0.f;
	wkm.r = 1.f, wkm.i = 0.f;
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
	    i__3 = j;
	    r_cnjg(&q__4, &t[k + j * t_dim1]);
	    q__3.r = wkm.r * q__4.r - wkm.i * q__4.i, q__3.i = wkm.r * q__4.i 
		    + wkm.i * q__4.r;
	    q__2.r = z__[i__3].r + q__3.r, q__2.i = z__[i__3].i + q__3.i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    sm += (r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(
		    r__2));
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &t[k + j * t_dim1]);
	    q__2.r = wk.r * q__3.r - wk.i * q__3.i, q__2.i = wk.r * q__3.i + 
		    wk.i * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = j;
	    s += (r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&z__[j]), 
		    dabs(r__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	q__1.r = wkm.r - wk.r, q__1.i = wkm.i - wk.i;
	w.r = q__1.r, w.i = q__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &t[k + j * t_dim1]);
	    q__2.r = w.r * q__3.r - w.i * q__3.i, q__2.i = w.r * q__3.i + w.i 
		    * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
/* L70: */
	}
L80:
L90:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L100: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     SOLVE T*Z = Y */

    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *n + 1 - kk;
	if (lower) {
	    k = kk;
	}
	i__2 = k;
	i__3 = k + k * t_dim1;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= (r__3 = t[i__3].r, dabs(r__3)) + (r__4 = r_imag(&t[
		k + k * t_dim1]), dabs(r__4))) {
	    goto L110;
	}
	i__2 = k + k * t_dim1;
	i__3 = k;
	s = ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * 
		t_dim1]), dabs(r__2))) / ((r__3 = z__[i__3].r, dabs(r__3)) + (
		r__4 = r_imag(&z__[k]), dabs(r__4)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L110:
	i__2 = k + k * t_dim1;
	if ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * t_dim1]
		), dabs(r__2)) != 0.f) {
	    i__3 = k;
	    c_div(&q__1, &z__[k], &t[k + k * t_dim1]);
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	}
	i__2 = k + k * t_dim1;
	if ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * t_dim1]
		), dabs(r__2)) == 0.f) {
	    i__3 = k;
	    z__[i__3].r = 1.f, z__[i__3].i = 0.f;
	}
	i1 = 1;
	if (lower) {
	    i1 = k + 1;
	}
	if (kk >= *n) {
	    goto L120;
	}
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	w.r = q__1.r, w.i = q__1.i;
	i__2 = *n - kk;
	caxpy_(&i__2, &w, &t[i1 + k * t_dim1], &c__1, &z__[i1], &c__1);
L120:
/* L130: */
	;
    }
/*     MAKE ZNORM = 1.0 */
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (tnorm != 0.f) {
	*rcond = ynorm / tnorm;
    }
    if (tnorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* ctrco_ */


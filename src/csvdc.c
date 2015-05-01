/* csvdc.f -- translated by f2c (version 12.02.01).
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
static complex c_b8 = {1.f,0.f};
static complex c_b53 = {-1.f,0.f};

/* DECK CSVDC */
/* Subroutine */ int csvdc_(complex *x, integer *ldx, integer *n, integer *p, 
	complex *s, complex *e, complex *u, integer *ldu, complex *v, integer 
	*ldv, complex *work, integer *job, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Local variables */
    static real b, c__, f, g;
    static integer i__, j, k, l, m;
    static complex r__, t;
    static real t1, el;
    static integer kk;
    static real cs;
    static integer ll, mm, ls;
    static real sl;
    static integer lu;
    static real sm, sn;
    static integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static real emm1, smm1;
    static integer kase, jobu, iter;
    static real test;
    static integer nctp1, nrtp1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    static real scale;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    static real shift;
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static integer maxit;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *), csrot_(integer *, complex *, 
	    integer *, complex *, integer *, real *, real *);
    static logical wantu, wantv;
    extern /* Subroutine */ int srotg_(real *, real *, real *, real *);
    static real ztest;
    extern doublereal scnrm2_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CSVDC */
/* ***PURPOSE  Perform the singular value decomposition of a rectangular */
/*            matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D6 */
/* ***TYPE      COMPLEX (SSVDC-S, DSVDC-D, CSVDC-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             SINGULAR VALUE DECOMPOSITION */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     CSVDC is a subroutine to reduce a complex NxP matrix X by */
/*     unitary transformations U and V to diagonal form.  The */
/*     diagonal elements S(I) are the singular values of X.  The */
/*     columns of U are the corresponding left singular vectors, */
/*     and the columns of V the right singular vectors. */

/*     On Entry */

/*         X         COMPLEX(LDX,P), where LDX .GE. N. */
/*                   X contains the matrix whose singular value */
/*                   decomposition is to be computed.  X is */
/*                   destroyed by CSVDC. */

/*         LDX       INTEGER. */
/*                   LDX is the leading dimension of the array X. */

/*         N         INTEGER. */
/*                   N is the number of rows of the matrix X. */

/*         P         INTEGER. */
/*                   P is the number of columns of the matrix X. */

/*         LDU       INTEGER. */
/*                   LDU is the leading dimension of the array U */
/*                   (see below). */

/*         LDV       INTEGER. */
/*                   LDV is the leading dimension of the array V */
/*                   (see below). */

/*         WORK      COMPLEX(N). */
/*                   WORK is a scratch array. */

/*         JOB       INTEGER. */
/*                   JOB controls the computation of the singular */
/*                   vectors.  It has the decimal expansion AB */
/*                   with the following meaning */

/*                        A .EQ. 0    Do not compute the left singular */
/*                                    vectors. */
/*                        A .EQ. 1    Return the N left singular vectors */
/*                                    in U. */
/*                        A .GE. 2    Return the first MIN(N,P) */
/*                                    left singular vectors in U. */
/*                        B .EQ. 0    Do not compute the right singular */
/*                                    vectors. */
/*                        B .EQ. 1    Return the right singular vectors */
/*                                    in V. */

/*     On Return */

/*         S         COMPLEX(MM), where MM = MIN(N+1,P). */
/*                   The first MIN(N,P) entries of S contain the */
/*                   singular values of X arranged in descending */
/*                   order of magnitude. */

/*         E         COMPLEX(P). */
/*                   E ordinarily contains zeros.  However see the */
/*                   discussion of INFO for exceptions. */

/*         U         COMPLEX(LDU,K), where LDU .GE. N.  If JOBA .EQ. 1 */
/*                                   then K .EQ. N.  If JOBA .GE. 2 then */
/*                                   K .EQ. MIN(N,P). */
/*                   U contains the matrix of right singular vectors. */
/*                   U is not referenced if JOBA .EQ. 0.  If N .LE. P */
/*                   or if JOBA .GT. 2, then U may be identified with X */
/*                   in the subroutine call. */

/*         V         COMPLEX(LDV,P), where LDV .GE. P. */
/*                   V contains the matrix of right singular vectors. */
/*                   V is not referenced if JOB .EQ. 0.  If P .LE. N, */
/*                   then V may be identified with X in the */
/*                   subroutine call. */

/*         INFO      INTEGER. */
/*                   The singular values (and their corresponding */
/*                   singular vectors) S(INFO+1),S(INFO+2),...,S(M) */
/*                   are correct (here M=MIN(N,P)).  Thus if */
/*                   INFO.EQ. 0, all the singular values and their */
/*                   vectors are correct.  In any event, the matrix */
/*                   B = CTRANS(U)*X*V is the bidiagonal matrix */
/*                   with the elements of S on its diagonal and the */
/*                   elements of E on its super-diagonal (CTRANS(U) */
/*                   is the conjugate-transpose of U).  Thus the */
/*                   singular values of X and B are the same. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTC, CSCAL, CSROT, CSWAP, SCNRM2, SROTG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790319  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSVDC */


/* ***FIRST EXECUTABLE STATEMENT  CSVDC */

/*     SET THE MAXIMUM NUMBER OF ITERATIONS. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --s;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;

    /* Function Body */
    maxit = 30;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = TRUE_;
    }
    if (*job % 10 != 0) {
	wantv = TRUE_;
    }

/*     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS */
/*     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND */
/*           PLACE THE L-TH DIAGONAL IN S(L). */

	i__2 = l;
	i__3 = *n - l + 1;
	r__1 = scnrm2_(&i__3, &x[l + l * x_dim1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = l;
	if ((r__1 = s[i__2].r, dabs(r__1)) + (r__2 = r_imag(&s[l]), dabs(r__2)
		) == 0.f) {
	    goto L10;
	}
	i__2 = l + l * x_dim1;
	if ((r__1 = x[i__2].r, dabs(r__1)) + (r__2 = r_imag(&x[l + l * x_dim1]
		), dabs(r__2)) != 0.f) {
	    i__3 = l;
	    r__3 = c_abs(&s[l]);
	    i__4 = l + l * x_dim1;
	    r__4 = c_abs(&x[l + l * x_dim1]);
	    q__2.r = x[i__4].r / r__4, q__2.i = x[i__4].i / r__4;
	    q__1.r = r__3 * q__2.r, q__1.i = r__3 * q__2.i;
	    s[i__3].r = q__1.r, s[i__3].i = q__1.i;
	}
	i__2 = *n - l + 1;
	c_div(&q__1, &c_b8, &s[l]);
	cscal_(&i__2, &q__1, &x[l + l * x_dim1], &c__1);
	i__2 = l + l * x_dim1;
	i__3 = l + l * x_dim1;
	q__1.r = x[i__3].r + 1.f, q__1.i = x[i__3].i + 0.f;
	x[i__2].r = q__1.r, x[i__2].i = q__1.i;
L10:
	i__2 = l;
	i__3 = l;
	q__1.r = -s[i__3].r, q__1.i = -s[i__3].i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    i__3 = l;
	    if ((r__1 = s[i__3].r, dabs(r__1)) + (r__2 = r_imag(&s[l]), dabs(
		    r__2)) == 0.f) {
		goto L30;
	    }

/*              APPLY THE TRANSFORMATION. */

	    i__3 = *n - l + 1;
	    cdotc_(&q__3, &i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1]
		    , &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &x[l + l * x_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           PLACE THE L-TH ROW OF X INTO  E FOR THE */
/*           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION. */

	    i__3 = j;
	    r_cnjg(&q__1, &x[l + j * x_dim1]);
	    e[i__3].r = q__1.r, e[i__3].i = q__1.i;
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK */
/*           MULTIPLICATION. */

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    i__4 = i__ + l * x_dim1;
	    u[i__3].r = x[i__4].r, u[i__3].i = x[i__4].i;
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE */
/*           L-TH SUPER-DIAGONAL IN E(L). */

	i__2 = l;
	i__3 = *p - l;
	r__1 = scnrm2_(&i__3, &e[lp1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	i__2 = l;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[l]), dabs(r__2)
		) == 0.f) {
	    goto L80;
	}
	i__2 = lp1;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[lp1]), dabs(
		r__2)) != 0.f) {
	    i__3 = l;
	    r__3 = c_abs(&e[l]);
	    i__4 = lp1;
	    r__4 = c_abs(&e[lp1]);
	    q__2.r = e[i__4].r / r__4, q__2.i = e[i__4].i / r__4;
	    q__1.r = r__3 * q__2.r, q__1.i = r__3 * q__2.i;
	    e[i__3].r = q__1.r, e[i__3].i = q__1.i;
	}
	i__2 = *p - l;
	c_div(&q__1, &c_b8, &e[l]);
	cscal_(&i__2, &q__1, &e[lp1], &c__1);
	i__2 = lp1;
	i__3 = lp1;
	q__1.r = e[i__3].r + 1.f, q__1.i = e[i__3].i + 0.f;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
L80:
	i__2 = l;
	r_cnjg(&q__2, &e[l]);
	q__1.r = -q__2.r, q__1.i = -q__2.i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	i__2 = l;
	if (lp1 > *n || (r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[l])
		, dabs(r__2)) == 0.f) {
	    goto L120;
	}

/*              APPLY THE TRANSFORMATION. */

	i__2 = *n;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    work[i__3].r = 0.f, work[i__3].i = 0.f;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    caxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    i__4 = j;
	    q__3.r = -e[i__4].r, q__3.i = -e[i__4].i;
	    c_div(&q__2, &q__3, &e[lp1]);
	    r_cnjg(&q__1, &q__2);
	    caxpy_(&i__3, &q__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT */
/*              BACK MULTIPLICATION. */

	i__2 = *p;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * v_dim1;
	    i__4 = i__;
	    v[i__3].r = e[i__4].r, v[i__3].i = e[i__4].i;
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	i__1 = nctp1;
	i__2 = nctp1 + nctp1 * x_dim1;
	s[i__1].r = x[i__2].r, s[i__1].i = x[i__2].i;
    }
    if (*n < m) {
	i__1 = m;
	s[i__1].r = 0.f, s[i__1].i = 0.f;
    }
    if (nrtp1 < m) {
	i__1 = nrtp1;
	i__2 = nrtp1 + m * x_dim1;
	e[i__1].r = x[i__2].r, e[i__1].i = x[i__2].i;
    }
    i__1 = m;
    e[i__1].r = 0.f, e[i__1].i = 0.f;

/*     IF REQUIRED, GENERATE U. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * u_dim1;
	    u[i__3].r = 0.f, u[i__3].i = 0.f;
/* L180: */
	}
	i__2 = j + j * u_dim1;
	u[i__2].r = 1.f, u[i__2].i = 0.f;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	i__2 = l;
	if ((r__1 = s[i__2].r, dabs(r__1)) + (r__2 = r_imag(&s[l]), dabs(r__2)
		) == 0.f) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    cdotc_(&q__3, &i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1]
		    , &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &u[l + l * u_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	cscal_(&i__2, &c_b53, &u[l + l * u_dim1], &c__1);
	i__2 = l + l * u_dim1;
	i__3 = l + l * u_dim1;
	q__1.r = u[i__3].r + 1.f, q__1.i = u[i__3].i + 0.f;
	u[i__2].r = q__1.r, u[i__2].i = q__1.i;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    u[i__3].r = 0.f, u[i__3].i = 0.f;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    u[i__3].r = 0.f, u[i__3].i = 0.f;
/* L260: */
	}
	i__2 = l + l * u_dim1;
	u[i__2].r = 1.f, u[i__2].i = 0.f;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     IF IT IS REQUIRED, GENERATE V. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	i__2 = l;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[l]), dabs(r__2)
		) == 0.f) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    cdotc_(&q__3, &i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &v[lp1 + l * v_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *p - l;
	    caxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * v_dim1;
	    v[i__3].r = 0.f, v[i__3].i = 0.f;
/* L330: */
	}
	i__2 = l + l * v_dim1;
	v[i__2].r = 1.f, v[i__2].i = 0.f;
/* L340: */
    }
L350:

/*     TRANSFORM S AND E SO THAT THEY ARE REAL. */

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((r__1 = s[i__2].r, dabs(r__1)) + (r__2 = r_imag(&s[i__]), dabs(
		r__2)) == 0.f) {
	    goto L360;
	}
	r__1 = c_abs(&s[i__]);
	q__1.r = r__1, q__1.i = 0.f;
	t.r = q__1.r, t.i = q__1.i;
	c_div(&q__1, &s[i__], &t);
	r__.r = q__1.r, r__.i = q__1.i;
	i__2 = i__;
	s[i__2].r = t.r, s[i__2].i = t.i;
	if (i__ < m) {
	    i__2 = i__;
	    c_div(&q__1, &e[i__], &r__);
	    e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	}
	if (wantu) {
	    cscal_(n, &r__, &u[i__ * u_dim1 + 1], &c__1);
	}
L360:
	if (i__ == m) {
	    goto L390;
	}
	i__2 = i__;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[i__]), dabs(
		r__2)) == 0.f) {
	    goto L370;
	}
	r__1 = c_abs(&e[i__]);
	q__1.r = r__1, q__1.i = 0.f;
	t.r = q__1.r, t.i = q__1.i;
	c_div(&q__1, &t, &e[i__]);
	r__.r = q__1.r, r__.i = q__1.i;
	i__2 = i__;
	e[i__2].r = t.r, e[i__2].i = t.i;
	i__2 = i__ + 1;
	i__3 = i__ + 1;
	q__1.r = s[i__3].r * r__.r - s[i__3].i * r__.i, q__1.i = s[i__3].r * 
		r__.i + s[i__3].i * r__.r;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	if (wantv) {
	    cscal_(p, &r__, &v[(i__ + 1) * v_dim1 + 1], &c__1);
	}
L370:
/* L380: */
	;
    }
L390:

/*     MAIN ITERATION LOOP FOR THE SINGULAR VALUES. */

    mm = m;
    iter = 0;
L400:

/*        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND. */

    if (m == 0) {
	goto L660;
    }

/*        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET */
/*        FLAG AND RETURN. */

    if (iter < maxit) {
	goto L410;
    }
    *info = m;
    goto L660;
L410:

/*        THIS SECTION OF THE PROGRAM INSPECTS FOR */
/*        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON */
/*        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS. */

/*           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M */
/*           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M */
/*           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND */
/*                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP). */
/*           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
	if (l == 0) {
	    goto L440;
	}
	test = c_abs(&s[l]) + c_abs(&s[l + 1]);
	ztest = test + c_abs(&e[l]);
	if (ztest != test) {
	    goto L420;
	}
	i__2 = l;
	e[i__2].r = 0.f, e[i__2].i = 0.f;
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (l != m - 1) {
	goto L450;
    }
    kase = 4;
    goto L520;
L450:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
	if (ls == l) {
	    goto L480;
	}
	test = 0.f;
	if (ls != m) {
	    test += c_abs(&e[ls]);
	}
	if (ls != l + 1) {
	    test += c_abs(&e[ls - 1]);
	}
	ztest = test + c_abs(&s[ls]);
	if (ztest != test) {
	    goto L460;
	}
	i__2 = ls;
	s[i__2].r = 0.f, s[i__2].i = 0.f;
	goto L480;
L460:
/* L470: */
	;
    }
L480:
    if (ls != l) {
	goto L490;
    }
    kase = 3;
    goto L510;
L490:
    if (ls != m) {
	goto L500;
    }
    kase = 1;
    goto L510;
L500:
    kase = 2;
    l = ls;
L510:
L520:
    ++l;

/*        PERFORM THE TASK INDICATED BY KASE. */

    switch (kase) {
	case 1:  goto L530;
	case 2:  goto L560;
	case 3:  goto L580;
	case 4:  goto L610;
    }

/*        DEFLATE NEGLIGIBLE S(M). */

L530:
    mm1 = m - 1;
    i__1 = m - 1;
    f = e[i__1].r;
    i__1 = m - 1;
    e[i__1].r = 0.f, e[i__1].i = 0.f;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	i__2 = k;
	t1 = s[i__2].r;
	srotg_(&t1, &f, &cs, &sn);
	i__2 = k;
	q__1.r = t1, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	if (k == l) {
	    goto L540;
	}
	i__2 = k - 1;
	f = -sn * e[i__2].r;
	i__2 = k - 1;
	i__3 = k - 1;
	q__1.r = cs * e[i__3].r, q__1.i = cs * e[i__3].i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
L540:
	if (wantv) {
	    csrot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L550: */
    }
    goto L650;

/*        SPLIT AT NEGLIGIBLE S(L). */

L560:
    i__1 = l - 1;
    f = e[i__1].r;
    i__1 = l - 1;
    e[i__1].r = 0.f, e[i__1].i = 0.f;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	i__2 = k;
	t1 = s[i__2].r;
	srotg_(&t1, &f, &cs, &sn);
	i__2 = k;
	q__1.r = t1, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = k;
	f = -sn * e[i__2].r;
	i__2 = k;
	i__3 = k;
	q__1.r = cs * e[i__3].r, q__1.i = cs * e[i__3].i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	if (wantu) {
	    csrot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L570: */
    }
    goto L650;

/*        PERFORM ONE QR STEP. */

L580:

/*           CALCULATE THE SHIFT. */

/* Computing MAX */
    r__1 = c_abs(&s[m]), r__2 = c_abs(&s[m - 1]), r__1 = max(r__1,r__2), r__2 
	    = c_abs(&e[m - 1]), r__1 = max(r__1,r__2), r__2 = c_abs(&s[l]), 
	    r__1 = max(r__1,r__2), r__2 = c_abs(&e[l]);
    scale = dmax(r__1,r__2);
    i__1 = m;
    sm = s[i__1].r / scale;
    i__1 = m - 1;
    smm1 = s[i__1].r / scale;
    i__1 = m - 1;
    emm1 = e[i__1].r / scale;
    i__1 = l;
    sl = s[i__1].r / scale;
    i__1 = l;
    el = e[i__1].r / scale;
/* Computing 2nd power */
    r__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + r__1 * r__1) / 2.f;
/* Computing 2nd power */
    r__1 = sm * emm1;
    c__ = r__1 * r__1;
    shift = 0.f;
    if (b == 0.f && c__ == 0.f) {
	goto L590;
    }
/* Computing 2nd power */
    r__1 = b;
    shift = sqrt(r__1 * r__1 + c__);
    if (b < 0.f) {
	shift = -shift;
    }
    shift = c__ / (b + shift);
L590:
    f = (sl + sm) * (sl - sm) - shift;
    g = sl * el;

/*           CHASE ZEROS. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	srotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    i__2 = k - 1;
	    q__1.r = f, q__1.i = 0.f;
	    e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	}
	i__2 = k;
	i__3 = k;
	f = cs * s[i__2].r + sn * e[i__3].r;
	i__2 = k;
	i__3 = k;
	q__2.r = cs * e[i__3].r, q__2.i = cs * e[i__3].i;
	i__4 = k;
	q__3.r = sn * s[i__4].r, q__3.i = sn * s[i__4].i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	i__2 = k + 1;
	g = sn * s[i__2].r;
	i__2 = k + 1;
	i__3 = k + 1;
	q__1.r = cs * s[i__3].r, q__1.i = cs * s[i__3].i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	if (wantv) {
	    csrot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	srotg_(&f, &g, &cs, &sn);
	i__2 = k;
	q__1.r = f, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = k;
	i__3 = k + 1;
	f = cs * e[i__2].r + sn * s[i__3].r;
	i__2 = k + 1;
	r__1 = -sn;
	i__3 = k;
	q__2.r = r__1 * e[i__3].r, q__2.i = r__1 * e[i__3].i;
	i__4 = k + 1;
	q__3.r = cs * s[i__4].r, q__3.i = cs * s[i__4].i;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = k + 1;
	g = sn * e[i__2].r;
	i__2 = k + 1;
	i__3 = k + 1;
	q__1.r = cs * e[i__3].r, q__1.i = cs * e[i__3].i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	if (wantu && k < *n) {
	    csrot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L600: */
    }
    i__1 = m - 1;
    q__1.r = f, q__1.i = 0.f;
    e[i__1].r = q__1.r, e[i__1].i = q__1.i;
    ++iter;
    goto L650;

/*        CONVERGENCE. */

L610:

/*           MAKE THE SINGULAR VALUE  POSITIVE */

    i__1 = l;
    if (s[i__1].r >= 0.f) {
	goto L620;
    }
    i__1 = l;
    i__2 = l;
    q__1.r = -s[i__2].r, q__1.i = -s[i__2].i;
    s[i__1].r = q__1.r, s[i__1].i = q__1.i;
    if (wantv) {
	cscal_(p, &c_b53, &v[l * v_dim1 + 1], &c__1);
    }
L620:

/*           ORDER THE SINGULAR VALUE. */

L630:
    if (l == mm) {
	goto L640;
    }
    i__1 = l;
    i__2 = l + 1;
    if (s[i__1].r >= s[i__2].r) {
	goto L640;
    }
    i__1 = l;
    t.r = s[i__1].r, t.i = s[i__1].i;
    i__1 = l;
    i__2 = l + 1;
    s[i__1].r = s[i__2].r, s[i__1].i = s[i__2].i;
    i__1 = l + 1;
    s[i__1].r = t.r, s[i__1].i = t.i;
    if (wantv && l < *p) {
	cswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	cswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L630;
L640:
    iter = 0;
    --m;
L650:
    goto L400;
L660:
    return 0;
} /* csvdc_ */


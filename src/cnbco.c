/* cnbco.f -- translated by f2c (version 12.02.01).
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

/* DECK CNBCO */
/* Subroutine */ int cnbco_(complex *abe, integer *lda, integer *n, integer *
	ml, integer *mu, integer *ipvt, real *rcond, complex *z__)
{
    /* System generated locals */
    integer abe_dim1, abe_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real s;
    static complex t;
    static integer kb;
    static complex ek;
    static integer lm, mm, nl, ju;
    static real sm;
    static complex wk;
    static integer nu, lz, ml1, kp1, ldb;
    static complex wkm;
    static integer info;
    extern /* Subroutine */ int cnbfa_(complex *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    static real anorm;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CNBCO */
/* ***PURPOSE  Factor a band matrix using Gaussian elimination and */
/*            estimate the condition number. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2C2 */
/* ***TYPE      COMPLEX (SNBCO-S, DNBCO-D, CNBCO-C) */
/* ***KEYWORDS  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION, */
/*             NONSYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*     CNBCO factors a complex band matrix by Gaussian */
/*     elimination and estimates the condition of the matrix. */

/*     If RCOND is not needed, CNBFA is slightly faster. */
/*     To solve  A*X = B , follow CNBCO by CNBSL. */
/*     To compute  INVERSE(A)*C , follow CNBCO by CNBSL. */
/*     To compute  DETERMINANT(A) , follow CNBCO by CNBDI. */

/*     On Entry */

/*        ABE     COMPLEX(LDA, NC) */
/*                contains the matrix in band storage.  The rows */
/*                of the original matrix are stored in the rows */
/*                of ABE and the diagonals of the original matrix */
/*                are stored in columns 1 through ML+MU+1 of ABE. */
/*                NC must be .GE. 2*ML+MU+1 . */
/*                See the comments below for details. */

/*        LDA     INTEGER */
/*                the leading dimension of the array ABE. */
/*                LDA must be .GE. N . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */
/*                0 .LE. MU .LT. N . */
/*                More efficient if ML .LE. MU . */

/*     On Return */

/*        ABE     an upper triangular matrix in band storage */
/*                and the multipliers which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        RCOND   REAL */
/*                an estimate of the reciprocal condition of  A . */
/*                For the system  A*X = B , relative perturbations */
/*                in  A  and  B  of size  EPSILON  may cause */
/*                relative perturbations in  X  of size  EPSILON/RCOND . */
/*                If  RCOND  is so small that the logical expression */
/*                         1.0 + RCOND .EQ. 1.0 */
/*                is true, then  A  may be singular to working */
/*                precision.  In particular,  RCOND  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        Z       COMPLEX(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  A  is close to a singular matrix, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     Band Storage */

/*           If  A  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ML = (band width below the diagonal) */
/*                   MU = (band width above the diagonal) */
/*                   DO 20 I = 1, N */
/*                      J1 = MAX(1, I-ML) */
/*                      J2 = MIN(N, I+MU) */
/*                      DO 10 J = J1, J2 */
/*                         K = J - I + ML + 1 */
/*                         ABE(I,K) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           This uses columns  1  through  ML+MU+1  of ABE . */
/*           Furthermore,  ML  additional columns are needed in */
/*           ABE  starting with column  ML+MU+2  for elements */
/*           generated during the triangularization.  The total */
/*           number of columns needed in  ABE  is  2*ML+MU+1 . */

/*     Example:  If the original matrix is */

/*           11 12 13  0  0  0 */
/*           21 22 23 24  0  0 */
/*            0 32 33 34 35  0 */
/*            0  0 43 44 45 46 */
/*            0  0  0 54 55 56 */
/*            0  0  0  0 65 66 */

/*      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain */

/*            * 11 12 13  +     , * = not used */
/*           21 22 23 24  +     , + = used for pivoting */
/*           32 33 34 35  + */
/*           43 44 45 46  + */
/*           54 55 56  *  + */
/*           65 66  *  *  + */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTC, CNBFA, CSSCAL, SCASUM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800730  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CNBCO */


/*     COMPUTE 1-NORM OF A */

/* ***FIRST EXECUTABLE STATEMENT  CNBCO */
    /* Parameter adjustments */
    abe_dim1 = *lda;
    abe_offset = 1 + abe_dim1;
    abe -= abe_offset;
    --ipvt;
    --z__;

    /* Function Body */
    ml1 = *ml + 1;
    ldb = *lda - 1;
    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = *mu, i__3 = j - 1;
	nu = min(i__2,i__3);
/* Computing MIN */
	i__2 = *ml, i__3 = *n - j;
	nl = min(i__2,i__3);
	l = nu + 1 + nl;
/* Computing MAX */
	r__1 = anorm, r__2 = scasum_(&l, &abe[j + nl + (ml1 - nl) * abe_dim1],
		 &ldb);
	anorm = dmax(r__1,r__2);
/* L10: */
    }

/*     FACTOR */

    cnbfa_(&abe[abe_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND CTRANS(A)*Y = E . */
/*     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF  W  WHERE CTRANS(U)*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE CTRANS(U)*W = E */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
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
	i__3 = k + ml1 * abe_dim1;
	if ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(r__2)) 
		<= (r__3 = abe[i__3].r, dabs(r__3)) + (r__4 = r_imag(&abe[k + 
		ml1 * abe_dim1]), dabs(r__4))) {
	    goto L30;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = k + ml1 * abe_dim1;
	s = ((r__1 = abe[i__3].r, dabs(r__1)) + (r__2 = r_imag(&abe[k + ml1 * 
		abe_dim1]), dabs(r__2))) / ((r__3 = q__1.r, dabs(r__3)) + (
		r__4 = r_imag(&q__1), dabs(r__4)));
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
	i__2 = k + ml1 * abe_dim1;
	if ((r__1 = abe[i__2].r, dabs(r__1)) + (r__2 = r_imag(&abe[k + ml1 * 
		abe_dim1]), dabs(r__2)) == 0.f) {
	    goto L40;
	}
	r_cnjg(&q__2, &abe[k + ml1 * abe_dim1]);
	c_div(&q__1, &wk, &q__2);
	wk.r = q__1.r, wk.i = q__1.i;
	r_cnjg(&q__2, &abe[k + ml1 * abe_dim1]);
	c_div(&q__1, &wkm, &q__2);
	wkm.r = q__1.r, wkm.i = q__1.i;
	goto L50;
L40:
	wk.r = 1.f, wk.i = 0.f;
	wkm.r = 1.f, wkm.i = 0.f;
L50:
	kp1 = k + 1;
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = ml1;
	if (kp1 > ju) {
	    goto L90;
	}
	i__2 = ju;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    ++mm;
	    i__3 = i__;
	    r_cnjg(&q__4, &abe[k + mm * abe_dim1]);
	    q__3.r = wkm.r * q__4.r - wkm.i * q__4.i, q__3.i = wkm.r * q__4.i 
		    + wkm.i * q__4.r;
	    q__2.r = z__[i__3].r + q__3.r, q__2.i = z__[i__3].i + q__3.i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    sm += (r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(
		    r__2));
	    i__3 = i__;
	    i__4 = i__;
	    r_cnjg(&q__3, &abe[k + mm * abe_dim1]);
	    q__2.r = wk.r * q__3.r - wk.i * q__3.i, q__2.i = wk.r * q__3.i + 
		    wk.i * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = i__;
	    s += (r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&z__[i__]),
		     dabs(r__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	q__1.r = wkm.r - wk.r, q__1.i = wkm.i - wk.i;
	t.r = q__1.r, t.i = q__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	mm = ml1;
	i__2 = ju;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    ++mm;
	    i__3 = i__;
	    i__4 = i__;
	    r_cnjg(&q__3, &abe[k + mm * abe_dim1]);
	    q__2.r = t.r * q__3.r - t.i * q__3.i, q__2.i = t.r * q__3.i + t.i 
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

/*     SOLVE CTRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	nl = min(i__2,i__3);
	if (k < *n) {
	    i__2 = k;
	    i__3 = k;
	    i__4 = -ldb;
	    cdotc_(&q__2, &nl, &abe[k + nl + (ml1 - nl) * abe_dim1], &i__4, &
		    z__[k + 1], &c__1);
	    q__1.r = z__[i__3].r + q__2.r, q__1.i = z__[i__3].i + q__2.i;
	    z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	}
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= 1.f) {
	    goto L110;
	}
	i__2 = k;
	s = 1.f / ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]),
		 dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* L120: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	nl = min(i__2,i__3);
	if (k < *n) {
	    i__2 = -ldb;
	    caxpy_(&nl, &t, &abe[k + nl + (ml1 - nl) * abe_dim1], &i__2, &z__[
		    k + 1], &c__1);
	}
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= 1.f) {
	    goto L130;
	}
	i__2 = k;
	s = 1.f / ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]),
		 dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k + ml1 * abe_dim1;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= (r__3 = abe[i__3].r, dabs(r__3)) + (r__4 = r_imag(&
		abe[k + ml1 * abe_dim1]), dabs(r__4))) {
	    goto L150;
	}
	i__2 = k + ml1 * abe_dim1;
	i__3 = k;
	s = ((r__1 = abe[i__2].r, dabs(r__1)) + (r__2 = r_imag(&abe[k + ml1 * 
		abe_dim1]), dabs(r__2))) / ((r__3 = z__[i__3].r, dabs(r__3)) 
		+ (r__4 = r_imag(&z__[k]), dabs(r__4)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	i__2 = k + ml1 * abe_dim1;
	if ((r__1 = abe[i__2].r, dabs(r__1)) + (r__2 = r_imag(&abe[k + ml1 * 
		abe_dim1]), dabs(r__2)) != 0.f) {
	    i__3 = k;
	    c_div(&q__1, &z__[k], &abe[k + ml1 * abe_dim1]);
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	}
	i__2 = k + ml1 * abe_dim1;
	if ((r__1 = abe[i__2].r, dabs(r__1)) + (r__2 = r_imag(&abe[k + ml1 * 
		abe_dim1]), dabs(r__2)) == 0.f) {
	    i__3 = k;
	    z__[i__3].r = 1.f, z__[i__3].i = 0.f;
	}
	lm = min(k,m) - 1;
	lz = k - lm;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = -ldb;
	caxpy_(&lm, &t, &abe[k - 1 + (*ml + 2) * abe_dim1], &i__2, &z__[lz], &
		c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0E0 */
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* cnbco_ */


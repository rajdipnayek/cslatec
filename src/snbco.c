/* snbco.f -- translated by f2c (version 12.02.01).
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

/* DECK SNBCO */
/* Subroutine */ int snbco_(real *abe, integer *lda, integer *n, integer *ml, 
	integer *mu, integer *ipvt, real *rcond, real *z__)
{
    /* System generated locals */
    integer abe_dim1, abe_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real s, t;
    static integer kb;
    static real ek;
    static integer lm, mm, nl, ju;
    static real sm, wk;
    static integer nu, lz, ml1, kp1, ldb;
    static real wkm;
    static integer info;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int snbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sscal_(integer *, real *, 
	    real *, integer *);
    static real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    static real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SNBCO */
/* ***PURPOSE  Factor a band matrix using Gaussian elimination and */
/*            estimate the condition number. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2A2 */
/* ***TYPE      SINGLE PRECISION (SNBCO-S, DNBCO-D, CNBCO-C) */
/* ***KEYWORDS  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION, */
/*             NONSYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*     SNBCO factors a real band matrix by Gaussian */
/*     elimination and estimates the condition of the matrix. */

/*     If RCOND is not needed, SNBFA is slightly faster. */
/*     To solve  A*X = B , follow SNBCO by SNBSL. */
/*     To compute  INVERSE(A)*C , follow SNBCO by SNBSL. */
/*     To compute  DETERMINANT(A) , follow SNBCO by SNBDI. */

/*     On Entry */

/*        ABE     REAL(LDA, NC) */
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
/*                The factorization can be written  A = L*U , where */
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

/*        Z       REAL(N) */
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
/* ***ROUTINES CALLED  SASUM, SAXPY, SDOT, SNBFA, SSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800723  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SNBCO */

/* ***FIRST EXECUTABLE STATEMENT  SNBCO */
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
	r__1 = anorm, r__2 = sasum_(&l, &abe[j + nl + (ml1 - nl) * abe_dim1], 
		&ldb);
	anorm = dmax(r__1,r__2);
/* L10: */
    }

/*     FACTOR */

    snbfa_(&abe[abe_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND TRANS(A)*Y = E . */
/*     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE */
/*     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF  W WHERE */
/*     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID */
/*     OVERFLOW. */

/*     SOLVE TRANS(U)*W = E */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= (r__2 = abe[k + ml1 * 
		abe_dim1], dabs(r__2))) {
	    goto L30;
	}
	s = (r__1 = abe[k + ml1 * abe_dim1], dabs(r__1)) / (r__2 = ek - z__[k]
		, dabs(r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	if (abe[k + ml1 * abe_dim1] == 0.f) {
	    goto L40;
	}
	wk /= abe[k + ml1 * abe_dim1];
	wkm /= abe[k + ml1 * abe_dim1];
	goto L50;
L40:
	wk = 1.f;
	wkm = 1.f;
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
	    sm += (r__1 = z__[i__] + wkm * abe[k + mm * abe_dim1], dabs(r__1))
		    ;
	    z__[i__] += wk * abe[k + mm * abe_dim1];
	    s += (r__1 = z__[i__], dabs(r__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	mm = ml1;
	i__2 = ju;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    ++mm;
	    z__[i__] += t * abe[k + mm * abe_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	nl = min(i__2,i__3);
	if (k < *n) {
	    i__2 = -ldb;
	    z__[k] += sdot_(&nl, &abe[k + nl + (ml1 - nl) * abe_dim1], &i__2, 
		    &z__[k + 1], &c__1);
	}
	if ((r__1 = z__[k], dabs(r__1)) <= 1.f) {
	    goto L110;
	}
	s = 1.f / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	nl = min(i__2,i__3);
	if (k < *n) {
	    i__2 = -ldb;
	    saxpy_(&nl, &t, &abe[k + nl + (ml1 - nl) * abe_dim1], &i__2, &z__[
		    k + 1], &c__1);
	}
	if ((r__1 = z__[k], dabs(r__1)) <= 1.f) {
	    goto L130;
	}
	s = 1.f / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = abe[k + ml1 * abe_dim1], 
		dabs(r__2))) {
	    goto L150;
	}
	s = (r__1 = abe[k + ml1 * abe_dim1], dabs(r__1)) / (r__2 = z__[k], 
		dabs(r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (abe[k + ml1 * abe_dim1] != 0.f) {
	    z__[k] /= abe[k + ml1 * abe_dim1];
	}
	if (abe[k + ml1 * abe_dim1] == 0.f) {
	    z__[k] = 1.f;
	}
	lm = min(k,m) - 1;
	lz = k - lm;
	t = -z__[k];
	i__2 = -ldb;
	saxpy_(&lm, &t, &abe[k - 1 + (*ml + 2) * abe_dim1], &i__2, &z__[lz], &
		c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0E0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* snbco_ */


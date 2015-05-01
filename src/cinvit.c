/* cinvit.f -- translated by f2c (version 12.02.01).
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

/* DECK CINVIT */
/* Subroutine */ int cinvit_(integer *nm, integer *n, real *ar, real *ai, 
	real *wr, real *wi, logical *select, integer *mm, integer *m, real *
	zr, real *zi, integer *ierr, real *rm1, real *rm2, real *rv1, real *
	rv2)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, rm1_dim1, rm1_offset, rm2_dim1, rm2_offset, 
	    i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, s;
    static real x, y;
    static integer ii, mp, uk, km1, ip1, its;
    static real eps3;
    extern /* Subroutine */ int cdiv_(real *, real *, real *, real *, real *, 
	    real *);
    static real norm, normv, ilambd, rlambd;
    extern doublereal pythag_(real *, real *);
    static real growto, ukroot;

/* ***BEGIN PROLOGUE  CINVIT */
/* ***PURPOSE  Compute the eigenvectors of a complex upper Hessenberg */
/*            associated with specified eigenvalues using inverse */
/*            iteration. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2B */
/* ***TYPE      COMPLEX (INVIT-S, CINVIT-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure CXINVIT */
/*     by Peters and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971). */

/*     This subroutine finds those eigenvectors of A COMPLEX UPPER */
/*     Hessenberg matrix corresponding to specified eigenvalues, */
/*     using inverse iteration. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR, AI, ZR and ZI, as declared in the */
/*          calling program dimension statement.  NM is an INTEGER */
/*          variable. */

/*        N is the order of the matrix A=(AR,AI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the complex upper Hessenberg matrix.  AR and AI are */
/*          two-dimensional REAL arrays, dimensioned AR(NM,N) */
/*          and AI(NM,N). */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues of the matrix.  The eigenvalues must be */
/*          stored in a manner identical to that of subroutine  COMLR, */
/*          which recognizes possible splitting of the matrix.  WR and */
/*          WI are one-dimensional REAL arrays, dimensioned WR(N) and */
/*          WI(N). */

/*        SELECT specifies the eigenvectors to be found.  The */
/*          eigenvector corresponding to the J-th eigenvalue is */
/*          specified by setting SELECT(J) to .TRUE.  SELECT is a */
/*          one-dimensional LOGICAL array, dimensioned SELECT(N). */

/*        MM should be set to an upper bound for the number of */
/*          eigenvectors to be found.  MM is an INTEGER variable. */

/*     On OUTPUT */

/*        AR, AI, WI, and SELECT are unaltered. */

/*        WR may have been altered since close eigenvalues are perturbed */
/*          slightly in searching for independent eigenvectors. */

/*        M is the number of eigenvectors actually found.  M is an */
/*          INTEGER variable. */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the eigenvectors corresponding to the flagged eigenvalues. */
/*          The eigenvectors are normalized so that the component of */
/*          largest magnitude is 1.  Any vector which fails the */
/*          acceptance test is set to zero.  ZR and ZI are */
/*          two-dimensional REAL arrays, dimensioned ZR(NM,MM) and */
/*          ZI(NM,MM). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          -(2*N+1)   if more than MM eigenvectors have been requested */
/*                     (the MM eigenvectors calculated to this point are */
/*                     in ZR and ZI), */
/*          -K         if the iteration corresponding to the K-th */
/*                     value fails (if this occurs more than once, K */
/*                     is the index of the last occurrence); the */
/*                     corresponding columns of ZR and ZI are set to */
/*                     zero vectors, */
/*          -(N+K)     if both error situations occur. */

/*        RV1 and RV2 are one-dimensional REAL arrays used for */
/*          temporary storage, dimensioned RV1(N) and RV2(N). */
/*          They hold the approximate eigenvectors during the inverse */
/*          iteration process. */

/*        RM1 and RM2 are two-dimensional REAL arrays used for */
/*          temporary storage, dimensioned RM1(N,N) and RM2(N,N). */
/*          These arrays hold the triangularized form of the upper */
/*          Hessenberg matrix used in the inverse iteration process. */

/*     The ALGOL procedure GUESSVEC appears in CINVIT in-line. */

/*     Calls PYTHAG(A,B) for sqrt(A**2 + B**2). */
/*     Calls CDIV for complex division. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  CDIV, PYTHAG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CINVIT */


/* ***FIRST EXECUTABLE STATEMENT  CINVIT */
    /* Parameter adjustments */
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1;
    zr -= zr_offset;
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
    rm2_dim1 = *n;
    rm2_offset = 1 + rm2_dim1;
    rm2 -= rm2_offset;
    rm1_dim1 = *n;
    rm1_offset = 1 + rm1_dim1;
    rm1 -= rm1_offset;
    --select;
    --wr;
    --wi;
    --rv1;
    --rv2;

    /* Function Body */
    *ierr = 0;
    uk = 0;
    s = 1;

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (! select[k]) {
	    goto L980;
	}
	if (s > *mm) {
	    goto L1000;
	}
	if (uk >= k) {
	    goto L200;
	}
/*     .......... CHECK FOR POSSIBLE SPLITTING .......... */
	i__2 = *n;
	for (uk = k; uk <= i__2; ++uk) {
	    if (uk == *n) {
		goto L140;
	    }
	    if (ar[uk + 1 + uk * ar_dim1] == 0.f && ai[uk + 1 + uk * ai_dim1] 
		    == 0.f) {
		goto L140;
	    }
/* L120: */
	}
/*     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK */
/*                (HESSENBERG) MATRIX .......... */
L140:
	norm = 0.f;
	mp = 1;

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = 0.f;

	    i__3 = uk;
	    for (j = mp; j <= i__3; ++j) {
/* L160: */
		x += pythag_(&ar[i__ + j * ar_dim1], &ai[i__ + j * ai_dim1]);
	    }

	    if (x > norm) {
		norm = x;
	    }
	    mp = i__;
/* L180: */
	}
/*     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION */
/*                AND CLOSE ROOTS ARE MODIFIED BY EPS3 .......... */
	if (norm == 0.f) {
	    norm = 1.f;
	}
	eps3 = norm;
L190:
	eps3 *= .5f;
	if (norm + eps3 > norm) {
	    goto L190;
	}
	eps3 *= 2.f;
/*     .......... GROWTO IS THE CRITERION FOR GROWTH .......... */
	ukroot = sqrt((real) uk);
	growto = .1f / ukroot;
L200:
	rlambd = wr[k];
	ilambd = wi[k];
	if (k == 1) {
	    goto L280;
	}
	km1 = k - 1;
	goto L240;
/*     .......... PERTURB EIGENVALUE IF IT IS CLOSE */
/*                TO ANY PREVIOUS EIGENVALUE .......... */
L220:
	rlambd += eps3;
/*     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- .......... */
L240:
	i__2 = km1;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = k - ii;
	    if (select[i__] && (r__1 = wr[i__] - rlambd, dabs(r__1)) < eps3 &&
		     (r__2 = wi[i__] - ilambd, dabs(r__2)) < eps3) {
		goto L220;
	    }
/* L260: */
	}

	wr[k] = rlambd;
/*     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I */
/*                AND INITIAL COMPLEX VECTOR .......... */
L280:
	mp = 1;

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {

	    i__3 = uk;
	    for (j = mp; j <= i__3; ++j) {
		rm1[i__ + j * rm1_dim1] = ar[i__ + j * ar_dim1];
		rm2[i__ + j * rm2_dim1] = ai[i__ + j * ai_dim1];
/* L300: */
	    }

	    rm1[i__ + i__ * rm1_dim1] -= rlambd;
	    rm2[i__ + i__ * rm2_dim1] -= ilambd;
	    mp = i__;
	    rv1[i__] = eps3;
/* L320: */
	}
/*     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES, */
/*                REPLACING ZERO PIVOTS BY EPS3 .......... */
	if (uk == 1) {
	    goto L420;
	}

	i__2 = uk;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    mp = i__ - 1;
	    if (pythag_(&rm1[i__ + mp * rm1_dim1], &rm2[i__ + mp * rm2_dim1]) 
		    <= pythag_(&rm1[mp + mp * rm1_dim1], &rm2[mp + mp * 
		    rm2_dim1])) {
		goto L360;
	    }

	    i__3 = uk;
	    for (j = mp; j <= i__3; ++j) {
		y = rm1[i__ + j * rm1_dim1];
		rm1[i__ + j * rm1_dim1] = rm1[mp + j * rm1_dim1];
		rm1[mp + j * rm1_dim1] = y;
		y = rm2[i__ + j * rm2_dim1];
		rm2[i__ + j * rm2_dim1] = rm2[mp + j * rm2_dim1];
		rm2[mp + j * rm2_dim1] = y;
/* L340: */
	    }

L360:
	    if (rm1[mp + mp * rm1_dim1] == 0.f && rm2[mp + mp * rm2_dim1] == 
		    0.f) {
		rm1[mp + mp * rm1_dim1] = eps3;
	    }
	    cdiv_(&rm1[i__ + mp * rm1_dim1], &rm2[i__ + mp * rm2_dim1], &rm1[
		    mp + mp * rm1_dim1], &rm2[mp + mp * rm2_dim1], &x, &y);
	    if (x == 0.f && y == 0.f) {
		goto L400;
	    }

	    i__3 = uk;
	    for (j = i__; j <= i__3; ++j) {
		rm1[i__ + j * rm1_dim1] = rm1[i__ + j * rm1_dim1] - x * rm1[
			mp + j * rm1_dim1] + y * rm2[mp + j * rm2_dim1];
		rm2[i__ + j * rm2_dim1] = rm2[i__ + j * rm2_dim1] - x * rm2[
			mp + j * rm2_dim1] - y * rm1[mp + j * rm1_dim1];
/* L380: */
	    }

L400:
	    ;
	}

L420:
	if (rm1[uk + uk * rm1_dim1] == 0.f && rm2[uk + uk * rm2_dim1] == 0.f) 
		{
	    rm1[uk + uk * rm1_dim1] = eps3;
	}
	its = 0;
/*     .......... BACK SUBSTITUTION */
/*                FOR I=UK STEP -1 UNTIL 1 DO -- .......... */
L660:
	i__2 = uk;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = uk + 1 - ii;
	    x = rv1[i__];
	    y = 0.f;
	    if (i__ == uk) {
		goto L700;
	    }
	    ip1 = i__ + 1;

	    i__3 = uk;
	    for (j = ip1; j <= i__3; ++j) {
		x = x - rm1[i__ + j * rm1_dim1] * rv1[j] + rm2[i__ + j * 
			rm2_dim1] * rv2[j];
		y = y - rm1[i__ + j * rm1_dim1] * rv2[j] - rm2[i__ + j * 
			rm2_dim1] * rv1[j];
/* L680: */
	    }

L700:
	    cdiv_(&x, &y, &rm1[i__ + i__ * rm1_dim1], &rm2[i__ + i__ * 
		    rm2_dim1], &rv1[i__], &rv2[i__]);
/* L720: */
	}
/*     .......... ACCEPTANCE TEST FOR EIGENVECTOR */
/*                AND NORMALIZATION .......... */
	++its;
	norm = 0.f;
	normv = 0.f;

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = pythag_(&rv1[i__], &rv2[i__]);
	    if (normv >= x) {
		goto L760;
	    }
	    normv = x;
	    j = i__;
L760:
	    norm += x;
/* L780: */
	}

	if (norm < growto) {
	    goto L840;
	}
/*     .......... ACCEPT VECTOR .......... */
	x = rv1[j];
	y = rv2[j];

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    cdiv_(&rv1[i__], &rv2[i__], &x, &y, &zr[i__ + s * zr_dim1], &zi[
		    i__ + s * zi_dim1]);
/* L820: */
	}

	if (uk == *n) {
	    goto L940;
	}
	j = uk + 1;
	goto L900;
/*     .......... IN-LINE PROCEDURE FOR CHOOSING */
/*                A NEW STARTING VECTOR .......... */
L840:
	if (its >= uk) {
	    goto L880;
	}
	x = ukroot;
	y = eps3 / (x + 1.f);
	rv1[1] = eps3;

	i__2 = uk;
	for (i__ = 2; i__ <= i__2; ++i__) {
/* L860: */
	    rv1[i__] = y;
	}

	j = uk - its + 1;
	rv1[j] -= eps3 * x;
	goto L660;
/*     .......... SET ERROR -- UNACCEPTED EIGENVECTOR .......... */
L880:
	j = 1;
	*ierr = -k;
/*     .......... SET REMAINING VECTOR COMPONENTS TO ZERO .......... */
L900:
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    zr[i__ + s * zr_dim1] = 0.f;
	    zi[i__ + s * zi_dim1] = 0.f;
/* L920: */
	}

L940:
	++s;
L980:
	;
    }

    goto L1001;
/*     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR */
/*                SPACE REQUIRED .......... */
L1000:
    if (*ierr != 0) {
	*ierr -= *n;
    }
    if (*ierr == 0) {
	*ierr = -((*n << 1) + 1);
    }
L1001:
    *m = s - 1;
    return 0;
} /* cinvit_ */


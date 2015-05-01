/* invit.f -- translated by f2c (version 12.02.01).
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

/* DECK INVIT */
/* Subroutine */ int invit_(integer *nm, integer *n, real *a, real *wr, real *
	wi, logical *select, integer *mm, integer *m, real *z__, integer *
	ierr, real *rm1, real *rv1, real *rv2)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, rm1_dim1, rm1_offset, i__1, 
	    i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, s;
    static real t, w, x, y;
    static integer n1, ii, ip, mp, ns, uk, km1, ip1, its;
    static real eps3;
    extern /* Subroutine */ int cdiv_(real *, real *, real *, real *, real *, 
	    real *);
    static real norm, normv, ilambd, rlambd;
    extern doublereal pythag_(real *, real *);
    static real growto, ukroot;

/* ***BEGIN PROLOGUE  INVIT */
/* ***PURPOSE  Compute the eigenvectors of a real upper Hessenberg */
/*            matrix associated with specified eigenvalues by inverse */
/*            iteration. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2B */
/* ***TYPE      SINGLE PRECISION (INVIT-S, CINVIT-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure INVIT */
/*     by Peters and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971). */

/*     This subroutine finds those eigenvectors of a REAL UPPER */
/*     Hessenberg matrix corresponding to specified eigenvalues, */
/*     using inverse iteration. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        A contains the upper Hessenberg matrix.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,N). */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues of the Hessenberg matrix.  The eigenvalues */
/*          must be stored in a manner identical to that output by */
/*          subroutine  HQR,  which recognizes possible splitting of the */
/*          matrix.  WR and WI are one-dimensional REAL arrays, */
/*          dimensioned WR(N) and WI(N). */

/*        SELECT specifies the eigenvectors to be found. The */
/*          eigenvector corresponding to the J-th eigenvalue is */
/*          specified by setting SELECT(J) to .TRUE.  SELECT is a */
/*          one-dimensional LOGICAL array, dimensioned SELECT(N). */

/*        MM should be set to an upper bound for the number of */
/*          columns required to store the eigenvectors to be found. */
/*          NOTE that two columns are required to store the */
/*          eigenvector corresponding to a complex eigenvalue.  One */
/*          column is required to store the eigenvector corresponding */
/*          to a real eigenvalue.  MM is an INTEGER variable. */

/*     On OUTPUT */

/*        A and WI are unaltered. */

/*        WR may have been altered since close eigenvalues are perturbed */
/*          slightly in searching for independent eigenvectors. */

/*        SELECT may have been altered.  If the elements corresponding */
/*          to a pair of conjugate complex eigenvalues were each */
/*          initially set to .TRUE., the program resets the second of */
/*          the two elements to .FALSE. */

/*        M is the number of columns actually used to store the */
/*          eigenvectors.  M is an INTEGER variable. */

/*        Z contains the real and imaginary parts of the eigenvectors. */
/*          The eigenvectors are packed into the columns of Z starting */
/*          at the first column.  If the next selected eigenvalue is */
/*          real, the next column of Z contains its eigenvector.  If the */
/*          eigenvalue is complex, the next two columns of Z contain the */
/*          real and imaginary parts of its eigenvector, with the real */
/*          part first.  The eigenvectors are normalized so that the */
/*          component of largest magnitude is 1. Any vector which fails */
/*          the acceptance test is set to zero.  Z is a two-dimensional */
/*          REAL array, dimensioned Z(NM,MM). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          -(2*N+1)   if more than MM columns of Z are necessary */
/*                     to store the eigenvectors corresponding to */
/*                     the specified eigenvalues (in this case, M is */
/*                     equal to the number of columns of Z containing */
/*                     eigenvectors already computed), */
/*          -K         if the iteration corresponding to the K-th */
/*                     value fails (if this occurs more than once, K */
/*                     is the index of the last occurrence); the */
/*                     corresponding columns of Z are set to zero */
/*                     vectors, */
/*          -(N+K)     if both error situations occur. */

/*        RM1 is a two-dimensional REAL array used for temporary storage. */
/*          This array holds the triangularized form of the upper */
/*          Hessenberg matrix used in the inverse iteration process. */
/*          RM1 is dimensioned RM1(N,N). */

/*        RV1 and RV2 are one-dimensional REAL arrays used for temporary */
/*          storage.  They hold the approximate eigenvectors during the */
/*          inverse iteration process.  RV1 and RV2 are dimensioned */
/*          RV1(N) and RV2(N). */

/*     The ALGOL procedure GUESSVEC appears in INVIT in-line. */

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
/* ***END PROLOGUE  INVIT */


/* ***FIRST EXECUTABLE STATEMENT  INVIT */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
/*     .......... IP = 0, REAL EIGENVALUE */
/*                     1, FIRST OF CONJUGATE COMPLEX PAIR */
/*                    -1, SECOND OF CONJUGATE COMPLEX PAIR .......... */
    ip = 0;
    n1 = *n - 1;

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (wi[k] == 0.f || ip < 0) {
	    goto L100;
	}
	ip = 1;
	if (select[k] && select[k + 1]) {
	    select[k + 1] = FALSE_;
	}
L100:
	if (! select[k]) {
	    goto L960;
	}
	if (wi[k] != 0.f) {
	    ++s;
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
	    if (a[uk + 1 + uk * a_dim1] == 0.f) {
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
		x += (r__1 = a[i__ + j * a_dim1], dabs(r__1));
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
/*     .......... GROWTO IS THE CRITERION FOR THE GROWTH .......... */
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
/*     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH .......... */
	ip1 = k + ip;
	wr[ip1] = rlambd;
/*     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED) */
/*                AND INITIAL REAL VECTOR .......... */
L280:
	mp = 1;

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {

	    i__3 = uk;
	    for (j = mp; j <= i__3; ++j) {
/* L300: */
		rm1[j + i__ * rm1_dim1] = a[i__ + j * a_dim1];
	    }

	    rm1[i__ + i__ * rm1_dim1] -= rlambd;
	    mp = i__;
	    rv1[i__] = eps3;
/* L320: */
	}

	its = 0;
	if (ilambd != 0.f) {
	    goto L520;
	}
/*     .......... REAL EIGENVALUE. */
/*                TRIANGULAR DECOMPOSITION WITH INTERCHANGES, */
/*                REPLACING ZERO PIVOTS BY EPS3 .......... */
	if (uk == 1) {
	    goto L420;
	}

	i__2 = uk;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    mp = i__ - 1;
	    if ((r__1 = rm1[mp + i__ * rm1_dim1], dabs(r__1)) <= (r__2 = rm1[
		    mp + mp * rm1_dim1], dabs(r__2))) {
		goto L360;
	    }

	    i__3 = uk;
	    for (j = mp; j <= i__3; ++j) {
		y = rm1[j + i__ * rm1_dim1];
		rm1[j + i__ * rm1_dim1] = rm1[j + mp * rm1_dim1];
		rm1[j + mp * rm1_dim1] = y;
/* L340: */
	    }

L360:
	    if (rm1[mp + mp * rm1_dim1] == 0.f) {
		rm1[mp + mp * rm1_dim1] = eps3;
	    }
	    x = rm1[mp + i__ * rm1_dim1] / rm1[mp + mp * rm1_dim1];
	    if (x == 0.f) {
		goto L400;
	    }

	    i__3 = uk;
	    for (j = i__; j <= i__3; ++j) {
/* L380: */
		rm1[j + i__ * rm1_dim1] -= x * rm1[j + mp * rm1_dim1];
	    }

L400:
	    ;
	}

L420:
	if (rm1[uk + uk * rm1_dim1] == 0.f) {
	    rm1[uk + uk * rm1_dim1] = eps3;
	}
/*     .......... BACK SUBSTITUTION FOR REAL VECTOR */
/*                FOR I=UK STEP -1 UNTIL 1 DO -- .......... */
L440:
	i__2 = uk;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = uk + 1 - ii;
	    y = rv1[i__];
	    if (i__ == uk) {
		goto L480;
	    }
	    ip1 = i__ + 1;

	    i__3 = uk;
	    for (j = ip1; j <= i__3; ++j) {
/* L460: */
		y -= rm1[j + i__ * rm1_dim1] * rv1[j];
	    }

L480:
	    rv1[i__] = y / rm1[i__ + i__ * rm1_dim1];
/* L500: */
	}

	goto L740;
/*     .......... COMPLEX EIGENVALUE. */
/*                TRIANGULAR DECOMPOSITION WITH INTERCHANGES, */
/*                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY */
/*                PARTS IN UPPER TRIANGLE STARTING AT (1,3) .......... */
L520:
	ns = *n - s;
	z__[(s - 1) * z_dim1 + 1] = -ilambd;
	z__[s * z_dim1 + 1] = 0.f;
	if (*n == 2) {
	    goto L550;
	}
	rm1[rm1_dim1 * 3 + 1] = -ilambd;
	z__[(s - 1) * z_dim1 + 1] = 0.f;
	if (*n == 3) {
	    goto L550;
	}

	i__2 = *n;
	for (i__ = 4; i__ <= i__2; ++i__) {
/* L540: */
	    rm1[i__ * rm1_dim1 + 1] = 0.f;
	}

L550:
	i__2 = uk;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    mp = i__ - 1;
	    w = rm1[mp + i__ * rm1_dim1];
	    if (i__ < *n) {
		t = rm1[mp + (i__ + 1) * rm1_dim1];
	    }
	    if (i__ == *n) {
		t = z__[mp + (s - 1) * z_dim1];
	    }
	    x = rm1[mp + mp * rm1_dim1] * rm1[mp + mp * rm1_dim1] + t * t;
	    if (w * w <= x) {
		goto L580;
	    }
	    x = rm1[mp + mp * rm1_dim1] / w;
	    y = t / w;
	    rm1[mp + mp * rm1_dim1] = w;
	    if (i__ < *n) {
		rm1[mp + (i__ + 1) * rm1_dim1] = 0.f;
	    }
	    if (i__ == *n) {
		z__[mp + (s - 1) * z_dim1] = 0.f;
	    }

	    i__3 = uk;
	    for (j = i__; j <= i__3; ++j) {
		w = rm1[j + i__ * rm1_dim1];
		rm1[j + i__ * rm1_dim1] = rm1[j + mp * rm1_dim1] - x * w;
		rm1[j + mp * rm1_dim1] = w;
		if (j < n1) {
		    goto L555;
		}
		l = j - ns;
		z__[i__ + l * z_dim1] = z__[mp + l * z_dim1] - y * w;
		z__[mp + l * z_dim1] = 0.f;
		goto L560;
L555:
		rm1[i__ + (j + 2) * rm1_dim1] = rm1[mp + (j + 2) * rm1_dim1] 
			- y * w;
		rm1[mp + (j + 2) * rm1_dim1] = 0.f;
L560:
		;
	    }

	    rm1[i__ + i__ * rm1_dim1] -= y * ilambd;
	    if (i__ < n1) {
		goto L570;
	    }
	    l = i__ - ns;
	    z__[mp + l * z_dim1] = -ilambd;
	    z__[i__ + l * z_dim1] += x * ilambd;
	    goto L640;
L570:
	    rm1[mp + (i__ + 2) * rm1_dim1] = -ilambd;
	    rm1[i__ + (i__ + 2) * rm1_dim1] += x * ilambd;
	    goto L640;
L580:
	    if (x != 0.f) {
		goto L600;
	    }
	    rm1[mp + mp * rm1_dim1] = eps3;
	    if (i__ < *n) {
		rm1[mp + (i__ + 1) * rm1_dim1] = 0.f;
	    }
	    if (i__ == *n) {
		z__[mp + (s - 1) * z_dim1] = 0.f;
	    }
	    t = 0.f;
	    x = eps3 * eps3;
L600:
	    w /= x;
	    x = rm1[mp + mp * rm1_dim1] * w;
	    y = -t * w;

	    i__3 = uk;
	    for (j = i__; j <= i__3; ++j) {
		if (j < n1) {
		    goto L610;
		}
		l = j - ns;
		t = z__[mp + l * z_dim1];
		z__[i__ + l * z_dim1] = -x * t - y * rm1[j + mp * rm1_dim1];
		goto L615;
L610:
		t = rm1[mp + (j + 2) * rm1_dim1];
		rm1[i__ + (j + 2) * rm1_dim1] = -x * t - y * rm1[j + mp * 
			rm1_dim1];
L615:
		rm1[j + i__ * rm1_dim1] = rm1[j + i__ * rm1_dim1] - x * rm1[j 
			+ mp * rm1_dim1] + y * t;
/* L620: */
	    }

	    if (i__ < n1) {
		goto L630;
	    }
	    l = i__ - ns;
	    z__[i__ + l * z_dim1] -= ilambd;
	    goto L640;
L630:
	    rm1[i__ + (i__ + 2) * rm1_dim1] -= ilambd;
L640:
	    ;
	}

	if (uk < n1) {
	    goto L650;
	}
	l = uk - ns;
	t = z__[uk + l * z_dim1];
	goto L655;
L650:
	t = rm1[uk + (uk + 2) * rm1_dim1];
L655:
	if (rm1[uk + uk * rm1_dim1] == 0.f && t == 0.f) {
	    rm1[uk + uk * rm1_dim1] = eps3;
	}
/*     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR */
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
		if (j < n1) {
		    goto L670;
		}
		l = j - ns;
		t = z__[i__ + l * z_dim1];
		goto L675;
L670:
		t = rm1[i__ + (j + 2) * rm1_dim1];
L675:
		x = x - rm1[j + i__ * rm1_dim1] * rv1[j] + t * rv2[j];
		y = y - rm1[j + i__ * rm1_dim1] * rv2[j] - t * rv1[j];
/* L680: */
	    }

L700:
	    if (i__ < n1) {
		goto L710;
	    }
	    l = i__ - ns;
	    t = z__[i__ + l * z_dim1];
	    goto L715;
L710:
	    t = rm1[i__ + (i__ + 2) * rm1_dim1];
L715:
	    cdiv_(&x, &y, &rm1[i__ + i__ * rm1_dim1], &t, &rv1[i__], &rv2[i__]
		    );
/* L720: */
	}
/*     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX */
/*                EIGENVECTOR AND NORMALIZATION .......... */
L740:
	++its;
	norm = 0.f;
	normv = 0.f;

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (ilambd == 0.f) {
		x = (r__1 = rv1[i__], dabs(r__1));
	    }
	    if (ilambd != 0.f) {
		x = pythag_(&rv1[i__], &rv2[i__]);
	    }
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
	if (ilambd == 0.f) {
	    x = 1.f / x;
	}
	if (ilambd != 0.f) {
	    y = rv2[j];
	}

	i__2 = uk;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (ilambd != 0.f) {
		goto L800;
	    }
	    z__[i__ + s * z_dim1] = rv1[i__] * x;
	    goto L820;
L800:
	    cdiv_(&rv1[i__], &rv2[i__], &x, &y, &z__[i__ + (s - 1) * z_dim1], 
		    &z__[i__ + s * z_dim1]);
L820:
	    ;
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
	if (ilambd == 0.f) {
	    goto L440;
	}
	goto L660;
/*     .......... SET ERROR -- UNACCEPTED EIGENVECTOR .......... */
L880:
	j = 1;
	*ierr = -k;
/*     .......... SET REMAINING VECTOR COMPONENTS TO ZERO .......... */
L900:
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    z__[i__ + s * z_dim1] = 0.f;
	    if (ilambd != 0.f) {
		z__[i__ + (s - 1) * z_dim1] = 0.f;
	    }
/* L920: */
	}

L940:
	++s;
L960:
	if (ip == -1) {
	    ip = 0;
	}
	if (ip == 1) {
	    ip = -1;
	}
/* L980: */
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
    *m = s - 1 - abs(ip);
    return 0;
} /* invit_ */


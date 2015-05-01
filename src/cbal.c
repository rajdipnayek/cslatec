/* cbal.f -- translated by f2c (version 12.02.01).
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

/* DECK CBAL */
/* Subroutine */ int cbal_(integer *nm, integer *n, real *ar, real *ai, 
	integer *low, integer *igh, real *scale)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real c__, f, g;
    static integer i__, j, k, l, m;
    static real r__, s, b2;
    static integer jj, iexc;
    static real radix;
    static logical noconv;

/* ***BEGIN PROLOGUE  CBAL */
/* ***PURPOSE  Balance a complex general matrix and isolate eigenvalues */
/*            whenever possible. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1A */
/* ***TYPE      COMPLEX (BALANC-S, CBAL-C) */
/* ***KEYWORDS  EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure */
/*     CBALANCE, which is a complex version of BALANCE, */
/*     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971). */

/*     This subroutine balances a COMPLEX matrix and isolates */
/*     eigenvalues whenever possible. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR and AI, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A=(AR,AI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        AR and AI contain the real and imaginary parts, */
/*          respectively, of the complex matrix to be balanced. */
/*          AR and AI are two-dimensional REAL arrays, dimensioned */
/*          AR(NM,N) and AI(NM,N). */

/*     On OUTPUT */

/*        AR and AI contain the real and imaginary parts, */
/*          respectively, of the balanced matrix. */

/*        LOW and IGH are two INTEGER variables such that AR(I,J) */
/*          and AI(I,J) are equal to zero if */
/*           (1) I is greater than J and */
/*           (2) J=1,...,LOW-1 or I=IGH+1,...,N. */

/*        SCALE contains information determining the permutations and */
/*          scaling factors used.  SCALE is a one-dimensional REAL array, */
/*          dimensioned SCALE(N). */

/*     Suppose that the principal submatrix in rows LOW through IGH */
/*     has been balanced, that P(J) denotes the index interchanged */
/*     with J during the permutation step, and that the elements */
/*     of the diagonal matrix used are denoted by D(I,J).  Then */
/*        SCALE(J) = P(J),    for J = 1,...,LOW-1 */
/*                 = D(J,J)       J = LOW,...,IGH */
/*                 = P(J)         J = IGH+1,...,N. */
/*     The order in which the interchanges are made is N to IGH+1, */
/*     then 1 to LOW-1. */

/*     Note that 1 is returned for IGH if IGH is zero formally. */

/*     The ALGOL procedure EXC contained in CBALANCE appears in */
/*     CBAL  in line.  (Note that the ALGOL roles of identifiers */
/*     K,L have been reversed.) */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CBAL */


/*     THE FOLLOWING PORTABLE VALUE OF RADIX WORKS WELL ENOUGH */
/*     FOR ALL MACHINES WHOSE BASE IS A POWER OF TWO. */

/* ***FIRST EXECUTABLE STATEMENT  CBAL */
    /* Parameter adjustments */
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
    --scale;

    /* Function Body */
    radix = 16.f;

    b2 = radix * radix;
    k = 1;
    l = *n;
    goto L100;
/*     .......... IN-LINE PROCEDURE FOR ROW AND */
/*                COLUMN EXCHANGE .......... */
L20:
    scale[m] = (real) j;
    if (j == m) {
	goto L50;
    }

    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f = ar[i__ + j * ar_dim1];
	ar[i__ + j * ar_dim1] = ar[i__ + m * ar_dim1];
	ar[i__ + m * ar_dim1] = f;
	f = ai[i__ + j * ai_dim1];
	ai[i__ + j * ai_dim1] = ai[i__ + m * ai_dim1];
	ai[i__ + m * ai_dim1] = f;
/* L30: */
    }

    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) {
	f = ar[j + i__ * ar_dim1];
	ar[j + i__ * ar_dim1] = ar[m + i__ * ar_dim1];
	ar[m + i__ * ar_dim1] = f;
	f = ai[j + i__ * ai_dim1];
	ai[j + i__ * ai_dim1] = ai[m + i__ * ai_dim1];
	ai[m + i__ * ai_dim1] = f;
/* L40: */
    }

L50:
    switch (iexc) {
	case 1:  goto L80;
	case 2:  goto L130;
    }
/*     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE */
/*                AND PUSH THEM DOWN .......... */
L80:
    if (l == 1) {
	goto L280;
    }
    --l;
/*     .......... FOR J=L STEP -1 UNTIL 1 DO -- .......... */
L100:
    i__1 = l;
    for (jj = 1; jj <= i__1; ++jj) {
	j = l + 1 - jj;

	i__2 = l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		goto L110;
	    }
	    if (ar[j + i__ * ar_dim1] != 0.f || ai[j + i__ * ai_dim1] != 0.f) 
		    {
		goto L120;
	    }
L110:
	    ;
	}

	m = l;
	iexc = 1;
	goto L20;
L120:
	;
    }

    goto L140;
/*     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE */
/*                AND PUSH THEM LEFT .......... */
L130:
    ++k;

L140:
    i__1 = l;
    for (j = k; j <= i__1; ++j) {

	i__2 = l;
	for (i__ = k; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		goto L150;
	    }
	    if (ar[i__ + j * ar_dim1] != 0.f || ai[i__ + j * ai_dim1] != 0.f) 
		    {
		goto L170;
	    }
L150:
	    ;
	}

	m = k;
	iexc = 2;
	goto L20;
L170:
	;
    }
/*     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L .......... */
    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
/* L180: */
	scale[i__] = 1.f;
    }
/*     .......... ITERATIVE LOOP FOR NORM REDUCTION .......... */
L190:
    noconv = FALSE_;

    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
	c__ = 0.f;
	r__ = 0.f;

	i__2 = l;
	for (j = k; j <= i__2; ++j) {
	    if (j == i__) {
		goto L200;
	    }
	    c__ = c__ + (r__1 = ar[j + i__ * ar_dim1], dabs(r__1)) + (r__2 = 
		    ai[j + i__ * ai_dim1], dabs(r__2));
	    r__ = r__ + (r__1 = ar[i__ + j * ar_dim1], dabs(r__1)) + (r__2 = 
		    ai[i__ + j * ai_dim1], dabs(r__2));
L200:
	    ;
	}
/*     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW .......... */
	if (c__ == 0.f || r__ == 0.f) {
	    goto L270;
	}
	g = r__ / radix;
	f = 1.f;
	s = c__ + r__;
L210:
	if (c__ >= g) {
	    goto L220;
	}
	f *= radix;
	c__ *= b2;
	goto L210;
L220:
	g = r__ * radix;
L230:
	if (c__ < g) {
	    goto L240;
	}
	f /= radix;
	c__ /= b2;
	goto L230;
/*     .......... NOW BALANCE .......... */
L240:
	if ((c__ + r__) / f >= s * .95f) {
	    goto L270;
	}
	g = 1.f / f;
	scale[i__] *= f;
	noconv = TRUE_;

	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
	    ar[i__ + j * ar_dim1] *= g;
	    ai[i__ + j * ai_dim1] *= g;
/* L250: */
	}

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    ar[j + i__ * ar_dim1] *= f;
	    ai[j + i__ * ai_dim1] *= f;
/* L260: */
	}

L270:
	;
    }

    if (noconv) {
	goto L190;
    }

L280:
    *low = k;
    *igh = l;
    return 0;
} /* cbal_ */


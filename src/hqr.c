/* hqr.f -- translated by f2c (version 12.02.01).
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

/* DECK HQR */
/* Subroutine */ int hqr_(integer *nm, integer *n, integer *low, integer *igh,
	 real *h__, real *wr, real *wi, integer *ierr)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real p, q, r__, s, t, w, x, y, s1, s2;
    static integer na, en, ll, mm;
    static real zz;
    static integer mp2, itn, its, enm2;
    static real norm;
    static logical notlas;

/* ***BEGIN PROLOGUE  HQR */
/* ***PURPOSE  Compute the eigenvalues of a real upper Hessenberg matrix */
/*            using the QR method. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2B */
/* ***TYPE      SINGLE PRECISION (HQR-S, COMQR-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure HQR, */
/*     NUM. MATH. 14, 219-231(1970) by Martin, Peters, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971). */

/*     This subroutine finds the eigenvalues of a REAL */
/*     UPPER Hessenberg matrix by the QR method. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, H, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix H.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  BALANC.  If  BALANC  has not been */
/*          used, set LOW=1 and IGH equal to the order of the matrix, N. */

/*        H contains the upper Hessenberg matrix.  Information about */
/*          the transformations used in the reduction to Hessenberg */
/*          form by  ELMHES  or  ORTHES, if performed, is stored */
/*          in the remaining triangle under the Hessenberg matrix. */
/*          H is a two-dimensional REAL array, dimensioned H(NM,N). */

/*     On OUTPUT */

/*        H has been destroyed.  Therefore, it must be saved before */
/*          calling  HQR  if subsequent calculation and back */
/*          transformation of eigenvectors is to be performed. */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues.  The eigenvalues are unordered except */
/*          that complex conjugate pairs of values appear consecutively */
/*          with the eigenvalue having the positive imaginary part first. */
/*          If an error exit is made, the eigenvalues should be correct */
/*          for indices IERR+1, IERR+2, ..., N.  WR and WI are one- */
/*          dimensional REAL arrays, dimensioned WR(N) and WI(N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30*N iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     IERR+1, IERR+2, ..., N. */

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
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HQR */


/* ***FIRST EXECUTABLE STATEMENT  HQR */
    /* Parameter adjustments */
    h_dim1 = *nm;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;

    /* Function Body */
    *ierr = 0;
    norm = 0.f;
    k = 1;
/*     .......... STORE ROOTS ISOLATED BY BALANC */
/*                AND COMPUTE MATRIX NORM .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
/* L40: */
	    norm += (r__1 = h__[i__ + j * h_dim1], dabs(r__1));
	}

	k = i__;
	if (i__ >= *low && i__ <= *igh) {
	    goto L50;
	}
	wr[i__] = h__[i__ + i__ * h_dim1];
	wi[i__] = 0.f;
L50:
	;
    }

    en = *igh;
    t = 0.f;
    itn = *n * 30;
/*     .......... SEARCH FOR NEXT EIGENVALUES .......... */
L60:
    if (en < *low) {
	goto L1001;
    }
    its = 0;
    na = en - 1;
    enm2 = na - 1;
/*     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT */
/*                FOR L=EN STEP -1 UNTIL LOW DO -- .......... */
L70:
    i__1 = en;
    for (ll = *low; ll <= i__1; ++ll) {
	l = en + *low - ll;
	if (l == *low) {
	    goto L100;
	}
	s = (r__1 = h__[l - 1 + (l - 1) * h_dim1], dabs(r__1)) + (r__2 = h__[
		l + l * h_dim1], dabs(r__2));
	if (s == 0.f) {
	    s = norm;
	}
	s2 = s + (r__1 = h__[l + (l - 1) * h_dim1], dabs(r__1));
	if (s2 == s) {
	    goto L100;
	}
/* L80: */
    }
/*     .......... FORM SHIFT .......... */
L100:
    x = h__[en + en * h_dim1];
    if (l == en) {
	goto L270;
    }
    y = h__[na + na * h_dim1];
    w = h__[en + na * h_dim1] * h__[na + en * h_dim1];
    if (l == na) {
	goto L280;
    }
    if (itn == 0) {
	goto L1000;
    }
    if (its != 10 && its != 20) {
	goto L130;
    }
/*     .......... FORM EXCEPTIONAL SHIFT .......... */
    t += x;

    i__1 = en;
    for (i__ = *low; i__ <= i__1; ++i__) {
/* L120: */
	h__[i__ + i__ * h_dim1] -= x;
    }

    s = (r__1 = h__[en + na * h_dim1], dabs(r__1)) + (r__2 = h__[na + enm2 * 
	    h_dim1], dabs(r__2));
    x = s * .75f;
    y = x;
    w = s * -.4375f * s;
L130:
    ++its;
    --itn;
/*     .......... LOOK FOR TWO CONSECUTIVE SMALL */
/*                SUB-DIAGONAL ELEMENTS. */
/*                FOR M=EN-2 STEP -1 UNTIL L DO -- .......... */
    i__1 = enm2;
    for (mm = l; mm <= i__1; ++mm) {
	m = enm2 + l - mm;
	zz = h__[m + m * h_dim1];
	r__ = x - zz;
	s = y - zz;
	p = (r__ * s - w) / h__[m + 1 + m * h_dim1] + h__[m + (m + 1) * 
		h_dim1];
	q = h__[m + 1 + (m + 1) * h_dim1] - zz - r__ - s;
	r__ = h__[m + 2 + (m + 1) * h_dim1];
	s = dabs(p) + dabs(q) + dabs(r__);
	p /= s;
	q /= s;
	r__ /= s;
	if (m == l) {
	    goto L150;
	}
	s1 = dabs(p) * ((r__1 = h__[m - 1 + (m - 1) * h_dim1], dabs(r__1)) + 
		dabs(zz) + (r__2 = h__[m + 1 + (m + 1) * h_dim1], dabs(r__2)))
		;
	s2 = s1 + (r__1 = h__[m + (m - 1) * h_dim1], dabs(r__1)) * (dabs(q) + 
		dabs(r__));
	if (s2 == s1) {
	    goto L150;
	}
/* L140: */
    }

L150:
    mp2 = m + 2;

    i__1 = en;
    for (i__ = mp2; i__ <= i__1; ++i__) {
	h__[i__ + (i__ - 2) * h_dim1] = 0.f;
	if (i__ == mp2) {
	    goto L160;
	}
	h__[i__ + (i__ - 3) * h_dim1] = 0.f;
L160:
	;
    }
/*     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND */
/*                COLUMNS M TO EN .......... */
    i__1 = na;
    for (k = m; k <= i__1; ++k) {
	notlas = k != na;
	if (k == m) {
	    goto L170;
	}
	p = h__[k + (k - 1) * h_dim1];
	q = h__[k + 1 + (k - 1) * h_dim1];
	r__ = 0.f;
	if (notlas) {
	    r__ = h__[k + 2 + (k - 1) * h_dim1];
	}
	x = dabs(p) + dabs(q) + dabs(r__);
	if (x == 0.f) {
	    goto L260;
	}
	p /= x;
	q /= x;
	r__ /= x;
L170:
	r__1 = sqrt(p * p + q * q + r__ * r__);
	s = r_sign(&r__1, &p);
	if (k == m) {
	    goto L180;
	}
	h__[k + (k - 1) * h_dim1] = -s * x;
	goto L190;
L180:
	if (l != m) {
	    h__[k + (k - 1) * h_dim1] = -h__[k + (k - 1) * h_dim1];
	}
L190:
	p += s;
	x = p / s;
	y = q / s;
	zz = r__ / s;
	q /= p;
	r__ /= p;
/*     .......... ROW MODIFICATION .......... */
	i__2 = en;
	for (j = k; j <= i__2; ++j) {
	    p = h__[k + j * h_dim1] + q * h__[k + 1 + j * h_dim1];
	    if (! notlas) {
		goto L200;
	    }
	    p += r__ * h__[k + 2 + j * h_dim1];
	    h__[k + 2 + j * h_dim1] -= p * zz;
L200:
	    h__[k + 1 + j * h_dim1] -= p * y;
	    h__[k + j * h_dim1] -= p * x;
/* L210: */
	}

/* Computing MIN */
	i__2 = en, i__3 = k + 3;
	j = min(i__2,i__3);
/*     .......... COLUMN MODIFICATION .......... */
	i__2 = j;
	for (i__ = l; i__ <= i__2; ++i__) {
	    p = x * h__[i__ + k * h_dim1] + y * h__[i__ + (k + 1) * h_dim1];
	    if (! notlas) {
		goto L220;
	    }
	    p += zz * h__[i__ + (k + 2) * h_dim1];
	    h__[i__ + (k + 2) * h_dim1] -= p * r__;
L220:
	    h__[i__ + (k + 1) * h_dim1] -= p * q;
	    h__[i__ + k * h_dim1] -= p;
/* L230: */
	}

L260:
	;
    }

    goto L70;
/*     .......... ONE ROOT FOUND .......... */
L270:
    wr[en] = x + t;
    wi[en] = 0.f;
    en = na;
    goto L60;
/*     .......... TWO ROOTS FOUND .......... */
L280:
    p = (y - x) / 2.f;
    q = p * p + w;
    zz = sqrt((dabs(q)));
    x += t;
    if (q < 0.f) {
	goto L320;
    }
/*     .......... REAL PAIR .......... */
    zz = p + r_sign(&zz, &p);
    wr[na] = x + zz;
    wr[en] = wr[na];
    if (zz != 0.f) {
	wr[en] = x - w / zz;
    }
    wi[na] = 0.f;
    wi[en] = 0.f;
    goto L330;
/*     .......... COMPLEX PAIR .......... */
L320:
    wr[na] = x + p;
    wr[en] = x + p;
    wi[na] = zz;
    wi[en] = -zz;
L330:
    en = enm2;
    goto L60;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30*N ITERATIONS .......... */
L1000:
    *ierr = en;
L1001:
    return 0;
} /* hqr_ */


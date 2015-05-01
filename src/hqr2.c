/* hqr2.f -- translated by f2c (version 12.02.01).
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

static real c_b46 = 0.f;

/* DECK HQR2 */
/* Subroutine */ int hqr2_(integer *nm, integer *n, integer *low, integer *
	igh, real *h__, real *wr, real *wi, real *z__, integer *ierr)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real p, q, r__, s, t, w, x, y, s1, s2;
    static integer na, ii, en, jj;
    static real ra, sa;
    static integer ll, mm, nn;
    static real vi, vr, zz;
    static integer mp2, itn, its, enm2;
    extern /* Subroutine */ int cdiv_(real *, real *, real *, real *, real *, 
	    real *);
    static real norm;
    static logical notlas;

/* ***BEGIN PROLOGUE  HQR2 */
/* ***PURPOSE  Compute the eigenvalues and eigenvectors of a real upper */
/*            Hessenberg matrix using QR method. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2B */
/* ***TYPE      SINGLE PRECISION (HQR2-S, COMQR2-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure HQR2, */
/*     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971). */

/*     This subroutine finds the eigenvalues and eigenvectors */
/*     of a REAL UPPER Hessenberg matrix by the QR method.  The */
/*     eigenvectors of a REAL GENERAL matrix can also be found */
/*     if  ELMHES  and  ELTRAN  or  ORTHES  and  ORTRAN  have */
/*     been used to reduce this general matrix to Hessenberg form */
/*     and to accumulate the similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, H and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix H.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  BALANC.  If  BALANC  has not been */
/*          used, set LOW=1 and IGH equal to the order of the matrix, N. */

/*        H contains the upper Hessenberg matrix.  H is a two-dimensional */
/*          REAL array, dimensioned H(NM,N). */

/*        Z contains the transformation matrix produced by  ELTRAN */
/*          after the reduction by  ELMHES, or by  ORTRAN  after the */
/*          reduction by  ORTHES, if performed.  If the eigenvectors */
/*          of the Hessenberg matrix are desired, Z must contain the */
/*          identity matrix.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,M). */

/*     On OUTPUT */

/*        H has been destroyed. */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues.  The eigenvalues are unordered except */
/*          that complex conjugate pairs of values appear consecutively */
/*          with the eigenvalue having the positive imaginary part first. */
/*          If an error exit is made, the eigenvalues should be correct */
/*          for indices IERR+1, IERR+2, ..., N.  WR and WI are one- */
/*          dimensional REAL arrays, dimensioned WR(N) and WI(N). */

/*        Z contains the real and imaginary parts of the eigenvectors. */
/*          If the J-th eigenvalue is real, the J-th column of Z */
/*          contains its eigenvector.  If the J-th eigenvalue is complex */
/*          with positive imaginary part, the J-th and (J+1)-th */
/*          columns of Z contain the real and imaginary parts of its */
/*          eigenvector.  The eigenvectors are unnormalized.  If an */
/*          error exit is made, none of the eigenvectors has been found. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30*N iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     IERR+1, IERR+2, ..., N, but no eigenvectors are */
/*                     computed. */

/*     Calls CDIV for complex division. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  CDIV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HQR2 */


/* ***FIRST EXECUTABLE STATEMENT  HQR2 */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
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
	goto L340;
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
	i__2 = *n;
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
	for (i__ = 1; i__ <= i__2; ++i__) {
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
/*     .......... ACCUMULATE TRANSFORMATIONS .......... */
	i__2 = *igh;
	for (i__ = *low; i__ <= i__2; ++i__) {
	    p = x * z__[i__ + k * z_dim1] + y * z__[i__ + (k + 1) * z_dim1];
	    if (! notlas) {
		goto L240;
	    }
	    p += zz * z__[i__ + (k + 2) * z_dim1];
	    z__[i__ + (k + 2) * z_dim1] -= p * r__;
L240:
	    z__[i__ + (k + 1) * z_dim1] -= p * q;
	    z__[i__ + k * z_dim1] -= p;
/* L250: */
	}

L260:
	;
    }

    goto L70;
/*     .......... ONE ROOT FOUND .......... */
L270:
    h__[en + en * h_dim1] = x + t;
    wr[en] = h__[en + en * h_dim1];
    wi[en] = 0.f;
    en = na;
    goto L60;
/*     .......... TWO ROOTS FOUND .......... */
L280:
    p = (y - x) / 2.f;
    q = p * p + w;
    zz = sqrt((dabs(q)));
    h__[en + en * h_dim1] = x + t;
    x = h__[en + en * h_dim1];
    h__[na + na * h_dim1] = y + t;
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
    x = h__[en + na * h_dim1];
    s = dabs(x) + dabs(zz);
    p = x / s;
    q = zz / s;
    r__ = sqrt(p * p + q * q);
    p /= r__;
    q /= r__;
/*     .......... ROW MODIFICATION .......... */
    i__1 = *n;
    for (j = na; j <= i__1; ++j) {
	zz = h__[na + j * h_dim1];
	h__[na + j * h_dim1] = q * zz + p * h__[en + j * h_dim1];
	h__[en + j * h_dim1] = q * h__[en + j * h_dim1] - p * zz;
/* L290: */
    }
/*     .......... COLUMN MODIFICATION .......... */
    i__1 = en;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zz = h__[i__ + na * h_dim1];
	h__[i__ + na * h_dim1] = q * zz + p * h__[i__ + en * h_dim1];
	h__[i__ + en * h_dim1] = q * h__[i__ + en * h_dim1] - p * zz;
/* L300: */
    }
/*     .......... ACCUMULATE TRANSFORMATIONS .......... */
    i__1 = *igh;
    for (i__ = *low; i__ <= i__1; ++i__) {
	zz = z__[i__ + na * z_dim1];
	z__[i__ + na * z_dim1] = q * zz + p * z__[i__ + en * z_dim1];
	z__[i__ + en * z_dim1] = q * z__[i__ + en * z_dim1] - p * zz;
/* L310: */
    }

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
/*     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND */
/*                VECTORS OF UPPER TRIANGULAR FORM .......... */
L340:
    if (norm == 0.f) {
	goto L1001;
    }
/*     .......... FOR EN=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (nn = 1; nn <= i__1; ++nn) {
	en = *n + 1 - nn;
	p = wr[en];
	q = wi[en];
	na = en - 1;
	if (q < 0.f) {
	    goto L710;
	} else if (q == 0) {
	    goto L600;
	} else {
	    goto L800;
	}
/*     .......... REAL VECTOR .......... */
L600:
	m = en;
	h__[en + en * h_dim1] = 1.f;
	if (na == 0) {
	    goto L800;
	}
/*     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- .......... */
	i__2 = na;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = en - ii;
	    w = h__[i__ + i__ * h_dim1] - p;
	    r__ = h__[i__ + en * h_dim1];
	    if (m > na) {
		goto L620;
	    }

	    i__3 = na;
	    for (j = m; j <= i__3; ++j) {
/* L610: */
		r__ += h__[i__ + j * h_dim1] * h__[j + en * h_dim1];
	    }

L620:
	    if (wi[i__] >= 0.f) {
		goto L630;
	    }
	    zz = w;
	    s = r__;
	    goto L700;
L630:
	    m = i__;
	    if (wi[i__] != 0.f) {
		goto L640;
	    }
	    t = w;
	    if (t != 0.f) {
		goto L635;
	    }
	    t = norm;
L632:
	    t *= .5f;
	    if (norm + t > norm) {
		goto L632;
	    }
	    t *= 2.f;
L635:
	    h__[i__ + en * h_dim1] = -r__ / t;
	    goto L700;
/*     .......... SOLVE REAL EQUATIONS .......... */
L640:
	    x = h__[i__ + (i__ + 1) * h_dim1];
	    y = h__[i__ + 1 + i__ * h_dim1];
	    q = (wr[i__] - p) * (wr[i__] - p) + wi[i__] * wi[i__];
	    t = (x * s - zz * r__) / q;
	    h__[i__ + en * h_dim1] = t;
	    if (dabs(x) <= dabs(zz)) {
		goto L650;
	    }
	    h__[i__ + 1 + en * h_dim1] = (-r__ - w * t) / x;
	    goto L700;
L650:
	    h__[i__ + 1 + en * h_dim1] = (-s - y * t) / zz;
L700:
	    ;
	}
/*     .......... END REAL VECTOR .......... */
	goto L800;
/*     .......... COMPLEX VECTOR .......... */
L710:
	m = na;
/*     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT */
/*                EIGENVECTOR MATRIX IS TRIANGULAR .......... */
	if ((r__1 = h__[en + na * h_dim1], dabs(r__1)) <= (r__2 = h__[na + en 
		* h_dim1], dabs(r__2))) {
	    goto L720;
	}
	h__[na + na * h_dim1] = q / h__[en + na * h_dim1];
	h__[na + en * h_dim1] = -(h__[en + en * h_dim1] - p) / h__[en + na * 
		h_dim1];
	goto L730;
L720:
	r__1 = -h__[na + en * h_dim1];
	r__2 = h__[na + na * h_dim1] - p;
	cdiv_(&c_b46, &r__1, &r__2, &q, &h__[na + na * h_dim1], &h__[na + en *
		 h_dim1]);
L730:
	h__[en + na * h_dim1] = 0.f;
	h__[en + en * h_dim1] = 1.f;
	enm2 = na - 1;
	if (enm2 == 0) {
	    goto L800;
	}
/*     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- .......... */
	i__2 = enm2;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = na - ii;
	    w = h__[i__ + i__ * h_dim1] - p;
	    ra = 0.f;
	    sa = h__[i__ + en * h_dim1];

	    i__3 = na;
	    for (j = m; j <= i__3; ++j) {
		ra += h__[i__ + j * h_dim1] * h__[j + na * h_dim1];
		sa += h__[i__ + j * h_dim1] * h__[j + en * h_dim1];
/* L760: */
	    }

	    if (wi[i__] >= 0.f) {
		goto L770;
	    }
	    zz = w;
	    r__ = ra;
	    s = sa;
	    goto L790;
L770:
	    m = i__;
	    if (wi[i__] != 0.f) {
		goto L780;
	    }
	    r__1 = -ra;
	    r__2 = -sa;
	    cdiv_(&r__1, &r__2, &w, &q, &h__[i__ + na * h_dim1], &h__[i__ + 
		    en * h_dim1]);
	    goto L790;
/*     .......... SOLVE COMPLEX EQUATIONS .......... */
L780:
	    x = h__[i__ + (i__ + 1) * h_dim1];
	    y = h__[i__ + 1 + i__ * h_dim1];
	    vr = (wr[i__] - p) * (wr[i__] - p) + wi[i__] * wi[i__] - q * q;
	    vi = (wr[i__] - p) * 2.f * q;
	    if (vr != 0.f || vi != 0.f) {
		goto L783;
	    }
	    s1 = norm * (dabs(w) + dabs(q) + dabs(x) + dabs(y) + dabs(zz));
	    vr = s1;
L782:
	    vr *= .5f;
	    if (s1 + vr > s1) {
		goto L782;
	    }
	    vr *= 2.f;
L783:
	    r__1 = x * r__ - zz * ra + q * sa;
	    r__2 = x * s - zz * sa - q * ra;
	    cdiv_(&r__1, &r__2, &vr, &vi, &h__[i__ + na * h_dim1], &h__[i__ + 
		    en * h_dim1]);
	    if (dabs(x) <= dabs(zz) + dabs(q)) {
		goto L785;
	    }
	    h__[i__ + 1 + na * h_dim1] = (-ra - w * h__[i__ + na * h_dim1] + 
		    q * h__[i__ + en * h_dim1]) / x;
	    h__[i__ + 1 + en * h_dim1] = (-sa - w * h__[i__ + en * h_dim1] - 
		    q * h__[i__ + na * h_dim1]) / x;
	    goto L790;
L785:
	    r__1 = -r__ - y * h__[i__ + na * h_dim1];
	    r__2 = -s - y * h__[i__ + en * h_dim1];
	    cdiv_(&r__1, &r__2, &zz, &q, &h__[i__ + 1 + na * h_dim1], &h__[
		    i__ + 1 + en * h_dim1]);
L790:
	    ;
	}
/*     .......... END COMPLEX VECTOR .......... */
L800:
	;
    }
/*     .......... END BACK SUBSTITUTION. */
/*                VECTORS OF ISOLATED ROOTS .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ >= *low && i__ <= *igh) {
	    goto L840;
	}

	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
/* L820: */
	    z__[i__ + j * z_dim1] = h__[i__ + j * h_dim1];
	}

L840:
	;
    }
/*     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE */
/*                VECTORS OF ORIGINAL FULL MATRIX. */
/*                FOR J=N STEP -1 UNTIL LOW DO -- .......... */
    i__1 = *n;
    for (jj = *low; jj <= i__1; ++jj) {
	j = *n + *low - jj;
	m = min(j,*igh);

	i__2 = *igh;
	for (i__ = *low; i__ <= i__2; ++i__) {
	    zz = 0.f;

	    i__3 = m;
	    for (k = *low; k <= i__3; ++k) {
/* L860: */
		zz += z__[i__ + k * z_dim1] * h__[k + j * h_dim1];
	    }

	    z__[i__ + j * z_dim1] = zz;
/* L880: */
	}
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30*N ITERATIONS .......... */
L1000:
    *ierr = en;
L1001:
    return 0;
} /* hqr2_ */


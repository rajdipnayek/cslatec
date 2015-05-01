/* comqr2.f -- translated by f2c (version 12.02.01).
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

/* DECK COMQR2 */
/* Subroutine */ int comqr2_(integer *nm, integer *n, integer *low, integer *
	igh, real *ortr, real *orti, real *hr, real *hi, real *wr, real *wi, 
	real *zr, real *zi, integer *ierr)
{
    /* System generated locals */
    integer hr_dim1, hr_offset, hi_dim1, hi_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real s1, s2;
    static integer ii, en, jj, ll, nn;
    static real si, ti, xi, yi, sr, tr, xr, yr;
    static integer ip1, lp1, itn, its;
    static real zzi, zzr;
    static integer enm1, iend;
    extern /* Subroutine */ int cdiv_(real *, real *, real *, real *, real *, 
	    real *);
    static real norm;
    extern doublereal pythag_(real *, real *);
    extern /* Subroutine */ int csroot_(real *, real *, real *, real *);

/* ***BEGIN PROLOGUE  COMQR2 */
/* ***PURPOSE  Compute the eigenvalues and eigenvectors of a complex upper */
/*            Hessenberg matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2B */
/* ***TYPE      COMPLEX (HQR2-S, COMQR2-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of a unitary analogue of the */
/*     ALGOL procedure  COMLR2, NUM. MATH. 16, 181-204(1970) by Peters */
/*     and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971). */
/*     The unitary analogue substitutes the QR algorithm of Francis */
/*     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm. */

/*     This subroutine finds the eigenvalues and eigenvectors */
/*     of a COMPLEX UPPER Hessenberg matrix by the QR */
/*     method.  The eigenvectors of a COMPLEX GENERAL matrix */
/*     can also be found if  CORTH  has been used to reduce */
/*     this general matrix to Hessenberg form. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, HR, HI, ZR, and ZI, as declared in the */
/*          calling program dimension statement.  NM is an INTEGER */
/*          variable. */

/*        N is the order of the matrix H=(HR,HI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  CBAL.  If  CBAL  has not been used, */
/*          set LOW=1 and IGH equal to the order of the matrix, N. */

/*        ORTR and ORTI contain information about the unitary trans- */
/*          formations used in the reduction by  CORTH, if performed. */
/*          Only elements LOW through IGH are used.  If the eigenvectors */
/*          of the Hessenberg matrix are desired, set ORTR(J) and */
/*          ORTI(J) to 0.0E0 for these elements.  ORTR and ORTI are */
/*          one-dimensional REAL arrays, dimensioned ORTR(IGH) and */
/*          ORTI(IGH). */

/*        HR and HI contain the real and imaginary parts, respectively, */
/*          of the complex upper Hessenberg matrix.  Their lower */
/*          triangles below the subdiagonal contain information about */
/*          the unitary transformations used in the reduction by  CORTH, */
/*          if performed.  If the eigenvectors of the Hessenberg matrix */
/*          are desired, these elements may be arbitrary.  HR and HI */
/*          are two-dimensional REAL arrays, dimensioned HR(NM,N) and */
/*          HI(NM,N). */

/*     On OUTPUT */

/*        ORTR, ORTI, and the upper Hessenberg portions of HR and HI */
/*          have been destroyed. */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues of the upper Hessenberg matrix.  If an */
/*          error exit is made, the eigenvalues should be correct for */
/*          indices IERR+1, IERR+2, ..., N.  WR and WI are one- */
/*          dimensional REAL arrays, dimensioned WR(N) and WI(N). */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the eigenvectors.  The eigenvectors are unnormalized. */
/*          If an error exit is made, none of the eigenvectors has been */
/*          found.  ZR and ZI are two-dimensional REAL arrays, */
/*          dimensioned ZR(NM,N) and ZI(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30*N iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     IERR+1, IERR+2, ..., N, but no eigenvectors are */
/*                     computed. */

/*     Calls CSROOT for complex square root. */
/*     Calls PYTHAG(A,B) for sqrt(A**2 + B**2). */
/*     Calls CDIV for complex division. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  CDIV, CSROOT, PYTHAG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COMQR2 */


/* ***FIRST EXECUTABLE STATEMENT  COMQR2 */
    /* Parameter adjustments */
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1;
    zr -= zr_offset;
    hi_dim1 = *nm;
    hi_offset = 1 + hi_dim1;
    hi -= hi_offset;
    hr_dim1 = *nm;
    hr_offset = 1 + hr_dim1;
    hr -= hr_offset;
    --ortr;
    --orti;
    --wr;
    --wi;

    /* Function Body */
    *ierr = 0;
/*     .......... INITIALIZE EIGENVECTOR MATRIX .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    zr[i__ + j * zr_dim1] = 0.f;
	    zi[i__ + j * zi_dim1] = 0.f;
	    if (i__ == j) {
		zr[i__ + j * zr_dim1] = 1.f;
	    }
/* L100: */
	}
    }
/*     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS */
/*                FROM THE INFORMATION LEFT BY CORTH .......... */
    iend = *igh - *low - 1;
    if (iend < 0) {
	goto L180;
    } else if (iend == 0) {
	goto L150;
    } else {
	goto L105;
    }
/*     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- .......... */
L105:
    i__2 = iend;
    for (ii = 1; ii <= i__2; ++ii) {
	i__ = *igh - ii;
	if (ortr[i__] == 0.f && orti[i__] == 0.f) {
	    goto L140;
	}
	if (hr[i__ + (i__ - 1) * hr_dim1] == 0.f && hi[i__ + (i__ - 1) * 
		hi_dim1] == 0.f) {
	    goto L140;
	}
/*     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH .......... */
	norm = hr[i__ + (i__ - 1) * hr_dim1] * ortr[i__] + hi[i__ + (i__ - 1) 
		* hi_dim1] * orti[i__];
	ip1 = i__ + 1;

	i__1 = *igh;
	for (k = ip1; k <= i__1; ++k) {
	    ortr[k] = hr[k + (i__ - 1) * hr_dim1];
	    orti[k] = hi[k + (i__ - 1) * hi_dim1];
/* L110: */
	}

	i__1 = *igh;
	for (j = i__; j <= i__1; ++j) {
	    sr = 0.f;
	    si = 0.f;

	    i__3 = *igh;
	    for (k = i__; k <= i__3; ++k) {
		sr = sr + ortr[k] * zr[k + j * zr_dim1] + orti[k] * zi[k + j *
			 zi_dim1];
		si = si + ortr[k] * zi[k + j * zi_dim1] - orti[k] * zr[k + j *
			 zr_dim1];
/* L115: */
	    }

	    sr /= norm;
	    si /= norm;

	    i__3 = *igh;
	    for (k = i__; k <= i__3; ++k) {
		zr[k + j * zr_dim1] = zr[k + j * zr_dim1] + sr * ortr[k] - si 
			* orti[k];
		zi[k + j * zi_dim1] = zi[k + j * zi_dim1] + sr * orti[k] + si 
			* ortr[k];
/* L120: */
	    }

/* L130: */
	}

L140:
	;
    }
/*     .......... CREATE REAL SUBDIAGONAL ELEMENTS .......... */
L150:
    l = *low + 1;

    i__2 = *igh;
    for (i__ = l; i__ <= i__2; ++i__) {
/* Computing MIN */
	i__1 = i__ + 1;
	ll = min(i__1,*igh);
	if (hi[i__ + (i__ - 1) * hi_dim1] == 0.f) {
	    goto L170;
	}
	norm = pythag_(&hr[i__ + (i__ - 1) * hr_dim1], &hi[i__ + (i__ - 1) * 
		hi_dim1]);
	yr = hr[i__ + (i__ - 1) * hr_dim1] / norm;
	yi = hi[i__ + (i__ - 1) * hi_dim1] / norm;
	hr[i__ + (i__ - 1) * hr_dim1] = norm;
	hi[i__ + (i__ - 1) * hi_dim1] = 0.f;

	i__1 = *n;
	for (j = i__; j <= i__1; ++j) {
	    si = yr * hi[i__ + j * hi_dim1] - yi * hr[i__ + j * hr_dim1];
	    hr[i__ + j * hr_dim1] = yr * hr[i__ + j * hr_dim1] + yi * hi[i__ 
		    + j * hi_dim1];
	    hi[i__ + j * hi_dim1] = si;
/* L155: */
	}

	i__1 = ll;
	for (j = 1; j <= i__1; ++j) {
	    si = yr * hi[j + i__ * hi_dim1] + yi * hr[j + i__ * hr_dim1];
	    hr[j + i__ * hr_dim1] = yr * hr[j + i__ * hr_dim1] - yi * hi[j + 
		    i__ * hi_dim1];
	    hi[j + i__ * hi_dim1] = si;
/* L160: */
	}

	i__1 = *igh;
	for (j = *low; j <= i__1; ++j) {
	    si = yr * zi[j + i__ * zi_dim1] + yi * zr[j + i__ * zr_dim1];
	    zr[j + i__ * zr_dim1] = yr * zr[j + i__ * zr_dim1] - yi * zi[j + 
		    i__ * zi_dim1];
	    zi[j + i__ * zi_dim1] = si;
/* L165: */
	}

L170:
	;
    }
/*     .......... STORE ROOTS ISOLATED BY CBAL .......... */
L180:
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (i__ >= *low && i__ <= *igh) {
	    goto L200;
	}
	wr[i__] = hr[i__ + i__ * hr_dim1];
	wi[i__] = hi[i__ + i__ * hi_dim1];
L200:
	;
    }

    en = *igh;
    tr = 0.f;
    ti = 0.f;
    itn = *n * 30;
/*     .......... SEARCH FOR NEXT EIGENVALUE .......... */
L220:
    if (en < *low) {
	goto L680;
    }
    its = 0;
    enm1 = en - 1;
/*     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT */
/*                FOR L=EN STEP -1 UNTIL LOW DO -- .......... */
L240:
    i__2 = en;
    for (ll = *low; ll <= i__2; ++ll) {
	l = en + *low - ll;
	if (l == *low) {
	    goto L300;
	}
	s1 = (r__1 = hr[l - 1 + (l - 1) * hr_dim1], dabs(r__1)) + (r__2 = hi[
		l - 1 + (l - 1) * hi_dim1], dabs(r__2)) + (r__3 = hr[l + l * 
		hr_dim1], dabs(r__3)) + (r__4 = hi[l + l * hi_dim1], dabs(
		r__4));
	s2 = s1 + (r__1 = hr[l + (l - 1) * hr_dim1], dabs(r__1));
	if (s2 == s1) {
	    goto L300;
	}
/* L260: */
    }
/*     .......... FORM SHIFT .......... */
L300:
    if (l == en) {
	goto L660;
    }
    if (itn == 0) {
	goto L1000;
    }
    if (its == 10 || its == 20) {
	goto L320;
    }
    sr = hr[en + en * hr_dim1];
    si = hi[en + en * hi_dim1];
    xr = hr[enm1 + en * hr_dim1] * hr[en + enm1 * hr_dim1];
    xi = hi[enm1 + en * hi_dim1] * hr[en + enm1 * hr_dim1];
    if (xr == 0.f && xi == 0.f) {
	goto L340;
    }
    yr = (hr[enm1 + enm1 * hr_dim1] - sr) / 2.f;
    yi = (hi[enm1 + enm1 * hi_dim1] - si) / 2.f;
/* Computing 2nd power */
    r__2 = yr;
/* Computing 2nd power */
    r__3 = yi;
    r__1 = r__2 * r__2 - r__3 * r__3 + xr;
    r__4 = yr * 2.f * yi + xi;
    csroot_(&r__1, &r__4, &zzr, &zzi);
    if (yr * zzr + yi * zzi >= 0.f) {
	goto L310;
    }
    zzr = -zzr;
    zzi = -zzi;
L310:
    r__1 = yr + zzr;
    r__2 = yi + zzi;
    cdiv_(&xr, &xi, &r__1, &r__2, &xr, &xi);
    sr -= xr;
    si -= xi;
    goto L340;
/*     .......... FORM EXCEPTIONAL SHIFT .......... */
L320:
    sr = (r__1 = hr[en + enm1 * hr_dim1], dabs(r__1)) + (r__2 = hr[enm1 + (en 
	    - 2) * hr_dim1], dabs(r__2));
    si = 0.f;

L340:
    i__2 = en;
    for (i__ = *low; i__ <= i__2; ++i__) {
	hr[i__ + i__ * hr_dim1] -= sr;
	hi[i__ + i__ * hi_dim1] -= si;
/* L360: */
    }

    tr += sr;
    ti += si;
    ++its;
    --itn;
/*     .......... REDUCE TO TRIANGLE (ROWS) .......... */
    lp1 = l + 1;

    i__2 = en;
    for (i__ = lp1; i__ <= i__2; ++i__) {
	sr = hr[i__ + (i__ - 1) * hr_dim1];
	hr[i__ + (i__ - 1) * hr_dim1] = 0.f;
	r__1 = pythag_(&hr[i__ - 1 + (i__ - 1) * hr_dim1], &hi[i__ - 1 + (i__ 
		- 1) * hi_dim1]);
	norm = pythag_(&r__1, &sr);
	xr = hr[i__ - 1 + (i__ - 1) * hr_dim1] / norm;
	wr[i__ - 1] = xr;
	xi = hi[i__ - 1 + (i__ - 1) * hi_dim1] / norm;
	wi[i__ - 1] = xi;
	hr[i__ - 1 + (i__ - 1) * hr_dim1] = norm;
	hi[i__ - 1 + (i__ - 1) * hi_dim1] = 0.f;
	hi[i__ + (i__ - 1) * hi_dim1] = sr / norm;

	i__1 = *n;
	for (j = i__; j <= i__1; ++j) {
	    yr = hr[i__ - 1 + j * hr_dim1];
	    yi = hi[i__ - 1 + j * hi_dim1];
	    zzr = hr[i__ + j * hr_dim1];
	    zzi = hi[i__ + j * hi_dim1];
	    hr[i__ - 1 + j * hr_dim1] = xr * yr + xi * yi + hi[i__ + (i__ - 1)
		     * hi_dim1] * zzr;
	    hi[i__ - 1 + j * hi_dim1] = xr * yi - xi * yr + hi[i__ + (i__ - 1)
		     * hi_dim1] * zzi;
	    hr[i__ + j * hr_dim1] = xr * zzr - xi * zzi - hi[i__ + (i__ - 1) *
		     hi_dim1] * yr;
	    hi[i__ + j * hi_dim1] = xr * zzi + xi * zzr - hi[i__ + (i__ - 1) *
		     hi_dim1] * yi;
/* L490: */
	}

/* L500: */
    }

    si = hi[en + en * hi_dim1];
    if (si == 0.f) {
	goto L540;
    }
    norm = pythag_(&hr[en + en * hr_dim1], &si);
    sr = hr[en + en * hr_dim1] / norm;
    si /= norm;
    hr[en + en * hr_dim1] = norm;
    hi[en + en * hi_dim1] = 0.f;
    if (en == *n) {
	goto L540;
    }
    ip1 = en + 1;

    i__2 = *n;
    for (j = ip1; j <= i__2; ++j) {
	yr = hr[en + j * hr_dim1];
	yi = hi[en + j * hi_dim1];
	hr[en + j * hr_dim1] = sr * yr + si * yi;
	hi[en + j * hi_dim1] = sr * yi - si * yr;
/* L520: */
    }
/*     .......... INVERSE OPERATION (COLUMNS) .......... */
L540:
    i__2 = en;
    for (j = lp1; j <= i__2; ++j) {
	xr = wr[j - 1];
	xi = wi[j - 1];

	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    yr = hr[i__ + (j - 1) * hr_dim1];
	    yi = 0.f;
	    zzr = hr[i__ + j * hr_dim1];
	    zzi = hi[i__ + j * hi_dim1];
	    if (i__ == j) {
		goto L560;
	    }
	    yi = hi[i__ + (j - 1) * hi_dim1];
	    hi[i__ + (j - 1) * hi_dim1] = xr * yi + xi * yr + hi[j + (j - 1) *
		     hi_dim1] * zzi;
L560:
	    hr[i__ + (j - 1) * hr_dim1] = xr * yr - xi * yi + hi[j + (j - 1) *
		     hi_dim1] * zzr;
	    hr[i__ + j * hr_dim1] = xr * zzr + xi * zzi - hi[j + (j - 1) * 
		    hi_dim1] * yr;
	    hi[i__ + j * hi_dim1] = xr * zzi - xi * zzr - hi[j + (j - 1) * 
		    hi_dim1] * yi;
/* L580: */
	}

	i__1 = *igh;
	for (i__ = *low; i__ <= i__1; ++i__) {
	    yr = zr[i__ + (j - 1) * zr_dim1];
	    yi = zi[i__ + (j - 1) * zi_dim1];
	    zzr = zr[i__ + j * zr_dim1];
	    zzi = zi[i__ + j * zi_dim1];
	    zr[i__ + (j - 1) * zr_dim1] = xr * yr - xi * yi + hi[j + (j - 1) *
		     hi_dim1] * zzr;
	    zi[i__ + (j - 1) * zi_dim1] = xr * yi + xi * yr + hi[j + (j - 1) *
		     hi_dim1] * zzi;
	    zr[i__ + j * zr_dim1] = xr * zzr + xi * zzi - hi[j + (j - 1) * 
		    hi_dim1] * yr;
	    zi[i__ + j * zi_dim1] = xr * zzi - xi * zzr - hi[j + (j - 1) * 
		    hi_dim1] * yi;
/* L590: */
	}

/* L600: */
    }

    if (si == 0.f) {
	goto L240;
    }

    i__2 = en;
    for (i__ = 1; i__ <= i__2; ++i__) {
	yr = hr[i__ + en * hr_dim1];
	yi = hi[i__ + en * hi_dim1];
	hr[i__ + en * hr_dim1] = sr * yr - si * yi;
	hi[i__ + en * hi_dim1] = sr * yi + si * yr;
/* L630: */
    }

    i__2 = *igh;
    for (i__ = *low; i__ <= i__2; ++i__) {
	yr = zr[i__ + en * zr_dim1];
	yi = zi[i__ + en * zi_dim1];
	zr[i__ + en * zr_dim1] = sr * yr - si * yi;
	zi[i__ + en * zi_dim1] = sr * yi + si * yr;
/* L640: */
    }

    goto L240;
/*     .......... A ROOT FOUND .......... */
L660:
    hr[en + en * hr_dim1] += tr;
    wr[en] = hr[en + en * hr_dim1];
    hi[en + en * hi_dim1] += ti;
    wi[en] = hi[en + en * hi_dim1];
    en = enm1;
    goto L220;
/*     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND */
/*                VECTORS OF UPPER TRIANGULAR FORM .......... */
L680:
    norm = 0.f;

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {

	i__1 = *n;
	for (j = i__; j <= i__1; ++j) {
	    norm = norm + (r__1 = hr[i__ + j * hr_dim1], dabs(r__1)) + (r__2 =
		     hi[i__ + j * hi_dim1], dabs(r__2));
/* L720: */
	}
    }

    if (*n == 1 || norm == 0.f) {
	goto L1001;
    }
/*     .......... FOR EN=N STEP -1 UNTIL 2 DO -- .......... */
    i__1 = *n;
    for (nn = 2; nn <= i__1; ++nn) {
	en = *n + 2 - nn;
	xr = wr[en];
	xi = wi[en];
	enm1 = en - 1;
/*     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- .......... */
	i__2 = enm1;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = en - ii;
	    zzr = hr[i__ + en * hr_dim1];
	    zzi = hi[i__ + en * hi_dim1];
	    if (i__ == enm1) {
		goto L760;
	    }
	    ip1 = i__ + 1;

	    i__3 = enm1;
	    for (j = ip1; j <= i__3; ++j) {
		zzr = zzr + hr[i__ + j * hr_dim1] * hr[j + en * hr_dim1] - hi[
			i__ + j * hi_dim1] * hi[j + en * hi_dim1];
		zzi = zzi + hr[i__ + j * hr_dim1] * hi[j + en * hi_dim1] + hi[
			i__ + j * hi_dim1] * hr[j + en * hr_dim1];
/* L740: */
	    }

L760:
	    yr = xr - wr[i__];
	    yi = xi - wi[i__];
	    if (yr != 0.f || yi != 0.f) {
		goto L775;
	    }
	    yr = norm;
L770:
	    yr *= .5f;
	    if (norm + yr > norm) {
		goto L770;
	    }
	    yr *= 2.f;
L775:
	    cdiv_(&zzr, &zzi, &yr, &yi, &hr[i__ + en * hr_dim1], &hi[i__ + en 
		    * hi_dim1]);
/* L780: */
	}

/* L800: */
    }
/*     .......... END BACKSUBSTITUTION .......... */
    enm1 = *n - 1;
/*     .......... VECTORS OF ISOLATED ROOTS .......... */
    i__1 = enm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ >= *low && i__ <= *igh) {
	    goto L840;
	}
	ip1 = i__ + 1;

	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    zr[i__ + j * zr_dim1] = hr[i__ + j * hr_dim1];
	    zi[i__ + j * zi_dim1] = hi[i__ + j * hi_dim1];
/* L820: */
	}

L840:
	;
    }
/*     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE */
/*                VECTORS OF ORIGINAL FULL MATRIX. */
/*                FOR J=N STEP -1 UNTIL LOW+1 DO -- .......... */
    i__1 = enm1;
    for (jj = *low; jj <= i__1; ++jj) {
	j = *n + *low - jj;
/* Computing MIN */
	i__2 = j - 1;
	m = min(i__2,*igh);

	i__2 = *igh;
	for (i__ = *low; i__ <= i__2; ++i__) {
	    zzr = zr[i__ + j * zr_dim1];
	    zzi = zi[i__ + j * zi_dim1];

	    i__3 = m;
	    for (k = *low; k <= i__3; ++k) {
		zzr = zzr + zr[i__ + k * zr_dim1] * hr[k + j * hr_dim1] - zi[
			i__ + k * zi_dim1] * hi[k + j * hi_dim1];
		zzi = zzi + zr[i__ + k * zr_dim1] * hi[k + j * hi_dim1] + zi[
			i__ + k * zi_dim1] * hr[k + j * hr_dim1];
/* L860: */
	    }

	    zr[i__ + j * zr_dim1] = zzr;
	    zi[i__ + j * zi_dim1] = zzi;
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
} /* comqr2_ */


/* comlr.f -- translated by f2c (version 12.02.01).
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

/* DECK COMLR */
/* Subroutine */ int comlr_(integer *nm, integer *n, integer *low, integer *
	igh, real *hr, real *hi, real *wr, real *wi, integer *ierr)
{
    /* System generated locals */
    integer hr_dim1, hr_offset, hi_dim1, hi_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, l, m;
    static real s1, s2;
    static integer en, ll, mm;
    static real si, ti, xi, yi, sr, tr, xr, yr;
    static integer im1, mp1, itn, its;
    static real zzi, zzr;
    static integer enm1;
    extern /* Subroutine */ int cdiv_(real *, real *, real *, real *, real *, 
	    real *), csroot_(real *, real *, real *, real *);

/* ***BEGIN PROLOGUE  COMLR */
/* ***PURPOSE  Compute the eigenvalues of a complex upper Hessenberg */
/*            matrix using the modified LR method. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2B */
/* ***TYPE      COMPLEX (COMLR-C) */
/* ***KEYWORDS  EIGENVALUES, EISPACK, LR METHOD */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure COMLR, */
/*     NUM. MATH. 12, 369-376(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971). */

/*     This subroutine finds the eigenvalues of a COMPLEX */
/*     UPPER Hessenberg matrix by the modified LR method. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, HR and HI, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix H=(HR,HI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  CBAL.  If  CBAL  has not been used, */
/*          set LOW=1 and IGH equal to the order of the matrix, N. */

/*        HR and HI contain the real and imaginary parts, respectively, */
/*          of the complex upper Hessenberg matrix.  Their lower */
/*          triangles below the subdiagonal contain the multipliers */
/*          which were used in the reduction by  COMHES, if performed. */
/*          HR and HI are two-dimensional REAL arrays, dimensioned */
/*          HR(NM,N) and HI(NM,N). */

/*     On OUTPUT */

/*        The upper Hessenberg portions of HR and HI have been */
/*          destroyed.  Therefore, they must be saved before calling */
/*          COMLR  if subsequent calculation of eigenvectors is to */
/*          be performed. */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues of the upper Hessenberg matrix.  If an */
/*          error exit is made, the eigenvalues should be correct for */
/*          indices IERR+1, IERR+2, ..., N.  WR and WI are one- */
/*          dimensional REAL arrays, dimensioned WR(N) and WI(N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30*N iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     IERR+1, IERR+2, ..., N. */

/*     Calls CSROOT for complex square root. */
/*     Calls CDIV for complex division. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  CDIV, CSROOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COMLR */


/* ***FIRST EXECUTABLE STATEMENT  COMLR */
    /* Parameter adjustments */
    hi_dim1 = *nm;
    hi_offset = 1 + hi_dim1;
    hi -= hi_offset;
    hr_dim1 = *nm;
    hr_offset = 1 + hr_dim1;
    hr -= hr_offset;
    --wr;
    --wi;

    /* Function Body */
    *ierr = 0;
/*     .......... STORE ROOTS ISOLATED BY CBAL .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
	goto L1001;
    }
    its = 0;
    enm1 = en - 1;
/*     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT */
/*                FOR L=EN STEP -1 UNTIL LOW E0 -- .......... */
L240:
    i__1 = en;
    for (ll = *low; ll <= i__1; ++ll) {
	l = en + *low - ll;
	if (l == *low) {
	    goto L300;
	}
	s1 = (r__1 = hr[l - 1 + (l - 1) * hr_dim1], dabs(r__1)) + (r__2 = hi[
		l - 1 + (l - 1) * hi_dim1], dabs(r__2)) + (r__3 = hr[l + l * 
		hr_dim1], dabs(r__3)) + (r__4 = hi[l + l * hi_dim1], dabs(
		r__4));
	s2 = s1 + (r__1 = hr[l + (l - 1) * hr_dim1], dabs(r__1)) + (r__2 = hi[
		l + (l - 1) * hi_dim1], dabs(r__2));
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
    xr = hr[enm1 + en * hr_dim1] * hr[en + enm1 * hr_dim1] - hi[enm1 + en * 
	    hi_dim1] * hi[en + enm1 * hi_dim1];
    xi = hr[enm1 + en * hr_dim1] * hi[en + enm1 * hi_dim1] + hi[enm1 + en * 
	    hi_dim1] * hr[en + enm1 * hr_dim1];
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
    si = (r__1 = hi[en + enm1 * hi_dim1], dabs(r__1)) + (r__2 = hi[enm1 + (en 
	    - 2) * hi_dim1], dabs(r__2));

L340:
    i__1 = en;
    for (i__ = *low; i__ <= i__1; ++i__) {
	hr[i__ + i__ * hr_dim1] -= sr;
	hi[i__ + i__ * hi_dim1] -= si;
/* L360: */
    }

    tr += sr;
    ti += si;
    ++its;
    --itn;
/*     .......... LOOK FOR TWO CONSECUTIVE SMALL */
/*                SUB-DIAGONAL ELEMENTS .......... */
    xr = (r__1 = hr[enm1 + enm1 * hr_dim1], dabs(r__1)) + (r__2 = hi[enm1 + 
	    enm1 * hi_dim1], dabs(r__2));
    yr = (r__1 = hr[en + enm1 * hr_dim1], dabs(r__1)) + (r__2 = hi[en + enm1 *
	     hi_dim1], dabs(r__2));
    zzr = (r__1 = hr[en + en * hr_dim1], dabs(r__1)) + (r__2 = hi[en + en * 
	    hi_dim1], dabs(r__2));
/*     .......... FOR M=EN-1 STEP -1 UNTIL L DO -- .......... */
    i__1 = enm1;
    for (mm = l; mm <= i__1; ++mm) {
	m = enm1 + l - mm;
	if (m == l) {
	    goto L420;
	}
	yi = yr;
	yr = (r__1 = hr[m + (m - 1) * hr_dim1], dabs(r__1)) + (r__2 = hi[m + (
		m - 1) * hi_dim1], dabs(r__2));
	xi = zzr;
	zzr = xr;
	xr = (r__1 = hr[m - 1 + (m - 1) * hr_dim1], dabs(r__1)) + (r__2 = hi[
		m - 1 + (m - 1) * hi_dim1], dabs(r__2));
	s1 = zzr / yi * (zzr + xr + xi);
	s2 = s1 + yr;
	if (s2 == s1) {
	    goto L420;
	}
/* L380: */
    }
/*     .......... TRIANGULAR DECOMPOSITION H=L*R .......... */
L420:
    mp1 = m + 1;

    i__1 = en;
    for (i__ = mp1; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	xr = hr[im1 + im1 * hr_dim1];
	xi = hi[im1 + im1 * hi_dim1];
	yr = hr[i__ + im1 * hr_dim1];
	yi = hi[i__ + im1 * hi_dim1];
	if (dabs(xr) + dabs(xi) >= dabs(yr) + dabs(yi)) {
	    goto L460;
	}
/*     .......... INTERCHANGE ROWS OF HR AND HI .......... */
	i__2 = en;
	for (j = im1; j <= i__2; ++j) {
	    zzr = hr[im1 + j * hr_dim1];
	    hr[im1 + j * hr_dim1] = hr[i__ + j * hr_dim1];
	    hr[i__ + j * hr_dim1] = zzr;
	    zzi = hi[im1 + j * hi_dim1];
	    hi[im1 + j * hi_dim1] = hi[i__ + j * hi_dim1];
	    hi[i__ + j * hi_dim1] = zzi;
/* L440: */
	}

	cdiv_(&xr, &xi, &yr, &yi, &zzr, &zzi);
	wr[i__] = 1.f;
	goto L480;
L460:
	cdiv_(&yr, &yi, &xr, &xi, &zzr, &zzi);
	wr[i__] = -1.f;
L480:
	hr[i__ + im1 * hr_dim1] = zzr;
	hi[i__ + im1 * hi_dim1] = zzi;

	i__2 = en;
	for (j = i__; j <= i__2; ++j) {
	    hr[i__ + j * hr_dim1] = hr[i__ + j * hr_dim1] - zzr * hr[im1 + j *
		     hr_dim1] + zzi * hi[im1 + j * hi_dim1];
	    hi[i__ + j * hi_dim1] = hi[i__ + j * hi_dim1] - zzr * hi[im1 + j *
		     hi_dim1] - zzi * hr[im1 + j * hr_dim1];
/* L500: */
	}

/* L520: */
    }
/*     .......... COMPOSITION R*L=H .......... */
    i__1 = en;
    for (j = mp1; j <= i__1; ++j) {
	xr = hr[j + (j - 1) * hr_dim1];
	xi = hi[j + (j - 1) * hi_dim1];
	hr[j + (j - 1) * hr_dim1] = 0.f;
	hi[j + (j - 1) * hi_dim1] = 0.f;
/*     .......... INTERCHANGE COLUMNS OF HR AND HI, */
/*                IF NECESSARY .......... */
	if (wr[j] <= 0.f) {
	    goto L580;
	}

	i__2 = j;
	for (i__ = l; i__ <= i__2; ++i__) {
	    zzr = hr[i__ + (j - 1) * hr_dim1];
	    hr[i__ + (j - 1) * hr_dim1] = hr[i__ + j * hr_dim1];
	    hr[i__ + j * hr_dim1] = zzr;
	    zzi = hi[i__ + (j - 1) * hi_dim1];
	    hi[i__ + (j - 1) * hi_dim1] = hi[i__ + j * hi_dim1];
	    hi[i__ + j * hi_dim1] = zzi;
/* L540: */
	}

L580:
	i__2 = j;
	for (i__ = l; i__ <= i__2; ++i__) {
	    hr[i__ + (j - 1) * hr_dim1] = hr[i__ + (j - 1) * hr_dim1] + xr * 
		    hr[i__ + j * hr_dim1] - xi * hi[i__ + j * hi_dim1];
	    hi[i__ + (j - 1) * hi_dim1] = hi[i__ + (j - 1) * hi_dim1] + xr * 
		    hi[i__ + j * hi_dim1] + xi * hr[i__ + j * hr_dim1];
/* L600: */
	}

/* L640: */
    }

    goto L240;
/*     .......... A ROOT FOUND .......... */
L660:
    wr[en] = hr[en + en * hr_dim1] + tr;
    wi[en] = hi[en + en * hi_dim1] + ti;
    en = enm1;
    goto L220;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30*N ITERATIONS .......... */
L1000:
    *ierr = en;
L1001:
    return 0;
} /* comlr_ */


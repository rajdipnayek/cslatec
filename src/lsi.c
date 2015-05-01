/* lsi.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__1 = 1;
static real c_b7 = 0.f;
static integer c__0 = 0;
static integer c__2 = 2;

/* DECK LSI */
/* Subroutine */ int lsi_(real *w, integer *mdw, integer *ma, integer *mg, 
	integer *n, real *prgopt, real *x, real *rnorm, integer *mode, real *
	ws, integer *ip)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, m, n1, n2, n3;
    extern /* Subroutine */ int h12_(integer *, integer *, integer *, integer 
	    *, real *, integer *, real *, real *, integer *, integer *, 
	    integer *);
    static real rb;
    static integer np1;
    static real fac, gam, tau;
    static integer key;
    static logical cov;
    static real tol;
    static integer map1, krm1, krp1;
    extern /* Subroutine */ int hfti_(real *, integer *, integer *, integer *,
	     real *, integer *, integer *, real *, integer *, real *, real *, 
	    real *, integer *);
    static integer link;
    extern /* Subroutine */ int lpdp_(real *, integer *, integer *, integer *,
	     integer *, real *, real *, real *, integer *, real *, integer *);
    static integer last;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer next;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer krank;
    static real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    );
    static real xnorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);
    extern doublereal r1mach_(integer *);
    static integer minman, mdlpdp;
    static logical sclcov;
    static real srelpr;

/* ***BEGIN PROLOGUE  LSI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to LSEI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (LSI-S, DLSI-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to LSEI.  The documentation for */
/*     LSEI has complete usage instructions. */

/*     Solve.. */
/*              AX = B,  A  MA by N  (least squares equations) */
/*     subject to.. */

/*              GX.GE.H, G  MG by N  (inequality constraints) */

/*     Input.. */

/*      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1. */
/*                       (G H) */

/*     MDW,MA,MG,N */
/*              contain (resp) var. dimension of W(*,*), */
/*              and matrix dimensions. */

/*     PRGOPT(*), */
/*              Program option vector. */

/*     OUTPUT.. */

/*      X(*),RNORM */

/*              Solution vector(unless MODE=2), length of AX-B. */

/*      MODE */
/*              =0   Inequality constraints are compatible. */
/*              =2   Inequality constraints contradictory. */

/*      WS(*), */
/*              Working storage of dimension K+N+(MG+2)*(N+7), */
/*              where K=MAX(MA+MG,N). */
/*      IP(MG+2*N+1) */
/*              Integer working storage */

/* ***ROUTINES CALLED  H12, HFTI, LPDP, R1MACH, SASUM, SAXPY, SCOPY, SDOT, */
/*                    SSCAL, SSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   920422  Changed CALL to HFTI to include variable MA.  (WRB) */
/* ***END PROLOGUE  LSI */



    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --ws;
    --ip;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  LSI */

/*     Set the nominal tolerance used in the code. */

    if (first) {
	srelpr = r1mach_(&c__4);
    }
    first = FALSE_;
    tol = sqrt(srelpr);

    *mode = 0;
    *rnorm = 0.f;
    m = *ma + *mg;
    np1 = *n + 1;
    krank = 0;
    if (*n <= 0 || m <= 0) {
	goto L370;
    }

/*     To process option vector. */

    cov = FALSE_;
    sclcov = TRUE_;
    last = 1;
    link = prgopt[1];

L100:
    if (link > 1) {
	key = prgopt[last + 1];
	if (key == 1) {
	    cov = prgopt[last + 2] != 0.f;
	}
	if (key == 10) {
	    sclcov = prgopt[last + 2] == 0.f;
	}
	if (key == 5) {
/* Computing MAX */
	    r__1 = srelpr, r__2 = prgopt[last + 2];
	    tol = dmax(r__1,r__2);
	}
	next = prgopt[link];
	last = link;
	link = next;
	goto L100;
    }

/*     Compute matrix norm of least squares equations. */

    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = anorm, r__2 = sasum_(ma, &w[j * w_dim1 + 1], &c__1);
	anorm = dmax(r__1,r__2);
/* L110: */
    }

/*     Set tolerance for HFTI( ) rank test. */

    tau = tol * anorm;

/*     Compute Householder orthogonal decomposition of matrix. */

    scopy_(n, &c_b7, &c__0, &ws[1], &c__1);
    scopy_(ma, &w[np1 * w_dim1 + 1], &c__1, &ws[1], &c__1);
    k = max(m,*n);
    minman = min(*ma,*n);
    n1 = k + 1;
    n2 = n1 + *n;
    hfti_(&w[w_offset], mdw, ma, n, &ws[1], ma, &c__1, &tau, &krank, rnorm, &
	    ws[n2], &ws[n1], &ip[1]);
    fac = 1.f;
    gam = (real) (*ma - krank);
    if (krank < *ma && sclcov) {
/* Computing 2nd power */
	r__1 = *rnorm;
	fac = r__1 * r__1 / gam;
    }

/*     Reduce to LPDP and solve. */

    map1 = *ma + 1;

/*     Compute inequality rt-hand side for LPDP. */

    if (*ma < m) {
	if (minman > 0) {
	    i__1 = m;
	    for (i__ = map1; i__ <= i__1; ++i__) {
		w[i__ + np1 * w_dim1] -= sdot_(n, &w[i__ + w_dim1], mdw, &ws[
			1], &c__1);
/* L120: */
	    }

/*           Apply permutations to col. of inequality constraint matrix. */

	    i__1 = minman;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sswap_(mg, &w[map1 + i__ * w_dim1], &c__1, &w[map1 + ip[i__] *
			 w_dim1], &c__1);
/* L130: */
	    }

/*           Apply Householder transformations to constraint matrix. */

	    if (krank > 0 && krank < *n) {
		for (i__ = krank; i__ >= 1; --i__) {
		    i__1 = krank + 1;
		    h12_(&c__2, &i__, &i__1, n, &w[i__ + w_dim1], mdw, &ws[n1 
			    + i__ - 1], &w[map1 + w_dim1], mdw, &c__1, mg);
/* L140: */
		}
	    }

/*           Compute permuted inequality constraint matrix times r-inv. */

	    i__1 = m;
	    for (i__ = map1; i__ <= i__1; ++i__) {
		i__2 = krank;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = j - 1;
		    w[i__ + j * w_dim1] = (w[i__ + j * w_dim1] - sdot_(&i__3, 
			    &w[j * w_dim1 + 1], &c__1, &w[i__ + w_dim1], mdw))
			     / w[j + j * w_dim1];
/* L150: */
		}
/* L160: */
	    }
	}

/*        Solve the reduced problem with LPDP algorithm, */
/*        the least projected distance problem. */

	i__1 = *n - krank;
	lpdp_(&w[map1 + w_dim1], mdw, mg, &krank, &i__1, &prgopt[1], &x[1], &
		xnorm, &mdlpdp, &ws[n2], &ip[*n + 1]);

/*        Compute solution in original coordinates. */

	if (mdlpdp == 1) {
	    for (i__ = krank; i__ >= 1; --i__) {
		i__1 = krank - i__;
		x[i__] = (x[i__] - sdot_(&i__1, &w[i__ + (i__ + 1) * w_dim1], 
			mdw, &x[i__ + 1], &c__1)) / w[i__ + i__ * w_dim1];
/* L170: */
	    }

/*           Apply Householder transformation to solution vector. */

	    if (krank < *n) {
		i__1 = krank;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = krank + 1;
		    h12_(&c__2, &i__, &i__2, n, &w[i__ + w_dim1], mdw, &ws[n1 
			    + i__ - 1], &x[1], &c__1, &c__1, &c__1);
/* L180: */
		}
	    }

/*           Repermute variables to their input order. */

	    if (minman > 0) {
		for (i__ = minman; i__ >= 1; --i__) {
		    sswap_(&c__1, &x[i__], &c__1, &x[ip[i__]], &c__1);
/* L190: */
		}

/*              Variables are now in original coordinates. */
/*              Add solution of unconstrained problem. */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    x[i__] += ws[i__];
/* L200: */
		}

/*              Compute the residual vector norm. */

/* Computing 2nd power */
		r__1 = *rnorm;
/* Computing 2nd power */
		r__2 = xnorm;
		*rnorm = sqrt(r__1 * r__1 + r__2 * r__2);
	    }
	} else {
	    *mode = 2;
	}
    } else {
	scopy_(n, &ws[1], &c__1, &x[1], &c__1);
    }

/*     Compute covariance matrix based on the orthogonal decomposition */
/*     from HFTI( ). */

    if (! cov || krank <= 0) {
	goto L370;
    }
    krm1 = krank - 1;
    krp1 = krank + 1;

/*     Copy diagonal terms to working array. */

    i__1 = *mdw + 1;
    scopy_(&krank, &w[w_offset], &i__1, &ws[n2], &c__1);

/*     Reciprocate diagonal terms. */

    i__1 = krank;
    for (j = 1; j <= i__1; ++j) {
	w[j + j * w_dim1] = 1.f / w[j + j * w_dim1];
/* L210: */
    }

/*     Invert the upper triangular QR factor on itself. */

    if (krank > 1) {
	i__1 = krm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = krank;
	    for (j = i__ + 1; j <= i__2; ++j) {
		i__3 = j - i__;
		w[i__ + j * w_dim1] = -sdot_(&i__3, &w[i__ + i__ * w_dim1], 
			mdw, &w[i__ + j * w_dim1], &c__1) * w[j + j * w_dim1];
/* L220: */
	    }
/* L230: */
	}
    }

/*     Compute the inverted factor times its transpose. */

    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = krank;
	for (j = i__; j <= i__2; ++j) {
	    i__3 = krank + 1 - j;
	    w[i__ + j * w_dim1] = sdot_(&i__3, &w[i__ + j * w_dim1], mdw, &w[
		    j + j * w_dim1], mdw);
/* L240: */
	}
/* L250: */
    }

/*     Zero out lower trapezoidal part. */
/*     Copy upper triangular to lower triangular part. */

    if (krank < *n) {
	i__1 = krank;
	for (j = 1; j <= i__1; ++j) {
	    scopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[j + w_dim1], mdw);
/* L260: */
	}

	i__1 = *n;
	for (i__ = krp1; i__ <= i__1; ++i__) {
	    scopy_(&i__, &c_b7, &c__0, &w[i__ + w_dim1], mdw);
/* L270: */
	}

/*        Apply right side transformations to lower triangle. */

	n3 = n2 + krp1;
	i__1 = krank;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    l = n1 + i__;
	    k = n2 + i__;
	    rb = ws[l - 1] * ws[k - 1];

/*           If RB.GE.0.E0, transformation can be regarded as zero. */

	    if (rb < 0.f) {
		rb = 1.f / rb;

/*              Store unscaled rank one Householder update in work array. */

		scopy_(n, &c_b7, &c__0, &ws[n3], &c__1);
		l = n1 + i__;
		k = n3 + i__;
		ws[k - 1] = ws[l - 1];

		i__2 = *n;
		for (j = krp1; j <= i__2; ++j) {
		    ws[n3 + j - 1] = w[i__ + j * w_dim1];
/* L280: */
		}

		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = j - i__;
		    i__4 = *n - j + 1;
		    ws[j] = rb * (sdot_(&i__3, &w[j + i__ * w_dim1], mdw, &ws[
			    n3 + i__ - 1], &c__1) + sdot_(&i__4, &w[j + j * 
			    w_dim1], &c__1, &ws[n3 + j - 1], &c__1));
/* L290: */
		}

		l = n3 + i__;
		i__2 = *n - i__ + 1;
		gam = rb * .5f * sdot_(&i__2, &ws[l - 1], &c__1, &ws[i__], &
			c__1);
		i__2 = *n - i__ + 1;
		saxpy_(&i__2, &gam, &ws[l - 1], &c__1, &ws[i__], &c__1);
		i__2 = *n;
		for (j = i__; j <= i__2; ++j) {
		    i__3 = i__ - 1;
		    for (l = 1; l <= i__3; ++l) {
			w[j + l * w_dim1] += ws[n3 + j - 1] * ws[l];
/* L300: */
		    }

		    i__3 = j;
		    for (l = i__; l <= i__3; ++l) {
			w[j + l * w_dim1] = w[j + l * w_dim1] + ws[j] * ws[n3 
				+ l - 1] + ws[l] * ws[n3 + j - 1];
/* L310: */
		    }
/* L320: */
		}
	    }
/* L330: */
	}

/*        Copy lower triangle to upper triangle to symmetrize the */
/*        covariance matrix. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scopy_(&i__, &w[i__ + w_dim1], mdw, &w[i__ * w_dim1 + 1], &c__1);
/* L340: */
	}
    }

/*     Repermute rows and columns. */

    for (i__ = minman; i__ >= 1; --i__) {
	k = ip[i__];
	if (i__ != k) {
	    sswap_(&c__1, &w[i__ + i__ * w_dim1], &c__1, &w[k + k * w_dim1], &
		    c__1);
	    i__1 = i__ - 1;
	    sswap_(&i__1, &w[i__ * w_dim1 + 1], &c__1, &w[k * w_dim1 + 1], &
		    c__1);
	    i__1 = k - i__ - 1;
	    sswap_(&i__1, &w[i__ + (i__ + 1) * w_dim1], mdw, &w[i__ + 1 + k * 
		    w_dim1], &c__1);
	    i__1 = *n - k;
	    sswap_(&i__1, &w[i__ + (k + 1) * w_dim1], mdw, &w[k + (k + 1) * 
		    w_dim1], mdw);
	}
/* L350: */
    }

/*     Put in normalized residual sum of squares scale factor */
/*     and symmetrize the resulting covariance matrix. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sscal_(&j, &fac, &w[j * w_dim1 + 1], &c__1);
	scopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[j + w_dim1], mdw);
/* L360: */
    }

L370:
    ip[1] = krank;
    ip[2] = *n + max(m,*n) + (*mg + 2) * (*n + 7);
    return 0;
} /* lsi_ */


/* dwnlsm.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b3 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;
static doublereal c_b34 = 0.;

/* DECK DWNLSM */
/* Subroutine */ int dwnlsm_(doublereal *w, integer *mdw, integer *mme, 
	integer *ma, integer *n, integer *l, doublereal *prgopt, doublereal *
	x, doublereal *rnorm, integer *mode, integer *ipivot, integer *itype, 
	doublereal *wd, doublereal *h__, doublereal *scale, doublereal *z__, 
	doublereal *temp, doublereal *d__)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, m;
    static doublereal t;
    static integer l1;
    static doublereal z2;
    static integer me, jp;
    static doublereal sm, zz, fac;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static integer key;
    static doublereal tau;
    static integer niv;
    static logical pos, done;
    static doublereal amax, dope[3];
    static integer jcon, link, imax;
    static doublereal alsq;
    static integer iter, last, isol;
    static doublereal wmax;
    static integer next, nopt;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer idope[3], krank;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer nlink;
    static doublereal bnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer itemp, itmax;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), drotm_(integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *);
    static integer iwmax, nsoln;
    extern doublereal d1mach_(integer *);
    static doublereal alamda;
    static logical feasbl;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal eanorm, sparam[5];
    static logical hitcon;
    static doublereal drelpr;
    extern /* Subroutine */ int drotmg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer ntimes;
    extern /* Subroutine */ int dwnlit_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, logical *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static doublereal blowup;

/* ***BEGIN PROLOGUE  DWNLSM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DWNNLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLSM-S, DWNLSM-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DWNNLS. */
/*     The documentation for DWNNLS has complete usage instructions. */

/*     In addition to the parameters discussed in the prologue to */
/*     subroutine DWNNLS, the following work arrays are used in */
/*     subroutine DWNLSM  (they are passed through the calling */
/*     sequence from DWNNLS for purposes of variable dimensioning). */
/*     Their contents will in general be of no interest to the user. */

/*     Variables of type REAL are DOUBLE PRECISION. */

/*         IPIVOT(*) */
/*            An array of length N.  Upon completion it contains the */
/*         pivoting information for the cols of W(*,*). */

/*         ITYPE(*) */
/*            An array of length M which is used to keep track */
/*         of the classification of the equations.  ITYPE(I)=0 */
/*         denotes equation I as an equality constraint. */
/*         ITYPE(I)=1 denotes equation I as a least squares */
/*         equation. */

/*         WD(*) */
/*            An array of length N.  Upon completion it contains the */
/*         dual solution vector. */

/*         H(*) */
/*            An array of length N.  Upon completion it contains the */
/*         pivot scalars of the Householder transformations performed */
/*         in the case KRANK.LT.L. */

/*         SCALE(*) */
/*            An array of length M which is used by the subroutine */
/*         to store the diagonal matrix of weights. */
/*         These are used to apply the modified Givens */
/*         transformations. */

/*         Z(*),TEMP(*) */
/*            Working arrays of length N. */

/*         D(*) */
/*            An array of length N that contains the */
/*         column scaling for the matrix (E). */
/*                                       (A) */

/* ***SEE ALSO  DWNNLS */
/* ***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM, */
/*                    DROTMG, DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and revised.  (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Fixed an error message.  (RWC) */
/*   900604  DP version created from SP version.  (RWC) */
/*   900911  Restriction on value of ALAMDA included.  (WRB) */
/* ***END PROLOGUE  DWNLSM */



    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --ipivot;
    --itype;
    --wd;
    --h__;
    --scale;
    --z__;
    --temp;
    --d__;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DWNLSM */

/*     Initialize variables. */
/*     DRELPR is the precision for the particular machine */
/*     being used.  This logic avoids resetting it every entry. */

    if (first) {
	drelpr = d1mach_(&c__4);
    }
    first = FALSE_;

/*     Set the nominal tolerance used in the code. */

    tau = sqrt(drelpr);

    m = *ma + *mme;
    me = *mme;
    *mode = 2;

/*     To process option vector */

    fac = 1e-4;

/*     Set the nominal blow up factor used in the code. */

    blowup = tau;

/*     The nominal column scaling used in the code is */
/*     the identity scaling. */

    dcopy_(n, &c_b3, &c__0, &d__[1], &c__1);

/*     Define bound for number of options to change. */

    nopt = 1000;

/*     Define bound for positive value of LINK. */

    nlink = 100000;
    ntimes = 0;
    last = 1;
    link = (integer) prgopt[1];
    if (link <= 0 || link > nlink) {
	xermsg_("SLATEC", "DWNLSM", "IN DWNNLS, THE OPTION VECTOR IS UNDEFIN"
		"ED", &c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)41);
	return 0;
    }

L100:
    if (link > 1) {
	++ntimes;
	if (ntimes > nopt) {
	    xermsg_("SLATEC", "DWNLSM", "IN DWNNLS, THE LINKS IN THE OPTION "
		    "VECTOR ARE CYCLING.", &c__3, &c__1, (ftnlen)6, (ftnlen)6, 
		    (ftnlen)54);
	    return 0;
	}

	key = (integer) prgopt[last + 1];
	if (key == 6 && prgopt[last + 2] != 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		t = dnrm2_(&m, &w[j * w_dim1 + 1], &c__1);
		if (t != 0.) {
		    t = 1. / t;
		}
		d__[j] = t;
/* L110: */
	    }
	}

	if (key == 7) {
	    dcopy_(n, &prgopt[last + 2], &c__1, &d__[1], &c__1);
	}
	if (key == 8) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    tau = max(d__1,d__2);
	}
	if (key == 9) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    blowup = max(d__1,d__2);
	}

	next = (integer) prgopt[link];
	if (next <= 0 || next > nlink) {
	    xermsg_("SLATEC", "DWNLSM", "IN DWNNLS, THE OPTION VECTOR IS UND"
		    "EFINED", &c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)41);
	    return 0;
	}

	last = link;
	link = next;
	goto L100;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dscal_(&m, &d__[j], &w[j * w_dim1 + 1], &c__1);
/* L120: */
    }

/*     Process option vector */

    done = FALSE_;
    iter = 0;
    itmax = (*n - *l) * 3;
    *mode = 0;
    nsoln = *l;
    l1 = min(m,*l);

/*     Compute scale factor to apply to equality constraint equations. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wd[j] = dasum_(&m, &w[j * w_dim1 + 1], &c__1);
/* L130: */
    }

    imax = idamax_(n, &wd[1], &c__1);
    eanorm = wd[imax];
    bnorm = dasum_(&m, &w[(*n + 1) * w_dim1 + 1], &c__1);
    alamda = eanorm / (drelpr * fac);

/*     On machines, such as the VAXes using D floating, with a very */
/*     limited exponent range for double precision values, the previously */
/*     computed value of ALAMDA may cause an overflow condition. */
/*     Therefore, this code further limits the value of ALAMDA. */

/* Computing MIN */
    d__1 = alamda, d__2 = sqrt(d1mach_(&c__2));
    alamda = min(d__1,d__2);

/*     Define scaling diagonal matrix for modified Givens usage and */
/*     classify equation types. */

/* Computing 2nd power */
    d__1 = alamda;
    alsq = d__1 * d__1;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        When equation I is heavily weighted ITYPE(I)=0, */
/*        else ITYPE(I)=1. */

	if (i__ <= me) {
	    t = alsq;
	    itemp = 0;
	} else {
	    t = 1.;
	    itemp = 1;
	}
	scale[i__] = t;
	itype[i__] = itemp;
/* L140: */
    }

/*     Set the solution vector X(*) to zero and the column interchange */
/*     matrix to the identity. */

    dcopy_(n, &c_b34, &c__0, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipivot[i__] = i__;
/* L150: */
    }

/*     Perform initial triangularization in the submatrix */
/*     corresponding to the unconstrained variables. */
/*     Set first L components of dual vector to zero because */
/*     these correspond to the unconstrained variables. */

    dcopy_(l, &c_b34, &c__0, &wd[1], &c__1);

/*     The arrays IDOPE(*) and DOPE(*) are used to pass */
/*     information to DWNLIT().  This was done to avoid */
/*     a long calling sequence or the use of COMMON. */

    idope[0] = me;
    idope[1] = nsoln;
    idope[2] = l1;

    dope[0] = alsq;
    dope[1] = eanorm;
    dope[2] = tau;
    dwnlit_(&w[w_offset], mdw, &m, n, l, &ipivot[1], &itype[1], &h__[1], &
	    scale[1], rnorm, idope, dope, &done);
    me = idope[0];
    krank = idope[1];
    niv = idope[2];

/*     Perform WNNLS algorithm using the following steps. */

/*     Until(DONE) */
/*        compute search direction and feasible point */
/*        when (HITCON) add constraints */
/*        else perform multiplier test and drop a constraint */
/*        fin */
/*     Compute-Final-Solution */

/*     To compute search direction and feasible point, */
/*     solve the triangular system of currently non-active */
/*     variables and store the solution in Z(*). */

/*     To solve system */
/*     Copy right hand side into TEMP vector to use overwriting method. */

L160:
    if (done) {
	goto L330;
    }
    isol = *l + 1;
    if (nsoln >= isol) {
	dcopy_(&niv, &w[(*n + 1) * w_dim1 + 1], &c__1, &temp[1], &c__1);
	i__1 = isol;
	for (j = nsoln; j >= i__1; --j) {
	    if (j > krank) {
		i__ = niv - nsoln + j;
	    } else {
		i__ = j;
	    }

	    if (j > krank && j <= *l) {
		z__[j] = 0.;
	    } else {
		z__[j] = temp[i__] / w[i__ + j * w_dim1];
		i__2 = i__ - 1;
		d__1 = -z__[j];
		daxpy_(&i__2, &d__1, &w[j * w_dim1 + 1], &c__1, &temp[1], &
			c__1);
	    }
/* L170: */
	}
    }

/*     Increment iteration counter and check against maximum number */
/*     of iterations. */

    ++iter;
    if (iter > itmax) {
	*mode = 1;
	done = TRUE_;
    }

/*     Check to see if any constraints have become active. */
/*     If so, calculate an interpolation factor so that all */
/*     active constraints are removed from the basis. */

    alpha = 2.;
    hitcon = FALSE_;
    i__1 = nsoln;
    for (j = *l + 1; j <= i__1; ++j) {
	zz = z__[j];
	if (zz <= 0.) {
	    t = x[j] / (x[j] - zz);
	    if (t < alpha) {
		alpha = t;
		jcon = j;
	    }
	    hitcon = TRUE_;
	}
/* L180: */
    }

/*     Compute search direction and feasible point */

    if (hitcon) {

/*        To add constraints, use computed ALPHA to interpolate between */
/*        last feasible solution X(*) and current unconstrained (and */
/*        infeasible) solution Z(*). */

	i__1 = nsoln;
	for (j = *l + 1; j <= i__1; ++j) {
	    x[j] += alpha * (z__[j] - x[j]);
/* L190: */
	}
	feasbl = FALSE_;

/*        Remove column JCON and shift columns JCON+1 through N to the */
/*        left.  Swap column JCON into the N th position.  This achieves */
/*        upper Hessenberg form for the nonactive constraints and */
/*        leaves an upper Hessenberg matrix to retriangularize. */

L200:
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = w[i__ + jcon * w_dim1];
	    i__2 = *n - jcon;
	    dcopy_(&i__2, &w[i__ + (jcon + 1) * w_dim1], mdw, &w[i__ + jcon * 
		    w_dim1], mdw);
	    w[i__ + *n * w_dim1] = t;
/* L210: */
	}

/*        Update permuted index vector to reflect this shift and swap. */

	itemp = ipivot[jcon];
	i__1 = *n - 1;
	for (i__ = jcon; i__ <= i__1; ++i__) {
	    ipivot[i__] = ipivot[i__ + 1];
/* L220: */
	}
	ipivot[*n] = itemp;

/*        Similarly permute X(*) vector. */

	i__1 = *n - jcon;
	dcopy_(&i__1, &x[jcon + 1], &c__1, &x[jcon], &c__1);
	x[*n] = 0.;
	--nsoln;
	--niv;

/*        Retriangularize upper Hessenberg matrix after adding */
/*        constraints. */

	i__ = krank + jcon - *l;
	i__1 = nsoln;
	for (j = jcon; j <= i__1; ++j) {
	    if (itype[i__] == 0 && itype[i__ + 1] == 0) {

/*              Zero IP1 to I in column J */

		if (w[i__ + 1 + j * w_dim1] != 0.) {
		    drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1]
			    , &w[i__ + 1 + j * w_dim1], sparam);
		    w[i__ + 1 + j * w_dim1] = 0.;
		    i__2 = *n + 1 - j;
		    drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 
			    + (j + 1) * w_dim1], mdw, sparam);
		}
	    } else if (itype[i__] == 1 && itype[i__ + 1] == 1) {

/*              Zero IP1 to I in column J */

		if (w[i__ + 1 + j * w_dim1] != 0.) {
		    drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1]
			    , &w[i__ + 1 + j * w_dim1], sparam);
		    w[i__ + 1 + j * w_dim1] = 0.;
		    i__2 = *n + 1 - j;
		    drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 
			    + (j + 1) * w_dim1], mdw, sparam);
		}
	    } else if (itype[i__] == 1 && itype[i__ + 1] == 0) {
		i__2 = *n + 1;
		dswap_(&i__2, &w[i__ + w_dim1], mdw, &w[i__ + 1 + w_dim1], 
			mdw);
		dswap_(&c__1, &scale[i__], &c__1, &scale[i__ + 1], &c__1);
		itemp = itype[i__ + 1];
		itype[i__ + 1] = itype[i__];
		itype[i__] = itemp;

/*              Swapped row was formerly a pivot element, so it will */
/*              be large enough to perform elimination. */
/*              Zero IP1 to I in column J. */

		if (w[i__ + 1 + j * w_dim1] != 0.) {
		    drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1]
			    , &w[i__ + 1 + j * w_dim1], sparam);
		    w[i__ + 1 + j * w_dim1] = 0.;
		    i__2 = *n + 1 - j;
		    drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 
			    + (j + 1) * w_dim1], mdw, sparam);
		}
	    } else if (itype[i__] == 0 && itype[i__ + 1] == 1) {
/* Computing 2nd power */
		d__1 = w[i__ + j * w_dim1];
/* Computing 2nd power */
		d__2 = tau * eanorm;
		if (scale[i__] * (d__1 * d__1) / alsq > d__2 * d__2) {

/*                 Zero IP1 to I in column J */

		    if (w[i__ + 1 + j * w_dim1] != 0.) {
			drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * 
				w_dim1], &w[i__ + 1 + j * w_dim1], sparam);
			w[i__ + 1 + j * w_dim1] = 0.;
			i__2 = *n + 1 - j;
			drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ 
				+ 1 + (j + 1) * w_dim1], mdw, sparam);
		    }
		} else {
		    i__2 = *n + 1;
		    dswap_(&i__2, &w[i__ + w_dim1], mdw, &w[i__ + 1 + w_dim1],
			     mdw);
		    dswap_(&c__1, &scale[i__], &c__1, &scale[i__ + 1], &c__1);
		    itemp = itype[i__ + 1];
		    itype[i__ + 1] = itype[i__];
		    itype[i__] = itemp;
		    w[i__ + 1 + j * w_dim1] = 0.;
		}
	    }
	    ++i__;
/* L230: */
	}

/*        See if the remaining coefficients in the solution set are */
/*        feasible.  They should be because of the way ALPHA was */
/*        determined.  If any are infeasible, it is due to roundoff */
/*        error.  Any that are non-positive will be set to zero and */
/*        removed from the solution set. */

	i__1 = nsoln;
	for (jcon = *l + 1; jcon <= i__1; ++jcon) {
	    if (x[jcon] <= 0.) {
		goto L250;
	    }
/* L240: */
	}
	feasbl = TRUE_;
L250:
	if (! feasbl) {
	    goto L200;
	}
    } else {

/*        To perform multiplier test and drop a constraint. */

	dcopy_(&nsoln, &z__[1], &c__1, &x[1], &c__1);
	if (nsoln < *n) {
	    i__1 = *n - nsoln;
	    dcopy_(&i__1, &c_b34, &c__0, &x[nsoln + 1], &c__1);
	}

/*        Reclassify least squares equations as equalities as necessary. */

	i__ = niv + 1;
L260:
	if (i__ <= me) {
	    if (itype[i__] == 0) {
		++i__;
	    } else {
		i__1 = *n + 1;
		dswap_(&i__1, &w[i__ + w_dim1], mdw, &w[me + w_dim1], mdw);
		dswap_(&c__1, &scale[i__], &c__1, &scale[me], &c__1);
		itemp = itype[i__];
		itype[i__] = itype[me];
		itype[me] = itemp;
		--me;
	    }
	    goto L260;
	}

/*        Form inner product vector WD(*) of dual coefficients. */

	i__1 = *n;
	for (j = nsoln + 1; j <= i__1; ++j) {
	    sm = 0.;
	    i__2 = m;
	    for (i__ = nsoln + 1; i__ <= i__2; ++i__) {
		sm += scale[i__] * w[i__ + j * w_dim1] * w[i__ + (*n + 1) * 
			w_dim1];
/* L270: */
	    }
	    wd[j] = sm;
/* L280: */
	}

/*        Find J such that WD(J)=WMAX is maximum.  This determines */
/*        that the incoming column J will reduce the residual vector */
/*        and be positive. */

L290:
	wmax = 0.;
	iwmax = nsoln + 1;
	i__1 = *n;
	for (j = nsoln + 1; j <= i__1; ++j) {
	    if (wd[j] > wmax) {
		wmax = wd[j];
		iwmax = j;
	    }
/* L300: */
	}
	if (wmax <= 0.) {
	    goto L330;
	}

/*        Set dual coefficients to zero for incoming column. */

	wd[iwmax] = 0.;

/*        WMAX .GT. 0.D0, so okay to move column IWMAX to solution set. */
/*        Perform transformation to retriangularize, and test for near */
/*        linear dependence. */

/*        Swap column IWMAX into NSOLN-th position to maintain upper */
/*        Hessenberg form of adjacent columns, and add new column to */
/*        triangular decomposition. */

	++nsoln;
	++niv;
	if (nsoln != iwmax) {
	    dswap_(&m, &w[nsoln * w_dim1 + 1], &c__1, &w[iwmax * w_dim1 + 1], 
		    &c__1);
	    wd[iwmax] = wd[nsoln];
	    wd[nsoln] = 0.;
	    itemp = ipivot[nsoln];
	    ipivot[nsoln] = ipivot[iwmax];
	    ipivot[iwmax] = itemp;
	}

/*        Reduce column NSOLN so that the matrix of nonactive constraints */
/*        variables is triangular. */

	i__1 = niv + 1;
	for (j = m; j >= i__1; --j) {
	    jp = j - 1;

/*           When operating near the ME line, test to see if the pivot */
/*           element is near zero.  If so, use the largest element above */
/*           it as the pivot.  This is to maintain the sharp interface */
/*           between weighted and non-weighted rows in all cases. */

	    if (j == me + 1) {
		imax = me;
/* Computing 2nd power */
		d__1 = w[me + nsoln * w_dim1];
		amax = scale[me] * (d__1 * d__1);
		i__2 = niv;
		for (jp = j - 1; jp >= i__2; --jp) {
/* Computing 2nd power */
		    d__1 = w[jp + nsoln * w_dim1];
		    t = scale[jp] * (d__1 * d__1);
		    if (t > amax) {
			imax = jp;
			amax = t;
		    }
/* L310: */
		}
		jp = imax;
	    }

	    if (w[j + nsoln * w_dim1] != 0.) {
		drotmg_(&scale[jp], &scale[j], &w[jp + nsoln * w_dim1], &w[j 
			+ nsoln * w_dim1], sparam);
		w[j + nsoln * w_dim1] = 0.;
		i__2 = *n + 1 - nsoln;
		drotm_(&i__2, &w[jp + (nsoln + 1) * w_dim1], mdw, &w[j + (
			nsoln + 1) * w_dim1], mdw, sparam);
	    }
/* L320: */
	}

/*        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if */
/*        this is nonpositive or too large.  If this was true or if the */
/*        pivot term was zero, reject the column as dependent. */

	if (w[niv + nsoln * w_dim1] != 0.) {
	    isol = niv;
	    z2 = w[isol + (*n + 1) * w_dim1] / w[isol + nsoln * w_dim1];
	    z__[nsoln] = z2;
	    pos = z2 > 0.;
	    if (z2 * eanorm >= bnorm && pos) {
		pos = ! (blowup * z2 * eanorm >= bnorm);
	    }

/*           Try to add row ME+1 as an additional equality constraint. */
/*           Check size of proposed new solution component. */
/*           Reject it if it is too large. */

	} else if (niv <= me && w[me + 1 + nsoln * w_dim1] != 0.) {
	    isol = me + 1;
	    if (pos) {

/*              Swap rows ME+1 and NIV, and scale factors for these rows. */

		i__1 = *n + 1;
		dswap_(&i__1, &w[me + 1 + w_dim1], mdw, &w[niv + w_dim1], mdw)
			;
		dswap_(&c__1, &scale[me + 1], &c__1, &scale[niv], &c__1);
		itemp = itype[me + 1];
		itype[me + 1] = itype[niv];
		itype[niv] = itemp;
		++me;
	    }
	} else {
	    pos = FALSE_;
	}

	if (! pos) {
	    --nsoln;
	    --niv;
	}
	if (! (pos || done)) {
	    goto L290;
	}
    }
    goto L160;

/*     Else perform multiplier test and drop a constraint.  To compute */
/*     final solution.  Solve system, store results in X(*). */

/*     Copy right hand side into TEMP vector to use overwriting method. */

L330:
    isol = 1;
    if (nsoln >= isol) {
	dcopy_(&niv, &w[(*n + 1) * w_dim1 + 1], &c__1, &temp[1], &c__1);
	i__1 = isol;
	for (j = nsoln; j >= i__1; --j) {
	    if (j > krank) {
		i__ = niv - nsoln + j;
	    } else {
		i__ = j;
	    }

	    if (j > krank && j <= *l) {
		z__[j] = 0.;
	    } else {
		z__[j] = temp[i__] / w[i__ + j * w_dim1];
		i__2 = i__ - 1;
		d__1 = -z__[j];
		daxpy_(&i__2, &d__1, &w[j * w_dim1 + 1], &c__1, &temp[1], &
			c__1);
	    }
/* L340: */
	}
    }

/*     Solve system. */

    dcopy_(&nsoln, &z__[1], &c__1, &x[1], &c__1);

/*     Apply Householder transformations to X(*) if KRANK.LT.L */

    if (krank < *l) {
	i__1 = krank;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = krank + 1;
	    dh12_(&c__2, &i__, &i__2, l, &w[i__ + w_dim1], mdw, &h__[i__], &x[
		    1], &c__1, &c__1, &c__1);
/* L350: */
	}
    }

/*     Fill in trailing zeroes for constrained variables not in solution. */

    if (nsoln < *n) {
	i__1 = *n - nsoln;
	dcopy_(&i__1, &c_b34, &c__0, &x[nsoln + 1], &c__1);
    }

/*     Permute solution vector to natural order. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__;
L360:
	if (ipivot[j] == i__) {
	    goto L370;
	}
	++j;
	goto L360;

L370:
	ipivot[j] = ipivot[i__];
	ipivot[i__] = j;
	dswap_(&c__1, &x[j], &c__1, &x[i__], &c__1);
/* L380: */
    }

/*     Rescale the solution using the column scaling. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] *= d__[j];
/* L390: */
    }

    i__1 = m;
    for (i__ = nsoln + 1; i__ <= i__1; ++i__) {
	t = w[i__ + (*n + 1) * w_dim1];
	if (i__ <= me) {
	    t /= alamda;
	}
	t = scale[i__] * t * t;
	*rnorm += t;
/* L400: */
    }

    *rnorm = sqrt(*rnorm);
    return 0;
} /* dwnlsm_ */


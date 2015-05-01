/* dwnlit.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b9 = 0.;
static integer c__0 = 0;

/* DECK DWNLIT */
/* Subroutine */ int dwnlit_(doublereal *w, integer *mdw, integer *m, integer 
	*n, integer *l, integer *ipivot, integer *itype, doublereal *h__, 
	doublereal *scale, doublereal *rnorm, integer *idope, doublereal *
	dope, logical *done)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer i1, j1, l1, lb, me, jj, jp, ir;
    static doublereal rn;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal tau;
    static integer niv;
    static doublereal hbar;
    static integer lend, mend;
    static doublereal amax;
    static integer imax;
    static doublereal alsq;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical indep;
    static integer krank;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), drotm_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer nsoln;
    extern /* Subroutine */ int dwnlt1_(integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern logical dwnlt2_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int dwnlt3_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);
    static logical recalc;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal factor, eanorm, sparam[5];
    extern /* Subroutine */ int drotmg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  DWNLIT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DWNNLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DWNNLS( ). */
/*     The documentation for DWNNLS( ) has complete usage instructions. */

/*     Note  The M by (N+1) matrix W( , ) contains the rt. hand side */
/*           B as the (N+1)st col. */

/*     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with */
/*     col interchanges. */

/* ***SEE ALSO  DWNNLS */
/* ***ROUTINES CALLED  DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1, */
/*                    DWNLT2, DWNLT3, IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and revised.  (WRB & RWC) */
/*   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900604  DP version created from SP version. .  (RWC) */
/* ***END PROLOGUE  DWNLIT */



/* ***FIRST EXECUTABLE STATEMENT  DWNLIT */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --ipivot;
    --itype;
    --h__;
    --scale;
    --idope;
    --dope;

    /* Function Body */
    me = idope[1];
    nsoln = idope[2];
    l1 = idope[3];

    alsq = dope[1];
    eanorm = dope[2];
    tau = dope[3];

/* Computing MIN */
    i__1 = *m - 1;
    lb = min(i__1,*l);
    recalc = TRUE_;
    *rnorm = 0.;
    krank = 0;

/*     We set FACTOR=1.0 so that the heavy weight ALAMDA will be */
/*     included in the test for column independence. */

    factor = 1.;
    lend = *l;
    i__1 = lb;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Set IR to point to the I-th row. */

	ir = i__;
	mend = *m;
	dwnlt1_(&i__, &lend, m, &ir, mdw, &recalc, &imax, &hbar, &h__[1], &
		scale[1], &w[w_offset]);

/*        Update column SS and find pivot column. */

	dwnlt3_(&i__, &imax, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);

/*        Perform column interchange. */
/*        Test independence of incoming column. */

L130:
	if (dwnlt2_(&me, &mend, &ir, &factor, &tau, &scale[1], &w[i__ * 
		w_dim1 + 1])) {

/*           Eliminate I-th column below diagonal using modified Givens */
/*           transformations applied to (A B). */

/*           When operating near the ME line, use the largest element */
/*           above it as the pivot. */

	    i__2 = i__ + 1;
	    for (j = *m; j >= i__2; --j) {
		jp = j - 1;
		if (j == me + 1) {
		    imax = me;
/* Computing 2nd power */
		    d__1 = w[me + i__ * w_dim1];
		    amax = scale[me] * (d__1 * d__1);
		    i__3 = i__;
		    for (jp = j - 1; jp >= i__3; --jp) {
/* Computing 2nd power */
			d__1 = w[jp + i__ * w_dim1];
			t = scale[jp] * (d__1 * d__1);
			if (t > amax) {
			    imax = jp;
			    amax = t;
			}
/* L150: */
		    }
		    jp = imax;
		}

		if (w[j + i__ * w_dim1] != 0.) {
		    drotmg_(&scale[jp], &scale[j], &w[jp + i__ * w_dim1], &w[
			    j + i__ * w_dim1], sparam);
		    w[j + i__ * w_dim1] = 0.;
		    i__3 = *n + 1 - i__;
		    drotm_(&i__3, &w[jp + (i__ + 1) * w_dim1], mdw, &w[j + (
			    i__ + 1) * w_dim1], mdw, sparam);
		}
/* L160: */
	    }
	} else if (lend > i__) {

/*           Column I is dependent.  Swap with column LEND. */
/*           Perform column interchange, */
/*           and find column in remaining set with largest SS. */

	    dwnlt3_(&i__, &lend, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);
	    --lend;
	    i__2 = lend - i__ + 1;
	    imax = idamax_(&i__2, &h__[i__], &c__1) + i__ - 1;
	    hbar = h__[imax];
	    goto L130;
	} else {
	    krank = i__ - 1;
	    goto L190;
	}
/* L180: */
    }
    krank = l1;

L190:
    if (krank < me) {
	factor = alsq;
	i__1 = me;
	for (i__ = krank + 1; i__ <= i__1; ++i__) {
	    dcopy_(l, &c_b9, &c__0, &w[i__ + w_dim1], mdw);
/* L200: */
	}

/*        Determine the rank of the remaining equality constraint */
/*        equations by eliminating within the block of constrained */
/*        variables.  Remove any redundant constraints. */

	recalc = TRUE_;
/* Computing MIN */
	i__1 = *l + me - krank;
	lb = min(i__1,*n);
	i__1 = lb;
	for (i__ = *l + 1; i__ <= i__1; ++i__) {
	    ir = krank + i__ - *l;
	    lend = *n;
	    mend = me;
	    dwnlt1_(&i__, &lend, &me, &ir, mdw, &recalc, &imax, &hbar, &h__[1]
		    , &scale[1], &w[w_offset]);

/*           Update col ss and find pivot col */

	    dwnlt3_(&i__, &imax, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);

/*           Perform column interchange */
/*           Eliminate elements in the I-th col. */

	    i__2 = ir + 1;
	    for (j = me; j >= i__2; --j) {
		if (w[j + i__ * w_dim1] != 0.) {
		    drotmg_(&scale[j - 1], &scale[j], &w[j - 1 + i__ * w_dim1]
			    , &w[j + i__ * w_dim1], sparam);
		    w[j + i__ * w_dim1] = 0.;
		    i__3 = *n + 1 - i__;
		    drotm_(&i__3, &w[j - 1 + (i__ + 1) * w_dim1], mdw, &w[j + 
			    (i__ + 1) * w_dim1], mdw, sparam);
		}
/* L240: */
	    }

/*           I=column being eliminated. */
/*           Test independence of incoming column. */
/*           Remove any redundant or dependent equality constraints. */

	    if (! dwnlt2_(&me, &mend, &ir, &factor, &tau, &scale[1], &w[i__ * 
		    w_dim1 + 1])) {
		jj = ir;
		i__2 = me;
		for (ir = jj; ir <= i__2; ++ir) {
		    dcopy_(n, &c_b9, &c__0, &w[ir + w_dim1], mdw);
		    *rnorm += scale[ir] * w[ir + (*n + 1) * w_dim1] / alsq * 
			    w[ir + (*n + 1) * w_dim1];
		    w[ir + (*n + 1) * w_dim1] = 0.;
		    scale[ir] = 1.;

/*                 Reclassify the zeroed row as a least squares equation. */

		    itype[ir] = 1;
/* L260: */
		}

/*              Reduce ME to reflect any discovered dependent equality */
/*              constraints. */

		me = jj - 1;
		goto L280;
	    }
/* L270: */
	}
    }

/*     Try to determine the variables KRANK+1 through L1 from the */
/*     least squares equations.  Continue the triangularization with */
/*     pivot element W(ME+1,I). */

L280:
    if (krank < l1) {
	recalc = TRUE_;

/*        Set FACTOR=ALSQ to remove effect of heavy weight from */
/*        test for column independence. */

	factor = alsq;
	i__1 = l1;
	for (i__ = krank + 1; i__ <= i__1; ++i__) {

/*           Set IR to point to the ME+1-st row. */

	    ir = me + 1;
	    lend = *l;
	    mend = *m;
	    dwnlt1_(&i__, l, m, &ir, mdw, &recalc, &imax, &hbar, &h__[1], &
		    scale[1], &w[w_offset]);

/*           Update column SS and find pivot column. */

	    dwnlt3_(&i__, &imax, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);

/*           Perform column interchange. */
/*           Eliminate I-th column below the IR-th element. */

	    i__2 = ir + 1;
	    for (j = *m; j >= i__2; --j) {
		if (w[j + i__ * w_dim1] != 0.) {
		    drotmg_(&scale[j - 1], &scale[j], &w[j - 1 + i__ * w_dim1]
			    , &w[j + i__ * w_dim1], sparam);
		    w[j + i__ * w_dim1] = 0.;
		    i__3 = *n + 1 - i__;
		    drotm_(&i__3, &w[j - 1 + (i__ + 1) * w_dim1], mdw, &w[j + 
			    (i__ + 1) * w_dim1], mdw, sparam);
		}
/* L320: */
	    }

/*           Test if new pivot element is near zero. */
/*           If so, the column is dependent. */
/*           Then check row norm test to be classified as independent. */

/* Computing 2nd power */
	    d__1 = w[ir + i__ * w_dim1];
	    t = scale[ir] * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = tau * eanorm;
	    indep = t > d__1 * d__1;
	    if (indep) {
		rn = 0.;
		i__2 = *m;
		for (i1 = ir; i1 <= i__2; ++i1) {
		    i__3 = *n;
		    for (j1 = i__ + 1; j1 <= i__3; ++j1) {
/* Computing MAX */
/* Computing 2nd power */
			d__3 = w[i1 + j1 * w_dim1];
			d__1 = rn, d__2 = scale[i1] * (d__3 * d__3);
			rn = max(d__1,d__2);
/* L330: */
		    }
/* L340: */
		}
/* Computing 2nd power */
		d__1 = tau;
		indep = t > rn * (d__1 * d__1);
	    }

/*           If independent, swap the IR-th and KRANK+1-th rows to */
/*           maintain the triangular form.  Update the rank indicator */
/*           KRANK and the equality constraint pointer ME. */

	    if (! indep) {
		goto L360;
	    }
	    i__2 = *n + 1;
	    dswap_(&i__2, &w[krank + 1 + w_dim1], mdw, &w[ir + w_dim1], mdw);
	    dswap_(&c__1, &scale[krank + 1], &c__1, &scale[ir], &c__1);

/*           Reclassify the least square equation as an equality */
/*           constraint and rescale it. */

	    itype[ir] = 0;
	    t = sqrt(scale[krank + 1]);
	    i__2 = *n + 1;
	    dscal_(&i__2, &t, &w[krank + 1 + w_dim1], mdw);
	    scale[krank + 1] = alsq;
	    ++me;
	    ++krank;
/* L350: */
	}
    }

/*     If pseudorank is less than L, apply Householder transformation. */
/*     from right. */

L360:
    if (krank < *l) {
	for (j = krank; j >= 1; --j) {
	    i__1 = krank + 1;
	    i__2 = j - 1;
	    dh12_(&c__1, &j, &i__1, l, &w[j + w_dim1], mdw, &h__[j], &w[
		    w_offset], mdw, &c__1, &i__2);
/* L370: */
	}
    }

    niv = krank + nsoln - *l;
    if (*l == *n) {
	*done = TRUE_;
    }

/*     End of initial triangularization. */

    idope[1] = me;
    idope[2] = krank;
    idope[3] = niv;
    return 0;
} /* dwnlit_ */


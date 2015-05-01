/* dexint.f -- translated by f2c (version 12.02.01).
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
static integer c__15 = 15;
static integer c__5 = 5;
static integer c__1 = 1;

/* DECK DEXINT */
/* Subroutine */ int dexint_(doublereal *x, integer *n, integer *kode, 
	integer *m, doublereal *tol, doublereal *en, integer *nz, integer *
	ierr)
{
    /* Initialized data */

    static doublereal xcut = 2.;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal a[99], b[99];
    static integer i__, k;
    static doublereal s, y[2], p1, p2, y1, y2, aa, cc, ah, ak, bk;
    static integer ic;
    static doublereal em;
    static integer ik;
    static doublereal at, bt, ct;
    static integer kk, kn, ml, nd, nm;
    static doublereal fx;
    static integer ks, ix, mu;
    static doublereal pt, tx, yt;
    static integer i1m, ind, ict;
    static doublereal fnm, emx, aams, etol;
    static integer jset;
    static doublereal xlim, xtol;
    static integer icase;
    static doublereal cnorm;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    extern doublereal dpsixn_(integer *);

/* ***BEGIN PROLOGUE  DEXINT */
/* ***PURPOSE  Compute an M member sequence of exponential integrals */
/*            E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C5 */
/* ***TYPE      DOUBLE PRECISION (EXINT-S, DEXINT-D) */
/* ***KEYWORDS  EXPONENTIAL INTEGRAL, SPECIAL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*         DEXINT computes M member sequences of exponential integrals */
/*         E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.  The */
/*         exponential integral is defined by */

/*         E(N,X)=integral on (1,infinity) of EXP(-XT)/T**N */

/*         where X=0.0 and N=1 cannot occur simultaneously.  Formulas */
/*         and notation are found in the NBS Handbook of Mathematical */
/*         Functions (ref. 1). */

/*         The power series is implemented for X .LE. XCUT and the */
/*         confluent hypergeometric representation */

/*                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X) */

/*         is computed for X .GT. XCUT.  Since sequences are computed in */
/*         a stable fashion by recurring away from X, A is selected as */
/*         the integer closest to X within the constraint N .LE. A .LE. */
/*         N+M-1.  For the U computation, A is further modified to be the */
/*         nearest even integer.  Indices are carried forward or */
/*         backward by the two term recursion relation */

/*                     K*E(K+1,X) + X*E(K,X) = EXP(-X) */

/*         once E(A,X) is computed.  The U function is computed by means */
/*         of the backward recursive Miller algorithm applied to the */
/*         three term contiguous relation for U(A+K,A,X), K=0,1,... */
/*         This produces accurate ratios and determines U(A+K,A,X), and */
/*         hence E(A,X), to within a multiplicative constant C. */
/*         Another contiguous relation applied to C*U(A,A,X) and */
/*         C*U(A+1,A,X) gets C*U(A+1,A+1,X), a quantity proportional to */
/*         E(A+1,X).  The normalizing constant C is obtained from the */
/*         two term recursion relation above with K=A. */

/*         The maximum number of significant digits obtainable */
/*         is the smaller of 14 and the number of digits carried in */
/*         double precision arithmetic. */

/*     Description of Arguments */

/*         Input     * X and TOL are double precision * */
/*           X       X .GT. 0.0 for N=1 and  X .GE. 0.0 for N .GE. 2 */
/*           N       order of the first member of the sequence, N .GE. 1 */
/*                   (X=0.0 and N=1 is an error) */
/*           KODE    a selection parameter for scaled values */
/*                   KODE=1   returns        E(N+K,X), K=0,1,...,M-1. */
/*                       =2   returns EXP(X)*E(N+K,X), K=0,1,...,M-1. */
/*           M       number of exponential integrals in the sequence, */
/*                   M .GE. 1 */
/*           TOL     relative accuracy wanted, ETOL .LE. TOL .LE. 0.1 */
/*                   ETOL is the larger of double precision unit */
/*                   roundoff = D1MACH(4) and 1.0D-18 */

/*         Output    * EN is a double precision vector * */
/*           EN      a vector of dimension at least M containing values */
/*                   EN(K) = E(N+K-1,X) or EXP(X)*E(N+K-1,X), K=1,M */
/*                   depending on KODE */
/*           NZ      underflow indicator */
/*                   NZ=0   a normal return */
/*                   NZ=M   X exceeds XLIM and an underflow occurs. */
/*                          EN(K)=0.0D0 , K=1,M returned on KODE=1 */
/*           IERR    error flag */
/*                   IERR=0, normal return, computation completed */
/*                   IERR=1, input error,   no computation */
/*                   IERR=2, error,         no computation */
/*                           algorithm termination condition not met */

/* ***REFERENCES  M. Abramowitz and I. A. Stegun, Handbook of */
/*                 Mathematical Functions, NBS AMS Series 55, U.S. Dept. */
/*                 of Commerce, 1955. */
/*               D. E. Amos, Computation of exponential integrals, ACM */
/*                 Transactions on Mathematical Software 6, (1980), */
/*                 pp. 365-377 and pp. 420-428. */
/* ***ROUTINES CALLED  D1MACH, DPSIXN, I1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   910408  Updated the REFERENCES section.  (WRB) */
/*   920207  Updated with code with a revision date of 880811 from */
/*           D. Amos.  Included correction of argument list.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DEXINT */
    /* Parameter adjustments */
    --en;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DEXINT */
    *ierr = 0;
    *nz = 0;
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    etol = max(d__1,5e-19);
    if (*x < 0.) {
	*ierr = 1;
    }
    if (*n < 1) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*m < 1) {
	*ierr = 1;
    }
    if (*tol < etol || *tol > .1) {
	*ierr = 1;
    }
    if (*x == 0. && *n == 1) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    i1m = -i1mach_(&c__15);
    pt = i1m * 2.3026 * d1mach_(&c__5);
    xlim = pt - 6.907755;
    bt = pt + (*n + *m - 1);
    if (bt > 1e3) {
	xlim = pt - log(bt);
    }

    if (*x > xcut) {
	goto L100;
    }
    if (*x == 0. && *n > 1) {
	goto L80;
    }
/* ----------------------------------------------------------------------- */
/*     SERIES FOR E(N,X) FOR X.LE.XCUT */
/* ----------------------------------------------------------------------- */
    tx = *x + .5;
    ix = (integer) tx;
/* ----------------------------------------------------------------------- */
/*     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1 */
/*     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N.GE.2 */
/* ----------------------------------------------------------------------- */
    icase = 2;
    if (ix > *n) {
	icase = 1;
    }
    nm = *n - icase + 1;
    nd = nm + 1;
    ind = 3 - icase;
    mu = *m - ind;
    ml = 1;
    ks = nd;
    fnm = (doublereal) nm;
    s = 0.;
    xtol = *tol * 3.;
    if (nd == 1) {
	goto L10;
    }
    xtol = *tol * .3333;
    s = 1. / fnm;
L10:
    aa = 1.;
    ak = 1.;
    ic = 35;
    if (*x < etol) {
	ic = 1;
    }
    i__1 = ic;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aa = -aa * *x / ak;
	if (i__ == nm) {
	    goto L30;
	}
	s -= aa / (ak - fnm);
	if (abs(aa) <= xtol * abs(s)) {
	    goto L20;
	}
	ak += 1.;
	goto L50;
L20:
	if (i__ < 2) {
	    goto L40;
	}
	if (nd - 2 > i__ || i__ > nd - 1) {
	    goto L60;
	}
	ak += 1.;
	goto L50;
L30:
	s += aa * (-log(*x) + dpsixn_(&nd));
	xtol = *tol * 3.;
L40:
	ak += 1.;
L50:
	;
    }
    if (ic != 1) {
	goto L340;
    }
L60:
    if (nd == 1) {
	s += -log(*x) + dpsixn_(&c__1);
    }
    if (*kode == 2) {
	s *= exp(*x);
    }
    en[1] = s;
    emx = 1.;
    if (*m == 1) {
	goto L70;
    }
    en[ind] = s;
    aa = (doublereal) ks;
    if (*kode == 1) {
	emx = exp(-(*x));
    }
    switch (icase) {
	case 1:  goto L220;
	case 2:  goto L240;
    }
L70:
    if (icase == 2) {
	return 0;
    }
    if (*kode == 1) {
	emx = exp(-(*x));
    }
    en[1] = (emx - s) / *x;
    return 0;
L80:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	en[i__] = 1. / (*n + i__ - 2);
/* L90: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECURSIVE MILLER ALGORITHM FOR */
/*              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X) */
/*     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X. */
/*     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION */
/* ----------------------------------------------------------------------- */
L100:
    emx = 1.;
    if (*kode == 2) {
	goto L130;
    }
    if (*x <= xlim) {
	goto L120;
    }
    *nz = *m;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	en[i__] = 0.;
/* L110: */
    }
    return 0;
L120:
    emx = exp(-(*x));
L130:
    tx = *x + .5;
    ix = (integer) tx;
    kn = *n + *m - 1;
    if (kn <= ix) {
	goto L140;
    }
    if (*n < ix && ix < kn) {
	goto L170;
    }
    if (*n >= ix) {
	goto L160;
    }
    goto L340;
L140:
    icase = 1;
    ks = kn;
    ml = *m - 1;
    mu = -1;
    ind = *m;
    if (kn > 1) {
	goto L180;
    }
L150:
    ks = 2;
    icase = 3;
    goto L180;
L160:
    icase = 2;
    ind = 1;
    ks = *n;
    mu = *m - 1;
    if (*n > 1) {
	goto L180;
    }
    if (kn == 1) {
	goto L150;
    }
    ix = 2;
L170:
    icase = 1;
    ks = ix;
    ml = ix - *n;
    ind = ml + 1;
    mu = kn - ix;
L180:
    ik = ks / 2;
    ah = (doublereal) ik;
    jset = ks + 1 - (ik + ik);
/* ----------------------------------------------------------------------- */
/*     START COMPUTATION FOR */
/*              EN(IND) = C*U( A , A ,X)    JSET=1 */
/*              EN(IND) = C*U(A+1,A+1,X)    JSET=2 */
/*     FOR AN EVEN INTEGER A. */
/* ----------------------------------------------------------------------- */
    ic = 0;
    aa = ah + ah;
    aams = aa - 1.;
    aams *= aams;
    tx = *x + *x;
    fx = tx + tx;
    ak = ah;
    xtol = *tol;
    if (*tol <= .001) {
	xtol = *tol * 20.;
    }
    ct = aams + fx * ah;
    em = (ah + 1.) / ((*x + aa) * xtol * sqrt(ct));
    bk = aa;
    cc = ah * ah;
/* ----------------------------------------------------------------------- */
/*     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD */
/*     RECURSION */
/* ----------------------------------------------------------------------- */
    p1 = 0.;
    p2 = 1.;
L190:
    if (ic == 99) {
	goto L340;
    }
    ++ic;
    ak += 1.;
    at = bk / (bk + ak + cc + ic);
    bk = bk + ak + ak;
    a[ic - 1] = at;
    bt = (ak + ak + *x) / (ak + 1.);
    b[ic - 1] = bt;
    pt = p2;
    p2 = bt * p2 - at * p1;
    p1 = pt;
    ct += fx;
    em = em * at * (1. - tx / ct);
    if (em * (ak + 1.) > p1 * p1) {
	goto L190;
    }
    ict = ic;
    kk = ic + 1;
    bt = tx / (ct + fx);
    y2 = bk / (bk + cc + kk) * (p1 / p2) * (1. - bt + bt * .375 * bt);
    y1 = 1.;
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECURRENCE FOR */
/*              Y1=             C*U( A ,A,X) */
/*              Y2= C*(A/(1+A/2))*U(A+1,A,X) */
/* ----------------------------------------------------------------------- */
    i__1 = ict;
    for (k = 1; k <= i__1; ++k) {
	--kk;
	yt = y1;
	y1 = (b[kk - 1] * y1 - y2) / a[kk - 1];
	y2 = yt;
/* L200: */
    }
/* ----------------------------------------------------------------------- */
/*     THE CONTIGUOUS RELATION */
/*              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X) */
/*     WITH  B=A+1 , C=A IS USED FOR */
/*              Y(2) = C * U(A+1,A+1,X) */
/*     X IS INCORPORATED INTO THE NORMALIZING RELATION */
/* ----------------------------------------------------------------------- */
    pt = y2 / y1;
    cnorm = 1.f - pt * (ah + 1.f) / aa;
    y[0] = 1.f / (cnorm * aa + *x);
    y[1] = cnorm * y[0];
    if (icase == 3) {
	goto L210;
    }
    en[ind] = emx * y[jset - 1];
    if (*m == 1) {
	return 0;
    }
    aa = (doublereal) ks;
    switch (icase) {
	case 1:  goto L220;
	case 2:  goto L240;
    }
/* ----------------------------------------------------------------------- */
/*     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX */
/* ----------------------------------------------------------------------- */
L210:
    en[1] = emx * (1.f - y[0]) / *x;
    return 0;
L220:
    k = ind - 1;
    i__1 = ml;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aa += -1.;
	en[k] = (emx - aa * en[k + 1]) / *x;
	--k;
/* L230: */
    }
    if (mu <= 0) {
	return 0;
    }
    aa = (doublereal) ks;
L240:
    k = ind;
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	en[k + 1] = (emx - *x * en[k]) / aa;
	aa += 1.;
	++k;
/* L250: */
    }
    return 0;
L340:
    *ierr = 2;
    return 0;
} /* dexint_ */


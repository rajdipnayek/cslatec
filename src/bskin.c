/* bskin.f -- translated by f2c (version 12.02.01).
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

static integer c__12 = 12;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__11 = 11;
static integer c__2 = 2;

/* DECK BSKIN */
/* Subroutine */ int bskin_(real *x, integer *n, integer *kode, integer *m, 
	real *y, integer *nz, integer *ierr)
{
    /* Initialized data */

    static real a[50] = { 1.f,.5f,.375f,.3125f,.2734375f,.24609375f,
	    .2255859375f,.20947265625f,.196380615234375f,.1854705810546875f,
	    .176197052001953125f,.168188095092773438f,.161180257797241211f,
	    .154981017112731934f,.149445980787277222f,.144464448094367981f,
	    .139949934091418982f,.135833759559318423f,.132060599571559578f,
	    .128585320635465905f,.125370687619579257f,.122385671247684513f,
	    .119604178719328047f,.117004087877603524f,.114566502713486784f,
	    .112275172659217048f,.110116034723462874f,.108076848895250599f,
	    .106146905164978267f,.104316786110409676f,.102578173008569515f,
	    .100923686347140974f,.0993467537479668965f,.0978414999033007314f,
	    .0964026543164874854f,.0950254735405376642f,.0937056752969190855f,
	    .09243938238750126f,.0912230747245078224f,.0900535481254756708f,
	    .0889278787739072249f,.0878433924473961612f,.0867976377754033498f,
	    .0857883629175498224f,.0848134951571231199f,.0838711229887106408f,
	    .0829594803475290034f,.0820769326842574183f,.0812219646354630702f,
	    .0803931690779583449f };
    static real hrtpi = .886226925452758014f;

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static real h__[31];
    static integer i__, k;
    static real w;
    static integer m3;
    static real t1, t2;
    static integer ne;
    static real fn;
    static integer il, kk;
    static real hn, gr;
    static integer nl, nn, np, ns, nt;
    static real ss, xp, ys[3];
    static integer i1m;
    static real exi[102], tol, yss[3];
    static integer nflg, nlim;
    static real xlim;
    static integer icase;
    extern /* Subroutine */ int bkias_(real *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, integer *);
    static real enlim;
    extern doublereal gamrn_(real *);
    extern /* Subroutine */ int bkisr_(real *, integer *, real *, integer *);
    static real xnlim;
    extern /* Subroutine */ int exint_(real *, integer *, integer *, integer *
	    , real *, real *, integer *, integer *);
    static integer ktrms;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  BSKIN */
/* ***PURPOSE  Compute repeated integrals of the K-zero Bessel function. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10F */
/* ***TYPE      SINGLE PRECISION (BSKIN-S, DBSKIN-D) */
/* ***KEYWORDS  BICKLEY FUNCTIONS, EXPONENTIAL INTEGRAL, */
/*             INTEGRALS OF BESSEL FUNCTIONS, K-ZERO BESSEL FUNCTION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*         The following definitions are used in BSKIN: */

/*   Definition 1 */
/*         KI(0,X) = K-zero Bessel function. */

/*   Definition 2 */
/*         KI(N,X) = Bickley Function */
/*                 =  integral from X to infinity of KI(N-1,t)dt */
/*                     for X .ge. 0 and N = 1,2,... */
/*   ____________________________________________________________________ */
/*      BSKIN computes sequences of Bickley functions (repeated integrals */
/*      of the K0 Bessel function); i.e. for fixed X and N and K=1,..., */
/*      BSKIN computes the M-member sequence */

/*                     Y(K) =        KI(N+K-1,X) for KODE=1 */
/*      or */
/*                     Y(K) = EXP(X)*KI(N+K-1,X) for KODE=2, */

/*      for N.ge.0 and X.ge.0 (N and X cannot be zero simultaneously). */

/*      INPUT */
/*        X      - Argument, X .ge. 0.0E0 */
/*        N      - Order of first member of the sequence N .ge. 0 */
/*        KODE   - Selection parameter */
/*                 KODE = 1 returns Y(K)=       KI(N+K-1,X), K=1,M */
/*                      = 2 returns Y(K)=EXP(X)*KI(N+K-1,X), K=1,M */
/*        M      - Number of members in the sequence, M.ge.1 */

/*      OUTPUT */
/*        Y      - A vector of dimension at least M containing the */
/*                 sequence selected by KODE. */
/*        NZ     - Underflow flag */
/*                 NZ = 0 means computation completed */
/*                    = M means an exponential underflow occurred on */
/*                        KODE=1.  Y(K)=0.0E0, K=1,...,M is returned */
/*        IERR   - Error flag */
/*                 IERR = 0, Normal return, computation completed. */
/*                      = 1, Input error,   no computation. */
/*                      = 2, Error,         no computation.  The */
/*                           termination condition was not met. */

/*      The nominal computational accuracy is the maximum of unit */
/*      roundoff (=R1MACH(4)) and 1.0e-18 since critical constants */
/*      are given to only 18 digits. */

/*      DBSKIN is the double precision version of BSKIN. */

/* *Long Description: */

/*         Numerical recurrence on */

/*      (L-1)*KI(L,X) = X(KI(L-3,X) - KI(L-1,X)) + (L-2)*KI(L-2,X) */

/*         is stable where recurrence is carried forward or backward */
/*         away from INT(X+0.5).  The power series for indices 0,1 and 2 */
/*         on 0.le.X.le. 2 starts a stable recurrence for indices */
/*         greater than 2.  If N is sufficiently large (N.gt.NLIM), the */
/*         uniform asymptotic expansion for N to INFINITY is more */
/*         economical.  On X.gt.2 the recursion is started by evaluating */
/*         the uniform expansion for the three members whose indices are */
/*         closest to INT(X+0.5) within the set N,...,N+M-1.  Forward */
/*         recurrence, backward recurrence or both, complete the */
/*         sequence depending on the relation of INT(X+0.5) to the */
/*         indices N,...,N+M-1. */

/* ***REFERENCES  D. E. Amos, Uniform asymptotic expansions for */
/*                 exponential integrals E(N,X) and Bickley functions */
/*                 KI(N,X), ACM Transactions on Mathematical Software, */
/*                 1983. */
/*               D. E. Amos, A portable Fortran subroutine for the */
/*                 Bickley functions KI(N,X), Algorithm 609, ACM */
/*                 Transactions on Mathematical Software, 1983. */
/* ***ROUTINES CALLED  BKIAS, BKISR, EXINT, GAMRN, I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced statement label.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BSKIN */
/* ----------------------------------------------------------------------- */
/*             COEFFICIENTS IN SERIES OF EXPONENTIAL INTEGRALS */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/*             SQRT(PI)/2 */
/* ----------------------------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  BSKIN */
    *ierr = 0;
    *nz = 0;
    if (*x < 0.f) {
	*ierr = 1;
    }
    if (*n < 0) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*m < 1) {
	*ierr = 1;
    }
    if (*x == 0.f && *n == 0) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    if (*x == 0.f) {
	goto L300;
    }
    i1m = -i1mach_(&c__12);
    t1 = r1mach_(&c__5) * 2.3026f * i1m;
    xlim = t1 - 3.228086f;
    t2 = t1 + *n + *m - 1;
    if (t2 > 1e3f) {
	xlim = t1 - (log(t2) - .451583f) * .5f;
    }
    if (*x > xlim && *kode == 1) {
	goto L320;
    }
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,1e-18f);
    i1m = i1mach_(&c__11);
/* ----------------------------------------------------------------------- */
/*     LN(NLIM) = 0.125*LN(EPS),   NLIM = 2*KTRMS+N */
/* ----------------------------------------------------------------------- */
    xnlim = (i1m - 1) * .287823f * r1mach_(&c__5);
    enlim = exp(xnlim);
    nlim = (integer) enlim + 2;
    nlim = min(100,nlim);
    nlim = max(20,nlim);
    m3 = min(*m,3);
    nl = *n + *m - 1;
    if (*x > 2.f) {
	goto L130;
    }
    if (*n > nlim) {
	goto L280;
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTATION BY SERIES FOR 0.LE.X.LE.2 */
/* ----------------------------------------------------------------------- */
    nflg = 0;
    nn = *n;
    if (nl <= 2) {
	goto L60;
    }
    m3 = 3;
    nn = 0;
    nflg = 1;
L60:
    xp = 1.f;
    if (*kode == 2) {
	xp = exp(*x);
    }
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bkisr_(x, &nn, &w, ierr);
	if (*ierr != 0) {
	    return 0;
	}
	w *= xp;
	if (nn < *n) {
	    goto L70;
	}
	kk = nn - *n + 1;
	y[kk] = w;
L70:
	ys[i__ - 1] = w;
	++nn;
/* L80: */
    }
    if (nflg == 0) {
	return 0;
    }
    ns = nn;
    xp = 1.f;
L90:
/* ----------------------------------------------------------------------- */
/*     FORWARD RECURSION SCALED BY EXP(X) ON ICASE=0,1,2 */
/* ----------------------------------------------------------------------- */
    fn = (real) (ns - 1);
    il = nl - ns + 1;
    if (il <= 0) {
	return 0;
    }
    i__1 = il;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1 = ys[1];
	t2 = ys[2];
	ys[2] = (*x * (ys[0] - ys[2]) + (fn - 1.f) * ys[1]) / fn;
	ys[1] = t2;
	ys[0] = t1;
	fn += 1.f;
	if (ns < *n) {
	    goto L100;
	}
	kk = ns - *n + 1;
	y[kk] = ys[2] * xp;
L100:
	++ns;
/* L110: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     COMPUTATION BY ASYMPTOTIC EXPANSION FOR X.GT.2 */
/* ----------------------------------------------------------------------- */
L130:
    w = *x + .5f;
    nt = (integer) w;
    if (nl > nt) {
	goto L270;
    }
/* ----------------------------------------------------------------------- */
/*     CASE NL.LE.NT, ICASE=0 */
/* ----------------------------------------------------------------------- */
    icase = 0;
    nn = nl;
/* Computing MIN */
    i__1 = *m - m3;
    nflg = min(i__1,1);
L140:
    kk = (nlim - nn) / 2;
    ktrms = max(0,kk);
    ns = nn + 1;
    np = nn - m3 + 1;
    xp = 1.f;
    if (*kode == 1) {
	xp = exp(-(*x));
    }
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kk = i__;
	bkias_(x, &np, &ktrms, a, &w, &kk, &ne, &gr, h__, ierr);
	if (*ierr != 0) {
	    return 0;
	}
	ys[i__ - 1] = w;
	++np;
/* L150: */
    }
/* ----------------------------------------------------------------------- */
/*     SUM SERIES OF EXPONENTIAL INTEGRALS BACKWARD */
/* ----------------------------------------------------------------------- */
    if (ktrms == 0) {
	goto L160;
    }
    ne = ktrms + ktrms + 1;
    np = nn - m3 + 2;
    exint_(x, &np, &c__2, &ne, &tol, exi, nz, ierr);
    if (*nz != 0) {
	goto L320;
    }
    if (*ierr == 2) {
	return 0;
    }
L160:
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ss = 0.f;
	if (ktrms == 0) {
	    goto L180;
	}
	kk = i__ + ktrms + ktrms - 2;
	il = ktrms;
	i__2 = ktrms;
	for (k = 1; k <= i__2; ++k) {
	    ss += a[il - 1] * exi[kk - 1];
	    kk += -2;
	    --il;
/* L170: */
	}
L180:
	ys[i__ - 1] += ss;
/* L190: */
    }
    if (icase == 1) {
	goto L200;
    }
    if (nflg != 0) {
	goto L220;
    }
L200:
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = ys[i__ - 1] * xp;
/* L210: */
    }
    if (icase == 1 && nflg == 1) {
	goto L90;
    }
    return 0;
L220:
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECURSION SCALED BY EXP(X) ICASE=0,2 */
/* ----------------------------------------------------------------------- */
    kk = nn - *n + 1;
    k = m3;
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[kk] = ys[k - 1] * xp;
	yss[i__ - 1] = ys[i__ - 1];
	--kk;
	--k;
/* L230: */
    }
    il = kk;
    if (il <= 0) {
	goto L250;
    }
    fn = (real) (nn - 3);
    i__1 = il;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1 = ys[1];
	t2 = ys[0];
	ys[0] = ys[1] + ((fn + 2.f) * ys[2] - (fn + 1.f) * ys[0]) / *x;
	ys[1] = t2;
	ys[2] = t1;
	y[kk] = ys[0] * xp;
	--kk;
	fn += -1.f;
/* L240: */
    }
L250:
    if (icase != 2) {
	return 0;
    }
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ys[i__ - 1] = yss[i__ - 1];
/* L260: */
    }
    goto L90;
L270:
    if (*n < nt) {
	goto L290;
    }
/* ----------------------------------------------------------------------- */
/*     ICASE=1, NT.LE.N.LE.NL WITH FORWARD RECURSION */
/* ----------------------------------------------------------------------- */
L280:
    nn = *n + m3 - 1;
/* Computing MIN */
    i__1 = *m - m3;
    nflg = min(i__1,1);
    icase = 1;
    goto L140;
/* ----------------------------------------------------------------------- */
/*     ICASE=2, N.LT.NT.LT.NL WITH BOTH FORWARD AND BACKWARD RECURSION */
/* ----------------------------------------------------------------------- */
L290:
    nn = nt + 1;
/* Computing MIN */
    i__1 = *m - m3;
    nflg = min(i__1,1);
    icase = 2;
    goto L140;
/* ----------------------------------------------------------------------- */
/*     X=0 CASE */
/* ----------------------------------------------------------------------- */
L300:
    fn = (real) (*n);
    hn = fn * .5f;
    gr = gamrn_(&hn);
    y[1] = hrtpi * gr;
    if (*m == 1) {
	return 0;
    }
    y[2] = hrtpi / (hn * gr);
    if (*m == 2) {
	return 0;
    }
    i__1 = *m;
    for (k = 3; k <= i__1; ++k) {
	y[k] = fn * y[k - 2] / (fn + 1.f);
	fn += 1.f;
/* L310: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     UNDERFLOW ON KODE=1, X.GT.XLIM */
/* ----------------------------------------------------------------------- */
L320:
    *nz = *m;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.f;
/* L330: */
    }
    return 0;
} /* bskin_ */


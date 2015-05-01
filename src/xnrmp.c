/* xnrmp.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static real c_b4 = 0.f;
static integer c__112 = 112;
static integer c__1 = 1;
static integer c__113 = 113;

/* DECK XNRMP */
/* Subroutine */ int xnrmp_(integer *nu, integer *mu1, integer *mu2, real *
	sarg, integer *mode, real *spn, integer *ipn, integer *isig, integer *
	ierror)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k;
    static real p, s, t, x, c1, c2, p1, p2, p3;
    static integer ip;
    static real rk;
    static integer mu;
    static real sx, tx;
    static integer ip1, ip2;
    extern /* Subroutine */ int xadd_(real *, integer *, real *, integer *, 
	    real *, integer *, integer *), xadj_(real *, integer *, integer *)
	    , xred_(real *, integer *, integer *), xset_(integer *, integer *,
	     real *, integer *, integer *), xermsg_(char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  XNRMP */
/* ***PURPOSE  Compute normalized Legendre polynomials. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      SINGLE PRECISION (XNRMP-S, DXNRMP-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */

/*        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS */
/*        (DXNRMP is double-precision version) */
/*        XNRMP calculates normalized Legendre polynomials of varying */
/*        order and fixed argument and degree. The order MU and degree */
/*        NU are non-negative integers and the argument is real. Because */
/*        the algorithm requires the use of numbers outside the normal */
/*        machine range, this subroutine employs a special arithmetic */
/*        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver, */
/*        and D.W. Lozier, Extended-Range Arithmetic and Normalized */
/*        Legendre Polynomials, ACM Transactions on Mathematical Soft- */
/*        ware, 93-105, March 1981, for a complete description of the */
/*        algorithm and special arithmetic. Also see program comments */
/*        in XSET. */

/*        The normalized Legendre polynomials are multiples of the */
/*        associated Legendre polynomials of the first kind where the */
/*        normalizing coefficients are chosen so as to make the integral */
/*        from -1 to 1 of the square of each function equal to 1. See */
/*        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions, */
/*        McGraw-Hill, New York, 1960, p. 121. */

/*        The input values to XNRMP are NU, MU1, MU2, SARG, and MODE. */
/*        These must satisfy */
/*          1. NU .GE. 0 specifies the degree of the normalized Legendre */
/*             polynomial that is wanted. */
/*          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre */
/*             polynomial that is wanted. */
/*          3. MU2 .GE. MU1 specifies the highest-order normalized Leg- */
/*             endre polynomial that is wanted. */
/*         4a. MODE = 1 and -1.0 .LE. SARG .LE. 1.0 specifies that */
/*             Normalized Legendre(NU, MU, SARG) is wanted for MU = MU1, */
/*             MU1 + 1, ..., MU2. */
/*         4b. MODE = 2 and -3.14159... .LT. SARG .LT. 3.14159... spec- */
/*             ifies that Normalized Legendre(NU, MU, COS(SARG)) is want- */
/*             ed for MU = MU1, MU1 + 1, ..., MU2. */

/*        The output of XNRMP consists of the two vectors SPN and IPN */
/*        and the error estimate ISIG. The computed values are stored as */
/*        extended-range numbers such that */
/*             (SPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,X) */
/*             (SPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,X) */
/*                . */
/*                . */
/*             (SPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,X) */
/*        where K = MU2 - MU1 + 1 and X = SARG or COS(SARG) according */
/*        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the */
/*        number of decimal digits lost through rounding errors in the */
/*        computation. For example if SARG is accurate to 12 significant */
/*        decimals, then the computed function values are accurate to */
/*        12 - ISIG significant decimals (except in neighborhoods of */
/*        zeros). */

/*        The interpretation of (SPN(I),IPN(I)) is SPN(I)*(IR**IPN(I)) */
/*        where IR is the internal radix of the computer arithmetic. When */
/*        IPN(I) = 0 the value of the normalized Legendre polynomial is */
/*        contained entirely in SPN(I) and subsequent single-precision */
/*        computations can be performed without further consideration of */
/*        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre- */
/*        sponding value of the normalized Legendre polynomial cannot be */
/*        represented in single-precision because of overflow or under- */
/*        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case */
/*        that IPN(I) is nonzero, the user should try using double pre- */
/*        cision if it has a wider exponent range. If double precision */
/*        fails, the user could rewrite his/her program to use extended- */
/*        range arithmetic. */

/*        The interpretation of (SPN(I),IPN(I)) can be changed to */
/*        SPN(I)*(10**IPN(I)) by calling the extended-range subroutine */
/*        XCON. This should be done before printing the computed values. */
/*        As an example of usage, the Fortran coding */
/*              J = K */
/*              DO 20 I = 1, K */
/*              CALL XCON(SPN(I), IPN(I),IERROR) */
/*              IF (IERROR.NE.0) RETURN */
/*              PRINT 10, SPN(I), IPN(I) */
/*           10 FORMAT(1X, E30.8 , I15) */
/*              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20 */
/*              J = I - 1 */
/*           20 CONTINUE */
/*        will print all computed values and determine the largest J */
/*        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the */
/*        change of representation caused by calling XCON, (SPN(I), */
/*        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent */
/*        extended-range computations. */

/*        IERROR is an error indicator. If no errors are detected, */
/*        IERROR=0 when control returns to the calling routine. If */
/*        an error is detected, IERROR is returned as nonzero. The */
/*        calling routine must check the value of IERROR. */

/*        If IERROR=112 or 113, invalid input was provided to XNRMP. */
/*        If IERROR=101,102,103, or 104, invalid input was provided */
/*        to XSET. */
/*        If IERROR=105 or 106, an internal consistency error occurred */
/*        in XSET (probably due to a software malfunction in the */
/*        library routine I1MACH). */
/*        If IERROR=107, an overflow or underflow of an extended-range */
/*        number was detected in XADJ. */
/*        If IERROR=108, an overflow or underflow of an extended-range */
/*        number was detected in XC210. */

/* ***SEE ALSO  XSET */
/* ***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and */
/*                 Normalized Legendre Polynomials, ACM Trans on Math */
/*                 Softw, v 7, n 1, March 1981, pp 93--105. */
/* ***ROUTINES CALLED  XADD, XADJ, XERMSG, XRED, XSET */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*           CALLs to XERROR changed to CALLs to XERMSG.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XNRMP */
/* CALL XSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE XSET */
/* LISTING FOR DETAILS) */
/* ***FIRST EXECUTABLE STATEMENT  XNRMP */
    /* Parameter adjustments */
    --ipn;
    --spn;

    /* Function Body */
    *ierror = 0;
    xset_(&c__0, &c__0, &c_b4, &c__0, ierror);
    if (*ierror != 0) {
	return 0;
    }

/*        TEST FOR PROPER INPUT VALUES. */

    if (*nu < 0) {
	goto L110;
    }
    if (*mu1 < 0) {
	goto L110;
    }
    if (*mu1 > *mu2) {
	goto L110;
    }
    if (*nu == 0) {
	goto L90;
    }
    if (*mode < 1 || *mode > 2) {
	goto L110;
    }
    switch (*mode) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    if (dabs(*sarg) > 1.f) {
	goto L120;
    }
    if (dabs(*sarg) == 1.f) {
	goto L90;
    }
    x = *sarg;
    sx = sqrt((dabs(x) + 1.f) * (.5f - dabs(x) + .5f));
    tx = x / sx;
/* Computing 2nd power */
    r__2 = tx;
    r__1 = *nu * 2.f * (r__2 * r__2 + 5.f);
    *isig = r_lg10(&r__1);
    goto L30;
L20:
    if (dabs(*sarg) > atan(1.f) * 4.f) {
	goto L120;
    }
    if (*sarg == 0.f) {
	goto L90;
    }
    x = cos(*sarg);
    sx = (r__1 = sin(*sarg), dabs(r__1));
    tx = x / sx;
    r__2 = *nu * 2.f * ((r__1 = *sarg * tx, dabs(r__1)) + 5.f);
    *isig = r_lg10(&r__2);

/*        BEGIN CALCULATION */

L30:
    mu = *mu2;
    i__ = *mu2 - *mu1 + 1;

/*        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0. */

L40:
    if (mu <= *nu) {
	goto L50;
    }
    spn[i__] = 0.f;
    ipn[i__] = 0;
    --i__;
    --mu;
    if (i__ > 0) {
	goto L40;
    }
    *isig = 0;
    goto L160;
L50:
    mu = *nu;

/*        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X) */

    p1 = 0.f;
    ip1 = 0;

/*        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X) */

    p2 = 1.f;
    ip2 = 0;
    p3 = .5f;
    rk = 2.f;
    i__1 = *nu;
    for (j = 1; j <= i__1; ++j) {
	p3 = (rk + 1.f) / rk * p3;
	p2 *= sx;
	xadj_(&p2, &ip2, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	rk += 2.f;
/* L60: */
    }
    p2 *= sqrt(p3);
    xadj_(&p2, &ip2, ierror);
    if (*ierror != 0) {
	return 0;
    }
    s = tx * 2.f;
    t = 1.f / *nu;
    if (*mu2 < *nu) {
	goto L70;
    }
    spn[i__] = p2;
    ipn[i__] = ip2;
    --i__;
    if (i__ == 0) {
	goto L140;
    }

/*        RECURRENCE PROCESS */

L70:
    p = mu * t;
    c1 = 1.f / sqrt((1.f - p + t) * (p + 1.f));
    c2 = s * p * c1 * p2;
    c1 = -sqrt((p + 1.f + t) * (1.f - p)) * c1 * p1;
    xadd_(&c2, &ip2, &c1, &ip1, &p, &ip, ierror);
    if (*ierror != 0) {
	return 0;
    }
    --mu;
    if (mu > *mu2) {
	goto L80;
    }

/*        STORE IN ARRAY SPN FOR RETURN TO CALLING ROUTINE. */

    spn[i__] = p;
    ipn[i__] = ip;
    --i__;
    if (i__ == 0) {
	goto L140;
    }
L80:
    p1 = p2;
    ip1 = ip2;
    p2 = p;
    ip2 = ip;
    if (mu <= *mu1) {
	goto L140;
    }
    goto L70;

/*        SPECIAL CASE WHEN X=-1 OR +1, OR NU=0. */

L90:
    k = *mu2 - *mu1 + 1;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	spn[i__] = 0.f;
	ipn[i__] = 0;
/* L100: */
    }
    *isig = 0;
    if (*mu1 > 0) {
	goto L160;
    }
    *isig = 1;
    spn[1] = sqrt(*nu + .5f);
    ipn[1] = 0;
    if (*nu % 2 == 0) {
	goto L160;
    }
    if (*mode == 1 && *sarg == 1.f) {
	goto L160;
    }
    if (*mode == 2) {
	goto L160;
    }
    spn[1] = -spn[1];
    goto L160;

/*          ERROR PRINTOUTS AND TERMINATION. */

L110:
    xermsg_("SLATEC", "XNRMP", "NU, MU1, MU2 or MODE not valid", &c__112, &
	    c__1, (ftnlen)6, (ftnlen)5, (ftnlen)30);
    *ierror = 112;
    return 0;
L120:
    xermsg_("SLATEC", "XNRMP", "SARG out of range", &c__113, &c__1, (ftnlen)6,
	     (ftnlen)5, (ftnlen)17);
    *ierror = 113;
    return 0;

/*        RETURN TO CALLING PROGRAM */

L140:
    k = *mu2 - *mu1 + 1;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xred_(&spn[i__], &ipn[i__], ierror);
	if (*ierror != 0) {
	    return 0;
	}
/* L150: */
    }
L160:
    return 0;
} /* xnrmp_ */


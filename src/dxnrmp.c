/* dxnrmp.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b4 = 0.;
static integer c__212 = 212;
static integer c__1 = 1;
static integer c__213 = 213;

/* DECK DXNRMP */
/* Subroutine */ int dxnrmp_(integer *nu, integer *mu1, integer *mu2, 
	doublereal *darg, integer *mode, doublereal *dpn, integer *ipn, 
	integer *isig, integer *ierror)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal p, s, t, x, c1, c2, p1, p2, p3, dk;
    static integer ip, mu;
    static doublereal sx, tx;
    static integer ip1, ip2;
    extern /* Subroutine */ int dxadd_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dxadj_(doublereal 
	    *, integer *, integer *), dxred_(doublereal *, integer *, integer 
	    *), dxset_(integer *, integer *, doublereal *, integer *, integer 
	    *), xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DXNRMP */
/* ***PURPOSE  Compute normalized Legendre polynomials. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      DOUBLE PRECISION (XNRMP-S, DXNRMP-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */

/*        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS */
/*        (XNRMP is single-precision version) */
/*        DXNRMP calculates normalized Legendre polynomials of varying */
/*        order and fixed argument and degree. The order MU and degree */
/*        NU are non-negative integers and the argument is real. Because */
/*        the algorithm requires the use of numbers outside the normal */
/*        machine range, this subroutine employs a special arithmetic */
/*        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver, */
/*        and D.W. Lozier, Extended-Range Arithmetic and Normalized */
/*        Legendre Polynomials, ACM Transactions on Mathematical Soft- */
/*        ware, 93-105, March 1981, for a complete description of the */
/*        algorithm and special arithmetic. Also see program comments */
/*        in DXSET. */

/*        The normalized Legendre polynomials are multiples of the */
/*        associated Legendre polynomials of the first kind where the */
/*        normalizing coefficients are chosen so as to make the integral */
/*        from -1 to 1 of the square of each function equal to 1. See */
/*        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions, */
/*        McGraw-Hill, New York, 1960, p. 121. */

/*        The input values to DXNRMP are NU, MU1, MU2, DARG, and MODE. */
/*        These must satisfy */
/*          1. NU .GE. 0 specifies the degree of the normalized Legendre */
/*             polynomial that is wanted. */
/*          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre */
/*             polynomial that is wanted. */
/*          3. MU2 .GE. MU1 specifies the highest-order normalized Leg- */
/*             endre polynomial that is wanted. */
/*         4a. MODE = 1 and -1.0D0 .LE. DARG .LE. 1.0D0 specifies that */
/*             Normalized Legendre(NU, MU, DARG) is wanted for MU = MU1, */
/*             MU1 + 1, ..., MU2. */
/*         4b. MODE = 2 and -3.14159... .LT. DARG .LT. 3.14159... spec- */
/*             ifies that Normalized Legendre(NU, MU, COS(DARG)) is */
/*             wanted for MU = MU1, MU1 + 1, ..., MU2. */

/*        The output of DXNRMP consists of the two vectors DPN and IPN */
/*        and the error estimate ISIG. The computed values are stored as */
/*        extended-range numbers such that */
/*             (DPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,DX) */
/*             (DPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,DX) */
/*                . */
/*                . */
/*             (DPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,DX) */
/*        where K = MU2 - MU1 + 1 and DX = DARG or COS(DARG) according */
/*        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the */
/*        number of decimal digits lost through rounding errors in the */
/*        computation. For example if DARG is accurate to 12 significant */
/*        decimals, then the computed function values are accurate to */
/*        12 - ISIG significant decimals (except in neighborhoods of */
/*        zeros). */

/*        The interpretation of (DPN(I),IPN(I)) is DPN(I)*(IR**IPN(I)) */
/*        where IR is the internal radix of the computer arithmetic. When */
/*        IPN(I) = 0 the value of the normalized Legendre polynomial is */
/*        contained entirely in DPN(I) and subsequent double-precision */
/*        computations can be performed without further consideration of */
/*        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre- */
/*        sponding value of the normalized Legendre polynomial cannot be */
/*        represented in double-precision because of overflow or under- */
/*        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case */
/*        that IPN(I) is nonzero, the user could rewrite his/her program */
/*        to use extended range arithmetic. */



/*        The interpretation of (DPN(I),IPN(I)) can be changed to */
/*        DPN(I)*(10**IPN(I)) by calling the extended-range subroutine */
/*        DXCON. This should be done before printing the computed values. */
/*        As an example of usage, the Fortran coding */
/*              J = K */
/*              DO 20 I = 1, K */
/*              CALL DXCON(DPN(I), IPN(I),IERROR) */
/*              IF (IERROR.NE.0) RETURN */
/*              PRINT 10, DPN(I), IPN(I) */
/*           10 FORMAT(1X, D30.18 , I15) */
/*              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20 */
/*              J = I - 1 */
/*           20 CONTINUE */
/*        will print all computed values and determine the largest J */
/*        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the */
/*        change of representation caused by calling DXCON, (DPN(I), */
/*        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent */
/*        extended-range computations. */

/*        IERROR is an error indicator. If no errors are detected, */
/*        IERROR=0 when control returns to the calling routine. If */
/*        an error is detected, IERROR is returned as nonzero. The */
/*        calling routine must check the value of IERROR. */

/*        If IERROR=212 or 213, invalid input was provided to DXNRMP. */
/*        If IERROR=201,202,203, or 204, invalid input was provided */
/*        to DXSET. */
/*        If IERROR=205 or 206, an internal consistency error occurred */
/*        in DXSET (probably due to a software malfunction in the */
/*        library routine I1MACH). */
/*        If IERROR=207, an overflow or underflow of an extended-range */
/*        number was detected in DXADJ. */
/*        If IERROR=208, an overflow or underflow of an extended-range */
/*        number was detected in DXC210. */

/* ***SEE ALSO  DXSET */
/* ***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and */
/*                 Normalized Legendre Polynomials, ACM Trans on Math */
/*                 Softw, v 7, n 1, March 1981, pp 93--105. */
/* ***ROUTINES CALLED  DXADD, DXADJ, DXRED, DXSET, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*           CALLs to XERROR changed to CALLs to XERMSG.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXNRMP */
/* CALL DXSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE DXSET */
/* LISTING FOR DETAILS) */
/* ***FIRST EXECUTABLE STATEMENT  DXNRMP */
    /* Parameter adjustments */
    --ipn;
    --dpn;

    /* Function Body */
    *ierror = 0;
    dxset_(&c__0, &c__0, &c_b4, &c__0, ierror);
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
    if (abs(*darg) > 1.) {
	goto L120;
    }
    if (abs(*darg) == 1.) {
	goto L90;
    }
    x = *darg;
    sx = sqrt((abs(x) + 1.) * (.5 - abs(x) + .5));
    tx = x / sx;
/* Computing 2nd power */
    d__2 = tx;
    d__1 = *nu * 2. * (d__2 * d__2 + 5.);
    *isig = (integer) d_lg10(&d__1);
    goto L30;
L20:
    if (abs(*darg) > atan(1.) * 4.) {
	goto L120;
    }
    if (*darg == 0.) {
	goto L90;
    }
    x = cos(*darg);
    sx = (d__1 = sin(*darg), abs(d__1));
    tx = x / sx;
    d__2 = *nu * 2. * ((d__1 = *darg * tx, abs(d__1)) + 5.);
    *isig = (integer) d_lg10(&d__2);

/*        BEGIN CALCULATION */

L30:
    mu = *mu2;
    i__ = *mu2 - *mu1 + 1;

/*        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0. */

L40:
    if (mu <= *nu) {
	goto L50;
    }
    dpn[i__] = 0.;
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

    p1 = 0.;
    ip1 = 0;

/*        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X) */

    p2 = 1.;
    ip2 = 0;
    p3 = .5;
    dk = 2.;
    i__1 = *nu;
    for (j = 1; j <= i__1; ++j) {
	p3 = (dk + 1.) / dk * p3;
	p2 *= sx;
	dxadj_(&p2, &ip2, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	dk += 2.;
/* L60: */
    }
    p2 *= sqrt(p3);
    dxadj_(&p2, &ip2, ierror);
    if (*ierror != 0) {
	return 0;
    }
    s = tx * 2.;
    t = 1. / *nu;
    if (*mu2 < *nu) {
	goto L70;
    }
    dpn[i__] = p2;
    ipn[i__] = ip2;
    --i__;
    if (i__ == 0) {
	goto L140;
    }

/*        RECURRENCE PROCESS */

L70:
    p = mu * t;
    c1 = 1. / sqrt((1. - p + t) * (p + 1.));
    c2 = s * p * c1 * p2;
    c1 = -sqrt((p + 1. + t) * (1. - p)) * c1 * p1;
    dxadd_(&c2, &ip2, &c1, &ip1, &p, &ip, ierror);
    if (*ierror != 0) {
	return 0;
    }
    --mu;
    if (mu > *mu2) {
	goto L80;
    }

/*        STORE IN ARRAY DPN FOR RETURN TO CALLING ROUTINE. */

    dpn[i__] = p;
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
	dpn[i__] = 0.;
	ipn[i__] = 0;
/* L100: */
    }
    *isig = 0;
    if (*mu1 > 0) {
	goto L160;
    }
    *isig = 1;
    dpn[1] = sqrt(*nu + .5);
    ipn[1] = 0;
    if (*nu % 2 == 0) {
	goto L160;
    }
    if (*mode == 1 && *darg == 1.) {
	goto L160;
    }
    if (*mode == 2) {
	goto L160;
    }
    dpn[1] = -dpn[1];
    goto L160;

/*          ERROR PRINTOUTS AND TERMINATION. */

L110:
    xermsg_("SLATEC", "DXNRMP", "NU, MU1, MU2 or MODE not valid", &c__212, &
	    c__1, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    *ierror = 212;
    return 0;
L120:
    xermsg_("SLATEC", "DXNRMP", "DARG out of range", &c__213, &c__1, (ftnlen)
	    6, (ftnlen)6, (ftnlen)17);
    *ierror = 213;
    return 0;

/*        RETURN TO CALLING PROGRAM */

L140:
    k = *mu2 - *mu1 + 1;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dxred_(&dpn[i__], &ipn[i__], ierror);
	if (*ierror != 0) {
	    return 0;
	}
/* L150: */
    }
L160:
    return 0;
} /* dxnrmp_ */


/* dxlegf.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b9 = 1.;
static integer c__210 = 210;
static integer c__1 = 1;
static integer c__211 = 211;

/* DECK DXLEGF */
/* Subroutine */ int dxlegf_(doublereal *dnu1, integer *nudiff, integer *mu1, 
	integer *mu2, doublereal *theta, integer *id, doublereal *pqa, 
	integer *ipqa, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, l;
    static doublereal x, sx, pi2, dnu2;
    extern /* Subroutine */ int dxred_(doublereal *, integer *, integer *), 
	    dxset_(integer *, integer *, doublereal *, integer *, integer *), 
	    dxpmu_(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *), dxqmu_(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dxqnu_(doublereal 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen), dxpnrm_(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *), dxpmup_(
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), dxpqnu_(doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);

/* ***BEGIN PROLOGUE  DXLEGF */
/* ***PURPOSE  Compute normalized Legendre polynomials and associated */
/*            Legendre functions. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      DOUBLE PRECISION (XLEGF-S, DXLEGF-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */

/*   DXLEGF: Extended-range Double-precision Legendre Functions */

/*   A feature of the DXLEGF subroutine for Legendre functions is */
/* the use of extended-range arithmetic, a software extension of */
/* ordinary floating-point arithmetic that greatly increases the */
/* exponent range of the representable numbers. This avoids the */
/* need for scaling the solutions to lie within the exponent range */
/* of the most restrictive manufacturer's hardware. The increased */
/* exponent range is achieved by allocating an integer storage */
/* location together with each floating-point storage location. */

/*   The interpretation of the pair (X,I) where X is floating-point */
/* and I is integer is X*(IR**I) where IR is the internal radix of */
/* the computer arithmetic. */

/*   This subroutine computes one of the following vectors: */

/* 1. Legendre function of the first kind of negative order, either */
/*    a. P(-MU1,NU,X), P(-MU1-1,NU,X), ..., P(-MU2,NU,X) or */
/*    b. P(-MU,NU1,X), P(-MU,NU1+1,X), ..., P(-MU,NU2,X) */
/* 2. Legendre function of the second kind, either */
/*    a. Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X) or */
/*    b. Q(MU,NU1,X), Q(MU,NU1+1,X), ..., Q(MU,NU2,X) */
/* 3. Legendre function of the first kind of positive order, either */
/*    a. P(MU1,NU,X), P(MU1+1,NU,X), ..., P(MU2,NU,X) or */
/*    b. P(MU,NU1,X), P(MU,NU1+1,X), ..., P(MU,NU2,X) */
/* 4. Normalized Legendre polynomials, either */
/*    a. PN(MU1,NU,X), PN(MU1+1,NU,X), ..., PN(MU2,NU,X) or */
/*    b. PN(MU,NU1,X), PN(MU,NU1+1,X), ..., PN(MU,NU2,X) */

/* where X = COS(THETA). */

/*   The input values to DXLEGF are DNU1, NUDIFF, MU1, MU2, THETA, */
/* and ID. These must satisfy */

/*    DNU1 is DOUBLE PRECISION and greater than or equal to -0.5; */
/*    NUDIFF is INTEGER and non-negative; */
/*    MU1 is INTEGER and non-negative; */
/*    MU2 is INTEGER and greater than or equal to MU1; */
/*    THETA is DOUBLE PRECISION and in the half-open interval (0,PI/2]; */
/*    ID is INTEGER and equal to 1, 2, 3 or 4; */

/* and  additionally either NUDIFF = 0 or MU2 = MU1. */

/*   If ID=1 and NUDIFF=0, a vector of type 1a above is computed */
/* with NU=DNU1. */

/*   If ID=1 and MU1=MU2, a vector of type 1b above is computed */
/* with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1. */

/*   If ID=2 and NUDIFF=0, a vector of type 2a above is computed */
/* with NU=DNU1. */

/*   If ID=2 and MU1=MU2, a vector of type 2b above is computed */
/* with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1. */

/*   If ID=3 and NUDIFF=0, a vector of type 3a above is computed */
/* with NU=DNU1. */

/*   If ID=3 and MU1=MU2, a vector of type 3b above is computed */
/* with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1. */

/*   If ID=4 and NUDIFF=0, a vector of type 4a above is computed */
/* with NU=DNU1. */

/*   If ID=4 and MU1=MU2, a vector of type 4b above is computed */
/* with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1. */

/*   In each case the vector of computed Legendre function values */
/* is returned in the extended-range vector (PQA(I),IPQA(I)). The */
/* length of this vector is either MU2-MU1+1 or NUDIFF+1. */

/*   Where possible, DXLEGF returns IPQA(I) as zero. In this case the */
/* value of the Legendre function is contained entirely in PQA(I), */
/* so it can be used in subsequent computations without further */
/* consideration of extended-range arithmetic. If IPQA(I) is nonzero, */
/* then the value of the Legendre function is not representable in */
/* floating-point because of underflow or overflow. The program that */
/* calls DXLEGF must test IPQA(I) to ensure correct usage. */

/*   IERROR is an error indicator. If no errors are detected, IERROR=0 */
/* when control returns to the calling routine. If an error is detected, */
/* IERROR is returned as nonzero. The calling routine must check the */
/* value of IERROR. */

/*   If IERROR=210 or 211, invalid input was provided to DXLEGF. */
/*   If IERROR=201,202,203, or 204, invalid input was provided to DXSET. */
/*   If IERROR=205 or 206, an internal consistency error occurred in */
/* DXSET (probably due to a software malfunction in the library routine */
/* I1MACH). */
/*   If IERROR=207, an overflow or underflow of an extended-range number */
/* was detected in DXADJ. */
/*   If IERROR=208, an overflow or underflow of an extended-range number */
/* was detected in DXC210. */

/* ***SEE ALSO  DXSET */
/* ***REFERENCES  Olver and Smith, Associated Legendre Functions on the */
/*                 Cut, J Comp Phys, v 51, n 3, Sept 1983, pp 502--518. */
/*               Smith, Olver and Lozier, Extended-Range Arithmetic and */
/*                 Normalized Legendre Polynomials, ACM Trans on Math */
/*                 Softw, v 7, n 1, March 1981, pp 93--105. */
/* ***ROUTINES CALLED  DXPMU, DXPMUP, DXPNRM, DXPQNU, DXQMU, DXQNU, DXRED, */
/*                    DXSET, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*           CALLs to XERROR changed to CALLs to XERMSG.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXLEGF */

/* ***FIRST EXECUTABLE STATEMENT  DXLEGF */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    dxset_(&c__0, &c__0, &c_b4, &c__0, ierror);
    if (*ierror != 0) {
	return 0;
    }
    pi2 = atan(1.) * 2.;

/*        ZERO OUTPUT ARRAYS */

    l = *mu2 - *mu1 + *nudiff + 1;
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pqa[i__] = 0.;
/* L290: */
	ipqa[i__] = 0;
    }

/*        CHECK FOR VALID INPUT VALUES */

    if (*nudiff < 0) {
	goto L400;
    }
    if (*dnu1 < -.5) {
	goto L400;
    }
    if (*mu2 < *mu1) {
	goto L400;
    }
    if (*mu1 < 0) {
	goto L400;
    }
    if (*theta <= 0. || *theta > pi2) {
	goto L420;
    }
    if (*id < 1 || *id > 4) {
	goto L400;
    }
    if (*mu1 != *mu2 && *nudiff > 0) {
	goto L400;
    }

/*        IF DNU1 IS NOT AN INTEGER, NORMALIZED P(MU,DNU,X) */
/*        CANNOT BE CALCULATED.  IF DNU1 IS AN INTEGER AND */
/*        MU1.GT.DNU2 THEN ALL VALUES OF P(+MU,DNU,X) AND */
/*        NORMALIZED P(MU,NU,X) WILL BE ZERO. */

    dnu2 = *dnu1 + *nudiff;
    if (*id == 3 && d_mod(dnu1, &c_b9) != 0.) {
	goto L295;
    }
    if (*id == 4 && d_mod(dnu1, &c_b9) != 0.) {
	goto L400;
    }
    if ((*id == 3 || *id == 4) && (doublereal) (*mu1) > dnu2) {
	return 0;
    }
L295:

    x = cos(*theta);
    sx = 1. / sin(*theta);
    if (*id == 2) {
	goto L300;
    }
    if (*mu2 - *mu1 <= 0) {
	goto L360;
    }

/*        FIXED NU, VARIABLE MU */
/*        CALL DXPMU TO CALCULATE P(-MU1,NU,X),....,P(-MU2,NU,X) */

    dxpmu_(dnu1, &dnu2, mu1, mu2, theta, &x, &sx, id, &pqa[1], &ipqa[1], 
	    ierror);
    if (*ierror != 0) {
	return 0;
    }
    goto L380;

L300:
    if (*mu2 == *mu1) {
	goto L320;
    }

/*        FIXED NU, VARIABLE MU */
/*        CALL DXQMU TO CALCULATE Q(MU1,NU,X),....,Q(MU2,NU,X) */

    dxqmu_(dnu1, &dnu2, mu1, mu2, theta, &x, &sx, id, &pqa[1], &ipqa[1], 
	    ierror);
    if (*ierror != 0) {
	return 0;
    }
    goto L390;

/*        FIXED MU, VARIABLE NU */
/*        CALL DXQNU TO CALCULATE Q(MU,DNU1,X),....,Q(MU,DNU2,X) */

L320:
    dxqnu_(dnu1, &dnu2, mu1, theta, &x, &sx, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    goto L390;

/*        FIXED MU, VARIABLE NU */
/*        CALL DXPQNU TO CALCULATE P(-MU,DNU1,X),....,P(-MU,DNU2,X) */

L360:
    dxpqnu_(dnu1, &dnu2, mu1, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }

/*        IF ID = 3, TRANSFORM P(-MU,NU,X) VECTOR INTO */
/*        P(MU,NU,X) VECTOR. */

L380:
    if (*id == 3) {
	dxpmup_(dnu1, &dnu2, mu1, mu2, &pqa[1], &ipqa[1], ierror);
    }
    if (*ierror != 0) {
	return 0;
    }

/*        IF ID = 4, TRANSFORM P(-MU,NU,X) VECTOR INTO */
/*        NORMALIZED P(MU,NU,X) VECTOR. */

    if (*id == 4) {
	dxpnrm_(dnu1, &dnu2, mu1, mu2, &pqa[1], &ipqa[1], ierror);
    }
    if (*ierror != 0) {
	return 0;
    }

/*        PLACE RESULTS IN REDUCED FORM IF POSSIBLE */
/*        AND RETURN TO MAIN PROGRAM. */

L390:
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dxred_(&pqa[i__], &ipqa[i__], ierror);
	if (*ierror != 0) {
	    return 0;
	}
/* L395: */
    }
    return 0;

/*        *****     ERROR TERMINATION     ***** */

L400:
    xermsg_("SLATEC", "DXLEGF", "DNU1, NUDIFF, MU1, MU2, or ID not valid", &
	    c__210, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)39);
    *ierror = 210;
    return 0;
L420:
    xermsg_("SLATEC", "DXLEGF", "THETA out of range", &c__211, &c__1, (ftnlen)
	    6, (ftnlen)6, (ftnlen)18);
    *ierror = 211;
    return 0;
} /* dxlegf_ */


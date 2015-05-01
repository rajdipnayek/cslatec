/* rj.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static doublereal c_b3 = .16666666666666666;
static integer c__1 = 1;
static doublereal c_b5 = .33333333333333331;
static integer c__2 = 2;
static integer c__6 = 6;
static integer c__10 = 10;
static integer c__9 = 9;

/* DECK RJ */
doublereal rj_(real *x, real *y, real *z__, real *p, integer *ier)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[6], a__2[10], a__3[9];
    integer i__1[6], i__2[10], i__3[9];
    real ret_val, r__1, r__2;
    doublereal d__1;
    char ch__1[88], ch__2[140], ch__3[130];

    /* Local variables */
    static real c1, c2, c3, c4, e2, e3, s1, s2, s3, ea, eb, ec;
    extern doublereal rc_(real *, real *, integer *);
    static real pn, mu, xn, yn, zn, alfa, beta;
    static char xern3[16], xern4[16], xern5[16], xern6[16], xern7[16];
    static real lamda, sigma, lolim, pndev, xndev, yndev, uplim, zndev;
    extern doublereal r1mach_(integer *);
    static real power4, epslon;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real errtol, xnroot, ynroot, znroot;

    /* Fortran I/O blocks */
    static icilist io___10 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___12 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___16 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___17 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___19 = { 0, xern6, 0, "(1PE15.6)", 16, 1 };
    static icilist io___21 = { 0, xern7, 0, "(1PE15.6)", 16, 1 };
    static icilist io___22 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___23 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___24 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___25 = { 0, xern6, 0, "(1PE15.6)", 16, 1 };
    static icilist io___26 = { 0, xern7, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  RJ */
/* ***PURPOSE  Compute the incomplete or complete (X or Y or Z is zero) */
/*            elliptic integral of the 3rd kind.  For X, Y, and Z non- */
/*            negative, at most one of them zero, and P positive, */
/*             RJ(X,Y,Z,P) = Integral from zero to infinity of */
/*                                  -1/2     -1/2     -1/2     -1 */
/*                        (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C14 */
/* ***TYPE      SINGLE PRECISION (RJ-S, DRJ-D) */
/* ***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM, */
/*             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND, */
/*             TAYLOR SERIES */
/* ***AUTHOR  Carlson, B. C. */
/*             Ames Laboratory-DOE */
/*             Iowa State University */
/*             Ames, IA  50011 */
/*           Notis, E. M. */
/*             Ames Laboratory-DOE */
/*             Iowa State University */
/*             Ames, IA  50011 */
/*           Pexton, R. L. */
/*             Lawrence Livermore National Laboratory */
/*             Livermore, CA  94550 */
/* ***DESCRIPTION */

/*   1.     RJ */
/*          Standard FORTRAN function routine */
/*          Single precision version */
/*          The routine calculates an approximation result to */
/*          RJ(X,Y,Z,P) = Integral from zero to infinity of */

/*                                -1/2     -1/2     -1/2     -1 */
/*                      (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt, */

/*          where X, Y, and Z are nonnegative, at most one of them is */
/*          zero, and P is positive.  If X or Y or Z is zero, the */
/*          integral is COMPLETE.  The duplication theorem is iterated */
/*          until the variables are nearly equal, and the function is */
/*          then expanded in Taylor series to fifth order. */


/*   2.     Calling Sequence */
/*          RJ( X, Y, Z, P, IER ) */

/*          Parameters On Entry */
/*          Values assigned by the calling routine */

/*          X      - Single precision, nonnegative variable */

/*          Y      - Single precision, nonnegative variable */

/*          Z      - Single precision, nonnegative variable */

/*          P      - Single precision, positive variable */


/*          On  Return     (values assigned by the RJ routine) */

/*          RJ     - Single precision approximation to the integral */

/*          IER    - Integer */

/*                   IER = 0 Normal and reliable termination of the */
/*                           routine.  It is assumed that the requested */
/*                           accuracy has been achieved. */

/*                   IER >  0 Abnormal termination of the routine */


/*          X, Y, Z, P are unaltered. */


/*   3.    Error Messages */

/*         Value of IER assigned by the RJ routine */

/*                  Value Assigned        Error Message Printed */
/*                  IER = 1               MIN(X,Y,Z) .LT. 0.0E0 */
/*                      = 2               MIN(X+Y,X+Z,Y+Z,P) .LT. LOLIM */
/*                      = 3               MAX(X,Y,Z,P) .GT. UPLIM */



/*   4.     Control Parameters */

/*                  Values of LOLIM, UPLIM, and ERRTOL are set by the */
/*                  routine. */


/*          LOLIM and UPLIM determine the valid range of X Y, Z, and P */

/*          LOLIM is not less than the cube root of the value */
/*          of LOLIM used in the routine for RC. */

/*          UPLIM is not greater than 0.3 times the cube root of */
/*          the value of UPLIM used in the routine for RC. */


/*                     Acceptable Values For:   LOLIM      UPLIM */
/*                     IBM 360/370 SERIES   :   2.0E-26     3.0E+24 */
/*                     CDC 6000/7000 SERIES :   5.0E-98     3.0E+106 */
/*                     UNIVAC 1100 SERIES   :   5.0E-13     6.0E+11 */
/*                     CRAY                 :   1.32E-822   1.4E+821 */
/*                     VAX 11 SERIES        :   2.5E-13     9.0E+11 */



/*          ERRTOL determines the accuracy of the answer */

/*                 The value assigned by the routine will result */
/*                 in solution precision within 1-2 decimals of */
/*                 "machine precision". */




/*          Relative error due to truncation of the series for RJ */
/*          is less than 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2. */



/*              The accuracy of the computed approximation to the inte- */
/*              gral can be controlled by choosing the value of ERRTOL. */
/*              Truncation of a Taylor series after terms of fifth order */
/*              Introduces an error less than the amount shown in the */
/*              second column of the following table for each value of */
/*              ERRTOL in the first column.  In addition to the trunca- */
/*              tion error there will be round-off error, but in prac- */
/*              tice the total error from both sources is usually less */
/*              than the amount given in the table. */



/*          Sample choices:  ERRTOL   Relative Truncation */
/*                                    error less than */
/*                           1.0E-3    4.0E-18 */
/*                           3.0E-3    3.0E-15 */
/*                           1.0E-2    4.0E-12 */
/*                           3.0E-2    3.0E-9 */
/*                           1.0E-1    4.0E-6 */

/*                    Decreasing ERRTOL by a factor of 10 yields six more */
/*                    decimal digits of accuracy at the expense of one or */
/*                    two more iterations of the duplication theorem. */

/* *Long Description: */

/*   RJ Special Comments */


/*          Check by addition theorem: RJ(X,X+Z,X+W,X+P) */
/*          + RJ(Y,Y+Z,Y+W,Y+P) + (A-B) * RJ(A,B,B,A) + 3 / SQRT(A) */
/*          = RJ(0,Z,W,P), where X,Y,Z,W,P are positive and X * Y */
/*          = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y), */
/*          and B - A = P * (P-Z) * (P-W).  The sum of the third and */
/*          fourth terms on the left side is 3 * RC(A,B). */


/*          On Input: */

/*          X, Y, Z, and P are the variables in the integral RJ(X,Y,Z,P). */


/*          On Output: */


/*          X, Y, Z, and P are unaltered. */

/*          ******************************************************** */

/*          Warning: Changes in the program may improve speed at the */
/*                   expense of robustness. */

/* ------------------------------------------------------------ */


/*   Special Functions via RJ and RF */


/*                  Legendre form of ELLIPTIC INTEGRAL of 3rd kind */
/*                  ---------------------------------------------- */


/*                               PHI         2         -1 */
/*                  P(PHI,K,N) = INT (1+N SIN (THETA) )   * */
/*                                0 */

/*                                      2    2         -1/2 */
/*                                 *(1-K  SIN (THETA) )     D THETA */


/*                                         2          2   2 */
/*                       = SIN (PHI) RF(COS (PHI), 1-K SIN (PHI),1) */

/*                                  3            2         2   2 */
/*                        -(N/3) SIN (PHI) RJ(COS (PHI),1-K SIN (PHI), */

/*                                 2 */
/*                        1,1+N SIN (PHI)) */



/*                  Bulirsch form of ELLIPTIC INTEGRAL of 3rd kind */
/*                  ---------------------------------------------- */


/*                                           2 2    2 */
/*                  EL3(X,KC,P) = X RF(1,1+KC X ,1+X ) + */

/*                                            3          2 2    2     2 */
/*                               +(1/3)(1-P) X  RJ(1,1+KC X ,1+X ,1+PX ) */


/*                                           2 */
/*                  CEL(KC,P,A,B) = A RF(0,KC ,1) + */

/*                                                     2 */
/*                                 +(1/3)(B-PA) RJ(0,KC ,1,P) */




/*                  Heuman's LAMBDA function */
/*                  ------------------------ */


/*                                 2                     2      2    1/2 */
/*                  L(A,B,P) = (COS(A)SIN(B)COS(B)/(1-COS (A)SIN (B))   ) */

/*                                           2         2       2 */
/*                            *(SIN(P) RF(COS (P),1-SIN (A) SIN (P),1) */

/*                                 2       3            2       2 */
/*                            +(SIN (A) SIN (P)/(3(1-COS (A) SIN (B)))) */

/*                                   2         2       2 */
/*                            *RJ(COS (P),1-SIN (A) SIN (P),1,1- */

/*                                2       2          2       2 */
/*                            -SIN (A) SIN (P)/(1-COS (A) SIN (B)))) */




/*                  (PI/2) LAMBDA0(A,B) =L(A,B,PI/2) = */


/*                    2                         2       2    -1/2 */
/*               = COS (A)  SIN(B) COS(B) (1-COS (A) SIN (B)) */

/*                           2                  2       2 */
/*                  *RF(0,COS (A),1) + (1/3) SIN (A) COS (A) */

/*                                       2       2    -3/2 */
/*                  *SIN(B) COS(B) (1-COS (A) SIN (B)) */

/*                           2         2       2          2       2 */
/*                  *RJ(0,COS (A),1,COS (A) COS (B)/(1-COS (A) SIN (B))) */



/*                  Jacobi ZETA function */
/*                  -------------------- */


/*                             2                     2   2    1/2 */
/*                  Z(B,K) = (K/3) SIN(B) COS(B) (1-K SIN (B)) */


/*                                      2      2   2                2 */
/*                             *RJ(0,1-K ,1,1-K SIN (B)) / RF (0,1-K ,1) */


/*    ------------------------------------------------------------------- */

/* ***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete */
/*                 elliptic integrals, ACM Transactions on Mathematical */
/*                 Software 7, 3 (September 1981), pp. 398-403. */
/*               B. C. Carlson, Computing elliptic integrals by */
/*                 duplication, Numerische Mathematik 33, (1979), */
/*                 pp. 1-16. */
/*               B. C. Carlson, Elliptic integrals of the first kind, */
/*                 SIAM Journal of Mathematical Analysis 8, (1977), */
/*                 pp. 231-242. */
/* ***ROUTINES CALLED  R1MACH, RC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced statement labels.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900510  Changed calls to XERMSG to standard form, and some */
/*           editorial changes.  (RWC)). */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RJ */

/* ***FIRST EXECUTABLE STATEMENT  RJ */
    if (first) {
	d__1 = (doublereal) (r1mach_(&c__3) / 3.f);
	errtol = pow_dd(&d__1, &c_b3);
	d__1 = (doublereal) (r1mach_(&c__1) * 5.f);
	lolim = pow_dd(&d__1, &c_b5);
	d__1 = (doublereal) (r1mach_(&c__2) / 5.f);
	uplim = pow_dd(&d__1, &c_b5) * .3f;

	c1 = .21428571428571427f;
	c2 = .33333333333333331f;
	c3 = .13636363636363635f;
	c4 = .11538461538461539f;
    }
    first = FALSE_;

/*         CALL ERROR HANDLER IF NECESSARY. */

    ret_val = 0.f;
/* Computing MIN */
    r__1 = min(*x,*y);
    if (dmin(r__1,*z__) < 0.f) {
	*ier = 1;
	s_wsfi(&io___10);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___12);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 26, a__1[0] = "MIN(X,Y,Z).LT.0 WHERE X = ";
	i__1[1] = 16, a__1[1] = xern3;
	i__1[2] = 5, a__1[2] = " Y = ";
	i__1[3] = 16, a__1[3] = xern4;
	i__1[4] = 9, a__1[4] = " AND Z = ";
	i__1[5] = 16, a__1[5] = xern5;
	s_cat(ch__1, a__1, i__1, &c__6, (ftnlen)88);
	xermsg_("SLATEC", "RJ", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)88);
	return ret_val;
    }

/* Computing MAX */
    r__1 = max(*x,*y), r__1 = max(r__1,*z__);
    if (dmax(r__1,*p) > uplim) {
	*ier = 3;
	s_wsfi(&io___15);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___16);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___17);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___19);
	do_fio(&c__1, (char *)&(*p), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___21);
	do_fio(&c__1, (char *)&uplim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 32, a__2[0] = "MAX(X,Y,Z,P).GT.UPLIM WHERE X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 5, a__2[2] = " Y = ";
	i__2[3] = 16, a__2[3] = xern4;
	i__2[4] = 5, a__2[4] = " Z = ";
	i__2[5] = 16, a__2[5] = xern5;
	i__2[6] = 5, a__2[6] = " P = ";
	i__2[7] = 16, a__2[7] = xern6;
	i__2[8] = 13, a__2[8] = " AND UPLIM = ";
	i__2[9] = 16, a__2[9] = xern7;
	s_cat(ch__2, a__2, i__2, &c__10, (ftnlen)140);
	xermsg_("SLATEC", "RJ", ch__2, &c__3, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)140);
	return ret_val;
    }

/* Computing MIN */
    r__1 = *x + *y, r__2 = *x + *z__, r__1 = min(r__1,r__2), r__2 = *y + *z__,
	     r__1 = min(r__1,r__2);
    if (dmin(r__1,*p) < lolim) {
	*ier = 2;
	s_wsfi(&io___22);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___23);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___24);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___25);
	do_fio(&c__1, (char *)&(*p), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___26);
	do_fio(&c__1, (char *)&lolim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__3[0] = 38, a__3[0] = "MIN(X+Y,X+Z,Y+Z,P).LT.LOLIM WHERE X = ";
	i__3[1] = 16, a__3[1] = xern3;
	i__3[2] = 5, a__3[2] = " Y = ";
	i__3[3] = 16, a__3[3] = xern4;
	i__3[4] = 5, a__3[4] = " Z = ";
	i__3[5] = 16, a__3[5] = xern5;
	i__3[6] = 5, a__3[6] = " P = ";
	i__3[7] = 16, a__3[7] = xern6;
	i__3[8] = 13, a__3[8] = " AND LOLIM = ";
	s_cat(ch__3, a__3, i__3, &c__9, (ftnlen)130);
	xermsg_("SLATEC", "RJ", ch__3, &c__2, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)130);
	return ret_val;
    }

    *ier = 0;
    xn = *x;
    yn = *y;
    zn = *z__;
    pn = *p;
    sigma = 0.f;
    power4 = 1.f;

L30:
    mu = (xn + yn + zn + pn + pn) * .2f;
    xndev = (mu - xn) / mu;
    yndev = (mu - yn) / mu;
    zndev = (mu - zn) / mu;
    pndev = (mu - pn) / mu;
/* Computing MAX */
    r__1 = dabs(xndev), r__2 = dabs(yndev), r__1 = max(r__1,r__2), r__2 = 
	    dabs(zndev), r__1 = max(r__1,r__2), r__2 = dabs(pndev);
    epslon = dmax(r__1,r__2);
    if (epslon < errtol) {
	goto L40;
    }
    xnroot = sqrt(xn);
    ynroot = sqrt(yn);
    znroot = sqrt(zn);
    lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
    alfa = pn * (xnroot + ynroot + znroot) + xnroot * ynroot * znroot;
    alfa *= alfa;
    beta = pn * (pn + lamda) * (pn + lamda);
    sigma += power4 * rc_(&alfa, &beta, ier);
    power4 *= .25f;
    xn = (xn + lamda) * .25f;
    yn = (yn + lamda) * .25f;
    zn = (zn + lamda) * .25f;
    pn = (pn + lamda) * .25f;
    goto L30;

L40:
    ea = xndev * (yndev + zndev) + yndev * zndev;
    eb = xndev * yndev * zndev;
    ec = pndev * pndev;
    e2 = ea - ec * 3.f;
    e3 = eb + pndev * 2.f * (ea - ec);
    s1 = e2 * (-c1 + c3 * .75f * e2 - c4 * 1.5f * e3) + 1.f;
    s2 = eb * (c2 * .5f + pndev * (-c3 - c3 + pndev * c4));
    s3 = pndev * ea * (c2 - pndev * c3) - c2 * pndev * ec;
    ret_val = sigma * 3.f + power4 * (s1 + s2 + s3) / (mu * sqrt(mu));
    return ret_val;
} /* rj_ */


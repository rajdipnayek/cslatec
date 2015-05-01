/* rd.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static doublereal c_b5 = .66666666666666663;
static integer c__1 = 1;
static doublereal c_b7 = .33333333333333331;
static doublereal c_b9 = 2.;
static integer c__4 = 4;
static integer c__8 = 8;

/* DECK RD */
doublereal rd_(real *x, real *y, real *z__, integer *ier)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[4], a__2[8];
    integer i__1[4], i__2[8];
    real ret_val, r__1, r__2;
    doublereal d__1;
    char ch__1[65], ch__2[117];

    /* Local variables */
    static real c1, c2, c3, c4, s1, s2, ea, eb, ec, ed, ef, mu, xn, yn, zn;
    static char xern3[16], xern4[16], xern5[16], xern6[16];
    static real lamda, sigma, lolim, xndev, yndev, uplim, zndev;
    extern doublereal r1mach_(integer *);
    static real power4, epslon;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real errtol, tuplim, xnroot, ynroot, znroot;

    /* Fortran I/O blocks */
    static icilist io___11 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___13 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___17 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___19 = { 0, xern6, 0, "(1PE15.6)", 16, 1 };
    static icilist io___20 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___21 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___22 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___23 = { 0, xern6, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  RD */
/* ***PURPOSE  Compute the incomplete or complete elliptic integral of the */
/*            2nd kind.  For X and Y nonnegative, X+Y and Z positive, */
/*             RD(X,Y,Z) = Integral from zero to infinity of */
/*                                -1/2     -1/2     -3/2 */
/*                      (3/2)(t+X)    (t+Y)    (t+Z)    dt. */
/*            If X or Y is zero, the integral is complete. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C14 */
/* ***TYPE      SINGLE PRECISION (RD-S, DRD-D) */
/* ***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM, */
/*             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND, */
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

/*   1.     RD */
/*          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL */
/*          of the second kind */
/*          Standard FORTRAN function routine */
/*          Single precision version */
/*          The routine calculates an approximation result to */
/*          RD(X,Y,Z) = Integral from zero to infinity of */
/*                              -1/2     -1/2     -3/2 */
/*                    (3/2)(t+X)    (t+Y)    (t+Z)    dt, */
/*          where X and Y are nonnegative, X + Y is positive, and Z is */
/*          positive.  If X or Y is zero, the integral is COMPLETE. */
/*          The duplication theorem is iterated until the variables are */
/*          nearly equal, and the function is then expanded in Taylor */
/*          series to fifth order. */

/*   2.     Calling Sequence */

/*          RD( X, Y, Z, IER ) */

/*          Parameters on Entry */
/*          Values assigned by the calling routine */

/*          X      - Single precision, nonnegative variable */

/*          Y      - Single precision, nonnegative variable */

/*                   X + Y is positive */

/*          Z      - Real, positive variable */



/*          On Return     (values assigned by the RD routine) */

/*          RD     - Real approximation to the integral */


/*          IER    - Integer */

/*                   IER = 0 Normal and reliable termination of the */
/*                           routine.  It is assumed that the requested */
/*                           accuracy has been achieved. */

/*                   IER >  0 Abnormal termination of the routine */


/*          X, Y, Z are unaltered. */

/*   3.    Error Messages */

/*         Value of IER assigned by the RD routine */

/*                  Value Assigned         Error Message Printed */
/*                  IER = 1                MIN(X,Y) .LT. 0.0E0 */
/*                      = 2                MIN(X + Y, Z ) .LT. LOLIM */
/*                      = 3                MAX(X,Y,Z) .GT. UPLIM */


/*   4.     Control Parameters */

/*                  Values of LOLIM, UPLIM, and ERRTOL are set by the */
/*                  routine. */

/*          LOLIM and UPLIM determine the valid range of X, Y, and Z */

/*          LOLIM  - Lower limit of valid arguments */

/*                    Not less  than 2 / (machine maximum) ** (2/3). */

/*          UPLIM  - Upper limit of valid arguments */

/*                    Not greater than (0.1E0 * ERRTOL / machine */
/*                    minimum) ** (2/3), where ERRTOL is described below. */
/*                    In the following table it is assumed that ERRTOL */
/*                    will never be chosen smaller than 1.0E-5. */


/*                    Acceptable Values For:   LOLIM      UPLIM */
/*                    IBM 360/370 SERIES   :   6.0E-51     1.0E+48 */
/*                    CDC 6000/7000 SERIES :   5.0E-215    2.0E+191 */
/*                    UNIVAC 1100 SERIES   :   1.0E-25     2.0E+21 */
/*                    CRAY                 :   3.0E-1644   1.69E+1640 */
/*                    VAX 11 SERIES        :   1.0E-25     4.5E+21 */


/*          ERRTOL determines the accuracy of the answer */

/*                 The value assigned by the routine will result */
/*                 in solution precision within 1-2 decimals of */
/*                 "machine precision". */

/*          ERRTOL    Relative error due to truncation is less than */
/*                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2. */



/*              The accuracy of the computed approximation to the inte- */
/*              gral can be controlled by choosing the value of ERRTOL. */
/*              Truncation of a Taylor series after terms of fifth order */
/*              introduces an error less than the amount shown in the */
/*              second column of the following table for each value of */
/*              ERRTOL in the first column.  In addition to the trunca- */
/*              tion error there will be round-off error, but in prac- */
/*              tice the total error from both sources is usually less */
/*              than the amount given in the table. */




/*          Sample Choices:  ERRTOL   Relative Truncation */
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

/*   RD Special Comments */



/*          Check: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) */
/*          = 3 /  SQRT(X * Y * Z), where X, Y, and Z are positive. */


/*          On Input: */

/*          X, Y, and Z are the variables in the integral RD(X,Y,Z). */


/*          On Output: */


/*          X, Y, and Z are unaltered. */



/*          ******************************************************** */

/*           WARNING: Changes in the program may improve speed at the */
/*                    expense of robustness. */



/*    ------------------------------------------------------------------- */


/*   Special Functions via RD and RF */


/*                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind */
/*                  ---------------------------------------------- */


/*                                            2         2   2 */
/*                  E(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1) - */

/*                     2      3            2         2   2 */
/*                  -(K/3) SIN (PHI) RD(COS (PHI),1-K SIN (PHI),1) */


/*                                 2        2           2 */
/*                  E(K) = RF(0,1-K ,1) - (K/3) RD(0,1-K ,1) */


/*                         PI/2     2   2      1/2 */
/*                       = INT  (1-K SIN (PHI) )  D PHI */
/*                          0 */



/*                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind */
/*                  ---------------------------------------------- */

/*                                              2 2    2 */
/*                  EL2(X,KC,A,B) = AX RF(1,1+KC X ,1+X ) + */

/*                                              3         2 2    2 */
/*                                 +(1/3)(B-A) X RD(1,1+KC X ,1+X ) */



/*                  Legendre form of alternative ELLIPTIC INTEGRAL of 2nd */
/*                  ----------------------------------------------------- */
/*                        kind */
/*                        ---- */

/*                            Q     2       2   2  -1/2 */
/*                  D(Q,K) = INT SIN P  (1-K SIN P)     DP */
/*                            0 */



/*                                   3          2     2   2 */
/*                  D(Q,K) =(1/3)(SIN Q)  RD(COS Q,1-K SIN Q,1) */





/*                  Lemniscate constant B */
/*                  --------------------- */



/*                       1    2    4 -1/2 */
/*                  B = INT  S (1-S )    DS */
/*                       0 */


/*                  B =(1/3)RD (0,2,1) */




/*                  Heuman's LAMBDA function */
/*                  ------------------------ */



/*                  (PI/2) LAMBDA0(A,B) = */

/*                                       2                2 */
/*                     = SIN(B) (RF(0,COS (A),1)-(1/3) SIN (A) * */

/*                               2              2         2       2 */
/*                      *RD(0,COS (A),1)) RF(COS (B),1-COS (A) SIN (B),1) */

/*                               2       3            2 */
/*                     -(1/3) COS (A) SIN (B) RF(0,COS (A),1) * */

/*                             2         2       2 */
/*                      *RD(COS (B),1-COS (A) SIN (B),1) */



/*                  Jacobi ZETA function */
/*                  -------------------- */


/*                             2                2       2   2 */
/*                  Z(B,K) = (K/3) SIN(B) RF(COS (B),1-K SIN (B),1) */


/*                                      2            2 */
/*                             *RD(0,1-K ,1)/RF(0,1-K ,1) */

/*                               2       3          2       2   2 */
/*                            -(K /3) SIN (B) RD(COS (B),1-K SIN (B),1) */


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
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900510  Modify calls to XERMSG to put in standard form.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RD */

/* ***FIRST EXECUTABLE STATEMENT  RD */
    if (first) {
	d__1 = (doublereal) (r1mach_(&c__3) / 3.f);
	errtol = pow_dd(&d__1, &c_b3);
	d__1 = (doublereal) r1mach_(&c__2);
	lolim = 2.f / pow_dd(&d__1, &c_b5);
	d__1 = (doublereal) r1mach_(&c__1);
	tuplim = pow_dd(&d__1, &c_b7);
	d__1 = (doublereal) (errtol * .1f);
	tuplim = pow_dd(&d__1, &c_b7) / tuplim;
	d__1 = (doublereal) tuplim;
	uplim = pow_dd(&d__1, &c_b9);

	c1 = .21428571428571427f;
	c2 = .16666666666666666f;
	c3 = .40909090909090912f;
	c4 = .11538461538461539f;
    }
    first = FALSE_;

/*         CALL ERROR HANDLER IF NECESSARY. */

    ret_val = 0.f;
    if (dmin(*x,*y) < 0.f) {
	*ier = 1;
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___13);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 24, a__1[0] = "MIN(X,Y).LT.0 WHERE X = ";
	i__1[1] = 16, a__1[1] = xern3;
	i__1[2] = 9, a__1[2] = " AND Y = ";
	i__1[3] = 16, a__1[3] = xern4;
	s_cat(ch__1, a__1, i__1, &c__4, (ftnlen)65);
	xermsg_("SLATEC", "RD", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)65);
	return ret_val;
    }

/* Computing MAX */
    r__1 = max(*x,*y);
    if (dmax(r__1,*z__) > uplim) {
	*ier = 3;
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___15);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___17);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___19);
	do_fio(&c__1, (char *)&uplim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 30, a__2[0] = "MAX(X,Y,Z).GT.UPLIM WHERE X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 5, a__2[2] = " Y = ";
	i__2[3] = 16, a__2[3] = xern4;
	i__2[4] = 5, a__2[4] = " Z = ";
	i__2[5] = 16, a__2[5] = xern5;
	i__2[6] = 13, a__2[6] = " AND UPLIM = ";
	i__2[7] = 16, a__2[7] = xern6;
	s_cat(ch__2, a__2, i__2, &c__8, (ftnlen)117);
	xermsg_("SLATEC", "RD", ch__2, &c__3, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)117);
	return ret_val;
    }

/* Computing MIN */
    r__1 = *x + *y;
    if (dmin(r__1,*z__) < lolim) {
	*ier = 2;
	s_wsfi(&io___20);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___21);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___22);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___23);
	do_fio(&c__1, (char *)&lolim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 30, a__2[0] = "MIN(X+Y,Z).LT.LOLIM WHERE X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 5, a__2[2] = " Y = ";
	i__2[3] = 16, a__2[3] = xern4;
	i__2[4] = 5, a__2[4] = " Z = ";
	i__2[5] = 16, a__2[5] = xern5;
	i__2[6] = 13, a__2[6] = " AND LOLIM = ";
	i__2[7] = 16, a__2[7] = xern6;
	s_cat(ch__2, a__2, i__2, &c__8, (ftnlen)117);
	xermsg_("SLATEC", "RD", ch__2, &c__2, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)117);
	return ret_val;
    }

    *ier = 0;
    xn = *x;
    yn = *y;
    zn = *z__;
    sigma = 0.f;
    power4 = 1.f;

L30:
    mu = (xn + yn + zn * 3.f) * .2f;
    xndev = (mu - xn) / mu;
    yndev = (mu - yn) / mu;
    zndev = (mu - zn) / mu;
/* Computing MAX */
    r__1 = dabs(xndev), r__2 = dabs(yndev), r__1 = max(r__1,r__2), r__2 = 
	    dabs(zndev);
    epslon = dmax(r__1,r__2);
    if (epslon < errtol) {
	goto L40;
    }
    xnroot = sqrt(xn);
    ynroot = sqrt(yn);
    znroot = sqrt(zn);
    lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
    sigma += power4 / (znroot * (zn + lamda));
    power4 *= .25f;
    xn = (xn + lamda) * .25f;
    yn = (yn + lamda) * .25f;
    zn = (zn + lamda) * .25f;
    goto L30;

L40:
    ea = xndev * yndev;
    eb = zndev * zndev;
    ec = ea - eb;
    ed = ea - eb * 6.f;
    ef = ed + ec + ec;
    s1 = ed * (-c1 + c3 * .25f * ed - c4 * 1.5f * zndev * ef);
    s2 = zndev * (c2 * ef + zndev * (-c3 * ec + zndev * c4 * ea));
    ret_val = sigma * 3.f + power4 * (s1 + 1.f + s2) / (mu * sqrt(mu));

    return ret_val;
} /* rd_ */


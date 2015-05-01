/* rf.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__6 = 6;
static integer c__8 = 8;

/* DECK RF */
doublereal rf_(real *x, real *y, real *z__, integer *ier)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[6], a__2[8];
    integer i__1[6], i__2[8];
    real ret_val, r__1, r__2;
    doublereal d__1;
    char ch__1[88], ch__2[117], ch__3[123];

    /* Local variables */
    static real s, c1, c2, c3, e2, e3, mu, xn, yn, zn;
    static char xern3[16], xern4[16], xern5[16], xern6[16];
    static real lamda, lolim, xndev, yndev, uplim, zndev;
    extern doublereal r1mach_(integer *);
    static real epslon;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real errtol, xnroot, ynroot, znroot;

    /* Fortran I/O blocks */
    static icilist io___9 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___11 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___13 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___16 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___18 = { 0, xern6, 0, "(1PE15.6)", 16, 1 };
    static icilist io___19 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___20 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___21 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___22 = { 0, xern6, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  RF */
/* ***PURPOSE  Compute the incomplete or complete elliptic integral of the */
/*            1st kind.  For X, Y, and Z non-negative and at most one of */
/*            them zero, RF(X,Y,Z) = Integral from zero to infinity of */
/*                                -1/2     -1/2     -1/2 */
/*                      (1/2)(t+X)    (t+Y)    (t+Z)    dt. */
/*            If X, Y or Z is zero, the integral is complete. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C14 */
/* ***TYPE      SINGLE PRECISION (RF-S, DRF-D) */
/* ***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM, */
/*             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND, */
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

/*   1.     RF */
/*          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL */
/*          of the first kind */
/*          Standard FORTRAN function routine */
/*          Single precision version */
/*          The routine calculates an approximation result to */
/*          RF(X,Y,Z) = Integral from zero to infinity of */

/*                               -1/2     -1/2     -1/2 */
/*                     (1/2)(t+X)    (t+Y)    (t+Z)    dt, */

/*          where X, Y, and Z are nonnegative and at most one of them */
/*          is zero.  If one of them is zero, the integral is COMPLETE. */
/*          The duplication theorem is iterated until the variables are */
/*          nearly equal, and the function is then expanded in Taylor */
/*          series to fifth order. */

/*   2.     Calling Sequence */
/*          RF( X, Y, Z, IER ) */

/*          Parameters on Entry */
/*          Values assigned by the calling routine */

/*          X      - Single precision, nonnegative variable */

/*          Y      - Single precision, nonnegative variable */

/*          Z      - Single precision, nonnegative variable */



/*          On Return     (values assigned by the RF routine) */

/*          RF     - Single precision approximation to the integral */

/*          IER    - Integer */

/*                   IER = 0 Normal and reliable termination of the */
/*                           routine.  It is assumed that the requested */
/*                           accuracy has been achieved. */

/*                   IER >  0 Abnormal termination of the routine */

/*          X, Y, Z are unaltered. */


/*   3.    Error Messages */

/*         Value of IER assigned by the RF routine */

/*                  Value assigned         Error Message Printed */
/*                  IER = 1                MIN(X,Y,Z) .LT. 0.0E0 */
/*                      = 2                MIN(X+Y,X+Z,Y+Z) .LT. LOLIM */
/*                      = 3                MAX(X,Y,Z) .GT. UPLIM */



/*   4.     Control Parameters */

/*                  Values of LOLIM, UPLIM, and ERRTOL are set by the */
/*                  routine. */

/*          LOLIM and UPLIM determine the valid range of X, Y and Z */

/*          LOLIM  - Lower limit of valid arguments */

/*                   Not less than 5 * (machine minimum). */

/*          UPLIM  - Upper limit of valid arguments */

/*                   Not greater than (machine maximum) / 5. */


/*                     Acceptable Values For:   LOLIM      UPLIM */
/*                     IBM 360/370 SERIES   :   3.0E-78     1.0E+75 */
/*                     CDC 6000/7000 SERIES :   1.0E-292    1.0E+321 */
/*                     UNIVAC 1100 SERIES   :   1.0E-37     1.0E+37 */
/*                     CRAY                 :   2.3E-2466   1.09E+2465 */
/*                     VAX 11 SERIES        :   1.5E-38     3.0E+37 */



/*          ERRTOL determines the accuracy of the answer */

/*                 The value assigned by the routine will result */
/*                 in solution precision within 1-2 decimals of */
/*                 "machine precision". */



/*          ERRTOL - Relative error due to truncation is less than */
/*                   ERRTOL ** 6 / (4 * (1-ERRTOL)  . */



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
/*                           1.0E-3    3.0E-19 */
/*                           3.0E-3    2.0E-16 */
/*                           1.0E-2    3.0E-13 */
/*                           3.0E-2    2.0E-10 */
/*                           1.0E-1    3.0E-7 */


/*                    Decreasing ERRTOL by a factor of 10 yields six more */
/*                    decimal digits of accuracy at the expense of one or */
/*                    two more iterations of the duplication theorem. */

/* *Long Description: */

/*   RF Special Comments */



/*          Check by addition theorem: RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W) */
/*          = RF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W. */


/*          On Input: */

/*          X, Y, and Z are the variables in the integral RF(X,Y,Z). */


/*          On Output: */


/*          X, Y, and Z are unaltered. */



/*          ******************************************************** */

/*          Warning: Changes in the program may improve speed at the */
/*                   expense of robustness. */



/*   Special Functions via RF */


/*                  Legendre form of ELLIPTIC INTEGRAL of 1st kind */
/*                  ---------------------------------------------- */


/*                                            2         2   2 */
/*                  F(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1) */


/*                                 2 */
/*                  K(K) = RF(0,1-K ,1) */

/*                         PI/2     2   2      -1/2 */
/*                       = INT  (1-K SIN (PHI) )   D PHI */
/*                          0 */





/*                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind */
/*                  ---------------------------------------------- */


/*                                         2 2    2 */
/*                  EL1(X,KC) = X RF(1,1+KC X ,1+X ) */




/*                  Lemniscate constant A */
/*                  --------------------- */


/*                       1      4 -1/2 */
/*                  A = INT (1-S )    DS = RF(0,1,2) = RF(0,2,1) */
/*                       0 */


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
/*   891009  Removed unreferenced statement labels.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900510  Changed calls to XERMSG to standard form, and some */
/*           editorial changes.  (RWC)) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RF */

/* ***FIRST EXECUTABLE STATEMENT  RF */

    if (first) {
	d__1 = (doublereal) (r1mach_(&c__3) * 4.f);
	errtol = pow_dd(&d__1, &c_b3);
	lolim = r1mach_(&c__1) * 5.f;
	uplim = r1mach_(&c__2) / 5.f;

	c1 = .041666666666666664f;
	c2 = .068181818181818177f;
	c3 = .071428571428571425f;
    }
    first = FALSE_;

/*         CALL ERROR HANDLER IF NECESSARY. */

    ret_val = 0.f;
/* Computing MIN */
    r__1 = min(*x,*y);
    if (dmin(r__1,*z__) < 0.f) {
	*ier = 1;
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___13);
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
	xermsg_("SLATEC", "RF", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)88);
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
	s_wsfi(&io___16);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___18);
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
	xermsg_("SLATEC", "RF", ch__2, &c__3, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)117);
	return ret_val;
    }

/* Computing MIN */
    r__1 = *x + *y, r__2 = *x + *z__, r__1 = min(r__1,r__2), r__2 = *y + *z__;
    if (dmin(r__1,r__2) < lolim) {
	*ier = 2;
	s_wsfi(&io___19);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___20);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___21);
	do_fio(&c__1, (char *)&(*z__), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___22);
	do_fio(&c__1, (char *)&lolim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 36, a__2[0] = "MIN(X+Y,X+Z,Y+Z).LT.LOLIM WHERE X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 5, a__2[2] = " Y = ";
	i__2[3] = 16, a__2[3] = xern4;
	i__2[4] = 5, a__2[4] = " Z = ";
	i__2[5] = 16, a__2[5] = xern5;
	i__2[6] = 13, a__2[6] = " AND LOLIM = ";
	i__2[7] = 16, a__2[7] = xern6;
	s_cat(ch__3, a__2, i__2, &c__8, (ftnlen)123);
	xermsg_("SLATEC", "RF", ch__3, &c__2, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)123);
	return ret_val;
    }

    *ier = 0;
    xn = *x;
    yn = *y;
    zn = *z__;

L30:
    mu = (xn + yn + zn) / 3.f;
    xndev = 2.f - (mu + xn) / mu;
    yndev = 2.f - (mu + yn) / mu;
    zndev = 2.f - (mu + zn) / mu;
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
    xn = (xn + lamda) * .25f;
    yn = (yn + lamda) * .25f;
    zn = (zn + lamda) * .25f;
    goto L30;

L40:
    e2 = xndev * yndev - zndev * zndev;
    e3 = xndev * yndev * zndev;
    s = (c1 * e2 - .1f - c2 * e3) * e2 + 1.f + c3 * e3;
    ret_val = s / sqrt(mu);

    return ret_val;
} /* rf_ */


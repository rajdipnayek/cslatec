/* rc.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;
static integer c__6 = 6;

/* DECK RC */
doublereal rc_(real *x, real *y, integer *ier)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[4], a__2[6];
    integer i__1[4], i__2[6];
    real ret_val;
    doublereal d__1;
    char ch__1[70], ch__2[94], ch__3[89];

    /* Local variables */
    static real s, c1, c2, sn, mu, xn, yn;
    static char xern3[16], xern4[16], xern5[16];
    static real lamda, lolim, uplim;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real errtol;

    /* Fortran I/O blocks */
    static icilist io___8 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___10 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___11 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___12 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___16 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___17 = { 0, xern5, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  RC */
/* ***PURPOSE  Calculate an approximation to */
/*             RC(X,Y) = Integral from zero to infinity of */
/*                              -1/2     -1 */
/*                    (1/2)(t+X)    (t+Y)  dt, */
/*            where X is nonnegative and Y is positive. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C14 */
/* ***TYPE      SINGLE PRECISION (RC-S, DRC-D) */
/* ***KEYWORDS  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS, */
/*             ELLIPTIC INTEGRAL, TAYLOR SERIES */
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

/*   1.     RC */
/*          Standard FORTRAN function routine */
/*          Single precision version */
/*          The routine calculates an approximation to */
/*           RC(X,Y) = Integral from zero to infinity of */

/*                              -1/2     -1 */
/*                    (1/2)(t+X)    (t+Y)  dt, */

/*          where X is nonnegative and Y is positive.  The duplication */
/*          theorem is iterated until the variables are nearly equal, */
/*          and the function is then expanded in Taylor series to fifth */
/*          order.  Logarithmic, inverse circular, and inverse hyper- */
/*          bolic functions can be expressed in terms of RC. */


/*   2.     Calling Sequence */
/*          RC( X, Y, IER ) */

/*          Parameters on Entry */
/*          Values assigned by the calling routine */

/*          X      - Single precision, nonnegative variable */

/*          Y      - Single precision, positive variable */



/*          On Return  (values assigned by the RC routine) */

/*          RC     - Single precision approximation to the integral */

/*          IER    - Integer to indicate normal or abnormal termination. */

/*                     IER = 0 Normal and reliable termination of the */
/*                             routine.  It is assumed that the requested */
/*                             accuracy has been achieved. */

/*                     IER > 0 Abnormal termination of the routine */

/*          X and Y are unaltered. */


/*   3.    Error Messages */

/*         Value of IER assigned by the RC routine */

/*                  Value Assigned         Error Message Printed */
/*                  IER = 1                X.LT.0.0E0.OR.Y.LE.0.0E0 */
/*                      = 2                X+Y.LT.LOLIM */
/*                      = 3                MAX(X,Y) .GT. UPLIM */


/*   4.     Control Parameters */

/*                  Values of LOLIM, UPLIM, and ERRTOL are set by the */
/*                  routine. */

/*          LOLIM and UPLIM determine the valid range of X and Y */

/*          LOLIM  - Lower limit of valid arguments */

/*                   Not less  than 5 * (machine minimum)  . */

/*          UPLIM  - Upper limit of valid arguments */

/*                   Not greater than (machine maximum) / 5 . */


/*                     Acceptable values for:   LOLIM       UPLIM */
/*                     IBM 360/370 SERIES   :   3.0E-78     1.0E+75 */
/*                     CDC 6000/7000 SERIES :   1.0E-292    1.0E+321 */
/*                     UNIVAC 1100 SERIES   :   1.0E-37     1.0E+37 */
/*                     CRAY                 :   2.3E-2466   1.09E+2465 */
/*                     VAX 11 SERIES        :   1.5E-38     3.0E+37 */

/*          ERRTOL determines the accuracy of the answer */

/*                 The value assigned by the routine will result */
/*                 in solution precision within 1-2 decimals of */
/*                 "machine precision". */


/*          ERRTOL  - Relative error due to truncation is less than */
/*                    16 * ERRTOL ** 6 / (1 - 2 * ERRTOL). */


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
/*                           1.0E-3    2.0E-17 */
/*                           3.0E-3    2.0E-14 */
/*                           1.0E-2    2.0E-11 */
/*                           3.0E-2    2.0E-8 */
/*                           1.0E-1    2.0E-5 */


/*                    Decreasing ERRTOL by a factor of 10 yields six more */
/*                    decimal digits of accuracy at the expense of one or */
/*                    two more iterations of the duplication theorem. */

/* *Long Description: */

/*   RC Special Comments */




/*                  Check: RC(X,X+Z) + RC(Y,Y+Z) = RC(0,Z) */

/*                  where X, Y, and Z are positive and X * Y = Z * Z */


/*          On Input: */

/*          X and Y are the variables in the integral RC(X,Y). */

/*          On Output: */

/*          X and Y are unaltered. */



/*                    RC(0,1/4)=RC(1/16,1/8)=PI=3.14159... */

/*                    RC(9/4,2)=LN(2) */



/*          ******************************************************** */

/*          Warning: Changes in the program may improve speed at the */
/*                   expense of robustness. */


/*   -------------------------------------------------------------------- */

/*   Special Functions via RC */



/*                  LN X                X .GT. 0 */

/*                                            2 */
/*                  LN(X) = (X-1) RC(((1+X)/2)  , X ) */


/*   -------------------------------------------------------------------- */

/*                  ARCSIN X            -1 .LE. X .LE. 1 */

/*                                      2 */
/*                  ARCSIN X = X RC (1-X  ,1 ) */

/*   -------------------------------------------------------------------- */

/*                  ARCCOS X            0 .LE. X .LE. 1 */


/*                                     2      2 */
/*                  ARCCOS X = SQRT(1-X ) RC(X  ,1 ) */

/*   -------------------------------------------------------------------- */

/*                  ARCTAN X            -INF .LT. X .LT. +INF */

/*                                       2 */
/*                  ARCTAN X = X RC(1,1+X  ) */

/*   -------------------------------------------------------------------- */

/*                  ARCCOT X            0 .LE. X .LT. INF */

/*                                 2   2 */
/*                  ARCCOT X = RC(X  ,X +1 ) */

/*   -------------------------------------------------------------------- */

/*                  ARCSINH X           -INF .LT. X .LT. +INF */

/*                                      2 */
/*                  ARCSINH X = X RC(1+X  ,1 ) */

/*   -------------------------------------------------------------------- */

/*                  ARCCOSH X           X .GE. 1 */

/*                                    2        2 */
/*                  ARCCOSH X = SQRT(X -1) RC(X  ,1 ) */

/*   -------------------------------------------------------------------- */

/*                  ARCTANH X           -1 .LT. X .LT. 1 */

/*                                        2 */
/*                  ARCTANH X = X RC(1,1-X  ) */

/*   -------------------------------------------------------------------- */

/*                  ARCCOTH X           X .GT. 1 */

/*                                  2   2 */
/*                  ARCCOTH X = RC(X  ,X -1 ) */

/*   -------------------------------------------------------------------- */

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
/* ***END PROLOGUE  RC */

/* ***FIRST EXECUTABLE STATEMENT  RC */
    if (first) {
	d__1 = (doublereal) (r1mach_(&c__3) / 16.f);
	errtol = pow_dd(&d__1, &c_b3);
	lolim = r1mach_(&c__1) * 5.f;
	uplim = r1mach_(&c__2) / 5.f;

	c1 = .14285714285714285f;
	c2 = .40909090909090912f;
    }
    first = FALSE_;

/*         CALL ERROR HANDLER IF NECESSARY. */

    ret_val = 0.f;
    if (*x < 0.f || *y <= 0.f) {
	*ier = 1;
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___10);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 29, a__1[0] = "X.LT.0 .OR. Y.LE.0 WHERE X = ";
	i__1[1] = 16, a__1[1] = xern3;
	i__1[2] = 9, a__1[2] = " AND Y = ";
	i__1[3] = 16, a__1[3] = xern4;
	s_cat(ch__1, a__1, i__1, &c__4, (ftnlen)70);
	xermsg_("SLATEC", "RC", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)70);
	return ret_val;
    }

    if (dmax(*x,*y) > uplim) {
	*ier = 3;
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___12);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&uplim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 28, a__2[0] = "MAX(X,Y).GT.UPLIM WHERE X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 5, a__2[2] = " Y = ";
	i__2[3] = 16, a__2[3] = xern4;
	i__2[4] = 13, a__2[4] = " AND UPLIM = ";
	i__2[5] = 16, a__2[5] = xern5;
	s_cat(ch__2, a__2, i__2, &c__6, (ftnlen)94);
	xermsg_("SLATEC", "RC", ch__2, &c__3, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)94);
	return ret_val;
    }

    if (*x + *y < lolim) {
	*ier = 2;
	s_wsfi(&io___15);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___16);
	do_fio(&c__1, (char *)&(*y), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___17);
	do_fio(&c__1, (char *)&lolim, (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 23, a__2[0] = "X+Y.LT.LOLIM WHERE X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 5, a__2[2] = " Y = ";
	i__2[3] = 16, a__2[3] = xern4;
	i__2[4] = 13, a__2[4] = " AND LOLIM = ";
	i__2[5] = 16, a__2[5] = xern5;
	s_cat(ch__3, a__2, i__2, &c__6, (ftnlen)89);
	xermsg_("SLATEC", "RC", ch__3, &c__2, &c__1, (ftnlen)6, (ftnlen)2, (
		ftnlen)89);
	return ret_val;
    }

    *ier = 0;
    xn = *x;
    yn = *y;

L30:
    mu = (xn + yn + yn) / 3.f;
    sn = (yn + mu) / mu - 2.f;
    if (dabs(sn) < errtol) {
	goto L40;
    }
    lamda = sqrt(xn) * 2.f * sqrt(yn) + yn;
    xn = (xn + lamda) * .25f;
    yn = (yn + lamda) * .25f;
    goto L30;

L40:
    s = sn * sn * (sn * (c1 + sn * (sn * c2 + .375f)) + .3f);
    ret_val = (s + 1.f) / sqrt(mu);
    return ret_val;
} /* rc_ */


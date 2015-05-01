/* dulsia.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__8 = 8;

/* DECK DULSIA */
/* Subroutine */ int dulsia_(doublereal *a, integer *mda, integer *m, integer 
	*n, doublereal *b, integer *mdb, integer *nb, doublereal *re, 
	doublereal *ae, integer *key, integer *mode, integer *np, integer *
	krank, integer *ksure, doublereal *rnorm, doublereal *w, integer *lw, 
	integer *iwork, integer *liw, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static integer i__, m1, m2, m3, m4, m5, it;
    static doublereal eps;
    extern /* Subroutine */ int du11us_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *), du12us_(doublereal *, integer *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *);
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DULSIA */
/* ***PURPOSE  Solve an underdetermined linear system of equations by */
/*            performing an LQ factorization of the matrix using */
/*            Householder transformations.  Emphasis is put on detecting */
/*            possible rank deficiency. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D9 */
/* ***TYPE      DOUBLE PRECISION (ULSIA-S, DULSIA-D) */
/* ***KEYWORDS  LINEAR LEAST SQUARES, LQ FACTORIZATION, */
/*             UNDERDETERMINED LINEAR SYSTEM */
/* ***AUTHOR  Manteuffel, T. A., (LANL) */
/* ***DESCRIPTION */

/*     DULSIA computes the minimal length solution(s) to the problem AX=B */
/*     where A is an M by N matrix with M.LE.N and B is the M by NB */
/*     matrix of right hand sides.  User input bounds on the uncertainty */
/*     in the elements of A are used to detect numerical rank deficiency. */
/*     The algorithm employs a row and column pivot strategy to */
/*     minimize the growth of uncertainty and round-off errors. */

/*     DULSIA requires (MDA+1)*N + (MDB+1)*NB + 6*M dimensioned space */

/*   ****************************************************************** */
/*   *                                                                * */
/*   *         WARNING - All input arrays are changed on exit.        * */
/*   *                                                                * */
/*   ****************************************************************** */

/*     Input.. All TYPE REAL variables are DOUBLE PRECISION */

/*     A(,)          Linear coefficient matrix of AX=B, with MDA the */
/*      MDA,M,N      actual first dimension of A in the calling program. */
/*                   M is the row dimension (no. of EQUATIONS of the */
/*                   problem) and N the col dimension (no. of UNKNOWNS). */
/*                   Must have MDA.GE.M and M.LE.N. */

/*     B(,)          Right hand side(s), with MDB the actual first */
/*      MDB,NB       dimension of B in the calling program. NB is the */
/*                   number of M by 1 right hand sides.  Since the */
/*                   solution is returned in B, must have MDB.GE.N.  If */
/*                   NB = 0, B is never accessed. */

/*   ****************************************************************** */
/*   *                                                                * */
/*   *         Note - Use of RE and AE are what make this             * */
/*   *                code significantly different from               * */
/*   *                other linear least squares solvers.             * */
/*   *                However, the inexperienced user is              * */
/*   *                advised to set RE=0.,AE=0.,KEY=0.               * */
/*   *                                                                * */
/*   ****************************************************************** */

/*     RE(),AE(),KEY */
/*     RE()          RE() is a vector of length N such that RE(I) is */
/*                   the maximum relative uncertainty in row I of */
/*                   the matrix A. The values of RE() must be between */
/*                   0 and 1. A minimum of 10*machine precision will */
/*                   be enforced. */

/*     AE()          AE() is a vector of length N such that AE(I) is */
/*                   the maximum absolute uncertainty in row I of */
/*                   the matrix A. The values of AE() must be greater */
/*                   than or equal to 0. */

/*     KEY           For ease of use, RE and AE may be input as either */
/*                   vectors or scalars. If a scalar is input, the algo- */
/*                   rithm will use that value for each column of A. */
/*                   The parameter KEY indicates whether scalars or */
/*                   vectors are being input. */
/*                        KEY=0     RE scalar  AE scalar */
/*                        KEY=1     RE vector  AE scalar */
/*                        KEY=2     RE scalar  AE vector */
/*                        KEY=3     RE vector  AE vector */


/*     MODE          The integer MODE indicates how the routine */
/*                   is to react if rank deficiency is detected. */
/*                   If MODE = 0 return immediately, no solution */
/*                             1 compute truncated solution */
/*                             2 compute minimal length least squares sol */
/*                   The inexperienced user is advised to set MODE=0 */

/*     NP            The first NP rows of A will not be interchanged */
/*                   with other rows even though the pivot strategy */
/*                   would suggest otherwise. */
/*                   The inexperienced user is advised to set NP=0. */

/*     WORK()        A real work array dimensioned 5*M.  However, if */
/*                   RE or AE have been specified as vectors, dimension */
/*                   WORK 4*M. If both RE and AE have been specified */
/*                   as vectors, dimension WORK 3*M. */

/*     LW            Actual dimension of WORK */

/*     IWORK()       Integer work array dimensioned at least N+M. */

/*     LIW           Actual dimension of IWORK. */


/*     INFO          Is a flag which provides for the efficient */
/*                   solution of subsequent problems involving the */
/*                   same A but different B. */
/*                   If INFO = 0 original call */
/*                      INFO = 1 subsequent calls */
/*                   On subsequent calls, the user must supply A, KRANK, */
/*                   LW, IWORK, LIW, and the first 2*M locations of WORK */
/*                   as output by the original call to DULSIA. MODE must */
/*                   be equal to the value of MODE in the original call. */
/*                   If MODE.LT.2, only the first N locations of WORK */
/*                   are accessed. AE, RE, KEY, and NP are not accessed. */




/*     Output..All TYPE REAL variables are DOUBLE PRECISION */

/*     A(,)          Contains the lower triangular part of the reduced */
/*                   matrix and the transformation information. It togeth */
/*                   with the first M elements of WORK (see below) */
/*                   completely specify the LQ factorization of A. */

/*     B(,)          Contains the N by NB solution matrix for X. */

/*     KRANK,KSURE   The numerical rank of A,  based upon the relative */
/*                   and absolute bounds on uncertainty, is bounded */
/*                   above by KRANK and below by KSURE. The algorithm */
/*                   returns a solution based on KRANK. KSURE provides */
/*                   an indication of the precision of the rank. */

/*     RNORM()       Contains the Euclidean length of the NB residual */
/*                   vectors  B(I)-AX(I), I=1,NB. If the matrix A is of */
/*                   full rank, then RNORM=0.0. */

/*     WORK()        The first M locations of WORK contain values */
/*                   necessary to reproduce the Householder */
/*                   transformation. */

/*     IWORK()       The first N locations contain the order in */
/*                   which the columns of A were used. The next */
/*                   M locations contain the order in which the */
/*                   rows of A were used. */

/*     INFO          Flag to indicate status of computation on completion */
/*                  -1   Parameter error(s) */
/*                   0 - Rank deficient, no solution */
/*                   1 - Rank deficient, truncated solution */
/*                   2 - Rank deficient, minimal length least squares sol */
/*                   3 - Numerical rank 0, zero solution */
/*                   4 - Rank .LT. NP */
/*                   5 - Full rank */

/* ***REFERENCES  T. Manteuffel, An interval analysis approach to rank */
/*                 determination in linear least squares problems, */
/*                 Report SAND80-0655, Sandia Laboratories, June 1980. */
/* ***ROUTINES CALLED  D1MACH, DU11US, DU12US, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810801  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Fixed an error message.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DULSIA */

/* ***FIRST EXECUTABLE STATEMENT  DULSIA */
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *mdb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --re;
    --ae;
    --rnorm;
    --w;
    --iwork;

    /* Function Body */
    if (*info < 0 || *info > 1) {
	goto L514;
    }
    it = *info;
    *info = -1;
    if (*nb == 0 && it == 1) {
	goto L501;
    }
    if (*m < 1) {
	goto L502;
    }
    if (*n < 1) {
	goto L503;
    }
    if (*n < *m) {
	goto L504;
    }
    if (*mda < *m) {
	goto L505;
    }
    if (*liw < *m + *n) {
	goto L506;
    }
    if (*mode < 0 || *mode > 3) {
	goto L515;
    }
    if (*nb == 0) {
	goto L4;
    }
    if (*nb < 0) {
	goto L507;
    }
    if (*mdb < *n) {
	goto L508;
    }
    if (it == 0) {
	goto L4;
    }
    goto L400;
L4:
    if (*key < 0 || *key > 3) {
	goto L509;
    }
    if (*key == 0 && *lw < *m * 5) {
	goto L510;
    }
    if (*key == 1 && *lw < *m << 2) {
	goto L510;
    }
    if (*key == 2 && *lw < *m << 2) {
	goto L510;
    }
    if (*key == 3 && *lw < *m * 3) {
	goto L510;
    }
    if (*np < 0 || *np > *m) {
	goto L516;
    }

    eps = d1mach_(&c__3) * 10.f;
    m1 = 1;
    m2 = m1 + *m;
    m3 = m2 + *m;
    m4 = m3 + *m;
    m5 = m4 + *m;

    if (*key == 1) {
	goto L100;
    }
    if (*key == 2) {
	goto L200;
    }
    if (*key == 3) {
	goto L300;
    }

    if (re[1] < 0.) {
	goto L511;
    }
    if (re[1] > 1.) {
	goto L512;
    }
    if (re[1] < eps) {
	re[1] = eps;
    }
    if (ae[1] < 0.) {
	goto L513;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[m4 - 1 + i__] = re[1];
	w[m5 - 1 + i__] = ae[1];
/* L20: */
    }
    du11us_(&a[a_offset], mda, m, n, &w[m4], &w[m5], mode, np, krank, ksure, &
	    w[m1], &w[m2], &w[m3], &iwork[m1], &iwork[m2]);
    goto L400;

L100:
    if (ae[1] < 0.) {
	goto L513;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (re[i__] < 0.) {
	    goto L511;
	}
	if (re[i__] > 1.) {
	    goto L512;
	}
	if (re[i__] < eps) {
	    re[i__] = eps;
	}
	w[m4 - 1 + i__] = ae[1];
/* L120: */
    }
    du11us_(&a[a_offset], mda, m, n, &re[1], &w[m4], mode, np, krank, ksure, &
	    w[m1], &w[m2], &w[m3], &iwork[m1], &iwork[m2]);
    goto L400;

L200:
    if (re[1] < 0.) {
	goto L511;
    }
    if (re[1] > 1.) {
	goto L512;
    }
    if (re[1] < eps) {
	re[1] = eps;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[m4 - 1 + i__] = re[1];
	if (ae[i__] < 0.) {
	    goto L513;
	}
/* L220: */
    }
    du11us_(&a[a_offset], mda, m, n, &w[m4], &ae[1], mode, np, krank, ksure, &
	    w[m1], &w[m2], &w[m3], &iwork[m1], &iwork[m2]);
    goto L400;

L300:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (re[i__] < 0.) {
	    goto L511;
	}
	if (re[i__] > 1.) {
	    goto L512;
	}
	if (re[i__] < eps) {
	    re[i__] = eps;
	}
	if (ae[i__] < 0.) {
	    goto L513;
	}
/* L320: */
    }
    du11us_(&a[a_offset], mda, m, n, &re[1], &ae[1], mode, np, krank, ksure, &
	    w[m1], &w[m2], &w[m3], &iwork[m1], &iwork[m2]);

/*     DETERMINE INFO */

L400:
    if (*krank != *m) {
	goto L402;
    }
    *info = 5;
    goto L410;
L402:
    if (*krank != 0) {
	goto L404;
    }
    *info = 3;
    goto L410;
L404:
    if (*krank >= *np) {
	goto L406;
    }
    *info = 4;
    return 0;
L406:
    *info = *mode;
    if (*mode == 0) {
	return 0;
    }
L410:
    if (*nb == 0) {
	return 0;
    }


/*     SOLUTION PHASE */

    m1 = 1;
    m2 = m1 + *m;
    m3 = m2 + *m;
    if (*info == 2) {
	goto L420;
    }
    if (*lw < m2 - 1) {
	goto L510;
    }
    du12us_(&a[a_offset], mda, m, n, &b[b_offset], mdb, nb, mode, krank, &
	    rnorm[1], &w[m1], &w[m1], &iwork[m1], &iwork[m2]);
    return 0;

L420:
    if (*lw < m3 - 1) {
	goto L510;
    }
    du12us_(&a[a_offset], mda, m, n, &b[b_offset], mdb, nb, mode, krank, &
	    rnorm[1], &w[m1], &w[m2], &iwork[m1], &iwork[m2]);
    return 0;

/*     ERROR MESSAGES */

L501:
    xermsg_("SLATEC", "DULSIA", "SOLUTION ONLY (INFO=1) BUT NO RIGHT HAND SI"
	    "DE (NB=0)", &c__1, &c__0, (ftnlen)6, (ftnlen)6, (ftnlen)52);
    return 0;
L502:
    xermsg_("SLATEC", "DULSIA", "M.LT.1", &c__2, &c__1, (ftnlen)6, (ftnlen)6, 
	    (ftnlen)6);
    return 0;
L503:
    xermsg_("SLATEC", "DULSIA", "N.LT.1", &c__2, &c__1, (ftnlen)6, (ftnlen)6, 
	    (ftnlen)6);
    return 0;
L504:
    xermsg_("SLATEC", "DULSIA", "N.LT.M", &c__2, &c__1, (ftnlen)6, (ftnlen)6, 
	    (ftnlen)6);
    return 0;
L505:
    xermsg_("SLATEC", "DULSIA", "MDA.LT.M", &c__2, &c__1, (ftnlen)6, (ftnlen)
	    6, (ftnlen)8);
    return 0;
L506:
    xermsg_("SLATEC", "DULSIA", "LIW.LT.M+N", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)10);
    return 0;
L507:
    xermsg_("SLATEC", "DULSIA", "NB.LT.0", &c__2, &c__1, (ftnlen)6, (ftnlen)6,
	     (ftnlen)7);
    return 0;
L508:
    xermsg_("SLATEC", "DULSIA", "MDB.LT.N", &c__2, &c__1, (ftnlen)6, (ftnlen)
	    6, (ftnlen)8);
    return 0;
L509:
    xermsg_("SLATEC", "DULSIA", "KEY OUT OF RANGE", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)16);
    return 0;
L510:
    xermsg_("SLATEC", "DULSIA", "INSUFFICIENT WORK SPACE", &c__8, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)23);
    *info = -1;
    return 0;
L511:
    xermsg_("SLATEC", "DULSIA", "RE(I) .LT. 0", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)12);
    return 0;
L512:
    xermsg_("SLATEC", "DULSIA", "RE(I) .GT. 1", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)12);
    return 0;
L513:
    xermsg_("SLATEC", "DULSIA", "AE(I) .LT. 0", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)12);
    return 0;
L514:
    xermsg_("SLATEC", "DULSIA", "INFO OUT OF RANGE", &c__2, &c__1, (ftnlen)6, 
	    (ftnlen)6, (ftnlen)17);
    return 0;
L515:
    xermsg_("SLATEC", "DULSIA", "MODE OUT OF RANGE", &c__2, &c__1, (ftnlen)6, 
	    (ftnlen)6, (ftnlen)17);
    return 0;
L516:
    xermsg_("SLATEC", "DULSIA", "NP OUT OF RANGE", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)15);
    return 0;
} /* dulsia_ */


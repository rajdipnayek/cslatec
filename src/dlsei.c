/* dlsei.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__8 = 8;
static doublereal c_b41 = 1.;
static integer c__0 = 0;
static doublereal c_b95 = 0.;

/* DECK DLSEI */
/* Subroutine */ int dlsei_(doublereal *w, integer *mdw, integer *me, integer 
	*ma, integer *mg, integer *n, doublereal *prgopt, doublereal *x, 
	doublereal *rnorme, doublereal *rnorml, integer *mode, doublereal *ws,
	 integer *ip)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[8], a__2[2];
    integer w_dim1, w_offset, i__1, i__2[8], i__3[2], i__4, i__5;
    doublereal d__1, d__2;
    char ch__1[131], ch__2[60], ch__3[61];

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t;
    static integer n1, n2;
    static doublereal rb, uj, rn, sn, vj, up;
    static integer jp1, np1;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal gam;
    static integer key;
    static doublereal tau;
    static logical cov;
    static integer mep1, lchk, mend;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dlsi_(doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer link, imax, last;
    static doublereal size;
    static integer next, nopt;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static char xern1[8], xern2[8], xern3[8], xern4[8];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer mdeqc;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer nlink;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal enorm, fnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal rnmax, snmax, xnrme;
    extern doublereal d1mach_(integer *);
    static doublereal xnorm;
    static integer mapke1, kranke;
    static doublereal drelpr;
    static integer ntimes;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern3, 0, "(I8)", 8, 1 };
    static icilist io___11 = { 0, xern4, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___14 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DLSEI */
/* ***PURPOSE  Solve a linearly constrained least squares problem with */
/*            equality and inequality constraints, and optionally compute */
/*            a covariance matrix. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1A2A, D9 */
/* ***TYPE      DOUBLE PRECISION (LSEI-S, DLSEI-D) */
/* ***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING, */
/*             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS, */
/*             QUADRATIC PROGRAMMING */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */

/*     This subprogram solves a linearly constrained least squares */
/*     problem with both equality and inequality constraints, and, if the */
/*     user requests, obtains a covariance matrix of the solution */
/*     parameters. */

/*     Suppose there are given matrices E, A and G of respective */
/*     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of */
/*     respective lengths ME, MA and MG.  This subroutine solves the */
/*     linearly constrained least squares problem */

/*                   EX = F, (E ME by N) (equations to be exactly */
/*                                       satisfied) */
/*                   AX = B, (A MA by N) (equations to be */
/*                                       approximately satisfied, */
/*                                       least squares sense) */
/*                   GX .GE. H,(G MG by N) (inequality constraints) */

/*     The inequalities GX .GE. H mean that every component of the */
/*     product GX must be .GE. the corresponding component of H. */

/*     In case the equality constraints cannot be satisfied, a */
/*     generalized inverse solution residual vector length is obtained */
/*     for F-EX.  This is the minimal length possible for F-EX. */

/*     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The */
/*     rank of the matrix E is estimated during the computation.  We call */
/*     this value KRANKE.  It is an output parameter in IP(1) defined */
/*     below.  Using a generalized inverse solution of EX=F, a reduced */
/*     least squares problem with inequality constraints is obtained. */
/*     The tolerances used in these tests for determining the rank */
/*     of E and the rank of the reduced least squares problem are */
/*     given in Sandia Tech. Rept. SAND-78-1290.  They can be */
/*     modified by the user if new values are provided in */
/*     the option list of the array PRGOPT(*). */

/*     The user must dimension all arrays appearing in the call list.. */
/*     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2) */
/*     where K=MAX(MA+MG,N).  This allows for a solution of a range of */
/*     problems in the given working space.  The dimension of WS(*) */
/*     given is a necessary overestimate.  Once a particular problem */
/*     has been run, the output parameter IP(3) gives the actual */
/*     dimension required for that problem. */

/*     The parameters for DLSEI( ) are */

/*     Input.. All TYPE REAL variables are DOUBLE PRECISION */

/*     W(*,*),MDW,   The array W(*,*) is doubly subscripted with */
/*     ME,MA,MG,N    first dimensioning parameter equal to MDW. */
/*                   For this discussion let us call M = ME+MA+MG.  Then */
/*                   MDW must satisfy MDW .GE. M.  The condition */
/*                   MDW .LT. M is an error. */

/*                   The array W(*,*) contains the matrices and vectors */

/*                                  (E  F) */
/*                                  (A  B) */
/*                                  (G  H) */

/*                   in rows and columns 1,...,M and 1,...,N+1 */
/*                   respectively. */

/*                   The integers ME, MA, and MG are the */
/*                   respective matrix row dimensions */
/*                   of E, A and G.  Each matrix has N columns. */

/*     PRGOPT(*)    This real-valued array is the option vector. */
/*                  If the user is satisfied with the nominal */
/*                  subprogram features set */

/*                  PRGOPT(1)=1 (or PRGOPT(1)=1.0) */

/*                  Otherwise PRGOPT(*) is a linked list consisting of */
/*                  groups of data of the following form */

/*                  LINK */
/*                  KEY */
/*                  DATA SET */

/*                  The parameters LINK and KEY are each one word. */
/*                  The DATA SET can be comprised of several words. */
/*                  The number of items depends on the value of KEY. */
/*                  The value of LINK points to the first */
/*                  entry of the next group of data within */
/*                  PRGOPT(*).  The exception is when there are */
/*                  no more options to change.  In that */
/*                  case, LINK=1 and the values KEY and DATA SET */
/*                  are not referenced.  The general layout of */
/*                  PRGOPT(*) is as follows. */

/*               ...PRGOPT(1) = LINK1 (link to first entry of next group) */
/*               .  PRGOPT(2) = KEY1 (key to the option change) */
/*               .  PRGOPT(3) = data value (data value for this change) */
/*               .       . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of */
/*               .                       next group) */
/*               .  PRGOPT(LINK1+1) = KEY2 (key to the option change) */
/*               .  PRGOPT(LINK1+2) = data value */
/*               ...     . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK) = 1 (no more options to change) */

/*                  Values of LINK that are nonpositive are errors. */
/*                  A value of LINK .GT. NLINK=100000 is also an error. */
/*                  This helps prevent using invalid but positive */
/*                  values of LINK that will probably extend */
/*                  beyond the program limits of PRGOPT(*). */
/*                  Unrecognized values of KEY are ignored.  The */
/*                  order of the options is arbitrary and any number */
/*                  of options can be changed with the following */
/*                  restriction.  To prevent cycling in the */
/*                  processing of the option array, a count of the */
/*                  number of options changed is maintained. */
/*                  Whenever this count exceeds NOPT=1000, an error */
/*                  message is printed and the subprogram returns. */

/*                  Options.. */

/*                  KEY=1 */
/*                         Compute in W(*,*) the N by N */
/*                  covariance matrix of the solution variables */
/*                  as an output parameter.  Nominally the */
/*                  covariance matrix will not be computed. */
/*                  (This requires no user input.) */
/*                  The data set for this option is a single value. */
/*                  It must be nonzero when the covariance matrix */
/*                  is desired.  If it is zero, the covariance */
/*                  matrix is not computed.  When the covariance matrix */
/*                  is computed, the first dimensioning parameter */
/*                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N). */

/*                  KEY=10 */
/*                         Suppress scaling of the inverse of the */
/*                  normal matrix by the scale factor RNORM**2/ */
/*                  MAX(1, no. of degrees of freedom).  This option */
/*                  only applies when the option for computing the */
/*                  covariance matrix (KEY=1) is used.  With KEY=1 and */
/*                  KEY=10 used as options the unscaled inverse of the */
/*                  normal matrix is returned in W(*,*). */
/*                  The data set for this option is a single value. */
/*                  When it is nonzero no scaling is done.  When it is */
/*                  zero scaling is done.  The nominal case is to do */
/*                  scaling so if option (KEY=1) is used alone, the */
/*                  matrix will be scaled on output. */

/*                  KEY=2 */
/*                         Scale the nonzero columns of the */
/*                         entire data matrix. */
/*                  (E) */
/*                  (A) */
/*                  (G) */

/*                  to have length one.  The data set for this */
/*                  option is a single value.  It must be */
/*                  nonzero if unit length column scaling */
/*                  is desired. */

/*                  KEY=3 */
/*                         Scale columns of the entire data matrix */
/*                  (E) */
/*                  (A) */
/*                  (G) */

/*                  with a user-provided diagonal matrix. */
/*                  The data set for this option consists */
/*                  of the N diagonal scaling factors, one for */
/*                  each matrix column. */

/*                  KEY=4 */
/*                         Change the rank determination tolerance for */
/*                  the equality constraint equations from */
/*                  the nominal value of SQRT(DRELPR).  This quantity can */
/*                  be no smaller than DRELPR, the arithmetic- */
/*                  storage precision.  The quantity DRELPR is the */
/*                  largest positive number such that T=1.+DRELPR */
/*                  satisfies T .EQ. 1.  The quantity used */
/*                  here is internally restricted to be at */
/*                  least DRELPR.  The data set for this option */
/*                  is the new tolerance. */

/*                  KEY=5 */
/*                         Change the rank determination tolerance for */
/*                  the reduced least squares equations from */
/*                  the nominal value of SQRT(DRELPR).  This quantity can */
/*                  be no smaller than DRELPR, the arithmetic- */
/*                  storage precision.  The quantity used */
/*                  here is internally restricted to be at */
/*                  least DRELPR.  The data set for this option */
/*                  is the new tolerance. */

/*                  For example, suppose we want to change */
/*                  the tolerance for the reduced least squares */
/*                  problem, compute the covariance matrix of */
/*                  the solution parameters, and provide */
/*                  column scaling for the data matrix.  For */
/*                  these options the dimension of PRGOPT(*) */
/*                  must be at least N+9.  The Fortran statements */
/*                  defining these options would be as follows: */

/*                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*)) */
/*                  PRGOPT(2)=1 (covariance matrix key) */
/*                  PRGOPT(3)=1 (covariance matrix wanted) */

/*                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*)) */
/*                  PRGOPT(5)=5 (least squares equas.  tolerance key) */
/*                  PRGOPT(6)=... (new value of the tolerance) */

/*                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*)) */
/*                  PRGOPT(8)=3 (user-provided column scaling key) */

/*                  CALL DCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N */
/*                    scaling factors from the user array D(*) */
/*                    to PRGOPT(9)-PRGOPT(N+8)) */

/*                  PRGOPT(N+9)=1 (no more options to change) */

/*                  The contents of PRGOPT(*) are not modified */
/*                  by the subprogram. */
/*                  The options for WNNLS( ) can also be included */
/*                  in this array.  The values of KEY recognized */
/*                  by WNNLS( ) are 6, 7 and 8.  Their functions */
/*                  are documented in the usage instructions for */
/*                  subroutine WNNLS( ).  Normally these options */
/*                  do not need to be modified when using DLSEI( ). */

/*     IP(1),       The amounts of working storage actually */
/*     IP(2)        allocated for the working arrays WS(*) and */
/*                  IP(*), respectively.  These quantities are */
/*                  compared with the actual amounts of storage */
/*                  needed by DLSEI( ).  Insufficient storage */
/*                  allocated for either WS(*) or IP(*) is an */
/*                  error.  This feature was included in DLSEI( ) */
/*                  because miscalculating the storage formulas */
/*                  for WS(*) and IP(*) might very well lead to */
/*                  subtle and hard-to-find execution errors. */

/*                  The length of WS(*) must be at least */

/*                  LW = 2*(ME+N)+K+(MG+2)*(N+7) */

/*                  where K = max(MA+MG,N) */
/*                  This test will not be made if IP(1).LE.0. */

/*                  The length of IP(*) must be at least */

/*                  LIP = MG+2*N+2 */
/*                  This test will not be made if IP(2).LE.0. */

/*     Output.. All TYPE REAL variables are DOUBLE PRECISION */

/*     X(*),RNORME,  The array X(*) contains the solution parameters */
/*     RNORML        if the integer output flag MODE = 0 or 1. */
/*                   The definition of MODE is given directly below. */
/*                   When MODE = 0 or 1, RNORME and RNORML */
/*                   respectively contain the residual vector */
/*                   Euclidean lengths of F - EX and B - AX.  When */
/*                   MODE=1 the equality constraint equations EX=F */
/*                   are contradictory, so RNORME .NE. 0.  The residual */
/*                   vector F-EX has minimal Euclidean length.  For */
/*                   MODE .GE. 2, none of these parameters is defined. */

/*     MODE          Integer flag that indicates the subprogram */
/*                   status after completion.  If MODE .GE. 2, no */
/*                   solution has been computed. */

/*                   MODE = */

/*                   0  Both equality and inequality constraints */
/*                      are compatible and have been satisfied. */

/*                   1  Equality constraints are contradictory. */
/*                      A generalized inverse solution of EX=F was used */
/*                      to minimize the residual vector length F-EX. */
/*                      In this sense, the solution is still meaningful. */

/*                   2  Inequality constraints are contradictory. */

/*                   3  Both equality and inequality constraints */
/*                      are contradictory. */

/*                   The following interpretation of */
/*                   MODE=1,2 or 3 must be made.  The */
/*                   sets consisting of all solutions */
/*                   of the equality constraints EX=F */
/*                   and all vectors satisfying GX .GE. H */
/*                   have no points in common.  (In */
/*                   particular this does not say that */
/*                   each individual set has no points */
/*                   at all, although this could be the */
/*                   case.) */

/*                   4  Usage error occurred.  The value */
/*                      of MDW is .LT. ME+MA+MG, MDW is */
/*                      .LT. N and a covariance matrix is */
/*                      requested, or the option vector */
/*                      PRGOPT(*) is not properly defined, */
/*                      or the lengths of the working arrays */
/*                      WS(*) and IP(*), when specified in */
/*                      IP(1) and IP(2) respectively, are not */
/*                      long enough. */

/*     W(*,*)        The array W(*,*) contains the N by N symmetric */
/*                   covariance matrix of the solution parameters, */
/*                   provided this was requested on input with */
/*                   the option vector PRGOPT(*) and the output */
/*                   flag is returned with MODE = 0 or 1. */

/*     IP(*)         The integer working array has three entries */
/*                   that provide rank and working array length */
/*                   information after completion. */

/*                      IP(1) = rank of equality constraint */
/*                              matrix.  Define this quantity */
/*                              as KRANKE. */

/*                      IP(2) = rank of reduced least squares */
/*                              problem. */

/*                      IP(3) = the amount of storage in the */
/*                              working array WS(*) that was */
/*                              actually used by the subprogram. */
/*                              The formula given above for the length */
/*                              of WS(*) is a necessary overestimate. */
/*                              If exactly the same problem matrices */
/*                              are used in subsequent executions, */
/*                              the declared dimension of WS(*) can */
/*                              be reduced to this output value. */
/*     User Designated */
/*     Working Arrays.. */

/*     WS(*),IP(*)              These are respectively type real */
/*                              and type integer working arrays. */
/*                              Their required minimal lengths are */
/*                              given above. */

/* ***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for */
/*                 linear least squares problems with equality and */
/*                 nonnegativity constraints, Report SAND77-0552, Sandia */
/*                 Laboratories, June 1978. */
/*               K. H. Haskell and R. J. Hanson, Selected algorithms for */
/*                 the linearly constrained least squares problem - a */
/*                 users guide, Report SAND78-1290, Sandia Laboratories, */
/*                 August 1979. */
/*               K. H. Haskell and R. J. Hanson, An algorithm for */
/*                 linear least squares problems with equality and */
/*                 nonnegativity constraints, Mathematical Programming */
/*                 21 (1981), pp. 98-118. */
/*               R. J. Hanson and K. H. Haskell, Two algorithms for the */
/*                 linearly constrained least squares problem, ACM */
/*                 Transactions on Mathematical Software, September 1982. */
/* ***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DDOT, DH12, DLSI, */
/*                    DNRM2, DSCAL, DSWAP, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   900604  DP version created from SP version.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DLSEI */



    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --ws;
    --ip;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DLSEI */

/*     Set the nominal tolerance used in the code for the equality */
/*     constraint equations. */

    if (first) {
	drelpr = d1mach_(&c__4);
    }
    first = FALSE_;
    tau = sqrt(drelpr);

/*     Check that enough storage was allocated in WS(*) and IP(*). */

    *mode = 4;
/* Computing MIN */
    i__1 = min(*n,*me), i__1 = min(i__1,*ma);
    if (min(i__1,*mg) < 0) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*me), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&(*ma), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&(*mg), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 78, a__1[0] = "ALL OF THE VARIABLES N, ME, MA, MG MUST BE "
		".GE. 0$$ENTERED ROUTINE WITH$$N  = ";
	i__2[1] = 8, a__1[1] = xern1;
	i__2[2] = 7, a__1[2] = "$$ME = ";
	i__2[3] = 8, a__1[3] = xern2;
	i__2[4] = 7, a__1[4] = "$$MA = ";
	i__2[5] = 8, a__1[5] = xern3;
	i__2[6] = 7, a__1[6] = "$$MG = ";
	i__2[7] = 8, a__1[7] = xern4;
	s_cat(ch__1, a__1, i__2, &c__8, (ftnlen)131);
	xermsg_("SLATEC", "LSEI", ch__1, &c__2, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)131);
	return 0;
    }

    if (ip[1] > 0) {
/* Computing MAX */
	i__1 = *ma + *mg;
	lchk = (*me + *n << 1) + max(i__1,*n) + (*mg + 2) * (*n + 7);
	if (ip[1] < lchk) {
	    s_wsfi(&io___13);
	    do_fio(&c__1, (char *)&lchk, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 52, a__2[0] = "INSUFFICIENT STORAGE ALLOCATED FOR WS(*"
		    "), NEED LW = ";
	    i__3[1] = 8, a__2[1] = xern1;
	    s_cat(ch__2, a__2, i__3, &c__2, (ftnlen)60);
	    xermsg_("SLATEC", "DLSEI", ch__2, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)60);
	    return 0;
	}
    }

    if (ip[2] > 0) {
	lchk = *mg + (*n << 1) + 2;
	if (ip[2] < lchk) {
	    s_wsfi(&io___14);
	    do_fio(&c__1, (char *)&lchk, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 53, a__2[0] = "INSUFFICIENT STORAGE ALLOCATED FOR IP(*"
		    "), NEED LIP = ";
	    i__3[1] = 8, a__2[1] = xern1;
	    s_cat(ch__3, a__2, i__3, &c__2, (ftnlen)61);
	    xermsg_("SLATEC", "DLSEI", ch__3, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)61);
	    return 0;
	}
    }

/*     Compute number of possible right multiplying Householder */
/*     transformations. */

    m = *me + *ma + *mg;
    if (*n <= 0 || m <= 0) {
	*mode = 0;
	*rnorme = 0.;
	*rnorml = 0.;
	return 0;
    }

    if (*mdw < m) {
	xermsg_("SLATEC", "DLSEI", "MDW.LT.ME+MA+MG IS AN ERROR", &c__2, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)27);
	return 0;
    }

    np1 = *n + 1;
    kranke = min(*me,*n);
    n1 = (kranke << 1) + 1;
    n2 = n1 + *n;

/*     Set nominal values. */

/*     The nominal column scaling used in the code is */
/*     the identity scaling. */

    dcopy_(n, &c_b41, &c__0, &ws[n1], &c__1);

/*     No covariance matrix is nominally computed. */

    cov = FALSE_;

/*     Process option vector. */
/*     Define bound for number of options to change. */

    nopt = 1000;
    ntimes = 0;

/*     Define bound for positive values of LINK. */

    nlink = 100000;
    last = 1;
    link = (integer) prgopt[1];
    if (link == 0 || link > nlink) {
	xermsg_("SLATEC", "DLSEI", "THE OPTION VECTOR IS UNDEFINED", &c__2, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)30);
	return 0;
    }

L100:
    if (link > 1) {
	++ntimes;
	if (ntimes > nopt) {
	    xermsg_("SLATEC", "DLSEI", "THE LINKS IN THE OPTION VECTOR ARE C"
		    "YCLING.", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)43);
	    return 0;
	}

	key = (integer) prgopt[last + 1];
	if (key == 1) {
	    cov = prgopt[last + 2] != 0.;
	} else if (key == 2 && prgopt[last + 2] != 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		t = dnrm2_(&m, &w[j * w_dim1 + 1], &c__1);
		if (t != 0.) {
		    t = 1. / t;
		}
		ws[j + n1 - 1] = t;
/* L110: */
	    }
	} else if (key == 3) {
	    dcopy_(n, &prgopt[last + 2], &c__1, &ws[n1], &c__1);
	} else if (key == 4) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    tau = max(d__1,d__2);
	}

	next = (integer) prgopt[link];
	if (next <= 0 || next > nlink) {
	    xermsg_("SLATEC", "DLSEI", "THE OPTION VECTOR IS UNDEFINED", &
		    c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)30);
	    return 0;
	}

	last = link;
	link = next;
	goto L100;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dscal_(&m, &ws[n1 + j - 1], &w[j * w_dim1 + 1], &c__1);
/* L120: */
    }

    if (cov && *mdw < *n) {
	xermsg_("SLATEC", "DLSEI", "MDW .LT. N WHEN COV MATRIX NEEDED, IS AN"
		" ERROR", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)46);
	return 0;
    }

/*     Problem definition and option vector OK. */

    *mode = 0;

/*     Compute norm of equality constraint matrix and right side. */

    enorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = enorm, d__2 = dasum_(me, &w[j * w_dim1 + 1], &c__1);
	enorm = max(d__1,d__2);
/* L130: */
    }

    fnorm = dasum_(me, &w[np1 * w_dim1 + 1], &c__1);
    snmax = 0.;
    rnmax = 0.;
    i__1 = kranke;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Compute maximum ratio of vector lengths. Partition is at */
/*        column I. */

	i__4 = *me;
	for (k = i__; k <= i__4; ++k) {
	    i__5 = *n - i__ + 1;
	    sn = ddot_(&i__5, &w[k + i__ * w_dim1], mdw, &w[k + i__ * w_dim1],
		     mdw);
	    i__5 = i__ - 1;
	    rn = ddot_(&i__5, &w[k + w_dim1], mdw, &w[k + w_dim1], mdw);
	    if (rn == 0. && sn > snmax) {
		snmax = sn;
		imax = k;
	    } else if (k == i__ || sn * rnmax > rn * snmax) {
		snmax = sn;
		rnmax = rn;
		imax = k;
	    }
/* L140: */
	}

/*        Interchange rows if necessary. */

	if (i__ != imax) {
	    dswap_(&np1, &w[i__ + w_dim1], mdw, &w[imax + w_dim1], mdw);
	}
/* Computing 2nd power */
	d__1 = tau;
	if (snmax > rnmax * (d__1 * d__1)) {

/*        Eliminate elements I+1,...,N in row I. */

	    i__4 = i__ + 1;
	    i__5 = m - i__;
	    dh12_(&c__1, &i__, &i__4, n, &w[i__ + w_dim1], mdw, &ws[i__], &w[
		    i__ + 1 + w_dim1], mdw, &c__1, &i__5);
	} else {
	    kranke = i__ - 1;
	    goto L160;
	}
/* L150: */
    }

/*     Save diagonal terms of lower trapezoidal matrix. */

L160:
    i__1 = *mdw + 1;
    dcopy_(&kranke, &w[w_offset], &i__1, &ws[kranke + 1], &c__1);

/*     Use Householder transformation from left to achieve */
/*     KRANKE by KRANKE upper triangular form. */

    if (kranke < *me) {
	for (k = kranke; k >= 1; --k) {

/*           Apply transformation to matrix cols. 1,...,K-1. */

	    i__1 = kranke + 1;
	    i__4 = k - 1;
	    dh12_(&c__1, &k, &i__1, me, &w[k * w_dim1 + 1], &c__1, &up, &w[
		    w_offset], &c__1, mdw, &i__4);

/*           Apply to rt side vector. */

	    i__1 = kranke + 1;
	    dh12_(&c__2, &k, &i__1, me, &w[k * w_dim1 + 1], &c__1, &up, &w[
		    np1 * w_dim1 + 1], &c__1, &c__1, &c__1);
/* L170: */
	}
    }

/*     Solve for variables 1,...,KRANKE in new coordinates. */

    dcopy_(&kranke, &w[np1 * w_dim1 + 1], &c__1, &x[1], &c__1);
    i__1 = kranke;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__4 = i__ - 1;
	x[i__] = (x[i__] - ddot_(&i__4, &w[i__ + w_dim1], mdw, &x[1], &c__1)) 
		/ w[i__ + i__ * w_dim1];
/* L180: */
    }

/*     Compute residuals for reduced problem. */

    mep1 = *me + 1;
    *rnorml = 0.;
    i__1 = m;
    for (i__ = mep1; i__ <= i__1; ++i__) {
	w[i__ + np1 * w_dim1] -= ddot_(&kranke, &w[i__ + w_dim1], mdw, &x[1], 
		&c__1);
	sn = ddot_(&kranke, &w[i__ + w_dim1], mdw, &w[i__ + w_dim1], mdw);
	i__4 = *n - kranke;
	rn = ddot_(&i__4, &w[i__ + (kranke + 1) * w_dim1], mdw, &w[i__ + (
		kranke + 1) * w_dim1], mdw);
/* Computing 2nd power */
	d__1 = tau;
	if (rn <= sn * (d__1 * d__1) && kranke < *n) {
	    i__4 = *n - kranke;
	    dcopy_(&i__4, &c_b95, &c__0, &w[i__ + (kranke + 1) * w_dim1], mdw)
		    ;
	}
/* L190: */
    }

/*     Compute equality constraint equations residual length. */

    i__1 = *me - kranke;
    *rnorme = dnrm2_(&i__1, &w[kranke + 1 + np1 * w_dim1], &c__1);

/*     Move reduced problem data upward if KRANKE.LT.ME. */

    if (kranke < *me) {
	i__1 = np1;
	for (j = 1; j <= i__1; ++j) {
	    i__4 = m - *me;
	    dcopy_(&i__4, &w[*me + 1 + j * w_dim1], &c__1, &w[kranke + 1 + j *
		     w_dim1], &c__1);
/* L200: */
	}
    }

/*     Compute solution of reduced problem. */

    i__1 = *n - kranke;
    dlsi_(&w[kranke + 1 + (kranke + 1) * w_dim1], mdw, ma, mg, &i__1, &prgopt[
	    1], &x[kranke + 1], rnorml, mode, &ws[n2], &ip[2]);

/*     Test for consistency of equality constraints. */

    if (*me > 0) {
	mdeqc = 0;
	xnrme = dasum_(&kranke, &w[np1 * w_dim1 + 1], &c__1);
	if (*rnorme > tau * (enorm * xnrme + fnorm)) {
	    mdeqc = 1;
	}
	*mode += mdeqc;

/*        Check if solution to equality constraints satisfies inequality */
/*        constraints when there are no degrees of freedom left. */

	if (kranke == *n && *mg > 0) {
	    xnorm = dasum_(n, &x[1], &c__1);
	    mapke1 = *ma + kranke + 1;
	    mend = *ma + kranke + *mg;
	    i__1 = mend;
	    for (i__ = mapke1; i__ <= i__1; ++i__) {
		size = dasum_(n, &w[i__ + w_dim1], mdw) * xnorm + (d__1 = w[
			i__ + np1 * w_dim1], abs(d__1));
		if (w[i__ + np1 * w_dim1] > tau * size) {
		    *mode += 2;
		    goto L290;
		}
/* L210: */
	    }
	}
    }

/*     Replace diagonal terms of lower trapezoidal matrix. */

    if (kranke > 0) {
	i__1 = *mdw + 1;
	dcopy_(&kranke, &ws[kranke + 1], &c__1, &w[w_offset], &i__1);

/*        Reapply transformation to put solution in original coordinates. */

	for (i__ = kranke; i__ >= 1; --i__) {
	    i__1 = i__ + 1;
	    dh12_(&c__2, &i__, &i__1, n, &w[i__ + w_dim1], mdw, &ws[i__], &x[
		    1], &c__1, &c__1, &c__1);
/* L220: */
	}

/*        Compute covariance matrix of equality constrained problem. */

	if (cov) {
/* Computing MIN */
	    i__1 = kranke, i__4 = *n - 1;
	    for (j = min(i__1,i__4); j >= 1; --j) {
		rb = ws[j] * w[j + j * w_dim1];
		if (rb != 0.) {
		    rb = 1. / rb;
		}
		jp1 = j + 1;
		i__1 = *n;
		for (i__ = jp1; i__ <= i__1; ++i__) {
		    i__4 = *n - j;
		    w[i__ + j * w_dim1] = rb * ddot_(&i__4, &w[i__ + jp1 * 
			    w_dim1], mdw, &w[j + jp1 * w_dim1], mdw);
/* L230: */
		}

		i__1 = *n - j;
		gam = rb * .5 * ddot_(&i__1, &w[jp1 + j * w_dim1], &c__1, &w[
			j + jp1 * w_dim1], mdw);
		i__1 = *n - j;
		daxpy_(&i__1, &gam, &w[j + jp1 * w_dim1], mdw, &w[jp1 + j * 
			w_dim1], &c__1);
		i__1 = *n;
		for (i__ = jp1; i__ <= i__1; ++i__) {
		    i__4 = *n;
		    for (k = i__; k <= i__4; ++k) {
			w[i__ + k * w_dim1] = w[i__ + k * w_dim1] + w[j + i__ 
				* w_dim1] * w[k + j * w_dim1] + w[i__ + j * 
				w_dim1] * w[j + k * w_dim1];
			w[k + i__ * w_dim1] = w[i__ + k * w_dim1];
/* L240: */
		    }
/* L250: */
		}
		uj = ws[j];
		vj = gam * uj;
		w[j + j * w_dim1] = uj * vj + uj * vj;
		i__1 = *n;
		for (i__ = jp1; i__ <= i__1; ++i__) {
		    w[j + i__ * w_dim1] = uj * w[i__ + j * w_dim1] + vj * w[j 
			    + i__ * w_dim1];
/* L260: */
		}
		i__1 = *n - j;
		dcopy_(&i__1, &w[j + jp1 * w_dim1], mdw, &w[jp1 + j * w_dim1],
			 &c__1);
/* L270: */
	    }
	}
    }

/*     Apply the scaling to the covariance matrix. */

    if (cov) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dscal_(n, &ws[i__ + n1 - 1], &w[i__ + w_dim1], mdw);
	    dscal_(n, &ws[i__ + n1 - 1], &w[i__ * w_dim1 + 1], &c__1);
/* L280: */
	}
    }

/*     Rescale solution vector. */

L290:
    if (*mode <= 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    x[j] *= ws[n1 + j - 1];
/* L300: */
	}
    }

    ip[1] = kranke;
    ip[3] = ip[3] + (kranke << 1) + *n;
    return 0;
} /* dlsei_ */


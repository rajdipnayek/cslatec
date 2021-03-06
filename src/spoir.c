/* spoir.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__4 = 4;
static integer c_n2 = -2;
static integer c__3 = 3;
static integer c_n3 = -3;
static integer c_n4 = -4;
static integer c_n10 = -10;
static integer c__0 = 0;

/* DECK SPOIR */
/* Subroutine */ int spoir_(real *a, integer *lda, integer *n, real *v, 
	integer *itask, integer *ind, real *work)
{
    /* System generated locals */
    address a__1[4], a__2[3];
    integer a_dim1, a_offset, work_dim1, work_offset, i__1[4], i__2[3], i__3, 
	    i__4, i__5;
    real r__1, r__2, r__3;
    char ch__1[40], ch__2[27], ch__3[31];

    /* Local variables */
    static integer j, info;
    static char xern1[8], xern2[8];
    extern /* Subroutine */ int spofa_(real *, integer *, integer *, integer *
	    );
    extern doublereal dsdot_(integer *, real *, integer *, real *, integer *);
    static real dnorm;
    extern doublereal sasum_(integer *, real *, integer *);
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sposl_(real *, integer *, integer *, real *);
    static real xnorm;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___4 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SPOIR */
/* ***PURPOSE  Solve a positive definite symmetric system of linear */
/*            equations.  Iterative refinement is used to obtain an error */
/*            estimate. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2B1B */
/* ***TYPE      SINGLE PRECISION (SPOIR-S, CPOIR-C) */
/* ***KEYWORDS  HERMITIAN, LINEAR EQUATIONS, POSITIVE DEFINITE, SYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*    Subroutine SPOIR solves a real positive definite symmetric */
/*    NxN system of single precision linear equations using LINPACK */
/*    subroutines SPOFA and SPOSL.  One pass of iterative refine- */
/*    ment is used only to obtain an estimate of the accuracy.  That */
/*    is, if A is an NxN real positive definite symmetric matrix */
/*    and if X and B are real N-vectors, then SPOIR solves the */
/*    equation */

/*                          A*X=B. */

/*    The matrix A is first factored into upper and lower */
/*    triangular matrices R and R-TRANSPOSE.  These */
/*    factors are used to calculate the solution, X. */
/*    Then the residual vector is found and used */
/*    to calculate an estimate of the relative error, IND. */
/*    IND estimates the accuracy of the solution only when the */
/*    input matrix and the right hand side are represented */
/*    exactly in the computer and does not take into account */
/*    any errors in the input data. */

/*    If the equation A*X=B is to be solved for more than one vector */
/*    B, the factoring of A does not need to be performed again and */
/*    the option to only solve (ITASK .GT. 1) will be faster for */
/*    the succeeding solutions.  In this case, the contents of A, */
/*    LDA, N, and WORK must not have been altered by the user */
/*    following factorization (ITASK=1).  IND will not be changed */
/*    by SPOIR in this case. */

/*  Argument Description *** */
/*    A      REAL(LDA,N) */
/*             the doubly subscripted array with dimension (LDA,N) */
/*             which contains the coefficient matrix.  Only the */
/*             upper triangle, including the diagonal, of the */
/*             coefficient matrix need be entered.  A is not */
/*             altered by the routine. */
/*    LDA    INTEGER */
/*             the leading dimension of the array A.  LDA must be great- */
/*             er than or equal to N.  (Terminal error message IND=-1) */
/*    N      INTEGER */
/*             the order of the matrix A.  N must be greater than */
/*             or equal to one.  (Terminal error message IND=-2) */
/*    V      REAL(N) */
/*             on entry, the singly subscripted array(vector) of di- */
/*               mension N which contains the right hand side B of a */
/*               system of simultaneous linear equations A*X=B. */
/*             on return, V contains the solution vector, X . */
/*    ITASK  INTEGER */
/*             If ITASK = 1, the matrix A is factored and then the */
/*               linear equation is solved. */
/*             If ITASK .GT. 1, the equation is solved using the existing */
/*               factored matrix A (stored in WORK). */
/*             If ITASK .LT. 1, then terminal terminal error IND=-3 is */
/*               printed. */
/*    IND    INTEGER */
/*             GT. 0  IND is a rough estimate of the number of digits */
/*                     of accuracy in the solution, X.  IND=75 means */
/*                     that the solution vector X is zero. */
/*             LT. 0  See error message corresponding to IND below. */
/*    WORK   REAL(N*(N+1)) */
/*             a singly subscripted array of dimension at least N*(N+1). */

/*  Error Messages Printed *** */

/*    IND=-1  terminal   N is greater than LDA. */
/*    IND=-2  terminal   N is less than one. */
/*    IND=-3  terminal   ITASK is less than one. */
/*    IND=-4  Terminal   The matrix A is computationally singular */
/*                         or is not positive definite. */
/*                         A solution has not been computed. */
/*    IND=-10 warning    The solution has no apparent significance. */
/*                         The solution may be inaccurate or the matrix */
/*                         A may be poorly scaled. */

/*               Note-  The above terminal(*fatal*) error messages are */
/*                      designed to be handled by XERMSG in which */
/*                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0 */
/*                      for warning error messages from XERMSG.  Unless */
/*                      the user provides otherwise, an error message */
/*                      will be printed followed by an abort. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DSDOT, R1MACH, SASUM, SCOPY, SPOFA, SPOSL, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800528  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SPOIR */

/* ***FIRST EXECUTABLE STATEMENT  SPOIR */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    work_dim1 = *n;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    --v;

    /* Function Body */
    if (*lda < *n) {
	*ind = -1;
	s_wsfi(&io___2);
	do_fio(&c__1, (char *)&(*lda), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___4);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 6, a__1[0] = "LDA = ";
	i__1[1] = 8, a__1[1] = xern1;
	i__1[2] = 18, a__1[2] = " IS LESS THAN N = ";
	i__1[3] = 8, a__1[3] = xern2;
	s_cat(ch__1, a__1, i__1, &c__4, (ftnlen)40);
	xermsg_("SLATEC", "SPOIR", ch__1, &c_n1, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)40);
	return 0;
    }

    if (*n <= 0) {
	*ind = -2;
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 4, a__2[0] = "N = ";
	i__2[1] = 8, a__2[1] = xern1;
	i__2[2] = 15, a__2[2] = " IS LESS THAN 1";
	s_cat(ch__2, a__2, i__2, &c__3, (ftnlen)27);
	xermsg_("SLATEC", "SPOIR", ch__2, &c_n2, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)27);
	return 0;
    }

    if (*itask < 1) {
	*ind = -3;
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&(*itask), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 8, a__2[0] = "ITASK = ";
	i__2[1] = 8, a__2[1] = xern1;
	i__2[2] = 15, a__2[2] = " IS LESS THAN 1";
	s_cat(ch__3, a__2, i__2, &c__3, (ftnlen)31);
	xermsg_("SLATEC", "SPOIR", ch__3, &c_n3, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)31);
	return 0;
    }

    if (*itask == 1) {

/*        MOVE MATRIX A TO WORK */

	i__3 = *n;
	for (j = 1; j <= i__3; ++j) {
	    scopy_(n, &a[j * a_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &
		    c__1);
/* L10: */
	}

/*        FACTOR MATRIX A INTO R */
	spofa_(&work[work_offset], n, n, &info);

/*        CHECK FOR  SINGULAR OR NOT POS.DEF. MATRIX */
	if (info != 0) {
	    *ind = -4;
	    xermsg_("SLATEC", "SPOIR", "SINGULAR OR NOT POSITIVE DEFINITE - "
		    "NO SOLUTION", &c_n4, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)
		    47);
	    return 0;
	}
    }

/*     SOLVE AFTER FACTORING */
/*     MOVE VECTOR B TO WORK */

    scopy_(n, &v[1], &c__1, &work[(*n + 1) * work_dim1 + 1], &c__1);
    sposl_(&work[work_offset], n, n, &v[1]);

/*     FORM NORM OF X0 */

    xnorm = sasum_(n, &v[1], &c__1);
    if (xnorm == 0.f) {
	*ind = 75;
	return 0;
    }

/*     COMPUTE  RESIDUAL */

    i__3 = *n;
    for (j = 1; j <= i__3; ++j) {
	i__4 = j - 1;
	i__5 = *n - j + 1;
	work[j + (*n + 1) * work_dim1] = -work[j + (*n + 1) * work_dim1] + 
		dsdot_(&i__4, &a[j * a_dim1 + 1], &c__1, &v[1], &c__1) + 
		dsdot_(&i__5, &a[j + j * a_dim1], lda, &v[j], &c__1);
/* L40: */
    }

/*     SOLVE A*DELTA=R */

    sposl_(&work[work_offset], n, n, &work[(*n + 1) * work_dim1 + 1]);

/*     FORM NORM OF DELTA */

    dnorm = sasum_(n, &work[(*n + 1) * work_dim1 + 1], &c__1);

/*     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS) */
/*     AND CHECK FOR IND GREATER THAN ZERO */

/* Computing MAX */
    r__2 = r1mach_(&c__4), r__3 = dnorm / xnorm;
    r__1 = dmax(r__2,r__3);
    *ind = -r_lg10(&r__1);
    if (*ind <= 0) {
	*ind = -10;
	xermsg_("SLATEC", "SPOIR", "SOLUTION MAY HAVE NO SIGNIFICANCE", &
		c_n10, &c__0, (ftnlen)6, (ftnlen)5, (ftnlen)33);
    }
    return 0;
} /* spoir_ */


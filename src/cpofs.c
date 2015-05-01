/* cpofs.f -- translated by f2c (version 12.02.01).
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

/* DECK CPOFS */
/* Subroutine */ int cpofs_(complex *a, integer *lda, integer *n, complex *v, 
	integer *itask, integer *ind, complex *work)
{
    /* System generated locals */
    address a__1[4], a__2[3];
    integer a_dim1, a_offset, i__1[4], i__2[3];
    real r__1;
    char ch__1[40], ch__2[27], ch__3[31];

    /* Local variables */
    static integer info;
    static char xern1[8], xern2[8];
    extern /* Subroutine */ int cpoco_(complex *, integer *, integer *, real *
	    , complex *, integer *);
    static real rcond;
    extern /* Subroutine */ int cposl_(complex *, integer *, integer *, 
	    complex *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___4 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  CPOFS */
/* ***PURPOSE  Solve a positive definite symmetric complex system of */
/*            linear equations. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2D1B */
/* ***TYPE      COMPLEX (SPOFS-S, DPOFS-D, CPOFS-C) */
/* ***KEYWORDS  HERMITIAN, LINEAR EQUATIONS, POSITIVE DEFINITE, SYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*    Subroutine CPOFS solves a  positive definite symmetric */
/*    NxN system of complex linear equations using LINPACK */
/*    subroutines CPOCO and CPOSL.  That is, if A is an NxN */
/*    complex positive definite symmetric matrix and if X and B */
/*    are complex N-vectors, then CPOFS solves the equation */

/*                          A*X=B. */

/*    Care should be taken not to use CPOFS with a non-Hermitian */
/*    matrix. */

/*    The matrix A is first factored into upper and lower tri- */
/*    angular matrices R and R-TRANSPOSE.  These factors are used to */
/*    find the solution vector X.  An approximate condition number is */
/*    calculated to provide a rough estimate of the number of */
/*    digits of accuracy in the computed solution. */

/*    If the equation A*X=B is to be solved for more than one vector */
/*    B, the factoring of a does not need to be performed again and */
/*    the option to only solve (ITASK .GT. 1) will be faster for */
/*    the succeeding solutions.  In this case, the contents of A, */
/*    LDA, and N must not have been altered by the user following */
/*    factorization (ITASK=1).  IND will not be changed by CPOFS */
/*    in this case. */

/*  Argument Description *** */

/*    A      COMPLEX(LDA,N) */
/*             on entry, the doubly subscripted array with dimension */
/*               (LDA,N) which contains the coefficient matrix.  Only */
/*               the upper triangle, including the diagonal, of the */
/*               coefficient matrix need be entered and will subse- */
/*               quently be referenced and changed by the routine. */
/*             on return, contains in its upper triangle an upper */
/*               triangular matrix R such that  A = (R-TRANSPOSE) * R . */
/*    LDA    INTEGER */
/*             the leading dimension of the array A.  LDA must be great- */
/*             er than or equal to N.  (terminal error message IND=-1) */
/*    N      INTEGER */
/*             the order of the matrix A.  N must be greater */
/*             than or equal to 1.  (terminal error message IND=-2) */
/*    V      COMPLEX(N) */
/*             on entry the singly subscripted array(vector) of di- */
/*               mension N which contains the right hand side B of a */
/*               system of simultaneous linear equations A*X=B. */
/*             on return, V contains the solution vector, X . */
/*    ITASK  INTEGER */
/*             if ITASK = 1, the matrix A is factored and then the */
/*               linear equation is solved. */
/*             if ITASK .GT. 1, the equation is solved using the existing */
/*               factored matrix A. */
/*             if ITASK .LT. 1, then terminal error message IND=-3 is */
/*               printed. */
/*    IND    INTEGER */
/*             GT. 0  IND is a rough estimate of the number of digits */
/*                     of accuracy in the solution, X. */
/*             LT. 0  see error message corresponding to IND below. */
/*    WORK   COMPLEX(N) */
/*             a singly subscripted array of dimension at least N. */

/*  Error Messages Printed *** */

/*    IND=-1  terminal   N is greater than LDA. */
/*    IND=-2  terminal   N is less than 1. */
/*    IND=-3  terminal   ITASK is less than 1. */
/*    IND=-4  terminal   The matrix A is computationally singular or */
/*                         is not positive definite.  A solution */
/*                         has not been computed. */
/*    IND=-10 warning    The solution has no apparent significance. */
/*                         The solution may be inaccurate or the */
/*                         matrix A may be poorly scaled. */

/*               NOTE-  The above terminal(*fatal*) error messages are */
/*                      designed to be handled by XERMSG in which */
/*                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0 */
/*                      for warning error messages from XERMSG.  Unless */
/*                      the user provides otherwise, an error message */
/*                      will be printed followed by an abort. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CPOCO, CPOSL, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800516  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to */
/*           IF-THEN-ELSE.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPOFS */

/* ***FIRST EXECUTABLE STATEMENT  CPOFS */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --v;
    --work;

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
	xermsg_("SLATEC", "CPOFS", ch__1, &c_n1, &c__1, (ftnlen)6, (ftnlen)5, 
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
	xermsg_("SLATEC", "CPOFS", ch__2, &c_n2, &c__1, (ftnlen)6, (ftnlen)5, 
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
	xermsg_("SLATEC", "CPOFS", ch__3, &c_n3, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)31);
	return 0;
    }

    if (*itask == 1) {

/*        FACTOR MATRIX A INTO R */

	cpoco_(&a[a_offset], lda, n, &rcond, &work[1], &info);

/*        CHECK FOR POSITIVE DEFINITE MATRIX */

	if (info != 0) {
	    *ind = -4;
	    xermsg_("SLATEC", "CPOFS", "SINGULAR OR NOT POSITIVE DEFINITE - "
		    "NO SOLUTION", &c_n4, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)
		    47);
	    return 0;
	}

/*        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS) */
/*        AND CHECK FOR IND GREATER THAN ZERO */

	r__1 = r1mach_(&c__4) / rcond;
	*ind = -r_lg10(&r__1);
	if (*ind <= 0) {
	    *ind = -10;
	    xermsg_("SLATEC", "CPOFS", "SOLUTION MAY HAVE NO SIGNIFICANCE", &
		    c_n10, &c__0, (ftnlen)6, (ftnlen)5, (ftnlen)33);
	}
    }

/*     SOLVE AFTER FACTORING */

    cposl_(&a[a_offset], lda, n, &v[1]);
    return 0;
} /* cpofs_ */


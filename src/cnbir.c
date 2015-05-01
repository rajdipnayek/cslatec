/* cnbir.f -- translated by f2c (version 12.02.01).
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
static integer c_n5 = -5;
static integer c_n6 = -6;
static integer c_n4 = -4;
static integer c__0 = 0;
static integer c_n10 = -10;

/* DECK CNBIR */
/* Subroutine */ int cnbir_(complex *abe, integer *lda, integer *n, integer *
	ml, integer *mu, complex *v, integer *itask, integer *ind, complex *
	work, integer *iwork)
{
    /* System generated locals */
    address a__1[4], a__2[3];
    integer abe_dim1, abe_offset, work_dim1, work_offset, i__1[4], i__2[3], 
	    i__3, i__4, i__5;
    real r__1, r__2, r__3;
    complex q__1, q__2;
    char ch__1[40], ch__2[27], ch__3[31], ch__4[29];

    /* Local variables */
    static integer j, k, l, m, nc, kk, info;
    static char xern1[8], xern2[8];
    extern /* Subroutine */ int cnbfa_(complex *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), cnbsl_(complex *, 
	    integer *, integer *, integer *, integer *, integer *, complex *, 
	    integer *), ccopy_(integer *, complex *, integer *, complex *, 
	    integer *);
    static real dnorm, xnorm;
    extern doublereal r1mach_(integer *);
    extern /* Complex */ void cdcdot_(complex *, integer *, complex *, 
	    complex *, integer *, complex *, integer *);
    extern doublereal scasum_(integer *, complex *, integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___4 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___8 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  CNBIR */
/* ***PURPOSE  Solve a general nonsymmetric banded system of linear */
/*            equations.  Iterative refinement is used to obtain an error */
/*            estimate. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2C2 */
/* ***TYPE      COMPLEX (SNBIR-S, CNBIR-C) */
/* ***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*    Subroutine CNBIR solves a general nonsymmetric banded NxN */
/*    system of single precision complex linear equations using */
/*    SLATEC subroutines CNBFA and CNBSL.  These are adaptations */
/*    of the LINPACK subroutines CGBFA and CGBSL which require */
/*    a different format for storing the matrix elements. */
/*    One pass of iterative refinement is used only to obtain an */
/*    estimate of the accuracy.  If  A  is an NxN complex banded */
/*    matrix and if  X  and  B  are complex N-vectors, then CNBIR */
/*    solves the equation */

/*                          A*X=B. */

/*    A band matrix is a matrix whose nonzero elements are all */
/*    fairly near the main diagonal, specifically  A(I,J) = 0 */
/*    if  I-J is greater than  ML  or  J-I  is greater than */
/*    MU .  The integers ML and MU are called the lower and upper */
/*    band widths and  M = ML+MU+1  is the total band width. */
/*    CNBIR uses less time and storage than the corresponding */
/*    program for general matrices (CGEIR) if 2*ML+MU .LT. N . */

/*    The matrix A is first factored into upper and lower tri- */
/*    angular matrices U and L using partial pivoting.  These */
/*    factors and the pivoting information are used to find the */
/*    solution vector X .  Then the residual vector is found and used */
/*    to calculate an estimate of the relative error, IND .  IND esti- */
/*    mates the accuracy of the solution only when the input matrix */
/*    and the right hand side are represented exactly in the computer */
/*    and does not take into account any errors in the input data. */

/*    If the equation A*X=B is to be solved for more than one vector */
/*    B, the factoring of A does not need to be performed again and */
/*    the option to only solve (ITASK .GT. 1) will be faster for */
/*    the succeeding solutions.  In this case, the contents of A, LDA, */
/*    N, WORK and IWORK must not have been altered by the user follow- */
/*    ing factorization (ITASK=1).  IND will not be changed by CNBIR */
/*    in this case. */


/*    Band Storage */

/*          If  A  is a band matrix, the following program segment */
/*          will set up the input. */

/*                  ML = (band width below the diagonal) */
/*                  MU = (band width above the diagonal) */
/*                  DO 20 I = 1, N */
/*                     J1 = MAX(1, I-ML) */
/*                     J2 = MIN(N, I+MU) */
/*                     DO 10 J = J1, J2 */
/*                        K = J - I + ML + 1 */
/*                        ABE(I,K) = A(I,J) */
/*               10    CONTINUE */
/*               20 CONTINUE */

/*          This uses columns  1  through  ML+MU+1  of ABE . */

/*    Example:  If the original matrix is */

/*          11 12 13  0  0  0 */
/*          21 22 23 24  0  0 */
/*           0 32 33 34 35  0 */
/*           0  0 43 44 45 46 */
/*           0  0  0 54 55 56 */
/*           0  0  0  0 65 66 */

/*     then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain */

/*           * 11 12 13        , * = not used */
/*          21 22 23 24 */
/*          32 33 34 35 */
/*          43 44 45 46 */
/*          54 55 56  * */
/*          65 66  *  * */


/*  Argument Description *** */

/*    ABE    COMPLEX(LDA,MM) */
/*             on entry, contains the matrix in band storage as */
/*               described above.  MM  must not be less than  M = */
/*               ML+MU+1 .  The user is cautioned to dimension  ABE */
/*               with care since MM is not an argument and cannot */
/*               be checked by CNBIR.  The rows of the original */
/*               matrix are stored in the rows of  ABE  and the */
/*               diagonals of the original matrix are stored in */
/*               columns  1  through  ML+MU+1  of  ABE .  ABE  is */
/*               not altered by the program. */
/*    LDA    INTEGER */
/*             the leading dimension of array ABE.  LDA must be great- */
/*             er than or equal to N.  (terminal error message IND=-1) */
/*    N      INTEGER */
/*             the order of the matrix A.  N must be greater */
/*             than or equal to 1 .  (terminal error message IND=-2) */
/*    ML     INTEGER */
/*             the number of diagonals below the main diagonal. */
/*             ML  must not be less than zero nor greater than or */
/*             equal to  N .  (terminal error message IND=-5) */
/*    MU     INTEGER */
/*             the number of diagonals above the main diagonal. */
/*             MU  must not be less than zero nor greater than or */
/*             equal to  N .  (terminal error message IND=-6) */
/*    V      COMPLEX(N) */
/*             on entry, the singly subscripted array(vector) of di- */
/*               mension N which contains the right hand side B of a */
/*               system of simultaneous linear equations A*X=B. */
/*             on return, V contains the solution vector, X . */
/*    ITASK  INTEGER */
/*             if ITASK=1, the matrix A is factored and then the */
/*               linear equation is solved. */
/*             if ITASK .GT. 1, the equation is solved using the existing */
/*               factored matrix A and IWORK. */
/*             if ITASK .LT. 1, then terminal error message IND=-3 is */
/*               printed. */
/*    IND    INTEGER */
/*             GT. 0  IND is a rough estimate of the number of digits */
/*                     of accuracy in the solution, X .  IND=75 means */
/*                     that the solution vector  X  is zero. */
/*             LT. 0  see error message corresponding to IND below. */
/*    WORK   COMPLEX(N*(NC+1)) */
/*             a singly subscripted array of dimension at least */
/*             N*(NC+1)  where  NC = 2*ML+MU+1 . */
/*    IWORK  INTEGER(N) */
/*             a singly subscripted array of dimension at least N. */

/*  Error Messages Printed *** */

/*    IND=-1  terminal   N is greater than LDA. */
/*    IND=-2  terminal   N is less than 1. */
/*    IND=-3  terminal   ITASK is less than 1. */
/*    IND=-4  terminal   The matrix A is computationally singular. */
/*                         A solution has not been computed. */
/*    IND=-5  terminal   ML is less than zero or is greater than */
/*                         or equal to N . */
/*    IND=-6  terminal   MU is less than zero or is greater than */
/*                         or equal to N . */
/*    IND=-10 warning    The solution has no apparent significance. */
/*                         The solution may be inaccurate or the matrix */
/*                         A may be poorly scaled. */

/*               NOTE-  The above terminal(*fatal*) error messages are */
/*                      designed to be handled by XERMSG in which */
/*                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0 */
/*                      for warning error messages from XERMSG.  Unless */
/*                      the user provides otherwise, an error message */
/*                      will be printed followed by an abort. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CCOPY, CDCDOT, CNBFA, CNBSL, R1MACH, SCASUM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800819  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to */
/*           IF-THEN-ELSE.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CNBIR */

/* ***FIRST EXECUTABLE STATEMENT  CNBIR */
    /* Parameter adjustments */
    abe_dim1 = *lda;
    abe_offset = 1 + abe_dim1;
    abe -= abe_offset;
    work_dim1 = *n;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    --v;
    --iwork;

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
	xermsg_("SLATEC", "CNBIR", ch__1, &c_n1, &c__1, (ftnlen)6, (ftnlen)5, 
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
	xermsg_("SLATEC", "CNBIR", ch__2, &c_n2, &c__1, (ftnlen)6, (ftnlen)5, 
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
	xermsg_("SLATEC", "CNBIR", ch__3, &c_n3, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)31);
	return 0;
    }

    if (*ml < 0 || *ml >= *n) {
	*ind = -5;
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*ml), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 5, a__2[0] = "ML = ";
	i__2[1] = 8, a__2[1] = xern1;
	i__2[2] = 16, a__2[2] = " IS OUT OF RANGE";
	s_cat(ch__4, a__2, i__2, &c__3, (ftnlen)29);
	xermsg_("SLATEC", "CNBIR", ch__4, &c_n5, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)29);
	return 0;
    }

    if (*mu < 0 || *mu >= *n) {
	*ind = -6;
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&(*mu), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 5, a__2[0] = "MU = ";
	i__2[1] = 8, a__2[1] = xern1;
	i__2[2] = 16, a__2[2] = " IS OUT OF RANGE";
	s_cat(ch__4, a__2, i__2, &c__3, (ftnlen)29);
	xermsg_("SLATEC", "CNBIR", ch__4, &c_n6, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)29);
	return 0;
    }

    nc = (*ml << 1) + *mu + 1;
    if (*itask == 1) {

/*        MOVE MATRIX ABE TO WORK */

	m = *ml + *mu + 1;
	i__3 = m;
	for (j = 1; j <= i__3; ++j) {
	    ccopy_(n, &abe[j * abe_dim1 + 1], &c__1, &work[j * work_dim1 + 1],
		     &c__1);
/* L10: */
	}

/*        FACTOR MATRIX A INTO LU */
	cnbfa_(&work[work_offset], n, n, ml, mu, &iwork[1], &info);

/*        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX */
	if (info != 0) {
	    *ind = -4;
	    xermsg_("SLATEC", "CNBIR", "SINGULAR MATRIX A - NO SOLUTION", &
		    c_n4, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)31);
	    return 0;
	}
    }

/*     SOLVE WHEN FACTORING COMPLETE */
/*     MOVE VECTOR B TO WORK */

    ccopy_(n, &v[1], &c__1, &work[(nc + 1) * work_dim1 + 1], &c__1);
    cnbsl_(&work[work_offset], n, n, ml, mu, &iwork[1], &v[1], &c__0);

/*     FORM NORM OF X0 */

    xnorm = scasum_(n, &v[1], &c__1);
    if (xnorm == 0.f) {
	*ind = 75;
	return 0;
    }

/*     COMPUTE  RESIDUAL */

    i__3 = *n;
    for (j = 1; j <= i__3; ++j) {
/* Computing MAX */
	i__4 = 1, i__5 = *ml + 2 - j;
	k = max(i__4,i__5);
/* Computing MAX */
	i__4 = 1, i__5 = j - *ml;
	kk = max(i__4,i__5);
/* Computing MIN */
	i__4 = j - 1;
/* Computing MIN */
	i__5 = *n - j;
	l = min(i__4,*ml) + min(i__5,*mu) + 1;
	i__4 = j + (nc + 1) * work_dim1;
	i__5 = j + (nc + 1) * work_dim1;
	q__2.r = -work[i__5].r, q__2.i = -work[i__5].i;
	cdcdot_(&q__1, &l, &q__2, &abe[j + k * abe_dim1], lda, &v[kk], &c__1);
	work[i__4].r = q__1.r, work[i__4].i = q__1.i;
/* L40: */
    }

/*     SOLVE A*DELTA=R */

    cnbsl_(&work[work_offset], n, n, ml, mu, &iwork[1], &work[(nc + 1) * 
	    work_dim1 + 1], &c__0);

/*     FORM NORM OF DELTA */

    dnorm = scasum_(n, &work[(nc + 1) * work_dim1 + 1], &c__1);

/*     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS) */
/*     AND CHECK FOR IND GREATER THAN ZERO */

/* Computing MAX */
    r__2 = r1mach_(&c__4), r__3 = dnorm / xnorm;
    r__1 = dmax(r__2,r__3);
    *ind = -r_lg10(&r__1);
    if (*ind <= 0) {
	*ind = -10;
	xermsg_("SLATEC", "CNBIR", "SOLUTION MAY HAVE NO SIGNIFICANCE", &
		c_n10, &c__0, (ftnlen)6, (ftnlen)5, (ftnlen)33);
    }
    return 0;
} /* cnbir_ */


/* dbndac.f -- translated by f2c (version 12.02.01).
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

/* DECK DBNDAC */
/* Subroutine */ int dbndac_(doublereal *g, integer *mdg, integer *nb, 
	integer *ip, integer *ir, integer *mt, integer *jt)
{
    /* System generated locals */
    integer g_dim1, g_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, l, ie, ig, jg, kh, mh, mu, ig1, ig2, lp1;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal rho;
    static integer nbp1, nerr, iopt;
    static doublereal zero;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBNDAC */
/* ***PURPOSE  Compute the LU factorization of a  banded matrices using */
/*            sequential accumulation of rows of the data matrix. */
/*            Exactly one right-hand side vector is permitted. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D9 */
/* ***TYPE      DOUBLE PRECISION (BNDACC-S, DBNDAC-D) */
/* ***KEYWORDS  BANDED MATRIX, CURVE FITTING, LEAST SQUARES */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*     These subroutines solve the least squares problem Ax = b for */
/*     banded matrices A using sequential accumulation of rows of the */
/*     data matrix.  Exactly one right-hand side vector is permitted. */

/*     These subroutines are intended for the type of least squares */
/*     systems that arise in applications such as curve or surface */
/*     fitting of data.  The least squares equations are accumulated and */
/*     processed using only part of the data.  This requires a certain */
/*     user interaction during the solution of Ax = b. */

/*     Specifically, suppose the data matrix (A B) is row partitioned */
/*     into Q submatrices.  Let (E F) be the T-th one of these */
/*     submatrices where E = (0 C 0).  Here the dimension of E is MT by N */
/*     and the dimension of C is MT by NB.  The value of NB is the */
/*     bandwidth of A.  The dimensions of the leading block of zeros in E */
/*     are MT by JT-1. */

/*     The user of the subroutine DBNDAC provides MT,JT,C and F for */
/*     T=1,...,Q.  Not all of this data must be supplied at once. */

/*     Following the processing of the various blocks (E F), the matrix */
/*     (A B) has been transformed to the form (R D) where R is upper */
/*     triangular and banded with bandwidth NB.  The least squares */
/*     system Rx = d is then easily solved using back substitution by */
/*     executing the statement CALL DBNDSL(1,...). The sequence of */
/*     values for JT must be nondecreasing.  This may require some */
/*     preliminary interchanges of rows and columns of the matrix A. */

/*     The primary reason for these subroutines is that the total */
/*     processing can take place in a working array of dimension MU by */
/*     NB+1.  An acceptable value for MU is */

/*                       MU = MAX(MT + N + 1), */

/*     where N is the number of unknowns. */

/*     Here the maximum is taken over all values of MT for T=1,...,Q. */
/*     Notice that MT can be taken to be a small as one, showing that */
/*     MU can be as small as N+2.  The subprogram DBNDAC processes the */
/*     rows more efficiently if MU is large enough so that each new */
/*     block (C F) has a distinct value of JT. */

/*     The four principle parts of these algorithms are obtained by the */
/*     following call statements */

/*     CALL DBNDAC(...)  Introduce new blocks of data. */

/*     CALL DBNDSL(1,...)Compute solution vector and length of */
/*                       residual vector. */

/*     CALL DBNDSL(2,...)Given any row vector H solve YR = H for the */
/*                       row vector Y. */

/*     CALL DBNDSL(3,...)Given any column vector W solve RZ = W for */
/*                       the column vector Z. */

/*     The dots in the above call statements indicate additional */
/*     arguments that will be specified in the following paragraphs. */

/*     The user must dimension the array appearing in the call list.. */
/*     G(MDG,NB+1) */

/*     Description of calling sequence for DBNDAC.. */

/*     The entire set of parameters for DBNDAC are */

/*     Input.. All Type REAL variables are DOUBLE PRECISION */

/*     G(*,*)            The working array into which the user will */
/*                       place the MT by NB+1 block (C F) in rows IR */
/*                       through IR+MT-1, columns 1 through NB+1. */
/*                       See descriptions of IR and MT below. */

/*     MDG               The number of rows in the working array */
/*                       G(*,*).  The value of MDG should be .GE. MU. */
/*                       The value of MU is defined in the abstract */
/*                       of these subprograms. */

/*     NB                The bandwidth of the data matrix A. */

/*     IP                Set by the user to the value 1 before the */
/*                       first call to DBNDAC.  Its subsequent value */
/*                       is controlled by DBNDAC to set up for the */
/*                       next call to DBNDAC. */

/*     IR                Index of the row of G(*,*) where the user is */
/*                       to place the new block of data (C F).  Set by */
/*                       the user to the value 1 before the first call */
/*                       to DBNDAC.  Its subsequent value is controlled */
/*                       by DBNDAC. A value of IR .GT. MDG is considered */
/*                       an error. */

/*     MT,JT             Set by the user to indicate respectively the */
/*                       number of new rows of data in the block and */
/*                       the index of the first nonzero column in that */
/*                       set of rows (E F) = (0 C 0 F) being processed. */

/*     Output.. All Type REAL variables are DOUBLE PRECISION */

/*     G(*,*)            The working array which will contain the */
/*                       processed rows of that part of the data */
/*                       matrix which has been passed to DBNDAC. */

/*     IP,IR             The values of these arguments are advanced by */
/*                       DBNDAC to be ready for storing and processing */
/*                       a new block of data in G(*,*). */

/*     Description of calling sequence for DBNDSL.. */

/*     The user must dimension the arrays appearing in the call list.. */

/*     G(MDG,NB+1), X(N) */

/*     The entire set of parameters for DBNDSL are */

/*     Input.. All Type REAL variables are DOUBLE PRECISION */

/*     MODE              Set by the user to one of the values 1, 2, or */
/*                       3.  These values respectively indicate that */
/*                       the solution of AX = B, YR = H or RZ = W is */
/*                       required. */

/*     G(*,*),MDG,       These arguments all have the same meaning and */
/*      NB,IP,IR         contents as following the last call to DBNDAC. */

/*     X(*)              With mode=2 or 3 this array contains, */
/*                       respectively, the right-side vectors H or W of */
/*                       the systems YR = H or RZ = W. */

/*     N                 The number of variables in the solution */
/*                       vector.  If any of the N diagonal terms are */
/*                       zero the subroutine DBNDSL prints an */
/*                       appropriate message.  This condition is */
/*                       considered an error. */

/*     Output.. All Type REAL variables are DOUBLE PRECISION */

/*     X(*)              This array contains the solution vectors X, */
/*                       Y or Z of the systems AX = B, YR = H or */
/*                       RZ = W depending on the value of MODE=1, */
/*                       2 or 3. */

/*     RNORM             If MODE=1 RNORM is the Euclidean length of the */
/*                       residual vector AX-B.  When MODE=2 or 3 RNORM */
/*                       is set to zero. */

/*     Remarks.. */

/*     To obtain the upper triangular matrix and transformed right-hand */
/*     side vector D so that the super diagonals of R form the columns */
/*     of G(*,*), execute the following Fortran statements. */

/*     NBP1=NB+1 */

/*     DO 10 J=1, NBP1 */

/*  10 G(IR,J) = 0.E0 */

/*     MT=1 */

/*     JT=N+1 */

/*     CALL DBNDAC(G,MDG,NB,IP,IR,MT,JT) */

/* ***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares */
/*                 Problems, Prentice-Hall, Inc., 1974, Chapter 27. */
/* ***ROUTINES CALLED  DH12, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBNDAC */
/* ***FIRST EXECUTABLE STATEMENT  DBNDAC */
    /* Parameter adjustments */
    g_dim1 = *mdg;
    g_offset = 1 + g_dim1;
    g -= g_offset;

    /* Function Body */
    zero = 0.;

/*              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE. */

    nbp1 = *nb + 1;
    if (*mt <= 0 || *nb <= 0) {
	return 0;
    }

    if (! (*mdg < *ir)) {
	goto L5;
    }
    nerr = 1;
    iopt = 2;
    xermsg_("SLATEC", "DBNDAC", "MDG.LT.IR, PROBABLE ERROR.", &nerr, &iopt, (
	    ftnlen)6, (ftnlen)6, (ftnlen)26);
    return 0;
L5:

/*                                             ALG. STEP 5 */
    if (*jt == *ip) {
	goto L70;
    }
/*                                             ALG. STEPS 6-7 */
    if (*jt <= *ir) {
	goto L30;
    }
/*                                             ALG. STEPS 8-9 */
    i__1 = *mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ig1 = *jt + *mt - i__;
	ig2 = *ir + *mt - i__;
	i__2 = nbp1;
	for (j = 1; j <= i__2; ++j) {
	    g[ig1 + j * g_dim1] = g[ig2 + j * g_dim1];
/* L10: */
	}
    }
/*                                             ALG. STEP 10 */
    ie = *jt - *ir;
    i__2 = ie;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ig = *ir + i__ - 1;
	i__1 = nbp1;
	for (j = 1; j <= i__1; ++j) {
	    g[ig + j * g_dim1] = zero;
/* L20: */
	}
    }
/*                                             ALG. STEP 11 */
    *ir = *jt;
/*                                             ALG. STEP 12 */
L30:
/* Computing MIN */
    i__1 = *nb - 1, i__2 = *ir - *ip - 1;
    mu = min(i__1,i__2);
    if (mu == 0) {
	goto L60;
    }
/*                                             ALG. STEP 13 */
    i__1 = mu;
    for (l = 1; l <= i__1; ++l) {
/*                                             ALG. STEP 14 */
/* Computing MIN */
	i__2 = l, i__3 = *jt - *ip;
	k = min(i__2,i__3);
/*                                             ALG. STEP 15 */
	lp1 = l + 1;
	ig = *ip + l;
	i__2 = *nb;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    jg = i__ - k;
	    g[ig + jg * g_dim1] = g[ig + i__ * g_dim1];
/* L40: */
	}
/*                                             ALG. STEP 16 */
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    jg = nbp1 - i__;
	    g[ig + jg * g_dim1] = zero;
/* L50: */
	}
    }
/*                                             ALG. STEP 17 */
L60:
    *ip = *jt;
/*                                             ALG. STEPS 18-19 */
L70:
    mh = *ir + *mt - *ip;
    kh = min(nbp1,mh);
/*                                             ALG. STEP 20 */
    i__2 = kh;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	i__3 = i__ + 1, i__4 = *ir - *ip + 1;
	i__1 = max(i__3,i__4);
	i__5 = nbp1 - i__;
	dh12_(&c__1, &i__, &i__1, &mh, &g[*ip + i__ * g_dim1], &c__1, &rho, &
		g[*ip + (i__ + 1) * g_dim1], &c__1, mdg, &i__5);
/* L80: */
    }
/*                                             ALG. STEP 21 */
    *ir = *ip + kh;
/*                                             ALG. STEP 22 */
    if (kh < nbp1) {
	goto L100;
    }
/*                                             ALG. STEP 23 */
    i__2 = *nb;
    for (i__ = 1; i__ <= i__2; ++i__) {
	g[*ir - 1 + i__ * g_dim1] = zero;
/* L90: */
    }
/*                                             ALG. STEP 24 */
L100:
/*                                             ALG. STEP 25 */
    return 0;
} /* dbndac_ */


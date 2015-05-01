/* dbndsl.f -- translated by f2c (version 12.02.01).
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

/* DECK DBNDSL */
/* Subroutine */ int dbndsl_(integer *mode, doublereal *g, integer *mdg, 
	integer *nb, integer *ip, integer *ir, doublereal *x, integer *n, 
	doublereal *rnorm)
{
    /* System generated locals */
    integer g_dim1, g_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l;
    static doublereal s;
    static integer i1, i2, ie, jg, ii, ix, np1;
    static doublereal rsq;
    static integer irm1, nerr, iopt;
    static doublereal zero;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBNDSL */
/* ***PURPOSE  Solve the least squares problem for a banded matrix using */
/*            sequential accumulation of rows of the data matrix. */
/*            Exactly one right-hand side vector is permitted. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D9 */
/* ***TYPE      DOUBLE PRECISION (BNDSOL-S, DBNDSL-D) */
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

/*     Input.. */

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

/*     Output.. */

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
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBNDSL */
/* ***FIRST EXECUTABLE STATEMENT  DBNDSL */
    /* Parameter adjustments */
    g_dim1 = *mdg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --x;

    /* Function Body */
    zero = 0.;

    *rnorm = zero;
    switch (*mode) {
	case 1:  goto L10;
	case 2:  goto L90;
	case 3:  goto L50;
    }
/*                                   ********************* MODE = 1 */
/*                                   ALG. STEP 26 */
L10:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = g[j + (*nb + 1) * g_dim1];
/* L20: */
    }
    rsq = zero;
    np1 = *n + 1;
    irm1 = *ir - 1;
    if (np1 > irm1) {
	goto L40;
    }
    i__1 = irm1;
    for (j = np1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = g[j + (*nb + 1) * g_dim1];
	rsq += d__1 * d__1;
/* L30: */
    }
    *rnorm = sqrt(rsq);
L40:
/*                                   ********************* MODE = 3 */
/*                                   ALG. STEP 27 */
L50:
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
/*                                   ALG. STEP 28 */
	s = zero;
/* Computing MAX */
	i__2 = 0, i__3 = i__ - *ip;
	l = max(i__2,i__3);
/*                                   ALG. STEP 29 */
	if (i__ == *n) {
	    goto L70;
	}
/*                                   ALG. STEP 30 */
/* Computing MIN */
	i__2 = *n + 1 - i__;
	ie = min(i__2,*nb);
	i__2 = ie;
	for (j = 2; j <= i__2; ++j) {
	    jg = j + l;
	    ix = i__ - 1 + j;
	    s += g[i__ + jg * g_dim1] * x[ix];
/* L60: */
	}
/*                                   ALG. STEP 31 */
L70:
	if (g[i__ + (l + 1) * g_dim1] != 0.) {
	    goto L80;
	} else {
	    goto L130;
	}
L80:
	x[i__] = (x[i__] - s) / g[i__ + (l + 1) * g_dim1];
    }
/*                                   ALG. STEP 32 */
    return 0;
/*                                   ********************* MODE = 2 */
L90:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s = zero;
	if (j == 1) {
	    goto L110;
	}
/* Computing MAX */
	i__2 = 1, i__3 = j - *nb + 1;
	i1 = max(i__2,i__3);
	i2 = j - 1;
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = 0, i__4 = i__ - *ip;
	    l = j - i__ + 1 + max(i__3,i__4);
	    s += x[i__] * g[i__ + l * g_dim1];
/* L100: */
	}
L110:
/* Computing MAX */
	i__2 = 0, i__3 = j - *ip;
	l = max(i__2,i__3);
	if (g[j + (l + 1) * g_dim1] != 0.) {
	    goto L120;
	} else {
	    goto L130;
	}
L120:
	x[j] = (x[j] - s) / g[j + (l + 1) * g_dim1];
    }
    return 0;

L130:
    nerr = 1;
    iopt = 2;
    xermsg_("SLATEC", "DBNDSL", "A ZERO DIAGONAL TERM IS IN THE N BY N UPPER"
	    " TRIANGULAR MATRIX.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)
	    62);
    return 0;
} /* dbndsl_ */


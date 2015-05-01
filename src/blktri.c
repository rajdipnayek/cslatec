/* blktri.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer npp, k;
    real eps, cnv;
    integer nm, ncmplx, ik;
} cblkt_;

#define cblkt_1 cblkt_

/* DECK BLKTRI */
/* Subroutine */ int blktri_(integer *iflg, integer *np, integer *n, real *an,
	 real *bn, real *cn, integer *mp, integer *m, real *am, real *bm, 
	real *cm, integer *idimy, real *y, integer *ierror, real *w)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer nh, nl, iw1, iw2, iw3, iwd, iwu, iww, iwah, iwbh;
    extern /* Subroutine */ int prod_();
    extern /* Subroutine */ int compb_(integer *, integer *, real *, real *, 
	    real *, real *, real *, real *);
    extern /* Subroutine */ int cprod_(), prodp_();
    extern /* Subroutine */ int blktr1_(integer *, real *, real *, real *, 
	    integer *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, U_fp, U_fp);
    extern /* Subroutine */ int cprodp_();

/* ***BEGIN PROLOGUE  BLKTRI */
/* ***PURPOSE  Solve a block tridiagonal system of linear equations */
/*            (usually resulting from the discretization of separable */
/*            two-dimensional elliptic equations). */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B4B */
/* ***TYPE      SINGLE PRECISION (BLKTRI-S, CBLKTR-C) */
/* ***KEYWORDS  ELLIPTIC PDE, FISHPACK, TRIDIAGONAL LINEAR SYSTEM */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine BLKTRI Solves a System of Linear Equations of the Form */

/*          AN(J)*X(I,J-1) + AM(I)*X(I-1,J) + (BN(J)+BM(I))*X(I,J) */

/*          + CN(J)*X(I,J+1) + CM(I)*X(I+1,J) = Y(I,J) */

/*               for I = 1,2,...,M  and  J = 1,2,...,N. */

/*     I+1 and I-1 are evaluated modulo M and J+1 and J-1 modulo N, i.e., */

/*          X(I,0) = X(I,N),  X(I,N+1) = X(I,1), */
/*          X(0,J) = X(M,J),  X(M+1,J) = X(1,J). */

/*     These equations usually result from the discretization of */
/*     separable elliptic equations.  Boundary conditions may be */
/*     Dirichlet, Neumann, or Periodic. */


/*     * * * * * * * * * *     ON INPUT     * * * * * * * * * * */

/*     IFLG */
/*       = 0  Initialization only.  Certain quantities that depend on NP, */
/*            N, AN, BN, and CN are computed and stored in the work */
/*            array  W. */
/*       = 1  The quantities that were computed in the initialization are */
/*            used to obtain the solution X(I,J). */

/*       NOTE   A call with IFLG=0 takes approximately one half the time */
/*              as a call with IFLG = 1  .  However, the */
/*              initialization does not have to be repeated unless NP, N, */
/*              AN, BN, or CN change. */

/*     NP */
/*       = 0  If AN(1) and CN(N) are not zero, which corresponds to */
/*            periodic boundary conditions. */
/*       = 1  If AN(1) and CN(N) are zero. */

/*     N */
/*       The number of unknowns in the J-direction. N must be greater */
/*       than 4. The operation count is proportional to MNlog2(N), hence */
/*       N should be selected less than or equal to M. */

/*     AN,BN,CN */
/*       One-dimensional arrays of length N that specify the coefficients */
/*       in the linear equations given above. */

/*     MP */
/*       = 0  If AM(1) and CM(M) are not zero, which corresponds to */
/*            periodic boundary conditions. */
/*       = 1  If AM(1) = CM(M) = 0  . */

/*     M */
/*       The number of unknowns in the I-direction. M must be greater */
/*       than 4. */

/*     AM,BM,CM */
/*       One-dimensional arrays of length M that specify the coefficients */
/*       in the linear equations given above. */

/*     IDIMY */
/*       The row (or first) dimension of the two-dimensional array Y as */
/*       it appears in the program calling BLKTRI.  This parameter is */
/*       used to specify the variable dimension of Y.  IDIMY must be at */
/*       least M. */

/*     Y */
/*       A two-dimensional array that specifies the values of the right */
/*       side of the linear system of equations given above.  Y must be */
/*       dimensioned at least M*N. */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space. */
/*             If NP=1 define K=INT(log2(N))+1 and set L=2**(K+1) then */
/*                     W must have dimension (K-2)*L+K+5+MAX(2N,6M) */

/*             If NP=0 define K=INT(log2(N-1))+1 and set L=2**(K+1) then */
/*                     W must have dimension (K-2)*L+K+5+2N+MAX(2N,6M) */

/*       **IMPORTANT** For purposes of checking, the required dimension */
/*                     of W is computed by BLKTRI and stored in W(1) */
/*                     in floating point format. */

/*     * * * * * * * * * *     On Output     * * * * * * * * * * */

/*     Y */
/*       Contains the solution X. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for number zero, a solution is not attempted. */

/*       = 0  No error. */
/*       = 1  M is less than 5. */
/*       = 2  N is less than 5. */
/*       = 3  IDIMY is less than M. */
/*       = 4  BLKTRI failed while computing results that depend on the */
/*            coefficient arrays AN, BN, CN.  Check these arrays. */
/*       = 5  AN(J)*CN(J-1) is less than 0 for some J. Possible reasons */
/*            for this condition are */
/*            1. The arrays AN and CN are not correct. */
/*            2. Too large a grid spacing was used in the discretization */
/*               of the elliptic equation. */
/*            3. The linear equations resulted from a partial */
/*               differential equation which was not elliptic. */

/*     W */
/*       Contains intermediate values that must not be destroyed if */
/*       BLKTRI will be called again with IFLG=1.  W(1) contains the */
/*       number of locations required by W in floating point format. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N) */
/*     Arguments      W(See argument list) */

/*     Latest         June 1979 */
/*     Revision */

/*     Required       BLKTRI,BLKTRI,PROD,PRODP,CPROD,CPRODP,COMPB,INDXA, */
/*     Subprograms    INDXB,INDXC,PPADD,PSGF,PPSGF,PPSPF,BSRH,TEVLS, */
/*                    R1MACH */

/*     Special        The Algorithm may fail if ABS(BM(I)+BN(J)) is less */
/*     Conditions     than ABS(AM(I))+ABS(AN(J))+ABS(CM(I))+ABS(CN(J)) */
/*                    for some I and J. The Algorithm will also fail if */
/*                    AN(J)*CN(J-1) is less than zero for some J. */
/*                    See the description of the output parameter IERROR. */

/*     Common         CBLKT */
/*     Blocks */

/*     I/O            None */

/*     Precision      Single */

/*     Specialist     Paul Swarztrauber */

/*     Language       FORTRAN */

/*     History        Version 1 September 1973 */
/*                    Version 2 April     1976 */
/*                    Version 3 June      1979 */

/*     Algorithm      Generalized Cyclic Reduction (See Reference below) */

/*     Space */
/*     Required       Control Data 7600 */

/*     Portability    American National Standards Institute Fortran. */
/*                    The machine accuracy is set using function R1MACH. */

/*     Required       None */
/*     Resident */
/*     Routines */

/*     References     Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN */
/*                    Subprograms For The Solution Of Elliptic Equations' */
/*                    NCAR TN/IA-109, July, 1975, 138 PP. */

/*                    Swarztrauber P. ,'A Direct Method For The Discrete */
/*                    Solution Of Separable Elliptic Equations', S.I.A.M. */
/*                    J. Numer. Anal.,11(1974) PP. 1136-1150. */

/* ***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran */
/*                 subprograms for the solution of elliptic equations, */
/*                 NCAR TN/IA-109, July 1975, 138 pp. */
/*               P. N. Swarztrauber, A direct method for the discrete */
/*                 solution of separable elliptic equations, SIAM Journal */
/*                 on Numerical Analysis 11, (1974), pp. 1136-1150. */
/* ***ROUTINES CALLED  BLKTR1, COMPB, CPROD, CPRODP, PROD, PRODP */
/* ***COMMON BLOCKS    CBLKT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BLKTRI */

/* ***FIRST EXECUTABLE STATEMENT  BLKTRI */
    /* Parameter adjustments */
    --an;
    --bn;
    --cn;
    --am;
    --bm;
    --cm;
    y_dim1 = *idimy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --w;

    /* Function Body */
    cblkt_1.nm = *n;
    *ierror = 0;
    if (*m - 5 >= 0) {
	goto L102;
    } else {
	goto L101;
    }
L101:
    *ierror = 1;
    goto L119;
L102:
    if (cblkt_1.nm - 3 >= 0) {
	goto L104;
    } else {
	goto L103;
    }
L103:
    *ierror = 2;
    goto L119;
L104:
    if (*idimy - *m >= 0) {
	goto L106;
    } else {
	goto L105;
    }
L105:
    *ierror = 3;
    goto L119;
L106:
    nh = *n;
    cblkt_1.npp = *np;
    if (cblkt_1.npp != 0) {
	goto L107;
    } else {
	goto L108;
    }
L107:
    ++nh;
L108:
    cblkt_1.ik = 2;
    cblkt_1.k = 1;
L109:
    cblkt_1.ik += cblkt_1.ik;
    ++cblkt_1.k;
    if (nh - cblkt_1.ik <= 0) {
	goto L110;
    } else {
	goto L109;
    }
L110:
    nl = cblkt_1.ik;
    cblkt_1.ik += cblkt_1.ik;
    --nl;
    iwah = (cblkt_1.k - 2) * cblkt_1.ik + cblkt_1.k + 6;
    if (cblkt_1.npp != 0) {
	goto L111;
    } else {
	goto L112;
    }

/*     DIVIDE W INTO WORKING SUB ARRAYS */

L111:
    iw1 = iwah;
    iwbh = iw1 + cblkt_1.nm;
/* Computing MAX */
    i__1 = cblkt_1.nm << 1, i__2 = *m * 6;
    w[1] = (real) (iw1 - 1 + max(i__1,i__2));
    goto L113;
L112:
    iwbh = iwah + cblkt_1.nm + cblkt_1.nm;
    iw1 = iwbh;
/* Computing MAX */
    i__1 = cblkt_1.nm << 1, i__2 = *m * 6;
    w[1] = (real) (iw1 - 1 + max(i__1,i__2));
    --cblkt_1.nm;

/* SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS */

L113:
    if (*ierror != 0) {
	goto L119;
    } else {
	goto L114;
    }
L114:
    iw2 = iw1 + *m;
    iw3 = iw2 + *m;
    iwd = iw3 + *m;
    iww = iwd + *m;
    iwu = iww + *m;
    if (*iflg != 0) {
	goto L116;
    } else {
	goto L115;
    }
L115:
    compb_(&nl, ierror, &an[1], &bn[1], &cn[1], &w[2], &w[iwah], &w[iwbh]);
    goto L119;
L116:
    if (*mp != 0) {
	goto L117;
    } else {
	goto L118;
    }

/* SUBROUTINE BLKTR1 SOLVES THE LINEAR SYSTEM */

L117:
    blktr1_(&nl, &an[1], &bn[1], &cn[1], m, &am[1], &bm[1], &cm[1], idimy, &y[
	    y_offset], &w[2], &w[iw1], &w[iw2], &w[iw3], &w[iwd], &w[iww], &w[
	    iwu], (U_fp)prod_, (U_fp)cprod_);
    goto L119;
L118:
    blktr1_(&nl, &an[1], &bn[1], &cn[1], m, &am[1], &bm[1], &cm[1], idimy, &y[
	    y_offset], &w[2], &w[iw1], &w[iw2], &w[iw3], &w[iwd], &w[iww], &w[
	    iwu], (U_fp)prodp_, (U_fp)cprodp_);
L119:
    return 0;
} /* blktri_ */


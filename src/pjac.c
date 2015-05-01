/* pjac.f -- translated by f2c (version 12.02.01).
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
    real rownd, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iownd[14], iowns[6], ier, jstart, kflag, l, meth, miter, maxord, 
	    n, nq, nst, nfe, nje, nqu;
} debdf1_;

#define debdf1_1 debdf1_

/* DECK PJAC */
/* Subroutine */ int pjac_(integer *neq, real *y, real *yh, integer *nyh, 
	real *ewt, real *ftem, real *savf, real *wm, integer *iwm, S_fp f, 
	S_fp jac, real *rpar, integer *ipar)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j;
    static real r__;
    static integer i1, i2, j1;
    static real r0, di;
    static integer ii, jj, ml, mu;
    static real yi, yj, hl0;
    static integer ml3;
    static real fac;
    static integer mba;
    static real con, yjj;
    static integer meb1, lenp;
    static real srur;
    static integer mband;
    extern /* Subroutine */ int sgbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sgefa_(real *, integer *, 
	    integer *, integer *, integer *);
    static integer meband;
    extern doublereal vnwrms_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  PJAC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PJAC-S, DPJAC-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   PJAC sets up the iteration matrix (involving the Jacobian) for the */
/*   integration package DEBDF. */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  SGBFA, SGEFA, VNWRMS */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920422  Changed DIMENSION statement.  (WRB) */
/* ***END PROLOGUE  PJAC */

/* LLL. OPTIMIZE */
/* ----------------------------------------------------------------------- */
/* PJAC IS CALLED BY STOD  TO COMPUTE AND PROCESS THE MATRIX */
/* P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN. */
/* HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE JAC IF */
/* MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5. */
/* IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED. */
/* J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN */
/* SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION */
/* OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE */
/* BY SGEFA IF MITER = 1 OR 2, AND BY SGBFA IF MITER = 4 OR 5. */

/* IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION */
/* WITH PJAC USES THE FOLLOWING.. */
/* Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY. */
/* FTEM = WORK ARRAY OF LENGTH N (ACOR IN STOD ). */
/* SAVF = ARRAY CONTAINING F EVALUATED AT PREDICTED Y. */
/* WM   = REAL WORK SPACE FOR MATRICES.  ON OUTPUT IT CONTAINS THE */
/*        INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU DECOMPOSITION */
/*        OF P IF MITER IS 1, 2 , 4, OR 5. */
/*        STORAGE OF MATRIX ELEMENTS STARTS AT WM(3). */
/*        WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA.. */
/*        WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS. */
/*        WM(2) = H*EL0, SAVED FOR LATER USE IF MITER = 3. */
/* IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT */
/*        IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS THE */
/*        BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS 4 OR 5. */
/* EL0  = EL(1) (INPUT). */
/* IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF */
/*        P MATRIX FOUND TO BE SINGULAR. */
/* THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND, */
/* MITER, N, NFE, AND NJE. */
/* ----------------------------------------------------------------------- */
/* ***FIRST EXECUTABLE STATEMENT  PJAC */
    /* Parameter adjustments */
    --y;
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --ewt;
    --ftem;
    --savf;
    --wm;
    --iwm;
    --rpar;
    --ipar;

    /* Function Body */
    ++debdf1_1.nje;
    hl0 = debdf1_1.h__ * debdf1_1.el0;
    switch (debdf1_1.miter) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
    }
/* IF MITER = 1, CALL JAC AND MULTIPLY BY SCALAR. ----------------------- */
L100:
    lenp = debdf1_1.n * debdf1_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	wm[i__ + 2] = 0.f;
    }
    (*jac)(&debdf1_1.tn, &y[1], &wm[3], &debdf1_1.n, &rpar[1], &ipar[1]);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	wm[i__ + 2] *= con;
    }
    goto L240;
/* IF MITER = 2, MAKE N CALLS TO F TO APPROXIMATE J. -------------------- */
L200:
    fac = vnwrms_(&debdf1_1.n, &savf[1], &ewt[1]);
    r0 = dabs(debdf1_1.h__) * 1e3f * debdf1_1.uround * debdf1_1.n * fac;
    if (r0 == 0.f) {
	r0 = 1.f;
    }
    srur = wm[1];
    j1 = 2;
    i__1 = debdf1_1.n;
    for (j = 1; j <= i__1; ++j) {
	yj = y[j];
/* Computing MAX */
	r__1 = srur * dabs(yj), r__2 = r0 * ewt[j];
	r__ = dmax(r__1,r__2);
	y[j] += r__;
	fac = -hl0 / r__;
	(*f)(&debdf1_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[1]);
	i__2 = debdf1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L220: */
	    wm[i__ + j1] = (ftem[i__] - savf[i__]) * fac;
	}
	y[j] = yj;
	j1 += debdf1_1.n;
/* L230: */
    }
    debdf1_1.nfe += debdf1_1.n;
/* ADD IDENTITY MATRIX. ------------------------------------------------- */
L240:
    j = 3;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[j] += 1.f;
/* L250: */
	j += debdf1_1.n + 1;
    }
/* DO LU DECOMPOSITION ON P. -------------------------------------------- */
    sgefa_(&wm[3], &debdf1_1.n, &debdf1_1.n, &iwm[21], &debdf1_1.ier);
    return 0;
/* IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND P. --------- */
L300:
    wm[2] = hl0;
    debdf1_1.ier = 0;
    r__ = debdf1_1.el0 * .1f;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L310: */
	y[i__] += r__ * (debdf1_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)]);
    }
    (*f)(&debdf1_1.tn, &y[1], &wm[3], &rpar[1], &ipar[1]);
    ++debdf1_1.nfe;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r0 = debdf1_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
	di = r0 * .1f - debdf1_1.h__ * (wm[i__ + 2] - savf[i__]);
	wm[i__ + 2] = 1.f;
	if (dabs(r0) < debdf1_1.uround * ewt[i__]) {
	    goto L320;
	}
	if (dabs(di) == 0.f) {
	    goto L330;
	}
	wm[i__ + 2] = r0 * .1f / di;
L320:
	;
    }
    return 0;
L330:
    debdf1_1.ier = -1;
    return 0;
/* IF MITER = 4, CALL JAC AND MULTIPLY BY SCALAR. ----------------------- */
L400:
    ml = iwm[1];
    mu = iwm[2];
    ml3 = 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * debdf1_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	wm[i__ + 2] = 0.f;
    }
    (*jac)(&debdf1_1.tn, &y[1], &wm[ml3], &meband, &rpar[1], &ipar[1]);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L420: */
	wm[i__ + 2] *= con;
    }
    goto L570;
/* IF MITER = 5, MAKE MBAND CALLS TO F TO APPROXIMATE J. ---------------- */
L500:
    ml = iwm[1];
    mu = iwm[2];
    mband = ml + mu + 1;
    mba = min(mband,debdf1_1.n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = wm[1];
    fac = vnwrms_(&debdf1_1.n, &savf[1], &ewt[1]);
    r0 = dabs(debdf1_1.h__) * 1e3f * debdf1_1.uround * debdf1_1.n * fac;
    if (r0 == 0.f) {
	r0 = 1.f;
    }
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = debdf1_1.n;
	i__3 = mband;
	for (i__ = j; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    yi = y[i__];
/* Computing MAX */
	    r__1 = srur * dabs(yi), r__2 = r0 * ewt[i__];
	    r__ = dmax(r__1,r__2);
/* L530: */
	    y[i__] += r__;
	}
	(*f)(&debdf1_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[1]);
	i__3 = debdf1_1.n;
	i__2 = mband;
	for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    y[jj] = yh[jj + yh_dim1];
	    yjj = y[jj];
/* Computing MAX */
	    r__1 = srur * dabs(yjj), r__2 = r0 * ewt[jj];
	    r__ = dmax(r__1,r__2);
	    fac = -hl0 / r__;
/* Computing MAX */
	    i__4 = jj - mu;
	    i1 = max(i__4,1);
/* Computing MIN */
	    i__4 = jj + ml;
	    i2 = min(i__4,debdf1_1.n);
	    ii = jj * meb1 - ml + 2;
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* L540: */
		wm[ii + i__] = (ftem[i__] - savf[i__]) * fac;
	    }
/* L550: */
	}
/* L560: */
    }
    debdf1_1.nfe += mba;
/* ADD IDENTITY MATRIX. ------------------------------------------------- */
L570:
    ii = mband + 2;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[ii] += 1.f;
/* L580: */
	ii += meband;
    }
/* DO LU DECOMPOSITION OF P. -------------------------------------------- */
    sgbfa_(&wm[3], &meband, &debdf1_1.n, &ml, &mu, &iwm[21], &debdf1_1.ier);
    return 0;
/* ----------------------- END OF SUBROUTINE PJAC ----------------------- */
} /* pjac_ */


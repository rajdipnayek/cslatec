/* dpjac.f -- translated by f2c (version 12.02.01).
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
    doublereal rownd, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iownd[14], iowns[6], ier, jstart, kflag, l, meth, miter, maxord, 
	    n, nq, nst, nfe, nje, nqu;
} ddebd1_;

#define ddebd1_1 ddebd1_

/* DECK DPJAC */
/* Subroutine */ int dpjac_(integer *neq, doublereal *y, doublereal *yh, 
	integer *nyh, doublereal *ewt, doublereal *ftem, doublereal *savf, 
	doublereal *wm, integer *iwm, S_fp df, S_fp djac, doublereal *rpar, 
	integer *ipar)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal r__;
    static integer i1, i2, j1;
    static doublereal r0, di;
    static integer ii, jj, ml, mu;
    static doublereal yi, yj, hl0;
    static integer ml3;
    static doublereal fac;
    static integer mba;
    static doublereal con, yjj;
    static integer meb1, lenp;
    static doublereal srur;
    extern /* Subroutine */ int dgbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), dgefa_(doublereal *, 
	    integer *, integer *, integer *, integer *);
    static integer mband, meband;
    extern doublereal dvnrms_(integer *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  DPJAC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (PJAC-S, DPJAC-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DPJAC sets up the iteration matrix (involving the Jacobian) for the */
/*   integration package DDEBDF. */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  DGBFA, DGEFA, DVNRMS */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920422  Changed DIMENSION statement.  (WRB) */
/* ***END PROLOGUE  DPJAC */

/*     ------------------------------------------------------------------ */
/*      DPJAC IS CALLED BY DSTOD  TO COMPUTE AND PROCESS THE MATRIX */
/*      P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN. */
/*      HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE DJAC IF */
/*      MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5. */
/*      IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED. */
/*      J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN */
/*      SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION */
/*      OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE */
/*      BY DGEFA IF MITER = 1 OR 2, AND BY DGBFA IF MITER = 4 OR 5. */

/*      IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION */
/*      WITH DPJAC USES THE FOLLOWING.. */
/*      Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY. */
/*      FTEM = WORK ARRAY OF LENGTH N (ACOR IN DSTOD ). */
/*      SAVF = ARRAY CONTAINING DF EVALUATED AT PREDICTED Y. */
/*      WM   = DOUBLE PRECISION WORK SPACE FOR MATRICES.  ON OUTPUT IT */
/*      CONTAINS THE */
/*             INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU */
/*             DECOMPOSITION OF P IF MITER IS 1, 2 , 4, OR 5. */
/*             STORAGE OF MATRIX ELEMENTS STARTS AT WM(3). */
/*             WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA.. */
/*             WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN */
/*             INCREMENTS.  WM(2) = H*EL0, SAVED FOR LATER USE IF MITER = */
/*             3. */
/*      IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING */
/*             AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS */
/*             THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER */
/*             IS 4 OR 5. */
/*      EL0  = EL(1) (INPUT). */
/*      IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF */
/*             P MATRIX FOUND TO BE SINGULAR. */
/*      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND, */
/*      MITER, N, NFE, AND NJE. */
/* ----------------------------------------------------------------------- */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 240 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 220 */
/*           BEGIN BLOCK PERMITTING ...EXITS TO 130 */
/*              BEGIN BLOCK PERMITTING ...EXITS TO 70 */
/* ***FIRST EXECUTABLE STATEMENT  DPJAC */
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
    ++ddebd1_1.nje;
    hl0 = ddebd1_1.h__ * ddebd1_1.el0;
    switch (ddebd1_1.miter) {
	case 1:  goto L10;
	case 2:  goto L40;
	case 3:  goto L90;
	case 4:  goto L140;
	case 5:  goto L170;
    }
/*                 IF MITER = 1, CALL DJAC AND MULTIPLY BY SCALAR. */
/*                 ----------------------- */
L10:
    lenp = ddebd1_1.n * ddebd1_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[i__ + 2] = 0.;
/* L20: */
    }
    (*djac)(&ddebd1_1.tn, &y[1], &wm[3], &ddebd1_1.n, &rpar[1], &ipar[1]);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[i__ + 2] *= con;
/* L30: */
    }
/*              ...EXIT */
    goto L70;
/*                 IF MITER = 2, MAKE N CALLS TO DF TO APPROXIMATE J. */
/*                 -------------------- */
L40:
    fac = dvnrms_(&ddebd1_1.n, &savf[1], &ewt[1]);
    r0 = abs(ddebd1_1.h__) * 1e3 * ddebd1_1.uround * ddebd1_1.n * fac;
    if (r0 == 0.) {
	r0 = 1.;
    }
    srur = wm[1];
    j1 = 2;
    i__1 = ddebd1_1.n;
    for (j = 1; j <= i__1; ++j) {
	yj = y[j];
/* Computing MAX */
	d__1 = srur * abs(yj), d__2 = r0 * ewt[j];
	r__ = max(d__1,d__2);
	y[j] += r__;
	fac = -hl0 / r__;
	(*df)(&ddebd1_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[1]);
	i__2 = ddebd1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wm[i__ + j1] = (ftem[i__] - savf[i__]) * fac;
/* L50: */
	}
	y[j] = yj;
	j1 += ddebd1_1.n;
/* L60: */
    }
    ddebd1_1.nfe += ddebd1_1.n;
L70:
/*              ADD IDENTITY MATRIX. */
/*              ------------------------------------------------- */
    j = 3;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[j] += 1.;
	j += ddebd1_1.n + 1;
/* L80: */
    }
/*              DO LU DECOMPOSITION ON P. */
/*              -------------------------------------------- */
    dgefa_(&wm[3], &ddebd1_1.n, &ddebd1_1.n, &iwm[21], &ddebd1_1.ier);
/*     .........EXIT */
    goto L240;
/*              IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND */
/*              P. --------- */
L90:
    wm[2] = hl0;
    ddebd1_1.ier = 0;
    r__ = ddebd1_1.el0 * .1;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] += r__ * (ddebd1_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)]);
/* L100: */
    }
    (*df)(&ddebd1_1.tn, &y[1], &wm[3], &rpar[1], &ipar[1]);
    ++ddebd1_1.nfe;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r0 = ddebd1_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
	di = r0 * .1 - ddebd1_1.h__ * (wm[i__ + 2] - savf[i__]);
	wm[i__ + 2] = 1.;
	if (abs(r0) < ddebd1_1.uround * ewt[i__]) {
	    goto L110;
	}
/*           .........EXIT */
	if (abs(di) == 0.) {
	    goto L130;
	}
	wm[i__ + 2] = r0 * .1 / di;
L110:
/* L120: */
	;
    }
/*     .........EXIT */
    goto L240;
L130:
    ddebd1_1.ier = -1;
/*     ......EXIT */
    goto L240;
/*           IF MITER = 4, CALL DJAC AND MULTIPLY BY SCALAR. */
/*           ----------------------- */
L140:
    ml = iwm[1];
    mu = iwm[2];
    ml3 = 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * ddebd1_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[i__ + 2] = 0.;
/* L150: */
    }
    (*djac)(&ddebd1_1.tn, &y[1], &wm[ml3], &meband, &rpar[1], &ipar[1]);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[i__ + 2] *= con;
/* L160: */
    }
/*        ...EXIT */
    goto L220;
/*           IF MITER = 5, MAKE MBAND CALLS TO DF TO APPROXIMATE J. */
/*           ---------------- */
L170:
    ml = iwm[1];
    mu = iwm[2];
    mband = ml + mu + 1;
    mba = min(mband,ddebd1_1.n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = wm[1];
    fac = dvnrms_(&ddebd1_1.n, &savf[1], &ewt[1]);
    r0 = abs(ddebd1_1.h__) * 1e3 * ddebd1_1.uround * ddebd1_1.n * fac;
    if (r0 == 0.) {
	r0 = 1.;
    }
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ddebd1_1.n;
	i__3 = mband;
	for (i__ = j; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    yi = y[i__];
/* Computing MAX */
	    d__1 = srur * abs(yi), d__2 = r0 * ewt[i__];
	    r__ = max(d__1,d__2);
	    y[i__] += r__;
/* L180: */
	}
	(*df)(&ddebd1_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[1]);
	i__3 = ddebd1_1.n;
	i__2 = mband;
	for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    y[jj] = yh[jj + yh_dim1];
	    yjj = y[jj];
/* Computing MAX */
	    d__1 = srur * abs(yjj), d__2 = r0 * ewt[jj];
	    r__ = max(d__1,d__2);
	    fac = -hl0 / r__;
/* Computing MAX */
	    i__4 = jj - mu;
	    i1 = max(i__4,1);
/* Computing MIN */
	    i__4 = jj + ml;
	    i2 = min(i__4,ddebd1_1.n);
	    ii = jj * meb1 - ml + 2;
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
		wm[ii + i__] = (ftem[i__] - savf[i__]) * fac;
/* L190: */
	    }
/* L200: */
	}
/* L210: */
    }
    ddebd1_1.nfe += mba;
L220:
/*        ADD IDENTITY MATRIX. */
/*        ------------------------------------------------- */
    ii = mband + 2;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[ii] += 1.;
	ii += meband;
/* L230: */
    }
/*        DO LU DECOMPOSITION OF P. */
/*        -------------------------------------------- */
    dgbfa_(&wm[3], &meband, &ddebd1_1.n, &ml, &mu, &iwm[21], &ddebd1_1.ier);
L240:
    return 0;
/*     ----------------------- END OF SUBROUTINE DPJAC */
/*     ----------------------- */
} /* dpjac_ */


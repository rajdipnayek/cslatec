/* dx.f -- translated by f2c (version 12.02.01).
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
    integer kswx, kswy, k, l;
    real ait, bit, cit, dit;
    integer mit, nit, is, ms, js, ns;
    real dlx, dly, tdlx3, tdly3, dlx4, dly4;
} splpcm_;

#define splpcm_1 splpcm_

/* DECK DX */
/* Subroutine */ int dx_(real *u, integer *idmn, integer *i__, integer *j, 
	real *uxxx, real *uxxxx)
{
    /* System generated locals */
    integer u_dim1, u_offset;

/* ***BEGIN PROLOGUE  DX */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPELI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DX-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This program computes second order finite difference */
/*     approximations to the third and fourth X */
/*     partial derivatives of U at the (I,J) mesh point. */

/* ***SEE ALSO  SEPELI */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPLPCM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DX */

/* ***FIRST EXECUTABLE STATEMENT  DX */
    /* Parameter adjustments */
    u_dim1 = *idmn;
    u_offset = 1 + u_dim1;
    u -= u_offset;

    /* Function Body */
    if (*i__ > 2 && *i__ < splpcm_1.k - 1) {
	goto L50;
    }
    if (*i__ == 1) {
	goto L10;
    }
    if (*i__ == 2) {
	goto L30;
    }
    if (*i__ == splpcm_1.k - 1) {
	goto L60;
    }
    if (*i__ == splpcm_1.k) {
	goto L80;
    }

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A */

L10:
    if (splpcm_1.kswx == 1) {
	goto L20;
    }
    *uxxx = (u[*j * u_dim1 + 1] * -5.f + u[*j * u_dim1 + 2] * 18.f - u[*j * 
	    u_dim1 + 3] * 24.f + u[*j * u_dim1 + 4] * 14.f - u[*j * u_dim1 + 
	    5] * 3.f) / splpcm_1.tdlx3;
    *uxxxx = (u[*j * u_dim1 + 1] * 3.f - u[*j * u_dim1 + 2] * 14.f + u[*j * 
	    u_dim1 + 3] * 26.f - u[*j * u_dim1 + 4] * 24.f + u[*j * u_dim1 + 
	    5] * 11.f - u[*j * u_dim1 + 6] * 2.f) / splpcm_1.dlx4;
    return 0;

/*     PERIODIC AT X=A */

L20:
    *uxxx = (-u[splpcm_1.k - 2 + *j * u_dim1] + u[splpcm_1.k - 1 + *j * 
	    u_dim1] * 2.f - u[*j * u_dim1 + 2] * 2.f + u[*j * u_dim1 + 3]) / 
	    splpcm_1.tdlx3;
    *uxxxx = (u[splpcm_1.k - 2 + *j * u_dim1] - u[splpcm_1.k - 1 + *j * 
	    u_dim1] * 4.f + u[*j * u_dim1 + 1] * 6.f - u[*j * u_dim1 + 2] * 
	    4.f + u[*j * u_dim1 + 3]) / splpcm_1.dlx4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX */

L30:
    if (splpcm_1.kswx == 1) {
	goto L40;
    }
    *uxxx = (u[*j * u_dim1 + 1] * -3.f + u[*j * u_dim1 + 2] * 10.f - u[*j * 
	    u_dim1 + 3] * 12.f + u[*j * u_dim1 + 4] * 6.f - u[*j * u_dim1 + 5]
	    ) / splpcm_1.tdlx3;
    *uxxxx = (u[*j * u_dim1 + 1] * 2.f - u[*j * u_dim1 + 2] * 9.f + u[*j * 
	    u_dim1 + 3] * 16.f - u[*j * u_dim1 + 4] * 14.f + u[*j * u_dim1 + 
	    5] * 6.f - u[*j * u_dim1 + 6]) / splpcm_1.dlx4;
    return 0;

/*     PERIODIC AT X=A+DLX */

L40:
    *uxxx = (-u[splpcm_1.k - 1 + *j * u_dim1] + u[*j * u_dim1 + 1] * 2.f - u[*
	    j * u_dim1 + 3] * 2.f + u[*j * u_dim1 + 4]) / splpcm_1.tdlx3;
    *uxxxx = (u[splpcm_1.k - 1 + *j * u_dim1] - u[*j * u_dim1 + 1] * 4.f + u[*
	    j * u_dim1 + 2] * 6.f - u[*j * u_dim1 + 3] * 4.f + u[*j * u_dim1 
	    + 4]) / splpcm_1.dlx4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR */

L50:
    *uxxx = (-u[*i__ - 2 + *j * u_dim1] + u[*i__ - 1 + *j * u_dim1] * 2.f - u[
	    *i__ + 1 + *j * u_dim1] * 2.f + u[*i__ + 2 + *j * u_dim1]) / 
	    splpcm_1.tdlx3;
    *uxxxx = (u[*i__ - 2 + *j * u_dim1] - u[*i__ - 1 + *j * u_dim1] * 4.f + u[
	    *i__ + *j * u_dim1] * 6.f - u[*i__ + 1 + *j * u_dim1] * 4.f + u[*
	    i__ + 2 + *j * u_dim1]) / splpcm_1.dlx4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX */

L60:
    if (splpcm_1.kswx == 1) {
	goto L70;
    }
    *uxxx = (u[splpcm_1.k - 4 + *j * u_dim1] - u[splpcm_1.k - 3 + *j * u_dim1]
	     * 6.f + u[splpcm_1.k - 2 + *j * u_dim1] * 12.f - u[splpcm_1.k - 
	    1 + *j * u_dim1] * 10.f + u[splpcm_1.k + *j * u_dim1] * 3.f) / 
	    splpcm_1.tdlx3;
    *uxxxx = (-u[splpcm_1.k - 5 + *j * u_dim1] + u[splpcm_1.k - 4 + *j * 
	    u_dim1] * 6.f - u[splpcm_1.k - 3 + *j * u_dim1] * 14.f + u[
	    splpcm_1.k - 2 + *j * u_dim1] * 16.f - u[splpcm_1.k - 1 + *j * 
	    u_dim1] * 9.f + u[splpcm_1.k + *j * u_dim1] * 2.f) / 
	    splpcm_1.dlx4;
    return 0;

/*     PERIODIC AT X=B-DLX */

L70:
    *uxxx = (-u[splpcm_1.k - 3 + *j * u_dim1] + u[splpcm_1.k - 2 + *j * 
	    u_dim1] * 2.f - u[*j * u_dim1 + 1] * 2.f + u[*j * u_dim1 + 2]) / 
	    splpcm_1.tdlx3;
    *uxxxx = (u[splpcm_1.k - 3 + *j * u_dim1] - u[splpcm_1.k - 2 + *j * 
	    u_dim1] * 4.f + u[splpcm_1.k - 1 + *j * u_dim1] * 6.f - u[*j * 
	    u_dim1 + 1] * 4.f + u[*j * u_dim1 + 2]) / splpcm_1.dlx4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B */

L80:
    *uxxx = -(u[splpcm_1.k - 4 + *j * u_dim1] * 3.f - u[splpcm_1.k - 3 + *j * 
	    u_dim1] * 14.f + u[splpcm_1.k - 2 + *j * u_dim1] * 24.f - u[
	    splpcm_1.k - 1 + *j * u_dim1] * 18.f + u[splpcm_1.k + *j * u_dim1]
	     * 5.f) / splpcm_1.tdlx3;
    *uxxxx = (u[splpcm_1.k - 5 + *j * u_dim1] * -2.f + u[splpcm_1.k - 4 + *j *
	     u_dim1] * 11.f - u[splpcm_1.k - 3 + *j * u_dim1] * 24.f + u[
	    splpcm_1.k - 2 + *j * u_dim1] * 26.f - u[splpcm_1.k - 1 + *j * 
	    u_dim1] * 14.f + u[splpcm_1.k + *j * u_dim1] * 3.f) / 
	    splpcm_1.dlx4;
    return 0;
} /* dx_ */


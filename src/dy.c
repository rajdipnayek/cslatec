/* dy.f -- translated by f2c (version 12.02.01).
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

/* DECK DY */
/* Subroutine */ int dy_(real *u, integer *idmn, integer *i__, integer *j, 
	real *uyyy, real *uyyyy)
{
    /* System generated locals */
    integer u_dim1, u_offset;

/* ***BEGIN PROLOGUE  DY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPELI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DY-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This program computes second order finite difference */
/*     approximations to the third and fourth Y */
/*     partial derivatives of U at the (I,J) mesh point. */

/* ***SEE ALSO  SEPELI */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPLPCM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DY */

/* ***FIRST EXECUTABLE STATEMENT  DY */
    /* Parameter adjustments */
    u_dim1 = *idmn;
    u_offset = 1 + u_dim1;
    u -= u_offset;

    /* Function Body */
    if (*j > 2 && *j < splpcm_1.l - 1) {
	goto L50;
    }
    if (*j == 1) {
	goto L10;
    }
    if (*j == 2) {
	goto L30;
    }
    if (*j == splpcm_1.l - 1) {
	goto L60;
    }
    if (*j == splpcm_1.l) {
	goto L80;
    }

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C */

L10:
    if (splpcm_1.kswy == 1) {
	goto L20;
    }
    *uyyy = (u[*i__ + u_dim1] * -5.f + u[*i__ + (u_dim1 << 1)] * 18.f - u[*
	    i__ + u_dim1 * 3] * 24.f + u[*i__ + (u_dim1 << 2)] * 14.f - u[*
	    i__ + u_dim1 * 5] * 3.f) / splpcm_1.tdly3;
    *uyyyy = (u[*i__ + u_dim1] * 3.f - u[*i__ + (u_dim1 << 1)] * 14.f + u[*
	    i__ + u_dim1 * 3] * 26.f - u[*i__ + (u_dim1 << 2)] * 24.f + u[*
	    i__ + u_dim1 * 5] * 11.f - u[*i__ + u_dim1 * 6] * 2.f) / 
	    splpcm_1.dly4;
    return 0;

/*     PERIODIC AT X=A */

L20:
    *uyyy = (-u[*i__ + (splpcm_1.l - 2) * u_dim1] + u[*i__ + (splpcm_1.l - 1) 
	    * u_dim1] * 2.f - u[*i__ + (u_dim1 << 1)] * 2.f + u[*i__ + u_dim1 
	    * 3]) / splpcm_1.tdly3;
    *uyyyy = (u[*i__ + (splpcm_1.l - 2) * u_dim1] - u[*i__ + (splpcm_1.l - 1) 
	    * u_dim1] * 4.f + u[*i__ + u_dim1] * 6.f - u[*i__ + (u_dim1 << 1)]
	     * 4.f + u[*i__ + u_dim1 * 3]) / splpcm_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY */

L30:
    if (splpcm_1.kswy == 1) {
	goto L40;
    }
    *uyyy = (u[*i__ + u_dim1] * -3.f + u[*i__ + (u_dim1 << 1)] * 10.f - u[*
	    i__ + u_dim1 * 3] * 12.f + u[*i__ + (u_dim1 << 2)] * 6.f - u[*i__ 
	    + u_dim1 * 5]) / splpcm_1.tdly3;
    *uyyyy = (u[*i__ + u_dim1] * 2.f - u[*i__ + (u_dim1 << 1)] * 9.f + u[*i__ 
	    + u_dim1 * 3] * 16.f - u[*i__ + (u_dim1 << 2)] * 14.f + u[*i__ + 
	    u_dim1 * 5] * 6.f - u[*i__ + u_dim1 * 6]) / splpcm_1.dly4;
    return 0;

/*     PERIODIC AT Y=C+DLY */

L40:
    *uyyy = (-u[*i__ + (splpcm_1.l - 1) * u_dim1] + u[*i__ + u_dim1] * 2.f - 
	    u[*i__ + u_dim1 * 3] * 2.f + u[*i__ + (u_dim1 << 2)]) / 
	    splpcm_1.tdly3;
    *uyyyy = (u[*i__ + (splpcm_1.l - 1) * u_dim1] - u[*i__ + u_dim1] * 4.f + 
	    u[*i__ + (u_dim1 << 1)] * 6.f - u[*i__ + u_dim1 * 3] * 4.f + u[*
	    i__ + (u_dim1 << 2)]) / splpcm_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR */

L50:
    *uyyy = (-u[*i__ + (*j - 2) * u_dim1] + u[*i__ + (*j - 1) * u_dim1] * 2.f 
	    - u[*i__ + (*j + 1) * u_dim1] * 2.f + u[*i__ + (*j + 2) * u_dim1])
	     / splpcm_1.tdly3;
    *uyyyy = (u[*i__ + (*j - 2) * u_dim1] - u[*i__ + (*j - 1) * u_dim1] * 4.f 
	    + u[*i__ + *j * u_dim1] * 6.f - u[*i__ + (*j + 1) * u_dim1] * 4.f 
	    + u[*i__ + (*j + 2) * u_dim1]) / splpcm_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY */

L60:
    if (splpcm_1.kswy == 1) {
	goto L70;
    }
    *uyyy = (u[*i__ + (splpcm_1.l - 4) * u_dim1] - u[*i__ + (splpcm_1.l - 3) *
	     u_dim1] * 6.f + u[*i__ + (splpcm_1.l - 2) * u_dim1] * 12.f - u[*
	    i__ + (splpcm_1.l - 1) * u_dim1] * 10.f + u[*i__ + splpcm_1.l * 
	    u_dim1] * 3.f) / splpcm_1.tdly3;
    *uyyyy = (-u[*i__ + (splpcm_1.l - 5) * u_dim1] + u[*i__ + (splpcm_1.l - 4)
	     * u_dim1] * 6.f - u[*i__ + (splpcm_1.l - 3) * u_dim1] * 14.f + u[
	    *i__ + (splpcm_1.l - 2) * u_dim1] * 16.f - u[*i__ + (splpcm_1.l - 
	    1) * u_dim1] * 9.f + u[*i__ + splpcm_1.l * u_dim1] * 2.f) / 
	    splpcm_1.dly4;
    return 0;

/*     PERIODIC AT Y=D-DLY */

L70:
    *uyyy = (-u[*i__ + (splpcm_1.l - 3) * u_dim1] + u[*i__ + (splpcm_1.l - 2) 
	    * u_dim1] * 2.f - u[*i__ + u_dim1] * 2.f + u[*i__ + (u_dim1 << 1)]
	    ) / splpcm_1.tdly3;
    *uyyyy = (u[*i__ + (splpcm_1.l - 3) * u_dim1] - u[*i__ + (splpcm_1.l - 2) 
	    * u_dim1] * 4.f + u[*i__ + (splpcm_1.l - 1) * u_dim1] * 6.f - u[*
	    i__ + u_dim1] * 4.f + u[*i__ + (u_dim1 << 1)]) / splpcm_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D */

L80:
    *uyyy = -(u[*i__ + (splpcm_1.l - 4) * u_dim1] * 3.f - u[*i__ + (
	    splpcm_1.l - 3) * u_dim1] * 14.f + u[*i__ + (splpcm_1.l - 2) * 
	    u_dim1] * 24.f - u[*i__ + (splpcm_1.l - 1) * u_dim1] * 18.f + u[*
	    i__ + splpcm_1.l * u_dim1] * 5.f) / splpcm_1.tdly3;
    *uyyyy = (u[*i__ + (splpcm_1.l - 5) * u_dim1] * -2.f + u[*i__ + (
	    splpcm_1.l - 4) * u_dim1] * 11.f - u[*i__ + (splpcm_1.l - 3) * 
	    u_dim1] * 24.f + u[*i__ + (splpcm_1.l - 2) * u_dim1] * 26.f - u[*
	    i__ + (splpcm_1.l - 1) * u_dim1] * 14.f + u[*i__ + splpcm_1.l * 
	    u_dim1] * 3.f) / splpcm_1.dly4;
    return 0;
} /* dy_ */


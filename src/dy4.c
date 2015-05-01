/* dy4.f -- translated by f2c (version 12.02.01).
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
} spl4_;

#define spl4_1 spl4_

/* DECK DY4 */
/* Subroutine */ int dy4_(real *u, integer *idmn, integer *i__, integer *j, 
	real *uyyy, real *uyyyy)
{
    /* System generated locals */
    integer u_dim1, u_offset;

/* ***BEGIN PROLOGUE  DY4 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPX4 */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DY4-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This program computes second order finite difference */
/*     approximations to the third and fourth Y */
/*     partial derivatives of U at the (I,J) mesh point. */

/* ***SEE ALSO  SEPX4 */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPL4 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DY4 */

/* ***FIRST EXECUTABLE STATEMENT  DY4 */
    /* Parameter adjustments */
    u_dim1 = *idmn;
    u_offset = 1 + u_dim1;
    u -= u_offset;

    /* Function Body */
    if (*j > 2 && *j < spl4_1.l - 1) {
	goto L50;
    }
    if (*j == 1) {
	goto L10;
    }
    if (*j == 2) {
	goto L30;
    }
    if (*j == spl4_1.l - 1) {
	goto L60;
    }
    if (*j == spl4_1.l) {
	goto L80;
    }

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C */

L10:
    if (spl4_1.kswy == 1) {
	goto L20;
    }
    *uyyy = (u[*i__ + u_dim1] * -5.f + u[*i__ + (u_dim1 << 1)] * 18.f - u[*
	    i__ + u_dim1 * 3] * 24.f + u[*i__ + (u_dim1 << 2)] * 14.f - u[*
	    i__ + u_dim1 * 5] * 3.f) / spl4_1.tdly3;
    *uyyyy = (u[*i__ + u_dim1] * 3.f - u[*i__ + (u_dim1 << 1)] * 14.f + u[*
	    i__ + u_dim1 * 3] * 26.f - u[*i__ + (u_dim1 << 2)] * 24.f + u[*
	    i__ + u_dim1 * 5] * 11.f - u[*i__ + u_dim1 * 6] * 2.f) / 
	    spl4_1.dly4;
    return 0;

/*     PERIODIC AT X=A */

L20:
    *uyyy = (-u[*i__ + (spl4_1.l - 2) * u_dim1] + u[*i__ + (spl4_1.l - 1) * 
	    u_dim1] * 2.f - u[*i__ + (u_dim1 << 1)] * 2.f + u[*i__ + u_dim1 * 
	    3]) / spl4_1.tdly3;
    *uyyyy = (u[*i__ + (spl4_1.l - 2) * u_dim1] - u[*i__ + (spl4_1.l - 1) * 
	    u_dim1] * 4.f + u[*i__ + u_dim1] * 6.f - u[*i__ + (u_dim1 << 1)] *
	     4.f + u[*i__ + u_dim1 * 3]) / spl4_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY */

L30:
    if (spl4_1.kswy == 1) {
	goto L40;
    }
    *uyyy = (u[*i__ + u_dim1] * -3.f + u[*i__ + (u_dim1 << 1)] * 10.f - u[*
	    i__ + u_dim1 * 3] * 12.f + u[*i__ + (u_dim1 << 2)] * 6.f - u[*i__ 
	    + u_dim1 * 5]) / spl4_1.tdly3;
    *uyyyy = (u[*i__ + u_dim1] * 2.f - u[*i__ + (u_dim1 << 1)] * 9.f + u[*i__ 
	    + u_dim1 * 3] * 16.f - u[*i__ + (u_dim1 << 2)] * 14.f + u[*i__ + 
	    u_dim1 * 5] * 6.f - u[*i__ + u_dim1 * 6]) / spl4_1.dly4;
    return 0;

/*     PERIODIC AT Y=C+DLY */

L40:
    *uyyy = (-u[*i__ + (spl4_1.l - 1) * u_dim1] + u[*i__ + u_dim1] * 2.f - u[*
	    i__ + u_dim1 * 3] * 2.f + u[*i__ + (u_dim1 << 2)]) / spl4_1.tdly3;
    *uyyyy = (u[*i__ + (spl4_1.l - 1) * u_dim1] - u[*i__ + u_dim1] * 4.f + u[*
	    i__ + (u_dim1 << 1)] * 6.f - u[*i__ + u_dim1 * 3] * 4.f + u[*i__ 
	    + (u_dim1 << 2)]) / spl4_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR */

L50:
    *uyyy = (-u[*i__ + (*j - 2) * u_dim1] + u[*i__ + (*j - 1) * u_dim1] * 2.f 
	    - u[*i__ + (*j + 1) * u_dim1] * 2.f + u[*i__ + (*j + 2) * u_dim1])
	     / spl4_1.tdly3;
    *uyyyy = (u[*i__ + (*j - 2) * u_dim1] - u[*i__ + (*j - 1) * u_dim1] * 4.f 
	    + u[*i__ + *j * u_dim1] * 6.f - u[*i__ + (*j + 1) * u_dim1] * 4.f 
	    + u[*i__ + (*j + 2) * u_dim1]) / spl4_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY */

L60:
    if (spl4_1.kswy == 1) {
	goto L70;
    }
    *uyyy = (u[*i__ + (spl4_1.l - 4) * u_dim1] - u[*i__ + (spl4_1.l - 3) * 
	    u_dim1] * 6.f + u[*i__ + (spl4_1.l - 2) * u_dim1] * 12.f - u[*i__ 
	    + (spl4_1.l - 1) * u_dim1] * 10.f + u[*i__ + spl4_1.l * u_dim1] * 
	    3.f) / spl4_1.tdly3;
    *uyyyy = (-u[*i__ + (spl4_1.l - 5) * u_dim1] + u[*i__ + (spl4_1.l - 4) * 
	    u_dim1] * 6.f - u[*i__ + (spl4_1.l - 3) * u_dim1] * 14.f + u[*i__ 
	    + (spl4_1.l - 2) * u_dim1] * 16.f - u[*i__ + (spl4_1.l - 1) * 
	    u_dim1] * 9.f + u[*i__ + spl4_1.l * u_dim1] * 2.f) / spl4_1.dly4;
    return 0;

/*     PERIODIC AT Y=D-DLY */

L70:
    *uyyy = (-u[*i__ + (spl4_1.l - 3) * u_dim1] + u[*i__ + (spl4_1.l - 2) * 
	    u_dim1] * 2.f - u[*i__ + u_dim1] * 2.f + u[*i__ + (u_dim1 << 1)]) 
	    / spl4_1.tdly3;
    *uyyyy = (u[*i__ + (spl4_1.l - 3) * u_dim1] - u[*i__ + (spl4_1.l - 2) * 
	    u_dim1] * 4.f + u[*i__ + (spl4_1.l - 1) * u_dim1] * 6.f - u[*i__ 
	    + u_dim1] * 4.f + u[*i__ + (u_dim1 << 1)]) / spl4_1.dly4;
    return 0;

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D */

L80:
    *uyyy = -(u[*i__ + (spl4_1.l - 4) * u_dim1] * 3.f - u[*i__ + (spl4_1.l - 
	    3) * u_dim1] * 14.f + u[*i__ + (spl4_1.l - 2) * u_dim1] * 24.f - 
	    u[*i__ + (spl4_1.l - 1) * u_dim1] * 18.f + u[*i__ + spl4_1.l * 
	    u_dim1] * 5.f) / spl4_1.tdly3;
    *uyyyy = (u[*i__ + (spl4_1.l - 5) * u_dim1] * -2.f + u[*i__ + (spl4_1.l - 
	    4) * u_dim1] * 11.f - u[*i__ + (spl4_1.l - 3) * u_dim1] * 24.f + 
	    u[*i__ + (spl4_1.l - 2) * u_dim1] * 26.f - u[*i__ + (spl4_1.l - 1)
	     * u_dim1] * 14.f + u[*i__ + spl4_1.l * u_dim1] * 3.f) / 
	    spl4_1.dly4;
    return 0;
} /* dy4_ */


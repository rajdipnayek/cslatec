/* passf5.f -- translated by f2c (version 12.02.01).
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

/* DECK PASSF5 */
/* Subroutine */ int passf5_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2, real *wa3, real *wa4)
{
    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static real pi, ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, 
	    cr4, ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5, 
	    ti11, ti12, tr11, tr12;

/* ***BEGIN PROLOGUE  PASSF5 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Calculate the fast Fourier transform of subvectors of */
/*            length five. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***TYPE      SINGLE PRECISION (PASSF5-S) */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           (a) changing dummy array size declarations (1) to (*), */
/*           (b) changing definition of variables PI, TI11, TI12, */
/*               TR11, TR12 by using FORTRAN intrinsic functions ATAN */
/*               and SIN instead of DATA statements. */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PASSF5 */
/* ***FIRST EXECUTABLE STATEMENT  PASSF5 */
    /* Parameter adjustments */
    cc_dim1 = *ido;
    cc_offset = 1 + cc_dim1 * 6;
    cc -= cc_offset;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
    ch -= ch_offset;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    /* Function Body */
    pi = atan(1.f) * 4.f;
    tr11 = sin(pi * .1f);
    ti11 = -sin(pi * .4f);
    tr12 = -sin(pi * .3f);
    ti12 = -sin(pi * .2f);
    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 2) * cc_dim1 + 2] - cc[(k * 5 + 5) * cc_dim1 + 2];
	ti2 = cc[(k * 5 + 2) * cc_dim1 + 2] + cc[(k * 5 + 5) * cc_dim1 + 2];
	ti4 = cc[(k * 5 + 3) * cc_dim1 + 2] - cc[(k * 5 + 4) * cc_dim1 + 2];
	ti3 = cc[(k * 5 + 3) * cc_dim1 + 2] + cc[(k * 5 + 4) * cc_dim1 + 2];
	tr5 = cc[(k * 5 + 2) * cc_dim1 + 1] - cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[(k * 5 + 2) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr4 = cc[(k * 5 + 3) * cc_dim1 + 1] - cc[(k * 5 + 4) * cc_dim1 + 1];
	tr3 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 4) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2 
		+ tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 5 + 1) * cc_dim1 + 2] + ti2 
		+ ti3;
	cr2 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2 + tr12 * tr3;
	ci2 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr11 * ti2 + tr12 * ti3;
	cr3 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2 + tr11 * tr3;
	ci3 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr12 * ti2 + tr11 * ti3;
	cr5 = ti11 * tr5 + ti12 * tr4;
	ci5 = ti11 * ti5 + ti12 * ti4;
	cr4 = ti12 * tr5 - ti11 * tr4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci3 + cr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ci3 - cr4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 2] = ci2 - cr5;
/* L101: */
    }
    return 0;
L102:
    if (*ido / 2 < *l1) {
	goto L105;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
/* DIR$ IVDEP */
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ti5 = cc[i__ + (k * 5 + 2) * cc_dim1] - cc[i__ + (k * 5 + 5) * 
		    cc_dim1];
	    ti2 = cc[i__ + (k * 5 + 2) * cc_dim1] + cc[i__ + (k * 5 + 5) * 
		    cc_dim1];
	    ti4 = cc[i__ + (k * 5 + 3) * cc_dim1] - cc[i__ + (k * 5 + 4) * 
		    cc_dim1];
	    ti3 = cc[i__ + (k * 5 + 3) * cc_dim1] + cc[i__ + (k * 5 + 4) * 
		    cc_dim1];
	    tr5 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] - cc[i__ - 1 + (k * 5 + 
		    5) * cc_dim1];
	    tr2 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] + cc[i__ - 1 + (k * 5 + 
		    5) * cc_dim1];
	    tr4 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] - cc[i__ - 1 + (k * 5 + 
		    4) * cc_dim1];
	    tr3 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] + cc[i__ - 1 + (k * 5 + 
		    4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 5 + 1) *
		     cc_dim1] + tr2 + tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 5 + 1) * 
		    cc_dim1] + ti2 + ti3;
	    cr2 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * 
		    tr3;
	    ci2 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * 
		    tr3;
	    ci3 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * dr2 
		    + wa1[i__] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * di2 - 
		    wa1[i__] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * dr3 + 
		    wa2[i__] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * di3 - wa2[
		    i__] * dr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * dr4 
		    + wa3[i__] * di4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * di4 - 
		    wa3[i__] * dr4;
	    ch[i__ - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * dr5 + 
		    wa4[i__] * di5;
	    ch[i__ + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * di5 - wa4[
		    i__] * dr5;
/* L103: */
	}
/* L104: */
    }
    return 0;
L105:
    i__1 = *ido;
    for (i__ = 2; i__ <= i__1; i__ += 2) {
/* DIR$ IVDEP */
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ti5 = cc[i__ + (k * 5 + 2) * cc_dim1] - cc[i__ + (k * 5 + 5) * 
		    cc_dim1];
	    ti2 = cc[i__ + (k * 5 + 2) * cc_dim1] + cc[i__ + (k * 5 + 5) * 
		    cc_dim1];
	    ti4 = cc[i__ + (k * 5 + 3) * cc_dim1] - cc[i__ + (k * 5 + 4) * 
		    cc_dim1];
	    ti3 = cc[i__ + (k * 5 + 3) * cc_dim1] + cc[i__ + (k * 5 + 4) * 
		    cc_dim1];
	    tr5 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] - cc[i__ - 1 + (k * 5 + 
		    5) * cc_dim1];
	    tr2 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] + cc[i__ - 1 + (k * 5 + 
		    5) * cc_dim1];
	    tr4 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] - cc[i__ - 1 + (k * 5 + 
		    4) * cc_dim1];
	    tr3 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] + cc[i__ - 1 + (k * 5 + 
		    4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 5 + 1) *
		     cc_dim1] + tr2 + tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 5 + 1) * 
		    cc_dim1] + ti2 + ti3;
	    cr2 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * 
		    tr3;
	    ci2 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * 
		    tr3;
	    ci3 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * dr2 
		    + wa1[i__] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * di2 - 
		    wa1[i__] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * dr3 + 
		    wa2[i__] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * di3 - wa2[
		    i__] * dr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * dr4 
		    + wa3[i__] * di4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * di4 - 
		    wa3[i__] * dr4;
	    ch[i__ - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * dr5 + 
		    wa4[i__] * di5;
	    ch[i__ + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * di5 - wa4[
		    i__] * dr5;
/* L106: */
	}
/* L107: */
    }
    return 0;
} /* passf5_ */


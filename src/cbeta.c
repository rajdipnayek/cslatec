/* cbeta.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK CBETA */
/* Complex */ void cbeta_(complex * ret_val, complex *a, complex *b)
{
    /* Initialized data */

    static real xmax = 0.f;

    /* System generated locals */
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Local variables */
    static real xmin, xmaxt;
    extern /* Complex */ void cgamma_(complex *, complex *), clbeta_(complex *
	    , complex *, complex *);
    extern /* Subroutine */ int gamlim_(real *, real *), xermsg_(char *, char 
	    *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CBETA */
/* ***PURPOSE  Compute the complete Beta function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7B */
/* ***TYPE      COMPLEX (BETA-S, DBETA-D, CBETA-C) */
/* ***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CBETA computes the complete beta function of complex parameters A */
/* and B. */
/* Input Parameters: */
/*       A   complex and the real part of A positive */
/*       B   complex and the real part of B positive */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CGAMMA, CLBETA, GAMLIM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  CBETA */
/* ***FIRST EXECUTABLE STATEMENT  CBETA */
    if (xmax == 0.f) {
	gamlim_(&xmin, &xmaxt);
	xmax = xmaxt;
    }

    if (a->r <= 0.f || b->r <= 0.f) {
	xermsg_("SLATEC", "CBETA", "REAL PART OF BOTH ARGUMENTS MUST BE GT 0",
		 &c__1, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)40);
    }

    if (a->r + b->r < xmax) {
	cgamma_(&q__2, a);
	cgamma_(&q__4, b);
	q__6.r = a->r + b->r, q__6.i = a->i + b->i;
	cgamma_(&q__5, &q__6);
	c_div(&q__3, &q__4, &q__5);
	q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i 
		+ q__2.i * q__3.r;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (a->r + b->r < xmax) {
	return ;
    }

    clbeta_(&q__2, a, b);
    c_exp(&q__1, &q__2);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* cbeta_ */


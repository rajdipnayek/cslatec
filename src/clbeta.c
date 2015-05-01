/* clbeta.f -- translated by f2c (version 12.02.01).
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

/* DECK CLBETA */
/* Complex */ void clbeta_(complex * ret_val, complex *a, complex *b)
{
    /* System generated locals */
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Local variables */
    extern /* Complex */ void clngam_(complex *, complex *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CLBETA */
/* ***PURPOSE  Compute the natural logarithm of the complete Beta */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7B */
/* ***TYPE      COMPLEX (ALBETA-S, DLBETA-D, CLBETA-C) */
/* ***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CLBETA computes the natural log of the complex valued complete beta */
/* function of complex parameters A and B.  This is a preliminary version */
/* which is not accurate. */

/* Input Parameters: */
/*       A   complex and the real part of A positive */
/*       B   complex and the real part of B positive */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CLNGAM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  CLBETA */
/* ***FIRST EXECUTABLE STATEMENT  CLBETA */
    if (a->r <= 0.f || b->r <= 0.f) {
	xermsg_("SLATEC", "CLBETA", "REAL PART OF BOTH ARGUMENTS MUST BE GT 0"
		, &c__1, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    }

    clngam_(&q__3, a);
    clngam_(&q__4, b);
    q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
    q__6.r = a->r + b->r, q__6.i = a->i + b->i;
    clngam_(&q__5, &q__6);
    q__1.r = q__2.r - q__5.r, q__1.i = q__2.i - q__5.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* clbeta_ */


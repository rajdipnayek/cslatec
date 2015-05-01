/* drotg.f -- translated by f2c (version 12.02.01).
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

/* DECK DROTG */
/* Subroutine */ int drotg_(doublereal *da, doublereal *db, doublereal *dc, 
	doublereal *ds)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal r__, u, v;

/* ***BEGIN PROLOGUE  DROTG */
/* ***PURPOSE  Construct a plane Givens rotation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B10 */
/* ***TYPE      DOUBLE PRECISION (SROTG-S, DROTG-D, CROTG-C) */
/* ***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION, */
/*             LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*       DA  double precision scalar */
/*       DB  double precision scalar */

/*     --Output-- */
/*       DA  double precision result R */
/*       DB  double precision result Z */
/*       DC  double precision result */
/*       DS  double precision result */

/*     Construct the Givens transformation */

/*         ( DC  DS ) */
/*     G = (        ) ,    DC**2 + DS**2 = 1 , */
/*         (-DS  DC ) */

/*     which zeros the second entry of the 2-vector  (DA,DB)**T . */

/*     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in */
/*     storage.  The value of DB is overwritten by a value Z which */
/*     allows DC and DS to be recovered by the following algorithm. */

/*           If Z=1  set  DC=0.0  and  DS=1.0 */
/*           If ABS(Z) .LT. 1  set  DC=SQRT(1-Z**2)  and  DS=Z */
/*           If ABS(Z) .GT. 1  set  DC=1/Z  and  DS=SQRT(1-DC**2) */

/*     Normally, the subprogram DROT(N,DX,INCX,DY,INCY,DC,DS) will */
/*     next be called to apply the transformation to a 2 by N matrix. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DROTG */
/* ***FIRST EXECUTABLE STATEMENT  DROTG */
    if (abs(*da) <= abs(*db)) {
	goto L10;
    }

/* *** HERE ABS(DA) .GT. ABS(DB) *** */

    u = *da + *da;
    v = *db / u;

/*     NOTE THAT U AND R HAVE THE SIGN OF DA */

/* Computing 2nd power */
    d__1 = v;
    r__ = sqrt(d__1 * d__1 + .25) * u;

/*     NOTE THAT DC IS POSITIVE */

    *dc = *da / r__;
    *ds = v * (*dc + *dc);
    *db = *ds;
    *da = r__;
    return 0;

/* *** HERE ABS(DA) .LE. ABS(DB) *** */

L10:
    if (*db == 0.) {
	goto L20;
    }
    u = *db + *db;
    v = *da / u;

/*     NOTE THAT U AND R HAVE THE SIGN OF DB */
/*     (R IS IMMEDIATELY STORED IN DA) */

/* Computing 2nd power */
    d__1 = v;
    *da = sqrt(d__1 * d__1 + .25) * u;

/*     NOTE THAT DS IS POSITIVE */

    *ds = *db / *da;
    *dc = v * (*ds + *ds);
    if (*dc == 0.) {
	goto L15;
    }
    *db = 1. / *dc;
    return 0;
L15:
    *db = 1.;
    return 0;

/* *** HERE DA = DB = 0.0 *** */

L20:
    *dc = 1.;
    *ds = 0.;
    return 0;

} /* drotg_ */


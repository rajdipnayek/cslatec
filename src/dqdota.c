/* dqdota.f -- translated by f2c (version 12.02.01).
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
    integer mpb, mpt, mpm, mplun, mpmxr, mpr[30];
} mpcom_;

#define mpcom_1 mpcom_

/* DECK DQDOTA */
doublereal dqdota_(integer *n, doublereal *db, integer *qc, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{
    /* Initialized data */

    static integer i1 = 0;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, ix, iy, qx[30], qy[30];
    extern /* Subroutine */ int mpadd_(integer *, integer *, integer *), 
	    mpcdm_(doublereal *, integer *), mpcmd_(integer *, doublereal *), 
	    mpmul_(integer *, integer *, integer *), mpblas_(integer *);

/* ***BEGIN PROLOGUE  DQDOTA */
/* ***PURPOSE  Compute the inner product of two vectors with extended */
/*            precision accumulation and result. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D1A4 */
/* ***TYPE      DOUBLE PRECISION (DQDOTA-D) */
/* ***KEYWORDS  DOT PRODUCT, INNER PRODUCT */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*               B L A S  Subprogram */
/*   Description of Parameters */

/*    --Input-- */
/*       N  number of elements in input vector(S) */
/*      DB  double precision scalar to be added to inner product */
/*      QC  extended precision scalar to be added to inner product */
/*      DX  double precision vector with N elements */
/*    INCX  storage spacing between elements of DX */
/*      DY  double precision vector with N elements */
/*    INCY  storage spacing between elements of DY */

/*    --Output-- */
/*  DQDOTA  double precision result */
/*      QC  extended precision result */

/*    D.P. dot product with extended precision accumulation (and result) */
/*    QC and DQDOTA are set = DB + QC + sum for I = 0 to N-1 of */
/*      DX(LX+I*INCX) * DY(LY+I*INCY),  where QC is an extended */
/*      precision result previously computed by DQDOTI or DQDOTA */
/*      and LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is */
/*      defined in a similar way using INCY.  The MP package by */
/*      Richard P. Brent is used for the extended precision arithmetic. */

/*    Fred T. Krogh,  JPL,  1977,  June 1 */

/*    The common block for the MP package is name MPCOM.  If local */
/*    variable I1 is zero, DQDOTA calls MPBLAS to initialize */
/*    the MP package and reset I1 to 1. */

/*    The argument QC(*) and the local variables QX and QY are INTEGER */
/*    arrays of size 30.  See the comments in the routine MPBLAS for the */
/*    reason for this choice. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  MPADD, MPBLAS, MPCDM, MPCMD, MPMUL */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   930124  Increased Array sizes for SUN -r8.  (RWC) */
/* ***END PROLOGUE  DQDOTA */
    /* Parameter adjustments */
    --dy;
    --dx;
    --qc;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DQDOTA */
    if (i1 == 0) {
	mpblas_(&i1);
    }
    if (*db == 0.) {
	goto L20;
    }
    mpcdm_(db, qx);
    mpadd_(&qc[1], qx, &qc[1]);
L20:
    if (*n == 0) {
	goto L40;
    }
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mpcdm_(&dx[ix], qx);
	mpcdm_(&dy[iy], qy);
	mpmul_(qx, qy, qx);
	mpadd_(&qc[1], qx, &qc[1]);
	ix += *incx;
	iy += *incy;
/* L30: */
    }
L40:
    mpcmd_(&qc[1], &ret_val);
    return ret_val;
} /* dqdota_ */


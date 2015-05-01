/* sods.f -- translated by f2c (version 12.02.01).
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

/* DECK SODS */
/* Subroutine */ int sods_(real *a, real *x, real *b, integer *neq, integer *
	nuk, integer *nrda, integer *iflag, real *work, integer *iwork)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer kc, kd, ip, is, ks, kt, kv, kz, iter;
    extern /* Subroutine */ int lssods_(real *, real *, real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *);

/* ***BEGIN PROLOGUE  SODS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SODS-S) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     SODS solves the overdetermined system of linear equations A X = B, */
/*     where A is NEQ by NUK and NEQ .GE. NUK. If rank A = NUK, */
/*     X is the UNIQUE least squares solution vector. That is, */
/*              R(1)**2 + ..... + R(NEQ)**2 = minimum */
/*     where R is the residual vector  R = B - A X. */
/*     If rank A .LT. NUK , the least squares solution of minimal */
/*     length can be provided. */
/*     SODS is an interfacing routine which calls subroutine LSSODS */
/*     for the solution. LSSODS in turn calls subroutine ORTHOL and */
/*     possibly subroutine OHTROR for the decomposition of A by */
/*     orthogonal transformations. In the process, ORTHOL calls upon */
/*     subroutine CSCALE for scaling. */

/* ********************************************************************** */
/*   Input */
/* ********************************************************************** */

/*     A -- Contains the matrix of NEQ equations in NUK unknowns and must */
/*          be dimensioned NRDA by NUK. The original A is destroyed */
/*     X -- Solution array of length at least NUK */
/*     B -- Given constant vector of length NEQ, B is destroyed */
/*     NEQ -- Number of equations, NEQ greater or equal to 1 */
/*     NUK -- Number of columns in the matrix (which is also the number */
/*            of unknowns), NUK not larger than NEQ */
/*     NRDA -- Row dimension of A, NRDA greater or equal to NEQ */
/*     IFLAG -- Status indicator */
/*            =0 For the first call (and for each new problem defined by */
/*               a new matrix A) when the matrix data is treated as exact */
/*           =-K For the first call (and for each new problem defined by */
/*               a new matrix A) when the matrix data is assumed to be */
/*               accurate to about K digits */
/*            =1 For subsequent calls whenever the matrix A has already */
/*               been decomposed (problems with new vectors B but */
/*               same matrix a can be handled efficiently) */
/*     WORK(*),IWORK(*) -- Arrays for storage of internal information, */
/*                     WORK must be dimensioned at least  2 + 5*NUK */
/*                     IWORK must be dimensioned at least NUK+2 */
/*     IWORK(2) -- Scaling indicator */
/*                 =-1 If the matrix A is to be pre-scaled by */
/*                 columns when appropriate */
/*                 If the scaling indicator is not equal to -1 */
/*                 no scaling will be attempted */
/*              For most problems scaling will probably not be necessary */

/* ********************************************************************** */
/*   OUTPUT */
/* ********************************************************************** */

/*     IFLAG -- Status indicator */
/*            =1 If solution was obtained */
/*            =2 If improper input is detected */
/*            =3 If rank of matrix is less than NUK */
/*               If the minimal length least squares solution is */
/*               desired, simply reset IFLAG=1 and call the code again */
/*     X -- Least squares solution of  A X = B */
/*     A -- Contains the strictly upper triangular part of the reduced */
/*           matrix and the transformation information */
/*     WORK(*),IWORK(*) -- Contains information needed on subsequent */
/*                         Calls (IFLAG=1 case on input) which must not */
/*                         be altered */
/*                         WORK(1) contains the Euclidean norm of */
/*                         the residual vector */
/*                         WORK(2) contains the Euclidean norm of */
/*                         the solution vector */
/*                         IWORK(1) contains the numerically determined */
/*                         rank of the matrix A */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***REFERENCES  G. Golub, Numerical methods for solving linear least */
/*                 squares problems, Numerische Mathematik 7, (1965), */
/*                 pp. 206-216. */
/*               P. Businger and G. Golub, Linear least squares */
/*                 solutions by Householder transformations, Numerische */
/*                 Mathematik  7, (1965), pp. 269-276. */
/*               H. A. Watts, Solving linear least squares problems */
/*                 using SODS/SUDS/CODS, Sandia Report SAND77-0683, */
/*                 Sandia Laboratories, 1977. */
/* ***ROUTINES CALLED  LSSODS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SODS */

/* ***FIRST EXECUTABLE STATEMENT  SODS */
    /* Parameter adjustments */
    --x;
    --b;
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --iwork;

    /* Function Body */
    iter = 0;
    is = 2;
    ip = 3;
    ks = 2;
    kd = 3;
    kz = kd + *nuk;
    kv = kz + *nuk;
    kt = kv + *nuk;
    kc = kt + *nuk;

    lssods_(&a[a_offset], &x[1], &b[1], neq, nuk, nrda, iflag, &iwork[1], &
	    iwork[is], &a[a_offset], &work[kd], &iwork[ip], &iter, &work[1], &
	    work[ks], &work[kz], &b[1], &work[kv], &work[kt], &work[kc]);

    return 0;
} /* sods_ */


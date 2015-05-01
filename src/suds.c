/* suds.f -- translated by f2c (version 12.02.01).
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

/* DECK SUDS */
/* Subroutine */ int suds_(real *a, real *x, real *b, integer *neq, integer *
	nuk, integer *nrda, integer *iflag, integer *mlso, real *work, 
	integer *iwork)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer il, ip, is, ks, kt, ku, kv;
    extern /* Subroutine */ int lssuds_(real *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, integer *, real *, real *, 
	    real *, integer *, real *);

/* ***BEGIN PROLOGUE  SUDS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SUDS-S, DSUDS-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     SUDS solves the underdetermined system of linear equations A Z = B */
/*     where A is NEQ by NUK and NEQ .LE. NUK. In particular, if rank A */
/*     equals IRA, a vector X and a matrix U are determined such that */
/*     X is the UNIQUE solution of smallest length, satisfying A X = B, */
/*     and the columns of U form an orthonormal basis for the null */
/*     space of A, satisfying A U = 0 . Then all solutions Z are */
/*     given by */
/*              Z = X + C(1)*U(1) + ..... + C(NUK-IRA)*U(NUK-IRA) */
/*     where U(J) represents the J-th column of U and the C(J) are */
/*     arbitrary constants. */
/*     If the system of equations are not compatible, only the least */
/*     squares solution of minimal length is computed. */
/*     SUDS is an interfacing routine which calls subroutine LSSUDS */
/*     for the solution. LSSUDS in turn calls subroutine ORTHOR and */
/*     possibly subroutine OHTROL for the decomposition of A by */
/*     orthogonal transformations. In the process, ORTHOR calls upon */
/*     subroutine CSCALE for scaling. */

/* ********************************************************************** */
/*   INPUT */
/* ********************************************************************** */

/*     A -- Contains the matrix of NEQ equations in NUK unknowns and must */
/*          be dimensioned NRDA by NUK. The original A is destroyed. */
/*     X -- Solution array of length at least NUK */
/*     B -- Given constant vector of length NEQ, B is destroyed */
/*     NEQ -- Number of equations, NEQ greater or equal to 1 */
/*     NUK -- Number of columns in the matrix (which is also the number */
/*            of unknowns), NUK not smaller than NEQ */
/*     NRDA -- Row dimension of A, NRDA greater or equal to NEQ */
/*     IFLAG -- Status indicator */
/*           =0  For the first call (and for each new problem defined by */
/*               a new matrix A) when the matrix data is treated as exact */
/*           =-K For the first call (and for each new problem defined by */
/*               a new matrix A) when the matrix data is assumed to be */
/*               accurate to about K digits */
/*           =1  For subsequent calls whenever the matrix A has already */
/*               been decomposed (problems with new vectors B but */
/*               same matrix A can be handled efficiently) */
/*     MLSO -- =0 If only the minimal length solution is wanted */
/*             =1 If the complete solution is wanted, includes the */
/*                linear space defined by the matrix U in the abstract */
/*     WORK(*),IWORK(*) -- Arrays for storage of internal information, */
/*                WORK must be dimensioned at least */
/*                       NUK + 3*NEQ + MLSO*NUK*(NUK-rank A) */
/*                where it is possible for   0 .LE. rank A .LE. NEQ */
/*                IWORK must be dimensioned at least   3 + NEQ */
/*     IWORK(2) -- Scaling indicator */
/*                 =-1 If the matrix is to be pre-scaled by */
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
/*            =3 If rank of matrix is less than NEQ */
/*               To continue simply reset IFLAG=1 and call SUDS again */
/*            =4 If the system of equations appears to be inconsistent. */
/*               However, the least squares solution of minimal length */
/*               was obtained. */
/*     X -- Minimal length least squares solution of  A X = B */
/*     A -- Contains the strictly upper triangular part of the reduced */
/*           matrix and transformation information */
/*     WORK(*),IWORK(*) -- Contains information needed on subsequent */
/*                         calls (IFLAG=1 case on input) which must not */
/*                         be altered. */
/*                         The matrix U described in the abstract is */
/*                         stored in the  NUK*(NUK-rank A) elements of */
/*                         the work array beginning at WORK(1+NUK+3*NEQ). */
/*                         However U is not defined when MLSO=0 or */
/*                         IFLAG=4. */
/*                         IWORK(1) Contains the numerically determined */
/*                         rank of the matrix A */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***REFERENCES  H. A. Watts, Solving linear least squares problems */
/*                 using SODS/SUDS/CODS, Sandia Report SAND77-0683, */
/*                 Sandia Laboratories, 1977. */
/* ***ROUTINES CALLED  LSSUDS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SUDS */

/* ***FIRST EXECUTABLE STATEMENT  SUDS */
    /* Parameter adjustments */
    --x;
    --b;
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --iwork;

    /* Function Body */
    is = 2;
    ip = 3;
    il = ip + *neq;
    kv = *neq + 1;
    kt = kv + *neq;
    ks = kt + *neq;
    ku = ks + *nuk;

    lssuds_(&a[a_offset], &x[1], &b[1], neq, nuk, nrda, &work[ku], nuk, iflag,
	     mlso, &iwork[1], &iwork[is], &a[a_offset], &work[1], &iwork[ip], 
	    &b[1], &work[kv], &work[kt], &iwork[il], &work[ks]);

    return 0;
} /* suds_ */


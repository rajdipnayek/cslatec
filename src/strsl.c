/* strsl.f -- translated by f2c (version 12.02.01).
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

/* DECK STRSL */
/* Subroutine */ int strsl_(real *t, integer *ldt, integer *n, real *b, 
	integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    static integer j, jj, case__;
    static real temp;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  STRSL */
/* ***PURPOSE  Solve a system of the form  T*X=B or TRANS(T)*X=B, where */
/*            T is a triangular matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2A3 */
/* ***TYPE      SINGLE PRECISION (STRSL-S, DTRSL-D, CTRSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, TRIANGULAR LINEAR SYSTEM, */
/*             TRIANGULAR MATRIX */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     STRSL solves systems of the form */

/*                   T * X = B */
/*     or */
/*                   TRANS(T) * X = B */

/*     where T is a triangular matrix of order N.  Here TRANS(T) */
/*     denotes the transpose of the matrix T. */

/*     On Entry */

/*         T         REAL(LDT,N) */
/*                   T contains the matrix of the system.  The zero */
/*                   elements of the matrix are not referenced, and */
/*                   the corresponding elements of the array can be */
/*                   used to store other information. */

/*         LDT       INTEGER */
/*                   LDT is the leading dimension of the array T. */

/*         N         INTEGER */
/*                   N is the order of the system. */

/*         B         REAL(N). */
/*                   B contains the right hand side of the system. */

/*         JOB       INTEGER */
/*                   JOB specifies what kind of system is to be solved. */
/*                   If JOB is */

/*                        00   solve T*X=B, T lower triangular, */
/*                        01   solve T*X=B, T upper triangular, */
/*                        10   solve TRANS(T)*X=B, T lower triangular, */
/*                        11   solve TRANS(T)*X=B, T upper triangular. */

/*     On Return */

/*         B         B contains the solution, if INFO .EQ. 0. */
/*                   Otherwise B is unaltered. */

/*         INFO      INTEGER */
/*                   INFO contains zero if the system is nonsingular. */
/*                   Otherwise INFO contains the index of */
/*                   the first zero diagonal element of T. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SAXPY, SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  STRSL */


/* ***FIRST EXECUTABLE STATEMENT  STRSL */

/*        CHECK FOR ZERO DIAGONAL ELEMENTS. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
	if (t[*info + *info * t_dim1] == 0.f) {
	    goto L150;
	}
/* L10: */
    }
    *info = 0;

/*        DETERMINE THE TASK AND GO TO IT. */

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch (case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

/*        SOLVE T*X=B FOR T LOWER TRIANGULAR */

L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	temp = -b[j - 1];
	i__2 = *n - j + 1;
	saxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L30: */
    }
L40:
    goto L140;

/*        SOLVE T*X=B FOR T UPPER TRIANGULAR. */

L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	temp = -b[j + 1];
	saxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L60: */
    }
L70:
    goto L140;

/*        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR. */

L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = jj - 1;
	b[j] -= sdot_(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L90: */
    }
L100:
    goto L140;

/*        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR. */

L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	b[j] -= sdot_(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* strsl_ */


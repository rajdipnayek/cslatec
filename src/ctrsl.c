/* ctrsl.f -- translated by f2c (version 12.02.01).
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

/* DECK CTRSL */
/* Subroutine */ int ctrsl_(complex *t, integer *ldt, integer *n, complex *b, 
	integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2;

    /* Local variables */
    static integer j, jj, case__;
    static complex temp;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CTRSL */
/* ***PURPOSE  Solve a system of the form  T*X=B or CTRANS(T)*X=B, where */
/*            T is a triangular matrix.  Here CTRANS(T) is the conjugate */
/*            transpose. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C3 */
/* ***TYPE      COMPLEX (STRSL-S, DTRSL-D, CTRSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, TRIANGULAR LINEAR SYSTEM, */
/*             TRIANGULAR MATRIX */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     CTRSL solves systems of the form */

/*                   T * X = B */
/*     or */
/*                   CTRANS(T) * X = B */

/*     where T is a triangular matrix of order N.  Here CTRANS(T) */
/*     denotes the conjugate transpose of the matrix T. */

/*     On Entry */

/*         T         COMPLEX(LDT,N) */
/*                   T contains the matrix of the system.  The zero */
/*                   elements of the matrix are not referenced, and */
/*                   the corresponding elements of the array can be */
/*                   used to store other information. */

/*         LDT       INTEGER */
/*                   LDT is the leading dimension of the array T. */

/*         N         INTEGER */
/*                   N is the order of the system. */

/*         B         COMPLEX(N). */
/*                   B contains the right hand side of the system. */

/*         JOB       INTEGER */
/*                   JOB specifies what kind of system is to be solved. */
/*                   If JOB is */

/*                        00   solve T*X = B, T lower triangular, */
/*                        01   solve T*X = B, T upper triangular, */
/*                        10   solve CTRANS(T)*X = B, T lower triangular, */
/*                        11   solve CTRANS(T)*X = B, T upper triangular. */

/*     On Return */

/*         B         B contains the solution, if INFO .EQ. 0. */
/*                   Otherwise B is unaltered. */

/*         INFO      INTEGER */
/*                   INFO contains zero if the system is nonsingular. */
/*                   Otherwise INFO contains the index of */
/*                   the first zero diagonal element of T. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CTRSL */


/* ***FIRST EXECUTABLE STATEMENT  CTRSL */

/*        CHECK FOR ZERO DIAGONAL ELEMENTS. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
	i__2 = *info + *info * t_dim1;
	if ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[*info + *info *
		 t_dim1]), dabs(r__2)) == 0.f) {
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
    c_div(&q__1, &b[1], &t[t_dim1 + 1]);
    b[1].r = q__1.r, b[1].i = q__1.i;
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	q__1.r = -b[i__2].r, q__1.i = -b[i__2].i;
	temp.r = q__1.r, temp.i = q__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	i__2 = j;
	c_div(&q__1, &b[j], &t[j + j * t_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L30: */
    }
L40:
    goto L140;

/*        SOLVE T*X=B FOR T UPPER TRIANGULAR. */

L50:
    i__1 = *n;
    c_div(&q__1, &b[*n], &t[*n + *n * t_dim1]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = j + 1;
	q__1.r = -b[i__2].r, q__1.i = -b[i__2].i;
	temp.r = q__1.r, temp.i = q__1.i;
	caxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	i__2 = j;
	c_div(&q__1, &b[j], &t[j + j * t_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L60: */
    }
L70:
    goto L140;

/*        SOLVE CTRANS(T)*X=B FOR T LOWER TRIANGULAR. */

L80:
    i__1 = *n;
    r_cnjg(&q__2, &t[*n + *n * t_dim1]);
    c_div(&q__1, &b[*n], &q__2);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = j;
	i__3 = j;
	i__4 = jj - 1;
	cdotc_(&q__2, &i__4, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	i__2 = j;
	r_cnjg(&q__2, &t[j + j * t_dim1]);
	c_div(&q__1, &b[j], &q__2);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L90: */
    }
L100:
    goto L140;

/*        SOLVE CTRANS(T)*X=B FOR T UPPER TRIANGULAR. */

L110:
    r_cnjg(&q__2, &t[t_dim1 + 1]);
    c_div(&q__1, &b[1], &q__2);
    b[1].r = q__1.r, b[1].i = q__1.i;
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j - 1;
	cdotc_(&q__2, &i__4, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	i__2 = j;
	r_cnjg(&q__2, &t[j + j * t_dim1]);
	c_div(&q__1, &b[j], &q__2);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* ctrsl_ */


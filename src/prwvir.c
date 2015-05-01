/* prwvir.f -- translated by f2c (version 12.02.01).
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

/* DECK PRWVIR */
/* Subroutine */ int prwvir_(integer *key, integer *ipage, integer *lpg, real 
	*sx, integer *ix)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;

    static integer iaddr;
    static logical first;
    static integer ipagef;
    extern /* Subroutine */ int sreadp_(integer *, integer *, real *, integer 
	    *, integer *), sopenm_(integer *, integer *);
    static integer istart;
    extern /* Subroutine */ int swritp_(integer *, integer *, real *, integer 
	    *, integer *);

/* ***BEGIN PROLOGUE  PRWVIR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PRWVIR-S, DPRWVR-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*     PRWVIR LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SPARSE MATRIX */
/*     STORAGE SCHEME.  THE PAGE STORAGE IS ON RANDOM ACCESS DISK. */
/*     PRWVIR IS PART OF THE SPARSE LP PACKAGE, SPLP. */

/*     KEY       IS A FLAG WHICH INDICATES WHETHER A READ OR WRITE */
/*               OPERATION IS TO BE PERFORMED. A VALUE OF KEY=1 INDICATES */
/*               A READ. A VALUE OF KEY=2 INDICATES A WRITE. */
/*     IPAGE     IS THE PAGE OF MATRIX MN WE ARE ACCESSING. */
/*     LPG       IS THE LENGTH OF THE PAGE. */
/*   SX(*),IX(*) IS THE MATRIX DATA. */

/*     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LRWVIR, */
/*     SANDIA LABS. REPT. SAND78-0785. */
/*     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  SOPENM, SREADP, SWRITP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   891009  Removed unreferenced variables.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB) */
/* ***END PROLOGUE  PRWVIR */
    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  PRWVIR */

/*     COMPUTE STARTING ADDRESS OF PAGE. */

    ipagef = sx[3];
    istart = ix[3] + 5;

/*     OPEN RANDOM ACCESS FILE NUMBER IPAGEF, IF FIRST PAGE WRITE. */

    first = sx[4] == zero;
    if (! first) {
	goto L20002;
    }
    sopenm_(&ipagef, lpg);
    sx[4] = one;

/*     PERFORM EITHER A READ OR A WRITE. */

L20002:
    iaddr = (*ipage << 1) - 1;
    if (! (*key == 1)) {
	goto L20005;
    }
    sreadp_(&ipagef, &ix[istart], &sx[istart], lpg, &iaddr);
    goto L20006;
L20005:
    if (! (*key == 2)) {
	goto L10001;
    }
    swritp_(&ipagef, &ix[istart], &sx[istart], lpg, &iaddr);
L10001:
L20006:
    return 0;
} /* prwvir_ */


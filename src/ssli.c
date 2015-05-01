/* ssli.f -- translated by f2c (version 12.02.01).
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

/* DECK SSLI */
/* Subroutine */ int ssli_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer nel;
    extern /* Subroutine */ int ssli2_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *);
    static integer locel, lociel, locjel;

/* ***BEGIN PROLOGUE  SSLI */
/* ***PURPOSE  SLAP MSOLVE for Lower Triangle Matrix. */
/*            This routine acts as an interface between the SLAP generic */
/*            MSOLVE calling convention and the routine that actually */
/*                      -1 */
/*            computes L  B = X. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A3 */
/* ***TYPE      SINGLE PRECISION (SSLI-S, DSLI-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */
/*       It is assumed that RWORK and IWORK have initialized with */
/*       the information required for SSLI2: */
/*          IWORK(1) = NEL */
/*          IWORK(2) = Starting location of IEL in IWORK. */
/*          IWORK(3) = Starting location of JEL in IWORK. */
/*          IWORK(4) = Starting location of EL in RWORK. */
/*       See the DESCRIPTION of SSLI2 for details. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  SSLI2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871119  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   921113  Corrected C***CATEGORY line.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  SSLI */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  SSLI */

    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    nel = iwork[1];
    lociel = iwork[2];
    locjel = iwork[3];
    locel = iwork[4];
    ssli2_(n, &b[1], &x[1], &nel, &iwork[lociel], &iwork[locjel], &rwork[
	    locel]);

    return 0;
/* ------------- LAST LINE OF SSLI FOLLOWS ---------------------------- */
} /* ssli_ */


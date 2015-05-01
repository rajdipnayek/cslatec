/* dsluti.f -- translated by f2c (version 12.02.01).
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

/* DECK DSLUTI */
/* Subroutine */ int dsluti_(integer *n, doublereal *b, doublereal *x, 
	integer *nelt, integer *ia, integer *ja, doublereal *a, integer *isym,
	 doublereal *rwork, integer *iwork)
{
    static integer locl, locu, locil, locjl, lociu, locju;
    extern /* Subroutine */ int dslui4_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *);
    static integer locdin;

/* ***BEGIN PROLOGUE  DSLUTI */
/* ***PURPOSE  SLAP MTSOLV for LDU Factorization. */
/*            This routine acts as an interface between the SLAP generic */
/*            MTSOLV calling convention and the routine that actually */
/*                           -T */
/*            computes  (LDU)  B = X. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2E */
/* ***TYPE      DOUBLE PRECISION (SSLUTI-S, DSLUTI-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */
/*       It is assumed that RWORK and IWORK have initialized with */
/*       the information required for DSLUI4: */
/*          IWORK(1) = Starting location of IL in IWORK. */
/*          IWORK(2) = Starting location of JL in IWORK. */
/*          IWORK(3) = Starting location of IU in IWORK. */
/*          IWORK(4) = Starting location of JU in IWORK. */
/*          IWORK(5) = Starting location of L in RWORK. */
/*          IWORK(6) = Starting location of DINV in RWORK. */
/*          IWORK(7) = Starting location of U in RWORK. */
/*       See the DESCRIPTION of DSLUI4 for details. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DSLUI4 */
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
/* ***END PROLOGUE  DSLUTI */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  DSLUTI */

/*         Pull out the pointers to the L, D and U matrices and call */
/*         the workhorse routine. */

    /* Parameter adjustments */
    --a;
    --x;
    --b;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    locil = iwork[1];
    locjl = iwork[2];
    lociu = iwork[3];
    locju = iwork[4];
    locl = iwork[5];
    locdin = iwork[6];
    locu = iwork[7];

    dslui4_(n, &b[1], &x[1], &iwork[locil], &iwork[locjl], &rwork[locl], &
	    rwork[locdin], &iwork[lociu], &iwork[locju], &rwork[locu]);

    return 0;
/* ------------- LAST LINE OF DSLUTI FOLLOWS ---------------------------- */
} /* dsluti_ */


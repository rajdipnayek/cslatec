/* ssllti.f -- translated by f2c (version 12.02.01).
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

/* DECK SSLLTI */
/* Subroutine */ int ssllti_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer nel, locel;
    extern /* Subroutine */ int sllti2_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, real *);
    static integer lociel, locdin, locjel;

/* ***BEGIN PROLOGUE  SSLLTI */
/* ***PURPOSE  SLAP MSOLVE for LDL' (IC) Factorization. */
/*            This routine acts as an interface between the SLAP generic */
/*            MSOLVE calling convention and the routine that actually */
/*                           -1 */
/*            computes (LDL')  B = X. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2E */
/* ***TYPE      SINGLE PRECISION (SSLLTI-S, DSLLTI-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */
/*       It is assumed that RWORK and IWORK have initialized with */
/*       the information required for SLLTI2: */
/*          IWORK(1) = NEL */
/*          IWORK(2) = Starting location of IEL in IWORK. */
/*          IWORK(3) = Starting location of JEL in IWORK. */
/*          IWORK(4) = Starting location of EL in RWORK. */
/*          IWORK(5) = Starting location of DINV in RWORK. */
/*       See the DESCRIPTION of SLLTI2 for details. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  SLLTI2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871119  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Corrected conversion error.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/*   921113  Corrected C***CATEGORY line.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  SSLLTI */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  SSLLTI */
    /* Parameter adjustments */
    --b;
    --x;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    nel = iwork[1];
    lociel = iwork[3];
    locjel = iwork[2];
    locel = iwork[4];
    locdin = iwork[5];
    sllti2_(n, &b[1], &x[1], &nel, &iwork[lociel], &iwork[locjel], &rwork[
	    locel], &rwork[locdin]);

    return 0;
/* ------------- LAST LINE OF SSLLTI FOLLOWS ---------------------------- */
} /* ssllti_ */


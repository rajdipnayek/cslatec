/* dbvsup.f -- translated by f2c (version 12.02.01).
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
    doublereal c__, xsav;
    integer igofxd, inhomo, ivp, ncompd, nfcd;
} dml8sz_;

#define dml8sz_1 dml8sz_

struct {
    doublereal aed, red, tol;
    integer nxptsd, nicd, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivd, numort, nfcc, icoco;
} dml18j_;

#define dml18j_1 dml18j_

struct {
    integer kkkzpw, needw, neediw, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, 
	    k11, l1, l2, kkkint, lllint;
} dml17b_;

#define dml17b_1 dml17b_

struct {
    doublereal px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} dml15t_;

#define dml15t_1 dml15t_

struct {
    doublereal uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} dml5mc_;

#define dml5mc_1 dml5mc_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__10 = 10;
static integer c__5 = 5;

/* DECK DBVSUP */
/* Subroutine */ int dbvsup_(doublereal *y, integer *nrowy, integer *ncomp, 
	doublereal *xpts, integer *nxpts, doublereal *a, integer *nrowa, 
	doublereal *alpha, integer *nic, doublereal *b, integer *nrowb, 
	doublereal *beta, integer *nfc, integer *igofx, doublereal *re, 
	doublereal *ae, integer *iflag, doublereal *work, integer *ndw, 
	integer *iwork, integer *ndiw, integer *neqivp)
{
    /* System generated locals */
    address a__1[10], a__2[5];
    integer a_dim1, a_offset, b_dim1, b_offset, y_dim1, y_offset, i__1, i__2[
	    10], i__3[5];
    char ch__1[194], ch__2[122];

    /* Local variables */
    static integer j, k, is, non, kkkg, ndeq, kkks, kkku, kkkv, kpts;
    static char xern1[8], xern2[8], xern3[8], xern4[8];
    static integer lllip, kkkws;
    extern /* Subroutine */ int dmacon_(void);
    static integer kkkcoe, kkkcof, lllcof;
    extern /* Subroutine */ int dexbvp_(doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer kkksud, kkksvc, lllsud, lllsvc, kkkyhp, nitemp, kkksto, 
	    llliws, mxnoni, nrtemp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer mxnonr, nxptsm;

    /* Fortran I/O blocks */
    static icilist io___29 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___31 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___33 = { 0, xern3, 0, "(I8)", 8, 1 };
    static icilist io___35 = { 0, xern4, 0, "(I8)", 8, 1 };
    static icilist io___36 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___37 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DBVSUP */
/* ***PURPOSE  Solve a linear two-point boundary value problem using */
/*            superposition coupled with an orthonormalization procedure */
/*            and a variable-step integration scheme. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  I1B1 */
/* ***TYPE      DOUBLE PRECISION (BVSUP-S, DBVSUP-D) */
/* ***KEYWORDS  ORTHONORMALIZATION, SHOOTING, */
/*             TWO-POINT BOUNDARY VALUE PROBLEM */
/* ***AUTHOR  Scott, M. R., (SNLA) */
/*           Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */

/*     Subroutine DBVSUP solves a linear two-point boundary-value problem */
/*     of the form */
/*                        DY/DX = MATRIX(X,U)*Y(X) + G(X,U) */
/*                A*Y(XINITIAL) = ALPHA ,  B*Y(XFINAL) = BETA */

/*     coupled with the solution of the initial value problem */

/*                        DU/DX = F(X,U) */
/*                      U(XINITIAL) = ETA */

/* ********************************************************************** */
/*     ABSTRACT */
/*        The method of solution uses superposition coupled with an */
/*     orthonormalization procedure and a variable-step integration */
/*     scheme.  Each time the superposition solutions start to */
/*     lose their numerical linear independence, the vectors are */
/*     reorthonormalized before integration proceeds.  The underlying */
/*     principle of the algorithm is then to piece together the */
/*     intermediate (orthogonalized) solutions, defined on the various */
/*     subintervals, to obtain the desired solutions. */

/* ********************************************************************** */
/*     INPUT to DBVSUP */
/* ********************************************************************** */

/*     NROWY = actual row dimension of Y in calling program. */
/*             NROWY must be .GE. NCOMP */

/*     NCOMP = number of components per solution vector. */
/*             NCOMP is equal to number of original differential */
/*             equations.  NCOMP = NIC + NFC. */

/*     XPTS = desired output points for solution. They must be monotonic. */
/*            XINITIAL = XPTS(1) */
/*            XFINAL = XPTS(NXPTS) */

/*     NXPTS = number of output points. */

/*     A(NROWA,NCOMP) = boundary condition matrix at XINITIAL */
/*                      must be contained in (NIC,NCOMP) sub-matrix. */

/*     NROWA = actual row dimension of A in calling program, */
/*             NROWA must be .GE. NIC. */

/*     ALPHA(NIC+NEQIVP) = boundary conditions at XINITIAL. */
/*                         If NEQIVP .GT. 0 (see below), the boundary */
/*                         conditions at XINITIAL for the initial value */
/*                         equations must be stored starting in */
/*                         position (NIC + 1) of ALPHA. */
/*                         Thus,  ALPHA(NIC+K) = ETA(K). */

/*     NIC = number of boundary conditions at XINITIAL. */

/*     B(NROWB,NCOMP) = boundary condition matrix at XFINAL. */
/*                      Must be contained in (NFC,NCOMP) sub-matrix. */

/*     NROWB = actual row dimension of B in calling program, */
/*             NROWB must be .GE. NFC. */

/*     BETA(NFC) = boundary conditions at XFINAL. */

/*     NFC = number of boundary conditions at XFINAL. */

/*     IGOFX =0 -- The inhomogeneous term G(X) is identically zero. */
/*           =1 -- The inhomogeneous term G(X) is not identically zero. */
/*                 (if IGOFX=1, then Subroutine DGVEC (or DUVEC) must be */
/*                  supplied). */

/*     RE = relative error tolerance used by the integrator. */
/*          (see one of the integrators) */

/*     AE = absolute error tolerance used by the integrator. */
/*          (see one of the integrators) */
/* **NOTE-  RE and AE should not both be zero. */

/*     IFLAG = a status parameter used principally for output. */
/*             However, for efficient solution of problems which */
/*             are originally defined as COMPLEX*16 valued (but */
/*             converted to double precision systems to use this code), */
/*             the user must set IFLAG=13 on input. See the comment */
/*             below for more information on solving such problems. */

/*     WORK(NDW) = floating point array used for internal storage. */

/*     NDW = actual dimension of work array allocated by user. */
/*           An estimate for NDW can be computed from the following */
/*            NDW = 130 + NCOMP**2 * (6 + NXPTS/2 + expected number of */
/*                                           orthonormalizations/8) */
/*           For the disk or tape storage mode, */
/*            NDW = 6 * NCOMP**2 + 10 * NCOMP + 130 */
/*  However, when the ADAMS integrator is to be used, the estimates are */
/*            NDW = 130 + NCOMP**2 * (13 + NXPTS/2 + expected number of */
/*                                           orthonormalizations/8) */
/*    and     NDW = 13 * NCOMP**2 + 22 * NCOMP + 130   , respectively. */

/*     IWORK(NDIW) = integer array used for internal storage. */

/*     NDIW = actual dimension of IWORK array allocated by user. */
/*            An estimate for NDIW can be computed from the following */
/*            NDIW = 68 + NCOMP * (1 + expected number of */
/*                                            orthonormalizations) */
/* **NOTE --  the amount of storage required is problem dependent and may */
/*            be difficult to predict in advance.  Experience has shown */
/*            that for most problems 20 or fewer orthonormalizations */
/*            should suffice. If the problem cannot be completed with the */
/*            allotted storage, then a message will be printed which */
/*            estimates the amount of storage necessary. In any case, the */
/*            user can examine the IWORK array for the actual storage */
/*            requirements, as described in the output information below. */

/*     NEQIVP = number of auxiliary initial value equations being added */
/*              to the boundary value problem. */
/* **NOTE -- Occasionally the coefficients  matrix  and/or  G  may be */
/*           functions which depend on the independent variable  X  and */
/*           on  U, the solution of an auxiliary initial value problem. */
/*           In order to avoid the difficulties associated with */
/*           interpolation, the auxiliary equations may be solved */
/*           simultaneously with the given boundary value problem. */
/*           This initial value problem may be linear or nonlinear. */
/*                 See SAND77-1328 for an example. */


/*     The user must supply subroutines DFMAT, DGVEC, DUIVP and DUVEC, */
/*     when needed (they must be so named), to evaluate the derivatives */
/*     as follows */

/*        A. DFMAT must be supplied. */

/*              SUBROUTINE DFMAT(X,Y,YP) */
/*              X = independent variable (input to DFMAT) */
/*              Y = dependent variable vector (input to DFMAT) */
/*              YP = DY/DX = derivative vector (output from DFMAT) */

/*            Compute the derivatives for the homogeneous problem */
/*              YP(I) = DY(I)/DX = MATRIX(X) * Y(I)  , I = 1,...,NCOMP */

/*            When (NEQIVP .GT. 0) and  matrix  is dependent on  U  as */
/*            well as on  X, the following common statement must be */
/*            included in DFMAT */
/*                    COMMON /DMLIVP/ NOFST */
/*            for convenience, the  U  vector is stored at the bottom */
/*            of the  Y  array.  Thus, during any call to DFMAT, */
/*            U(I) is referenced by  Y(NOFST + I). */


/*            Subroutine DBVDER calls DFMAT NFC times to evaluate the */
/*            homogeneous equations and, if necessary, it calls DFMAT */
/*            once in evaluating the particular solution. since X remains */
/*            unchanged in this sequence of calls it is possible to */
/*            realize considerable computational savings for complicated */
/*            and expensive evaluations of the matrix entries. To do this */
/*            the user merely passes a variable, say XS, via common where */
/*            XS is defined in the main program to be any value except */
/*            the initial X. Then the non-constant elements of matrix(x) */
/*            appearing in the differential equations need only be */
/*            computed if X is unequal to XS, whereupon XS is reset to X. */


/*        B. If  NEQIVP .GT. 0 ,  DUIVP must also be supplied. */

/*              SUBROUTINE DUIVP(X,U,UP) */
/*              X = independent variable (input to DUIVP) */
/*              U = dependent variable vector (input to DUIVP) */
/*              UP = DU/DX = derivative vector (output from DUIVP) */

/*            Compute the derivatives for the auxiliary initial value eqs */
/*              UP(I) = DU(I)/DX, I = 1,...,NEQIVP. */

/*            Subroutine DBVDER calls DUIVP once to evaluate the */
/*            derivatives for the auxiliary initial value equations. */


/*        C. If  NEQIVP = 0  and  IGOFX = 1 ,  DGVEC must be supplied. */

/*              SUBROUTINE DGVEC(X,G) */
/*              X = independent variable (input to DGVEC) */
/*              G = vector of inhomogeneous terms G(X) (output from */
/*              DGVEC) */

/*            Compute the inhomogeneous terms G(X) */
/*                G(I) = G(X) values for I = 1,...,NCOMP. */

/*            Subroutine DBVDER calls DGVEC in evaluating the particular */
/*            solution provided G(X) is not identically zero. Thus, when */
/*            IGOFX=0, the user need not write a DGVEC subroutine. Also, */
/*            the user does not have to bother with the computational */
/*            savings scheme for DGVEC as this is automatically achieved */
/*            via the DBVDER subroutine. */


/*        D. If  NEQIVP .GT. 0  and  IGOFX = 1 ,  DUVEC must be supplied. */

/*             SUBROUTINE DUVEC(X,U,G) */
/*             X = independent variable (input to DUVEC) */
/*             U = dependent variable vector from the auxiliary initial */
/*                 value problem    (input to DUVEC) */
/*             G = array of inhomogeneous terms G(X,U)(output from DUVEC) */

/*            Compute the inhomogeneous terms G(X,U) */
/*                G(I) = G(X,U) values for I = 1,...,NCOMP. */

/*            Subroutine DBVDER calls DUVEC in evaluating the particular */
/*            solution provided G(X,U) is not identically zero.  Thus, */
/*            when IGOFX=0, the user need not write a DUVEC subroutine. */



/*     The following is optional input to DBVSUP to give user more */
/*     flexibility in use of code.  See SAND75-0198, SAND77-1328, */
/*     SAND77-1690, SAND78-0522, and SAND78-1501 for more information. */

/* ****CAUTION -- The user must zero out IWORK(1),...,IWORK(15) */
/*                prior to calling DBVSUP. These locations define */
/*                optional input and must be zero unless set to special */
/*                values by the user as described below. */

/*     IWORK(1) -- number of orthonormalization points. */
/*                 A value need be set only if IWORK(11) = 1 */

/*     IWORK(9) -- integrator and orthonormalization parameter */
/*                 (default value is 1) */
/*                 1 = RUNGE-KUTTA-FEHLBERG code using GRAM-SCHMIDT test. */
/*                 2 = ADAMS code using GRAM-SCHMIDT test. */

/*     IWORK(11) -- orthonormalization points parameter */
/*                  (default value is 0) */
/*                  0 - orthonormalization points not pre-assigned. */
/*                  1 - orthonormalization points pre-assigned in */
/*                      the first IWORK(1) positions of work. */

/*     IWORK(12) -- storage parameter */
/*                  (default value is 0) */
/*                  0 - all storage in core. */
/*                  LUN - homogeneous and inhomogeneous solutions at */
/*                      output points and orthonormalization information */
/*                      are stored on disk.  The logical unit number to */
/*                      be used for disk I/O (NTAPE) is set to IWORK(12). */

/*     WORK(1),... -- pre-assigned orthonormalization points, stored */
/*                    monotonically, corresponding to the direction */
/*                    of integration. */



/*                 ****************************************************** */
/*                 *** COMPLEX*16 VALUED PROBLEM *** */
/*                 ****************************************************** */
/* **NOTE*** */
/*       Suppose the original boundary value problem is NC equations */
/*     of the form */
/*                   DW/DX = MAT(X,U)*W(X) + H(X,U) */
/*                 R*W(XINITIAL)=GAMMA , S*W(XFINAL)=DELTA */
/*     where all variables are COMPLEX*16 valued. The DBVSUP code can be */
/*     used by converting to a double precision system of size 2*NC. To */
/*     solve the larger dimensioned problem efficiently, the user must */
/*     initialize IFLAG=13 on input and order the vector components */
/*     according to Y(1)=DOUBLE PRECISION(W(1)),...,Y(NC)=DOUBLE */
/*     PRECISION(W(NC)),Y(NC+1)=IMAG(W(1)),...., Y(2*NC)=IMAG(W(NC)). */
/*     Then define */
/*                        ............................................... */
/*                        . DOUBLE PRECISION(MAT)    -IMAG(MAT) . */
/*            MATRIX  =   .                         . */
/*                        . IMAG(MAT)     DOUBLE PRECISION(MAT) . */
/*                        ............................................... */

/*     The matrices A,B and vectors G,ALPHA,BETA must be defined */
/*     similarly. Further details can be found in SAND78-1501. */


/* ********************************************************************** */
/*     OUTPUT from DBVSUP */
/* ********************************************************************** */

/*     Y(NROWY,NXPTS) = solution at specified output points. */

/*     IFLAG Output Values */
/*            =-5 algorithm ,for obtaining starting vectors for the */
/*                special COMPLEX*16 problem structure, was unable to */
/*                obtain the initial vectors satisfying the necessary */
/*                independence criteria. */
/*            =-4 rank of boundary condition matrix A is less than NIC, */
/*                as determined by DLSSUD. */
/*            =-2 invalid input parameters. */
/*            =-1 insufficient number of storage locations allocated for */
/*                WORK or IWORK. */

/*            =0 indicates successful solution. */

/*            =1 a computed solution is returned but uniqueness of the */
/*               solution of the boundary-value problem is questionable. */
/*               For an eigenvalue problem, this should be treated as a */
/*               successful execution since this is the expected mode */
/*               of return. */
/*            =2 a computed solution is returned but the existence of the */
/*               solution to the boundary-value problem is questionable. */
/*            =3 a nontrivial solution approximation is returned although */
/*               the boundary condition matrix B*Y(XFINAL) is found to be */
/*               nonsingular (to the desired accuracy level) while the */
/*               right hand side vector is zero. To eliminate this type */
/*               of return, the accuracy of the eigenvalue parameter */
/*               must be improved. */
/*            ***NOTE-We attempt to diagnose the correct problem behavior */
/*               and report possible difficulties by the appropriate */
/*               error flag.  However, the user should probably resolve */
/*               the problem using smaller error tolerances and/or */
/*               perturbations in the boundary conditions or other */
/*               parameters. This will often reveal the correct */
/*               interpretation for the problem posed. */

/*            =13 maximum number of orthonormalizations attained before */
/*                reaching XFINAL. */
/*            =20-flag from integrator (DDERKF or DDEABM) values can */
/*                range from 21 to 25. */
/*            =30 solution vectors form a dependent set. */

/*     WORK(1),...,WORK(IWORK(1)) = orthonormalization points */
/*                                  determined by DBVPOR. */

/*     IWORK(1) = number of orthonormalizations performed by DBVPOR. */

/*     IWORK(2) = maximum number of orthonormalizations allowed as */
/*                calculated from storage allocated by user. */

/*     IWORK(3),IWORK(4),IWORK(5),IWORK(6)   give information about */
/*                actual storage requirements for WORK and IWORK */
/*                arrays.  In particular, */
/*                       required storage for  work array is */
/*        IWORK(3) + IWORK(4)*(expected number of orthonormalizations) */

/*                       required storage for IWORK array is */
/*        IWORK(5) + IWORK(6)*(expected number of orthonormalizations) */

/*     IWORK(8) = final value of exponent parameter used in tolerance */
/*                test for orthonormalization. */

/*     IWORK(16) = number of independent vectors returned from DMGSBV. */
/*                It is only of interest when IFLAG=30 is obtained. */

/*     IWORK(17) = numerically estimated rank of the boundary */
/*                 condition matrix defined from B*Y(XFINAL) */

/* ********************************************************************** */

/*     Necessary machine constants are defined in the Function */
/*     Routine D1MACH. The user must make sure that the values */
/*     set in D1MACH are relevant to the computer being used. */

/* ********************************************************************** */
/* ********************************************************************** */

/* ***REFERENCES  M. R. Scott and H. A. Watts, SUPORT - a computer code */
/*                 for two-point boundary-value problems via */
/*                 orthonormalization, SIAM Journal of Numerical */
/*                 Analysis 14, (1977), pp. 40-70. */
/*               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications */
/*                 of SUPORT, a linear boundary value problem solver */
/*                 Part I - pre-assigning orthonormalization points, */
/*                 auxiliary initial value problem, disk or tape storage, */
/*                 Report SAND77-1328, Sandia Laboratories, Albuquerque, */
/*                 New Mexico, 1977. */
/*               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications */
/*                 of SUPORT, a linear boundary value problem solver */
/*                 Part II - inclusion of an Adams integrator, Report */
/*                 SAND77-1690, Sandia Laboratories, Albuquerque, */
/*                 New Mexico, 1977. */
/*               M. E. Lord and H. A. Watts, Modifications of SUPORT, */
/*                 a linear boundary value problem solver Part III - */
/*                 orthonormalization improvements, Report SAND78-0522, */
/*                 Sandia Laboratories, Albuquerque, New Mexico, 1978. */
/*               H. A. Watts, M. R. Scott and M. E. Lord, Computational */
/*                 solution of complex*16 valued boundary problems, */
/*                 Report SAND78-1501, Sandia Laboratories, */
/*                 Albuquerque, New Mexico, 1978. */
/* ***ROUTINES CALLED  DEXBVP, DMACON, XERMSG */
/* ***COMMON BLOCKS    DML15T, DML17B, DML18J, DML5MC, DML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   890921  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900510  Convert XERRWV calls to XERMSG calls, remove some extraneous */
/*           comments.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBVSUP */
/* ********************************************************************** */


/*     ****************************************************************** */
/*         THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH SUBROUTINE */
/*         DBVDER.  THE USER SHOULD NOT ALTER OR USE THIS COMMON BLOCK IN */
/*         THE CALLING PROGRAM. */


/*     ****************************************************************** */
/*         THESE COMMON BLOCKS AID IN REDUCING THE NUMBER OF SUBROUTINE */
/*         ARGUMENTS PREVALENT IN THIS MODULAR STRUCTURE */


/*     ****************************************************************** */
/*         THIS COMMON BLOCK IS USED IN SUBROUTINES DBVSUP,DBVPOR,DRKFAB, */
/*         DREORT, AND DSTWAY. IT CONTAINS INFORMATION NECESSARY */
/*         FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND A BACKUP */
/*         RESTARTING CAPABILITY. */


/*     ****************************************************************** */
/*         THIS COMMON BLOCK CONTAINS THE MACHINE DEPENDENT PARAMETERS */
/*         USED BY THE CODE */


/*      ***************************************************************** */
/*          SET UP MACHINE DEPENDENT CONSTANTS. */

/* ***FIRST EXECUTABLE STATEMENT  DBVSUP */
    /* Parameter adjustments */
    y_dim1 = *nrowy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --xpts;
    a_dim1 = *nrowa;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --alpha;
    b_dim1 = *nrowb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --beta;
    --work;
    --iwork;

    /* Function Body */
    dmacon_();

/*                       ************************************************ */
/*                           TEST FOR INVALID INPUT */

    if (*nrowy < *ncomp) {
	goto L80;
    }
    if (*ncomp != *nic + *nfc) {
	goto L80;
    }
    if (*nxpts < 2) {
	goto L80;
    }
    if (*nic <= 0) {
	goto L80;
    }
    if (*nrowa < *nic) {
	goto L80;
    }
    if (*nfc <= 0) {
	goto L80;
    }
    if (*nrowb < *nfc) {
	goto L80;
    }
    if (*igofx < 0 || *igofx > 1) {
	goto L80;
    }
    if (*re < 0.) {
	goto L80;
    }
    if (*ae < 0.) {
	goto L80;
    }
    if (*re == 0. && *ae == 0.) {
	goto L80;
    }
/*                          BEGIN BLOCK PERMITTING ...EXITS TO 70 */
    is = 1;
    if (xpts[*nxpts] < xpts[1]) {
	is = 2;
    }
    nxptsm = *nxpts - 1;
    i__1 = nxptsm;
    for (k = 1; k <= i__1; ++k) {
	if (is == 2) {
	    goto L10;
	}
/*                          .........EXIT */
	if (xpts[k + 1] <= xpts[k]) {
	    goto L70;
	}
	goto L20;
L10:
/*                          .........EXIT */
	if (xpts[k] <= xpts[k + 1]) {
	    goto L70;
	}
L20:
/* L30: */
	;
    }

/*                             ****************************************** */
/*                                 CHECK FOR DISK STORAGE */

    kpts = *nxpts;
    dml18j_1.ndisk = 0;
    if (iwork[12] == 0) {
	goto L40;
    }
    dml18j_1.ntape = iwork[12];
    kpts = 1;
    dml18j_1.ndisk = 1;
L40:

/*                             ****************************************** */
/*                                 SET INTEG PARAMETER ACCORDING TO */
/*                                 CHOICE OF INTEGRATOR. */

    dml18j_1.integ = 1;
    if (iwork[9] == 2) {
	dml18j_1.integ = 2;
    }

/*                             ****************************************** */
/*                                 COMPUTE INHOMO */

/*                 ............EXIT */
    if (*igofx == 1) {
	goto L100;
    }
    i__1 = *nic;
    for (j = 1; j <= i__1; ++j) {
/*                 ...............EXIT */
	if (alpha[j] != 0.) {
	    goto L100;
	}
/* L50: */
    }
    i__1 = *nfc;
    for (j = 1; j <= i__1; ++j) {
/*                    ............EXIT */
	if (beta[j] != 0.) {
	    goto L90;
	}
/* L60: */
    }
    dml8sz_1.inhomo = 3;
/*              ...............EXIT */
    goto L110;
L70:
L80:
    *iflag = -2;
/*     ..................EXIT */
    goto L220;
L90:
    dml8sz_1.inhomo = 2;
/*              ......EXIT */
    goto L110;
L100:
    dml8sz_1.inhomo = 1;
L110:

/*              ********************************************************* */
/*                  TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN */
/*                  SOLVING A COMPLEX*16 VALUED PROBLEM,WE INTRODUCE */
/*                  NFCC=NFC WHILE CHANGING THE INTERNAL VALUE OF NFC */

    dml18j_1.nfcc = *nfc;
    if (*iflag == 13) {
	*nfc /= 2;
    }

/*              ********************************************************* */
/*                  DETERMINE NECESSARY STORAGE REQUIREMENTS */

/*              FOR BASIC ARRAYS IN DBVPOR */
    kkkyhp = *ncomp * (*nfc + 1) + *neqivp;
    kkku = *ncomp * *nfc * kpts;
    kkkv = *ncomp * kpts;
    kkkcoe = dml18j_1.nfcc;
    kkks = *nfc + 1;
    kkksto = *ncomp * (*nfc + 1) + *neqivp + 1;
    kkkg = *ncomp;

/*              FOR ORTHONORMALIZATION RELATED MATTERS */
    dml18j_1.ntp = dml18j_1.nfcc * (dml18j_1.nfcc + 1) / 2;
    dml17b_1.kkkzpw = dml18j_1.ntp + 1 + dml18j_1.nfcc;
    lllip = dml18j_1.nfcc;

/*              FOR ADDITIONAL REQUIRED WORK SPACE */
/*                (DLSSUD) */
    kkksud = (*nic << 2) + (*nrowa + 1) * *ncomp;
    lllsud = *nic;
/*              (DVECS) */
/* Computing 2nd power */
    i__1 = dml18j_1.nfcc;
    kkksvc = (dml18j_1.nfcc << 2) + 1 + (i__1 * i__1 << 1);
    lllsvc = dml18j_1.nfcc << 1;

    ndeq = *ncomp * *nfc + *neqivp;
    if (dml8sz_1.inhomo == 1) {
	ndeq += *ncomp;
    }
    switch (dml18j_1.integ) {
	case 1:  goto L120;
	case 2:  goto L130;
    }
/*              (DDERKF) */
L120:
    dml17b_1.kkkint = ndeq * 7 + 33;
    dml17b_1.lllint = 34;
    goto L140;
/*              (DDEABM) */
L130:
    dml17b_1.kkkint = ndeq * 21 + 130;
    dml17b_1.lllint = 51;
L140:

/*              (COEF) */
/* Computing 2nd power */
    i__1 = dml18j_1.nfcc;
    kkkcof = dml18j_1.nfcc * 5 + i__1 * i__1;
    lllcof = dml18j_1.nfcc + 3;

/* Computing MAX */
    i__1 = max(kkksud,kkksvc), i__1 = max(i__1,dml17b_1.kkkint);
    kkkws = max(i__1,kkkcof);
/* Computing MAX */
    i__1 = max(lllsud,lllsvc), i__1 = max(i__1,dml17b_1.lllint);
    llliws = max(i__1,lllcof);

    dml17b_1.needw = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto + kkkg + 
	    dml17b_1.kkkzpw + kkkws;
    dml17b_1.neediw = lllip + 17 + llliws;
/*              ********************************************************* */
/*                  COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS */
/*                  WITH THE ALLOTTED STORAGE */

    iwork[3] = dml17b_1.needw;
    iwork[4] = dml17b_1.kkkzpw;
    iwork[5] = dml17b_1.neediw;
    iwork[6] = lllip;
    nrtemp = *ndw - dml17b_1.needw;
    nitemp = *ndiw - dml17b_1.neediw;
/*           ...EXIT */
    if (nrtemp < 0) {
	goto L180;
    }
/*           ...EXIT */
    if (nitemp < 0) {
	goto L180;
    }

    if (dml18j_1.ndisk == 0) {
	goto L150;
    }
    non = 0;
    dml18j_1.mxnon = nrtemp;
    goto L160;
L150:

    mxnonr = nrtemp / dml17b_1.kkkzpw;
    mxnoni = nitemp / lllip;
    dml18j_1.mxnon = min(mxnonr,mxnoni);
    non = dml18j_1.mxnon;
L160:

    iwork[2] = dml18j_1.mxnon;

/*              ********************************************************* */
/*                  CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS */

    dml18j_1.nopg = 0;
/*        ......EXIT */
    if (iwork[11] != 1) {
	goto L210;
    }
    if (dml18j_1.mxnon < iwork[1]) {
	goto L170;
    }
    dml18j_1.nopg = 1;
    dml18j_1.mxnon = iwork[1];
    work[dml18j_1.mxnon + 1] = xpts[*nxpts] * 2. - xpts[1];
/*        .........EXIT */
    goto L210;
L170:
L180:

    *iflag = -1;
    if (dml18j_1.ndisk != 1) {
	s_wsfi(&io___29);
	do_fio(&c__1, (char *)&dml17b_1.needw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___31);
	do_fio(&c__1, (char *)&dml17b_1.kkkzpw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___33);
	do_fio(&c__1, (char *)&dml17b_1.neediw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___35);
	do_fio(&c__1, (char *)&lllip, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 35, a__1[0] = "REQUIRED STORAGE FOR WORK ARRAY IS ";
	i__2[1] = 8, a__1[1] = xern1;
	i__2[2] = 3, a__1[2] = " + ";
	i__2[3] = 8, a__1[3] = xern2;
	i__2[4] = 44, a__1[4] = "*(EXPECTED NUMBER OF ORTHONORMALIZATIONS) $$"
		;
	i__2[5] = 36, a__1[5] = "REQUIRED STORAGE FOR IWORK ARRAY IS ";
	i__2[6] = 8, a__1[6] = xern3;
	i__2[7] = 3, a__1[7] = " + ";
	i__2[8] = 8, a__1[8] = xern4;
	i__2[9] = 41, a__1[9] = "*(EXPECTED NUMBER OF ORTHONORMALIZATIONS)";
	s_cat(ch__1, a__1, i__2, &c__10, (ftnlen)194);
	xermsg_("SLATEC", "DBVSUP", ch__1, &c__1, &c__0, (ftnlen)6, (ftnlen)6,
		 (ftnlen)194);
    } else {
	s_wsfi(&io___36);
	do_fio(&c__1, (char *)&dml17b_1.needw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___37);
	do_fio(&c__1, (char *)&dml17b_1.neediw, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__3[0] = 35, a__2[0] = "REQUIRED STORAGE FOR WORK ARRAY IS ";
	i__3[1] = 8, a__2[1] = xern1;
	i__3[2] = 35, a__2[2] = " + NUMBER OF ORTHONOMALIZATIONS. $$";
	i__3[3] = 36, a__2[3] = "REQUIRED STORAGE FOR IWORK ARRAY IS ";
	i__3[4] = 8, a__2[4] = xern2;
	s_cat(ch__2, a__2, i__3, &c__5, (ftnlen)122);
	xermsg_("SLATEC", "DBVSUP", ch__2, &c__1, &c__0, (ftnlen)6, (ftnlen)6,
		 (ftnlen)122);
    }
    return 0;

/*        *************************************************************** */
/*            ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS */

/*         (Z) */
L210:
    dml17b_1.k1 = dml18j_1.mxnon + 2;
/*        (P) */
    dml17b_1.k2 = dml17b_1.k1 + dml18j_1.ntp * (non + 1);
/*        (W) */
    dml17b_1.k3 = dml17b_1.k2 + dml18j_1.nfcc * (non + 1);
/*        (YHP) */
    dml17b_1.k4 = dml17b_1.k3 + kkkyhp;
/*        (U) */
    dml17b_1.k5 = dml17b_1.k4 + kkku;
/*        (V) */
    dml17b_1.k6 = dml17b_1.k5 + kkkv;
/*        (COEF) */
    dml17b_1.k7 = dml17b_1.k6 + kkkcoe;
/*        (S) */
    dml17b_1.k8 = dml17b_1.k7 + kkks;
/*        (STOWA) */
    dml17b_1.k9 = dml17b_1.k8 + kkksto;
/*        (G) */
    dml17b_1.k10 = dml17b_1.k9 + kkkg;
    dml17b_1.k11 = dml17b_1.k10 + kkkws;
/*                  REQUIRED ADDITIONAL DOUBLE PRECISION WORK SPACE */
/*                  STARTS AT WORK(K10) AND EXTENDS TO WORK(K11-1) */

/*           FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL */
/*           INPUT AND OUTPUT ITEMS */
/*        (IP) */
    dml17b_1.l1 = dml18j_1.nfcc * (non + 1) + 18;
    dml17b_1.l2 = dml17b_1.l1 + llliws;
/*                   REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1) */
/*                   AND EXTENDS TO IWORK(L2-1) */

/*        *************************************************************** */
/*            SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION */

    dml18j_1.nps = 0;
    if (iwork[10] == 1) {
	dml18j_1.nps = 1;
    }

/*        *************************************************************** */
/*            SET PIVOTING PARAMETER */

    dml18j_1.indpvt = 0;
    if (iwork[15] == 1) {
	dml18j_1.indpvt = 1;
    }

/*        *************************************************************** */
/*            SET OTHER COMMON BLOCK PARAMETERS */

    dml8sz_1.nfcd = *nfc;
    dml8sz_1.ncompd = *ncomp;
    dml8sz_1.igofxd = *igofx;
    dml18j_1.nxptsd = *nxpts;
    dml18j_1.nicd = *nic;
    dml18j_1.red = *re;
    dml18j_1.aed = *ae;
    dml18j_1.neqivd = *neqivp;
    dml15t_1.mnswot = 20;
    if (iwork[13] == -1) {
	dml15t_1.mnswot = max(1,iwork[14]);
    }
    dml15t_1.xbeg = xpts[1];
    dml15t_1.xend = xpts[*nxpts];
    dml8sz_1.xsav = dml15t_1.xend;
    dml18j_1.icoco = 1;
    if (dml8sz_1.inhomo == 3 && dml18j_1.nopg == 1) {
	work[dml18j_1.mxnon + 1] = dml15t_1.xend;
    }

/*        *************************************************************** */

    dexbvp_(&y[y_offset], nrowy, &xpts[1], &a[a_offset], nrowa, &alpha[1], &b[
	    b_offset], nrowb, &beta[1], iflag, &work[1], &iwork[1]);
    *nfc = dml18j_1.nfcc;
    iwork[17] = iwork[dml17b_1.l1];
L220:
    return 0;
} /* dbvsup_ */


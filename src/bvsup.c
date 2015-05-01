/* bvsup.f -- translated by f2c (version 12.02.01).
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
    real c__, xsav;
    integer igofxd, inhomo, ivp, ncompd, nfcd;
} ml8sz_;

#define ml8sz_1 ml8sz_

struct {
    real aed, red, tol;
    integer nxptsd, nicd, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivd, numort, nfcc, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

struct {
    integer kkkzpw, needw, neediw, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, 
	    k11, l1, l2, kkkint, lllint;
} ml17bw_;

#define ml17bw_1 ml17bw_

struct {
    real px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} ml15to_;

#define ml15to_1 ml15to_

struct {
    real uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} ml5mco_;

#define ml5mco_1 ml5mco_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__10 = 10;
static integer c__5 = 5;

/* DECK BVSUP */
/* Subroutine */ int bvsup_(real *y, integer *nrowy, integer *ncomp, real *
	xpts, integer *nxpts, real *a, integer *nrowa, real *alpha, integer *
	nic, real *b, integer *nrowb, real *beta, integer *nfc, integer *
	igofx, real *re, real *ae, integer *iflag, real *work, integer *ndw, 
	integer *iwork, integer *ndiw, integer *neqivp)
{
    /* System generated locals */
    address a__1[10], a__2[5];
    integer y_dim1, y_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2[
	    10], i__3[5];
    char ch__1[194], ch__2[122];

    /* Local variables */
    static integer j, k, is, non, kkkg, ndeq, kkks, kkku, kkkv, kpts;
    static char xern1[8], xern2[8], xern3[8], xern4[8];
    extern /* Subroutine */ int macon_(void);
    static integer lllip;
    extern /* Subroutine */ int exbvp_(real *, integer *, real *, real *, 
	    integer *, real *, real *, integer *, real *, integer *, real *, 
	    integer *);
    static integer kkkws, kkkcoe, kkkcof, lllcof, kkksud, kkksvc, lllsud, 
	    lllsvc, kkkyhp, nitemp, nrtemp, kkksto, llliws;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer mxnoni, mxnonr, nxptsm;

    /* Fortran I/O blocks */
    static icilist io___26 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___28 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___30 = { 0, xern3, 0, "(I8)", 8, 1 };
    static icilist io___32 = { 0, xern4, 0, "(I8)", 8, 1 };
    static icilist io___33 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___34 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  BVSUP */
/* ***PURPOSE  Solve a linear two-point boundary value problem using */
/*            superposition coupled with an orthonormalization procedure */
/*            and a variable-step integration scheme. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  I1B1 */
/* ***TYPE      SINGLE PRECISION (BVSUP-S, DBVSUP-D) */
/* ***KEYWORDS  ORTHONORMALIZATION, SHOOTING, */
/*             TWO-POINT BOUNDARY VALUE PROBLEM */
/* ***AUTHOR  Scott, M. R., (SNLA) */
/*           Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*     Subroutine BVSUP solves a LINEAR two-point boundary-value problem */
/*     of the form */
/*                        dY/dX = MATRIX(X,U)*Y(X) + G(X,U) */
/*                A*Y(Xinitial) = ALPHA ,  B*Y(Xfinal) = BETA */

/*     Coupled with the solution of the initial value problem */

/*                        dU/dX = F(X,U) */
/*                      U(Xinitial) = ETA */

/* ********************************************************************** */
/*     Abstract */
/*        The method of solution uses superposition coupled with an */
/*     orthonormalization procedure and a variable-step integration */
/*     scheme.  Each time the superposition solutions start to */
/*     lose their numerical linear independence, the vectors are */
/*     reorthonormalized before integration proceeds.  The underlying */
/*     principle of the algorithm is then to piece together the */
/*     intermediate (orthogonalized) solutions, defined on the various */
/*     subintervals, to obtain the desired solutions. */

/* ********************************************************************** */
/*     INPUT to BVSUP */
/* ********************************************************************** */

/*     NROWY = Actual row dimension of Y in calling program. */
/*             NROWY must be .GE. NCOMP */

/*     NCOMP = Number of components per solution vector. */
/*             NCOMP is equal to number of original differential */
/*             equations.  NCOMP = NIC + NFC. */

/*     XPTS = Desired output points for solution. They must be monotonic. */
/*            Xinitial = XPTS(1) */
/*            Xfinal = XPTS(NXPTS) */

/*     NXPTS = Number of output points */

/*     A(NROWA,NCOMP) = Boundary condition matrix at Xinitial, */
/*                      must be contained in (NIC,NCOMP) sub-matrix. */

/*     NROWA = Actual row dimension of A in calling program, */
/*             NROWA must be .GE. NIC. */

/*     ALPHA(NIC+NEQIVP) = Boundary conditions at Xinitial. */
/*                         If NEQIVP .GT. 0 (see below), the boundary */
/*                         conditions at Xinitial for the initial value */
/*                         equations must be stored starting in */
/*                         position (NIC + 1) of ALPHA. */
/*                         Thus,  ALPHA(NIC+K) = ETA(K). */

/*     NIC = Number of boundary conditions at Xinitial. */

/*     B(NROWB,NCOMP) = Boundary condition matrix at Xfinal, */
/*                      must be contained in (NFC,NCOMP) sub-matrix. */

/*     NROWB = Actual row dimension of B in calling program, */
/*             NROWB must be .GE. NFC. */

/*     BETA(NFC) = Boundary conditions at Xfinal. */

/*     NFC = Number of boundary conditions at Xfinal */

/*     IGOFX =0 -- The inhomogeneous term G(X) is identically zero. */
/*           =1 -- The inhomogeneous term G(X) is not identically zero. */
/*                 (if IGOFX=1, then subroutine GVEC (or UVEC) must be */
/*                  supplied). */

/*     RE = Relative error tolerance used by the integrator */
/*          (see one of the integrators) */

/*     AE = Absolute error tolerance used by the integrator */
/*          (see one of the integrators) */
/* **NOTE-  RE and AE should not both be zero. */

/*     IFLAG = A status parameter used principally for output. */
/*             However, for efficient solution of problems which */
/*             are originally defined as complex valued (but */
/*             converted to real systems to use this code), the */
/*             user must set IFLAG=13 on input. See the comment below */
/*             for more information on solving such problems. */

/*     WORK(NDW) = Floating point array used for internal storage. */

/*     NDW = Actual dimension of WORK array allocated by user. */
/*           An estimate for NDW can be computed from the following */
/*            NDW = 130 + NCOMP**2 * (6 + NXPTS/2 + expected number of */
/*                                                orthonormalizations/8) */
/*             For the DISK or TAPE storage mode, */
/*            NDW = 6 * NCOMP**2 + 10 * NCOMP + 130 */
/*  However, when the ADAMS integrator is to be used, the estimates are */
/*            NDW = 130 + NCOMP**2 * (13 + NXPTS/2 + expected number of */
/*                                                orthonormalizations/8) */
/*    and     NDW = 13 * NCOMP**2 + 22 * NCOMP + 130   , respectively. */

/*     IWORK(NDIW) = Integer array used for internal storage. */

/*     NDIW = Actual dimension of IWORK array allocated by user. */
/*            An estimate for NDIW can be computed from the following */
/*            NDIW = 68 + NCOMP * (1 + expected number of */
/*                                        orthonormalizations) */
/* **NOTE --  The amount of storage required is problem dependent and may */
/*            be difficult to predict in advance. Experience has shown */
/*            that for most problems 20 or fewer orthonormalizations */
/*            should suffice. If the problem cannot be completed with the */
/*            allotted storage, then a message will be printed which */
/*            estimates the amount of storage necessary. In any case, the */
/*            user can examine the IWORK array for the actual storage */
/*            requirements, as described in the output information below. */

/*     NEQIVP = Number of auxiliary initial value equations being added */
/*              to the boundary value problem. */
/* **NOTE -- Occasionally the coefficients  MATRIX  and/or  G  may be */
/*           functions which depend on the independent variable  X  and */
/*           on  U, the solution of an auxiliary initial value problem. */
/*           In order to avoid the difficulties associated with */
/*           interpolation, the auxiliary equations may be solved */
/*           simultaneously with the given boundary value problem. */
/*           This initial value problem may be LINEAR or NONLINEAR. */
/*                 See SAND77-1328 for an example. */


/*     The user must supply subroutines FMAT, GVEC, UIVP and UVEC, when */
/*     needed (they MUST be so named), to evaluate the derivatives */
/*     as follows */

/*        A. FMAT must be supplied. */

/*              SUBROUTINE FMAT(X,Y,YP) */
/*              X = Independent variable (input to FMAT) */
/*              Y = Dependent variable vector (input to FMAT) */
/*              YP = dY/dX = Derivative vector (output from FMAT) */

/*            Compute the derivatives for the HOMOGENEOUS problem */
/*              YP(I) = dY(I)/dX = MATRIX(X) * Y(I)  , I = 1,...,NCOMP */

/*            When (NEQIVP .GT. 0) and  MATRIX  is dependent on  U  as */
/*            well as on  X, the following common statement must be */
/*            included in FMAT */
/*                    COMMON /MLIVP/ NOFST */
/*            For convenience, the  U  vector is stored at the bottom */
/*            of the  Y  array.  Thus, during any call to FMAT, */
/*            U(I) is referenced by  Y(NOFST + I). */


/*            Subroutine BVDER calls FMAT NFC times to evaluate the */
/*            homogeneous equations and, if necessary, it calls FMAT once */
/*            in evaluating the particular solution. Since X remains */
/*            unchanged in this sequence of calls it is possible to */
/*            realize considerable computational savings for complicated */
/*            and expensive evaluations of the MATRIX entries. To do this */
/*            the user merely passes a variable, say XS, via COMMON where */
/*            XS is defined in the main program to be any value except */
/*            the initial X. Then the non-constant elements of MATRIX(X) */
/*            appearing in the differential equations need only be */
/*            computed if X is unequal to XS, whereupon XS is reset to X. */


/*        B. If  NEQIVP .GT. 0 ,  UIVP must also be supplied. */

/*              SUBROUTINE UIVP(X,U,UP) */
/*              X = Independent variable (input to UIVP) */
/*              U = Dependent variable vector (input to UIVP) */
/*              UP = dU/dX = Derivative vector (output from UIVP) */

/*            Compute the derivatives for the auxiliary initial value eqs */
/*              UP(I) = dU(I)/dX, I = 1,...,NEQIVP. */

/*            Subroutine BVDER calls UIVP once to evaluate the */
/*            derivatives for the auxiliary initial value equations. */


/*        C. If  NEQIVP = 0  and  IGOFX = 1 ,  GVEC must be supplied. */

/*              SUBROUTINE GVEC(X,G) */
/*              X = Independent variable (input to GVEC) */
/*              G = Vector of inhomogeneous terms G(X) (output from GVEC) */

/*            Compute the inhomogeneous terms G(X) */
/*                G(I) = G(X) values for I = 1,...,NCOMP. */

/*            Subroutine BVDER calls GVEC in evaluating the particular */
/*            solution provided G(X) is NOT identically zero. Thus, when */
/*            IGOFX=0, the user need NOT write a GVEC subroutine. Also, */
/*            the user does not have to bother with the computational */
/*            savings scheme for GVEC as this is automatically achieved */
/*            via the BVDER subroutine. */


/*        D. If  NEQIVP .GT. 0  and  IGOFX = 1 ,  UVEC must be supplied. */

/*              SUBROUTINE UVEC(X,U,G) */
/*              X = Independent variable (input to UVEC) */
/*              U = Dependent variable vector from the auxiliary initial */
/*                  value problem    (input to UVEC) */
/*              G = Array of inhomogeneous terms G(X,U)(output from UVEC) */

/*            Compute the inhomogeneous terms G(X,U) */
/*                G(I) = G(X,U) values for I = 1,...,NCOMP. */

/*            Subroutine BVDER calls UVEC in evaluating the particular */
/*            solution provided G(X,U) is NOT identically zero.  Thus, */
/*            when IGOFX=0, the user need NOT write a UVEC subroutine. */



/*     The following is optional input to BVSUP to give the user more */
/*     flexibility in use of the code.  See SAND75-0198 , SAND77-1328 , */
/*     SAND77-1690,SAND78-0522, and SAND78-1501 for more information. */

/* ****CAUTION -- The user MUST zero out IWORK(1),...,IWORK(15) */
/*                prior to calling BVSUP. These locations define optional */
/*                input and MUST be zero UNLESS set to special values by */
/*                the user as described below. */

/*     IWORK(1) -- Number of orthonormalization points. */
/*                 A value need be set only if IWORK(11) = 1 */

/*     IWORK(9) -- Integrator and orthonormalization parameter */
/*                 (default value is 1) */
/*                 1 = RUNGE-KUTTA-FEHLBERG code using GRAM-SCHMIDT test. */
/*                 2 = ADAMS code using GRAM-SCHMIDT TEST. */

/*     IWORK(11) -- Orthonormalization points parameter */
/*                  (default value is 0) */
/*                  0 - Orthonormalization points not pre-assigned. */
/*                  1 - Orthonormalization points pre-assigned in */
/*                      the first IWORK(1) positions of WORK. */

/*     IWORK(12) -- Storage parameter */
/*                  (default value is 0) */
/*                  0 - All storage IN CORE */
/*                LUN - Homogeneous and inhomogeneous solutions at */
/*                     output points and orthonormalization information */
/*                     are stored on DISK.  The logical unit number to be */
/*                     used for DISK I/O (NTAPE) is set to IWORK(12). */

/*     WORK(1),... -- Pre-assigned orthonormalization points, stored */
/*                    monotonically, corresponding to the direction */
/*                    of integration. */



/*                 ****************************** */
/*                 *** COMPLEX VALUED PROBLEM *** */
/*                 ****************************** */
/* **NOTE*** */
/*       Suppose the original boundary value problem is NC equations */
/*     of the form */
/*                   dW/dX = MAT(X,U)*W(X) + H(X,U) */
/*                 R*W(Xinitial)=GAMMA , S*W(Xfinal)=DELTA */

/*     where all variables are complex valued. The BVSUP code can be */
/*     used by converting to a real system of size 2*NC. To solve the */
/*     larger dimensioned problem efficiently,  the user must initialize */
/*     IFLAG=13 on input and order the vector components according to */
/*     Y(1)=real(W(1)),...,Y(NC)=real(W(NC)),Y(NC+1)=imag(W(1)),...., */
/*     Y(2*NC)=imag(W(NC)). Then define */
/*                        ........................... */
/*                        . real(MAT)    -imag(MAT) . */
/*            MATRIX  =   .                         . */
/*                        . imag(MAT)     real(MAT) . */
/*                        ........................... */

/*     The matrices A,B and vectors G,ALPHA,BETA must be defined */
/*     similarly. Further details can be found in SAND78-1501. */


/* ********************************************************************** */
/*     OUTPUT from BVSUP */
/* ********************************************************************** */

/*     Y(NROWY,NXPTS) = Solution at specified output points. */

/*     IFLAG output values */
/*            =-5 Algorithm ,for obtaining starting vectors for the */
/*                special complex problem structure, was unable to obtain */
/*                the initial vectors satisfying the necessary */
/*                independence criteria. */
/*            =-4 Rank of boundary condition matrix A is less than NIC, */
/*                as determined by LSSUDS. */
/*            =-2 Invalid input parameters. */
/*            =-1 Insufficient number of storage locations allocated for */
/*                WORK or IWORK. */

/*            =0 Indicates successful solution */

/*            =1 A computed solution is returned but UNIQUENESS of the */
/*               solution of the boundary-value problem is questionable. */
/*               For an eigenvalue problem, this should be treated as a */
/*               successful execution since this is the expected mode */
/*               of return. */
/*            =2 A computed solution is returned but the EXISTENCE of the */
/*               solution to the boundary-value problem is questionable. */
/*            =3 A nontrivial solution approximation is returned although */
/*               the boundary condition matrix B*Y(Xfinal) is found to be */
/*               nonsingular (to the desired accuracy level) while the */
/*               right hand side vector is zero. To eliminate this type */
/*               of return, the accuracy of the eigenvalue parameter */
/*               must be improved. */
/*           ***NOTE- We attempt to diagnose the correct problem behavior */
/*               and report possible difficulties by the appropriate */
/*               error flag.  However, the user should probably resolve */
/*               the problem using smaller error tolerances and/or */
/*               perturbations in the boundary conditions or other */
/*               parameters. This will often reveal the correct */
/*               interpretation for the problem posed. */

/*            =13 Maximum number of orthonormalizations attained before */
/*                reaching Xfinal. */
/*            =20-flag from integrator (DERKF or DEABM) values can range */
/*                from 21 to 25. */
/*            =30 Solution vectors form a dependent set. */

/*     WORK(1),...,WORK(IWORK(1)) = Orthonormalization points */
/*                                  determined by BVPOR. */

/*     IWORK(1) = Number of orthonormalizations performed by BVPOR. */

/*     IWORK(2) = Maximum number of orthonormalizations allowed as */
/*                calculated from storage allocated by user. */

/*     IWORK(3),IWORK(4),IWORK(5),IWORK(6)   Give information about */
/*                actual storage requirements for WORK and IWORK */
/*                arrays.  In particular, */
/*                       required storage for  WORK array is */
/*        IWORK(3) + IWORK(4)*(expected number of orthonormalizations) */

/*                       required storage for IWORK array is */
/*        IWORK(5) + IWORK(6)*(expected number of orthonormalizations) */

/*     IWORK(8) = Final value of exponent parameter used in tolerance */
/*                test for orthonormalization. */

/*     IWORK(16) = Number of independent vectors returned from MGSBV. */
/*                 It is only of interest when IFLAG=30 is obtained. */

/*     IWORK(17) = Numerically estimated rank of the boundary */
/*                 condition matrix defined from B*Y(Xfinal) */

/* ********************************************************************** */

/*     Necessary machine constants are defined in the function */
/*     routine R1MACH. The user must make sure that the values */
/*     set in R1MACH are relevant to the computer being used. */

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
/* ***ROUTINES CALLED  EXBVP, MACON, XERMSG */
/* ***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   890921  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BVSUP */
/* ********************************************************************** */



/* ********************************************************************** */
/*     THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH SUBROUTINE */
/*     BVDER.  THE USER SHOULD NOT ALTER OR USE THIS COMMON BLOCK IN THE */
/*     CALLING PROGRAM. */


/* ********************************************************************** */
/*     THESE COMMON BLOCKS AID IN REDUCING THE NUMBER OF SUBROUTINE */
/*     ARGUMENTS PREVALENT IN THIS MODULAR STRUCTURE */


/* ********************************************************************** */
/*     THIS COMMON BLOCK IS USED IN SUBROUTINES BVSUP,BVPOR,RKFAB, */
/*     REORT, AND STWAY. IT CONTAINS INFORMATION NECESSARY */
/*     FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND A BACKUP */
/*     RESTARTING CAPABILITY. */


/* ********************************************************************** */
/*     THIS COMMON BLOCK CONTAINS THE MACHINE DEPENDENT PARAMETERS */
/*     USED BY THE CODE */


/* ********************************************************************** */
/*     SET UP MACHINE DEPENDENT CONSTANTS. */

/* ***FIRST EXECUTABLE STATEMENT  BVSUP */
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
    macon_();

/* ********************************************************************** */
/*     TEST FOR INVALID INPUT */

    if (*nrowy < *ncomp) {
	goto L20;
    }
    if (*ncomp != *nic + *nfc) {
	goto L20;
    }
    if (*nxpts < 2) {
	goto L20;
    }
    if (*nic <= 0) {
	goto L20;
    }
    if (*nrowa < *nic) {
	goto L20;
    }
    if (*nfc <= 0) {
	goto L20;
    }
    if (*nrowb < *nfc) {
	goto L20;
    }
    if (*igofx < 0 || *igofx > 1) {
	goto L20;
    }
    if (*re < 0.f) {
	goto L20;
    }
    if (*ae < 0.f) {
	goto L20;
    }
    if (*re == 0.f && *ae == 0.f) {
	goto L20;
    }
    is = 1;
    if (xpts[*nxpts] < xpts[1]) {
	is = 2;
    }
    nxptsm = *nxpts - 1;
    i__1 = nxptsm;
    for (k = 1; k <= i__1; ++k) {
	if (is == 2) {
	    goto L12;
	}
	if (xpts[k + 1] <= xpts[k]) {
	    goto L20;
	}
	goto L13;
L12:
	if (xpts[k] <= xpts[k + 1]) {
	    goto L20;
	}
L13:
	;
    }
    goto L30;
L20:
    *iflag = -2;
    return 0;
L30:

/* ********************************************************************** */
/*     CHECK FOR DISK STORAGE */

    kpts = *nxpts;
    ml18jr_1.ndisk = 0;
    if (iwork[12] == 0) {
	goto L35;
    }
    ml18jr_1.ntape = iwork[12];
    kpts = 1;
    ml18jr_1.ndisk = 1;
L35:

/* ********************************************************************** */
/*     SET INTEG PARAMETER ACCORDING TO CHOICE OF INTEGRATOR. */

    ml18jr_1.integ = 1;
    if (iwork[9] == 2) {
	ml18jr_1.integ = 2;
    }

/* ********************************************************************** */
/*     COMPUTE INHOMO */

    if (*igofx == 1) {
	goto L43;
    }
    i__1 = *nic;
    for (j = 1; j <= i__1; ++j) {
	if (alpha[j] != 0.f) {
	    goto L43;
	}
/* L40: */
    }
    i__1 = *nfc;
    for (j = 1; j <= i__1; ++j) {
	if (beta[j] != 0.f) {
	    goto L42;
	}
/* L41: */
    }
    ml8sz_1.inhomo = 3;
    goto L45;
L42:
    ml8sz_1.inhomo = 2;
    goto L45;
L43:
    ml8sz_1.inhomo = 1;
L45:

/* ********************************************************************** */
/*     TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN SOLVING A */
/*     COMPLEX VALUED PROBLEM,WE INTRODUCE NFCC=NFC WHILE CHANGING */
/*     THE INTERNAL VALUE OF NFC */

    ml18jr_1.nfcc = *nfc;
    if (*iflag == 13) {
	*nfc /= 2;
    }

/* ********************************************************************** */
/*     DETERMINE NECESSARY STORAGE REQUIREMENTS */

/* FOR BASIC ARRAYS IN BVPOR */
    kkkyhp = *ncomp * (*nfc + 1) + *neqivp;
    kkku = *ncomp * *nfc * kpts;
    kkkv = *ncomp * kpts;
    kkkcoe = ml18jr_1.nfcc;
    kkks = *nfc + 1;
    kkksto = *ncomp * (*nfc + 1) + *neqivp + 1;
    kkkg = *ncomp;

/* FOR ORTHONORMALIZATION RELATED MATTERS */
    ml18jr_1.ntp = ml18jr_1.nfcc * (ml18jr_1.nfcc + 1) / 2;
    ml17bw_1.kkkzpw = ml18jr_1.ntp + 1 + ml18jr_1.nfcc;
    lllip = ml18jr_1.nfcc;

/* FOR ADDITIONAL REQUIRED WORK SPACE */
/*   (LSSUDS) */
    kkksud = (*nic << 2) + (*nrowa + 1) * *ncomp;
    lllsud = *nic;
/*   (SVECS) */
/* Computing 2nd power */
    i__1 = ml18jr_1.nfcc;
    kkksvc = (ml18jr_1.nfcc << 2) + 1 + (i__1 * i__1 << 1);
    lllsvc = ml18jr_1.nfcc << 1;

    ndeq = *ncomp * *nfc + *neqivp;
    if (ml8sz_1.inhomo == 1) {
	ndeq += *ncomp;
    }
    switch (ml18jr_1.integ) {
	case 1:  goto L51;
	case 2:  goto L52;
    }
/*   (DERKF) */
L51:
    ml17bw_1.kkkint = ndeq * 7 + 33;
    ml17bw_1.lllint = 34;
    goto L55;
/*   (DEABM) */
L52:
    ml17bw_1.kkkint = ndeq * 21 + 130;
    ml17bw_1.lllint = 51;

/*   (COEF) */
L55:
/* Computing 2nd power */
    i__1 = ml18jr_1.nfcc;
    kkkcof = ml18jr_1.nfcc * 5 + i__1 * i__1;
    lllcof = ml18jr_1.nfcc + 3;

/* Computing MAX */
    i__1 = max(kkksud,kkksvc), i__1 = max(i__1,ml17bw_1.kkkint);
    kkkws = max(i__1,kkkcof);
/* Computing MAX */
    i__1 = max(lllsud,lllsvc), i__1 = max(i__1,ml17bw_1.lllint);
    llliws = max(i__1,lllcof);

    ml17bw_1.needw = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto + kkkg + 
	    ml17bw_1.kkkzpw + kkkws;
    ml17bw_1.neediw = lllip + 17 + llliws;
/* ********************************************************************** */
/*     COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS WITH THE */
/*     ALLOTTED STORAGE */

    iwork[3] = ml17bw_1.needw;
    iwork[4] = ml17bw_1.kkkzpw;
    iwork[5] = ml17bw_1.neediw;
    iwork[6] = lllip;
    nrtemp = *ndw - ml17bw_1.needw;
    nitemp = *ndiw - ml17bw_1.neediw;
    if (nrtemp < 0) {
	goto L70;
    }
    if (nitemp >= 0) {
	goto L75;
    }

L70:
    *iflag = -1;
    if (ml18jr_1.ndisk != 1) {
	s_wsfi(&io___26);
	do_fio(&c__1, (char *)&ml17bw_1.needw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___28);
	do_fio(&c__1, (char *)&ml17bw_1.kkkzpw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___30);
	do_fio(&c__1, (char *)&ml17bw_1.neediw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___32);
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
	xermsg_("SLATEC", "BVSUP", ch__1, &c__1, &c__0, (ftnlen)6, (ftnlen)5, 
		(ftnlen)194);
    } else {
	s_wsfi(&io___33);
	do_fio(&c__1, (char *)&ml17bw_1.needw, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___34);
	do_fio(&c__1, (char *)&ml17bw_1.neediw, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__3[0] = 35, a__2[0] = "REQUIRED STORAGE FOR WORK ARRAY IS ";
	i__3[1] = 8, a__2[1] = xern1;
	i__3[2] = 35, a__2[2] = " + NUMBER OF ORTHONOMALIZATIONS. $$";
	i__3[3] = 36, a__2[3] = "REQUIRED STORAGE FOR IWORK ARRAY IS ";
	i__3[4] = 8, a__2[4] = xern2;
	s_cat(ch__2, a__2, i__3, &c__5, (ftnlen)122);
	xermsg_("SLATEC", "BVSUP", ch__2, &c__1, &c__0, (ftnlen)6, (ftnlen)5, 
		(ftnlen)122);
    }
    return 0;

L75:
    if (ml18jr_1.ndisk == 0) {
	goto L77;
    }
    non = 0;
    ml18jr_1.mxnon = nrtemp;
    goto L78;

L77:
    mxnonr = nrtemp / ml17bw_1.kkkzpw;
    mxnoni = nitemp / lllip;
    ml18jr_1.mxnon = min(mxnonr,mxnoni);
    non = ml18jr_1.mxnon;

L78:
    iwork[2] = ml18jr_1.mxnon;

/* ********************************************************************** */
/*     CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS */

    ml18jr_1.nopg = 0;
    if (iwork[11] != 1) {
	goto L85;
    }
    if (ml18jr_1.mxnon < iwork[1]) {
	goto L70;
    }
    ml18jr_1.nopg = 1;
    ml18jr_1.mxnon = iwork[1];
    work[ml18jr_1.mxnon + 1] = xpts[*nxpts] * 2.f - xpts[1];
L85:

/* ********************************************************************** */
/*     ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS */

/*  (Z) */
    ml17bw_1.k1 = ml18jr_1.mxnon + 2;
/*  (P) */
    ml17bw_1.k2 = ml17bw_1.k1 + ml18jr_1.ntp * (non + 1);
/*  (W) */
    ml17bw_1.k3 = ml17bw_1.k2 + ml18jr_1.nfcc * (non + 1);
/*  (YHP) */
    ml17bw_1.k4 = ml17bw_1.k3 + kkkyhp;
/*  (U) */
    ml17bw_1.k5 = ml17bw_1.k4 + kkku;
/*  (V) */
    ml17bw_1.k6 = ml17bw_1.k5 + kkkv;
/*  (COEF) */
    ml17bw_1.k7 = ml17bw_1.k6 + kkkcoe;
/*  (S) */
    ml17bw_1.k8 = ml17bw_1.k7 + kkks;
/*  (STOWA) */
    ml17bw_1.k9 = ml17bw_1.k8 + kkksto;
/*  (G) */
    ml17bw_1.k10 = ml17bw_1.k9 + kkkg;
    ml17bw_1.k11 = ml17bw_1.k10 + kkkws;
/*            REQUIRED ADDITIONAL REAL WORK SPACE STARTS AT WORK(K10) */
/*            AND EXTENDS TO WORK(K11-1) */

/*     FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL */
/*     INPUT AND OUTPUT ITEMS */
/*  (IP) */
    ml17bw_1.l1 = ml18jr_1.nfcc * (non + 1) + 18;
    ml17bw_1.l2 = ml17bw_1.l1 + llliws;
/*            REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1) */
/*            AND EXTENDS TO IWORK(L2-1) */

/* ********************************************************************** */
/*     SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION */

    ml18jr_1.nps = 0;
    if (iwork[10] == 1) {
	ml18jr_1.nps = 1;
    }

/* ********************************************************************** */
/*     SET PIVOTING PARAMETER */

    ml18jr_1.indpvt = 0;
    if (iwork[15] == 1) {
	ml18jr_1.indpvt = 1;
    }

/* ********************************************************************** */
/*     SET OTHER COMMON BLOCK PARAMETERS */

    ml8sz_1.nfcd = *nfc;
    ml8sz_1.ncompd = *ncomp;
    ml8sz_1.igofxd = *igofx;
    ml18jr_1.nxptsd = *nxpts;
    ml18jr_1.nicd = *nic;
    ml18jr_1.red = *re;
    ml18jr_1.aed = *ae;
    ml18jr_1.neqivd = *neqivp;
    ml15to_1.mnswot = 20;
    if (iwork[13] == -1) {
	ml15to_1.mnswot = max(1,iwork[14]);
    }
    ml15to_1.xbeg = xpts[1];
    ml15to_1.xend = xpts[*nxpts];
    ml8sz_1.xsav = ml15to_1.xend;
    ml18jr_1.icoco = 1;
    if (ml8sz_1.inhomo == 3 && ml18jr_1.nopg == 1) {
	work[ml18jr_1.mxnon + 1] = ml15to_1.xend;
    }

/* ********************************************************************** */

    exbvp_(&y[y_offset], nrowy, &xpts[1], &a[a_offset], nrowa, &alpha[1], &b[
	    b_offset], nrowb, &beta[1], iflag, &work[1], &iwork[1]);
    *nfc = ml18jr_1.nfcc;
    iwork[17] = iwork[ml17bw_1.l1];
    return 0;
} /* bvsup_ */


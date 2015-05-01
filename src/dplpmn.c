/* dplpmn.f -- translated by f2c (version 12.02.01).
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
    doublereal small;
    integer lp, lenl, lenu, ncp, lrow, lcol;
} la05dd_;

#define la05dd_1 la05dd_

/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static integer c__5 = 5;

/* DECK DPLPMN */
/* Subroutine */ int dplpmn_(U_fp dusrmt, integer *mrelas, integer *nvars__, 
	doublereal *costs, doublereal *prgopt, doublereal *dattrv, doublereal 
	*bl, doublereal *bu, integer *ind, integer *info, doublereal *primal, 
	doublereal *duals, doublereal *amat, doublereal *csc, doublereal *
	colnrm, doublereal *erd, doublereal *erp, doublereal *basmat, 
	doublereal *wr, doublereal *rz, doublereal *rg, doublereal *rprim, 
	doublereal *rhs, doublereal *ww, integer *lmx, integer *lbm, integer *
	ibasis, integer *ibb, integer *imat, integer *ibrc, integer *ipr, 
	integer *iwr)
{
    /* Format strings */
    static char fmt_20013[] = "";
    static char fmt_20018[] = "";
    static char fmt_20019[] = "";
    static char fmt_20020[] = "";
    static char fmt_20024[] = "";
    static char fmt_20029[] = "";
    static char fmt_20030[] = "";
    static char fmt_20031[] = "";
    static char fmt_20032[] = "";
    static char fmt_20036[] = "";
    static char fmt_20037[] = "";
    static char fmt_20038[] = "";
    static char fmt_20039[] = "";
    static char fmt_20040[] = "";
    static char fmt_20044[] = "";
    static char fmt_20050[] = "";
    static char fmt_20051[] = "";
    static char fmt_20096[] = "";
    static char fmt_20134[] = "";
    static char fmt_20135[] = "";
    static char fmt_20141[] = "";
    static char fmt_20154[] = "";
    static char fmt_20204[] = "";
    static char fmt_20242[] = "";
    static char fmt_20243[] = "";
    static char fmt_20246[] = "";
    static char fmt_20233[] = "";
    static char fmt_20237[] = "";
    static char fmt_20244[] = "";
    static char fmt_20245[] = "";
    static char fmt_20251[] = "";
    static char fmt_20275[] = "";
    static char fmt_20267[] = "";

    /* System generated locals */
    address a__1[5];
    integer ibrc_dim1, ibrc_offset, i__1, i__2, i__3[5];
    doublereal d__1;
    char ch__1[118];
    alist al__1, al__2;
    static logical equiv_7[8];
    static integer equiv_15[8];
    static doublereal equiv_22[7];

    /* Local variables */
    static integer i__, j, k;
    static doublereal gg;
    static integer np;
    static doublereal uu, aij;
#define idg (equiv_15)
    static doublereal one;
    static integer lpg;
#define eps (equiv_22)
    static integer key;
#define npp (equiv_15 + 6)
    static integer lpr;
    static doublereal rzj;
    static integer n20080, n20206, n20046, n20119, n20172, n20058, n20247, 
	    n20252, n20271, n20098, n20276, n20283, n20290, lpr1;
#define abig (equiv_22 + 2)
    static logical feas;
    static integer ibas;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer idum[1], nx0091, nx0106, nx0066;
#define lprg (equiv_15 + 7)
    static integer nerr;
    static doublereal rdum[1];
    static integer itlp;
    static doublereal size;
#define tune (equiv_22 + 5)
    static doublereal xval;
    static integer iopt;
#define lopt (equiv_7)
    static doublereal zero;
#define ropt (equiv_22)
    static integer npr010, npr011, npr012, npr004, npr005, npr006, npr007, 
	    npr008, npr009, npr013, npr014, npr015;
    static char xern1[8], xern2[8];
    extern /* Subroutine */ int la05bd_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, logical *);
    static integer ipage, nredc;
#define itbrc (equiv_15 + 5)
    static doublereal scalr, theta;
    static logical unbnd;
#define isave (equiv_15 + 2)
    static doublereal upbnd;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static doublereal anorm;
    static logical found;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nparm;
    extern /* Subroutine */ int dpopt_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *);
    static logical trans;
#define tolls (equiv_22 + 4)
    extern /* Subroutine */ int dvout_(integer *, doublereal *, char *, 
	    integer *, ftnlen);
    static integer jstrt;
    extern /* Subroutine */ int ivout_(integer *, integer *, char *, integer *
	    , ftnlen);
#define ipagef (equiv_15 + 1)
    static integer iplace;
    static logical redbas;
    static integer ileave;
    static doublereal xlamda;
    extern /* Subroutine */ int dplpce_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, logical *);
#define asmall (equiv_22 + 1)
    extern /* Subroutine */ int dplpfe_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *);
    static doublereal factor;
    static logical finite;
    extern /* Subroutine */ int dplpdm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, logical *), 
	    dplpfl_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, logical *);
#define tolabs (equiv_22 + 6)
#define colscp (equiv_7 + 4)
    static doublereal erdnrm;
#define savedt (equiv_7 + 3)
#define minprb (equiv_7 + 6)
    static doublereal dirnrm;
#define contin (equiv_7)
#define costsc (equiv_22 + 3)
#define cstscp (equiv_7 + 5)
    static logical singlr;
    static doublereal dulnrm;
#define stpedg (equiv_7 + 7)
#define usrbas (equiv_7 + 1)
    extern /* Subroutine */ int dpintm_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal resnrm;
#define kprint (equiv_15 + 4)
    extern /* Subroutine */ int dplpup_(U_fp, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, logical *, doublereal *, 
	    doublereal *);
    static doublereal rhsnrm;
    extern /* Subroutine */ int dpinit_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, logical *), xermsg_(char *, char *, char *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen);
#define intopt (equiv_15)
    static doublereal scosts, rprnrm;
#define sizeup (equiv_7 + 2)
    static logical zerolv;
#define mxitlp (equiv_15 + 3)
    extern /* Subroutine */ int dpnnzr_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    static integer ienter;
    extern /* Subroutine */ int dprwpg_(integer *, integer *, integer *, 
	    doublereal *, integer *);
    static integer ntries;
    extern /* Subroutine */ int dpincw_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *), dplpmu_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, logical *, logical *, 
	    logical *), sclosm_(integer *);

    /* Fortran I/O blocks */
    static cilist io___74 = { 0, 0, 0, 0, 0 };
    static cilist io___79 = { 0, 0, 0, 0, 0 };
    static cilist io___81 = { 0, 0, 0, 0, 0 };
    static cilist io___82 = { 0, 0, 0, 0, 0 };
    static cilist io___83 = { 0, 0, 0, 0, 0 };
    static cilist io___84 = { 0, 0, 0, 0, 0 };
    static icilist io___100 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___102 = { 0, xern2, 0, "(I8)", 8, 1 };


    /* Assigned format variables */
    static char *npr004_fmt, *npr005_fmt, *npr006_fmt, *npr007_fmt, *
	    npr008_fmt, *npr009_fmt, *npr010_fmt, *npr011_fmt, *npr012_fmt, *
	    npr013_fmt, *npr014_fmt, *npr015_fmt;

/* ***BEGIN PROLOGUE  DPLPMN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPLPMN-S, DPLPMN-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     MARVEL OPTION(S).. OUTPUT=YES/NO TO ELIMINATE PRINTED OUTPUT. */
/*     THIS DOES NOT APPLY TO THE CALLS TO THE ERROR PROCESSOR. */

/*     MAIN SUBROUTINE FOR DSPLP PACKAGE. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DASUM, DCOPY, DDOT, DPINCW, DPINIT, DPINTM, DPLPCE, */
/*                    DPLPDM, DPLPFE, DPLPFL, DPLPMU, DPLPUP, DPNNZR, */
/*                    DPOPT, DPRWPG, DVOUT, IVOUT, LA05BD, SCLOSM, XERMSG */
/* ***COMMON BLOCKS    LA05DD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  DPLPMN */


/*     ARRAY LOCAL VARIABLES */
/*     NAME(LENGTH)          DESCRIPTION */

/*     COSTS(NVARS)          COST COEFFICIENTS */
/*     PRGOPT( )             OPTION VECTOR */
/*     DATTRV( )             DATA TRANSFER VECTOR */
/*     PRIMAL(NVARS+MRELAS)  AS OUTPUT IT IS PRIMAL SOLUTION OF LP. */
/*                           INTERNALLY, THE FIRST NVARS POSITIONS HOLD */
/*                           THE COLUMN CHECK SUMS.  THE NEXT MRELAS */
/*                           POSITIONS HOLD THE CLASSIFICATION FOR THE */
/*                           BASIC VARIABLES  -1 VIOLATES LOWER */
/*                           BOUND, 0 FEASIBLE, +1 VIOLATES UPPER BOUND */
/*     DUALS(MRELAS+NVARS)   DUAL SOLUTION. INTERNALLY HOLDS R.H. SIDE */
/*                           AS FIRST MRELAS ENTRIES. */
/*     AMAT(LMX)             SPARSE FORM OF DATA MATRIX */
/*     IMAT(LMX)             SPARSE FORM OF DATA MATRIX */
/*     BL(NVARS+MRELAS)      LOWER BOUNDS FOR VARIABLES */
/*     BU(NVARS+MRELAS)      UPPER BOUNDS FOR VARIABLES */
/*     IND(NVARS+MRELAS)     INDICATOR FOR VARIABLES */
/*     CSC(NVARS)            COLUMN SCALING */
/*     IBASIS(NVARS+MRELAS)  COLS. 1-MRELAS ARE BASIC, REST ARE NON-BASIC */
/*     IBB(NVARS+MRELAS)     INDICATOR FOR NON-BASIC VARS., POLARITY OF */
/*                           VARS., AND POTENTIALLY INFINITE VARS. */
/*                           IF IBB(J).LT.0, VARIABLE J IS BASIC */
/*                           IF IBB(J).GT.0, VARIABLE J IS NON-BASIC */
/*                           IF IBB(J).EQ.0, VARIABLE J HAS TO BE IGNORED */
/*                           BECAUSE IT WOULD CAUSE UNBOUNDED SOLN. */
/*                           WHEN MOD(IBB(J),2).EQ.0, VARIABLE IS AT ITS */
/*                           UPPER BOUND, OTHERWISE IT IS AT ITS LOWER */
/*                           BOUND */
/*     COLNRM(NVARS)         NORM OF COLUMNS */
/*     ERD(MRELAS)           ERRORS IN DUAL VARIABLES */
/*     ERP(MRELAS)           ERRORS IN PRIMAL VARIABLES */
/*     BASMAT(LBM)           BASIS MATRIX FOR HARWELL SPARSE CODE */
/*     IBRC(LBM,2)           ROW AND COLUMN POINTERS FOR BASMAT(*) */
/*     IPR(2*MRELAS)         WORK ARRAY FOR HARWELL SPARSE CODE */
/*     IWR(8*MRELAS)         WORK ARRAY FOR HARWELL SPARSE CODE */
/*     WR(MRELAS)            WORK ARRAY FOR HARWELL SPARSE CODE */
/*     RZ(NVARS+MRELAS)      REDUCED COSTS */
/*     RPRIM(MRELAS)         INTERNAL PRIMAL SOLUTION */
/*     RG(NVARS+MRELAS)      COLUMN WEIGHTS */
/*     WW(MRELAS)            WORK ARRAY */
/*     RHS(MRELAS)           HOLDS TRANSLATED RIGHT HAND SIDE */

/*     SCALAR LOCAL VARIABLES */
/*     NAME       TYPE         DESCRIPTION */

/*     LMX        INTEGER      LENGTH OF AMAT(*) */
/*     LPG        INTEGER      LENGTH OF PAGE FOR AMAT(*) */
/*     EPS        DOUBLE       MACHINE PRECISION */
/*     TUNE       DOUBLE       PARAMETER TO SCALE ERROR ESTIMATES */
/*     TOLLS      DOUBLE       RELATIVE TOLERANCE FOR SMALL RESIDUALS */
/*     TOLABS     DOUBLE       ABSOLUTE TOLERANCE FOR SMALL RESIDUALS. */
/*                             USED IF RELATIVE ERROR TEST FAILS. */
/*                             IN CONSTRAINT EQUATIONS */
/*     FACTOR     DOUBLE      .01--DETERMINES IF BASIS IS SINGULAR */
/*                             OR COMPONENT IS FEASIBLE.  MAY NEED TO */
/*                             BE INCREASED TO 1.D0 ON SHORT WORD */
/*                             LENGTH MACHINES. */
/*     ASMALL     DOUBLE       LOWER BOUND FOR NON-ZERO MAGN. IN AMAT(*) */
/*     ABIG       DOUBLE       UPPER BOUND FOR NON-ZERO MAGN. IN AMAT(*) */
/*     MXITLP     INTEGER      MAXIMUM NUMBER OF ITERATIONS FOR LP */
/*     ITLP       INTEGER      ITERATION COUNTER FOR TOTAL LP ITERS */
/*     COSTSC     DOUBLE       COSTS(*) SCALING */
/*     SCOSTS     DOUBLE       TEMP LOC. FOR COSTSC. */
/*     XLAMDA     DOUBLE       WEIGHT PARAMETER FOR PEN. METHOD. */
/*     ANORM      DOUBLE       NORM OF DATA MATRIX AMAT(*) */
/*     RPRNRM     DOUBLE       NORM OF THE SOLUTION */
/*     DULNRM     DOUBLE       NORM OF THE DUALS */
/*     ERDNRM     DOUBLE       NORM OF ERROR IN DUAL VARIABLES */
/*     DIRNRM     DOUBLE       NORM OF THE DIRECTION VECTOR */
/*     RHSNRM     DOUBLE       NORM OF TRANSLATED RIGHT HAND SIDE VECTOR */
/*     RESNRM     DOUBLE       NORM OF RESIDUAL VECTOR FOR CHECKING */
/*                             FEASIBILITY */
/*     NZBM       INTEGER      NUMBER OF NON-ZEROS IN BASMAT(*) */
/*     LBM        INTEGER      LENGTH OF BASMAT(*) */
/*     SMALL      DOUBLE       EPS*ANORM  USED IN HARWELL SPARSE CODE */
/*     LP         INTEGER      USED IN HARWELL LA05*() PACK AS OUTPUT */
/*                             FILE NUMBER. SET=I1MACH(4) NOW. */
/*     UU         DOUBLE       0.1--USED IN HARWELL SPARSE CODE */
/*                             FOR RELATIVE PIVOTING TOLERANCE. */
/*     GG         DOUBLE       OUTPUT INFO FLAG IN HARWELL SPARSE CODE */
/*     IPLACE     INTEGER      INTEGER USED BY SPARSE MATRIX CODES */
/*     IENTER     INTEGER      NEXT COLUMN TO ENTER BASIS */
/*     NREDC      INTEGER      NO. OF FULL REDECOMPOSITIONS */
/*     KPRINT     INTEGER      LEVEL OF OUTPUT, =0-3 */
/*     IDG        INTEGER      FORMAT AND PRECISION OF OUTPUT */
/*     ITBRC      INTEGER      NO. OF ITERS. BETWEEN RECALCULATING */
/*                             THE ERROR IN THE PRIMAL SOLUTION. */
/*     NPP        INTEGER      NO. OF NEGATIVE REDUCED COSTS REQUIRED */
/*                             IN PARTIAL PRICING */
/*     JSTRT      INTEGER      STARTING PLACE FOR PARTIAL PRICING. */


/*     COMMON BLOCK USED BY LA05 () PACKAGE.. */

/*     SET LP=0 SO NO ERROR MESSAGES WILL PRINT WITHIN LA05 () PACKAGE. */
/* ***FIRST EXECUTABLE STATEMENT  DPLPMN */
    /* Parameter adjustments */
    --costs;
    --prgopt;
    --dattrv;
    --bl;
    --bu;
    --ind;
    --primal;
    --duals;
    --amat;
    --csc;
    --colnrm;
    --erd;
    --erp;
    --basmat;
    --wr;
    --rz;
    --rg;
    --rprim;
    --rhs;
    --ww;
    ibrc_dim1 = *lbm;
    ibrc_offset = 1 + ibrc_dim1;
    ibrc -= ibrc_offset;
    --ibasis;
    --ibb;
    --imat;
    --ipr;
    --iwr;

    /* Function Body */
    la05dd_1.lp = 0;

/*     THE VALUES ZERO AND ONE. */
    zero = 0.;
    one = 1.;
    factor = .01;
    lpg = *lmx - (*nvars__ + 4);
    iopt = 1;
    *info = 0;
    unbnd = FALSE_;
    jstrt = 1;

/*     PROCESS USER OPTIONS IN PRGOPT(*). */
/*     CHECK THAT ANY USER-GIVEN CHANGES ARE WELL-DEFINED. */
    dpopt_(&prgopt[1], mrelas, nvars__, info, &csc[1], &ibasis[1], ropt, 
	    intopt, lopt);
    if (! (*info < 0)) {
	goto L20002;
    }
    goto L30001;
L20002:
    if (! (*contin)) {
	goto L20003;
    }
    goto L30002;
L20006:
    goto L20004;

/*     INITIALIZE SPARSE DATA MATRIX, AMAT(*) AND IMAT(*). */
L20003:
    dpintm_(mrelas, nvars__, &amat[1], &imat[1], lmx, ipagef);

/*     UPDATE MATRIX DATA AND CHECK BOUNDS FOR CONSISTENCY. */
L20004:
    dplpup_((U_fp)dusrmt, mrelas, nvars__, &prgopt[1], &dattrv[1], &bl[1], &
	    bu[1], &ind[1], info, &amat[1], &imat[1], sizeup, asmall, abig);
    if (! (*info < 0)) {
	goto L20007;
    }
    goto L30001;

/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
L20007:
    if (! (*kprint >= 1)) {
	goto L20008;
    }
    goto L30003;
L20011:
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */

/*     INITIALIZATION. SCALE DATA, NORMALIZE BOUNDS, FORM COLUMN */
/*     CHECK SUMS, AND FORM INITIAL BASIS MATRIX. */
L20008:
    dpinit_(mrelas, nvars__, &costs[1], &bl[1], &bu[1], &ind[1], &primal[1], 
	    info, &amat[1], &csc[1], costsc, &colnrm[1], &xlamda, &anorm, &
	    rhs[1], &rhsnrm, &ibasis[1], &ibb[1], &imat[1], lopt);
    if (! (*info < 0)) {
	goto L20012;
    }
    goto L30001;

L20012:
    nredc = 0;
    npr004 = 0;
    npr004_fmt = fmt_20013;
    goto L30004;
L20013:
    if (! singlr) {
	goto L20014;
    }
    nerr = 23;
    xermsg_("SLATEC", "DPLPMN", "IN DSPLP,  A SINGULAR INITIAL BASIS WAS ENC"
	    "OUNTERED.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)52);
    *info = -nerr;
    goto L30001;
L20014:
    npr005 = 0;
    npr005_fmt = fmt_20018;
    goto L30005;
L20018:
    npr006 = 0;
    npr006_fmt = fmt_20019;
    goto L30006;
L20019:
    npr007 = 0;
    npr007_fmt = fmt_20020;
    goto L30007;
L20020:
    if (! (*usrbas)) {
	goto L20021;
    }
    npr008 = 0;
    npr008_fmt = fmt_20024;
    goto L30008;
L20024:
    if (feas) {
	goto L20025;
    }
    nerr = 24;
    xermsg_("SLATEC", "DPLPMN", "IN DSPLP, AN INFEASIBLE INITIAL BASIS WAS E"
	    "NCOUNTERED.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)54);
    *info = -nerr;
    goto L30001;
L20025:
L20021:
    itlp = 0;

/*     LAMDA HAS BEEN SET TO A CONSTANT, PERFORM PENALTY METHOD. */
    npr009 = 0;
    npr009_fmt = fmt_20029;
    goto L30009;
L20029:
    npr010 = 0;
    npr010_fmt = fmt_20030;
    goto L30010;
L20030:
    npr006 = 1;
    npr006_fmt = fmt_20031;
    goto L30006;
L20031:
    npr008 = 1;
    npr008_fmt = fmt_20032;
    goto L30008;
L20032:
    if (feas) {
	goto L20033;
    }

/*     SET LAMDA TO INFINITY BY SETTING COSTSC TO ZERO (SAVE THE VALUE OF */
/*     COSTSC) AND PERFORM STANDARD PHASE-1. */
    if (*kprint >= 2) {
	ivout_(&c__0, idum, "(' ENTER STANDARD PHASE-1')", idg, (ftnlen)27);
    }
    scosts = *costsc;
    *costsc = zero;
    npr007 = 1;
    npr007_fmt = fmt_20036;
    goto L30007;
L20036:
    npr009 = 1;
    npr009_fmt = fmt_20037;
    goto L30009;
L20037:
    npr010 = 1;
    npr010_fmt = fmt_20038;
    goto L30010;
L20038:
    npr006 = 2;
    npr006_fmt = fmt_20039;
    goto L30006;
L20039:
    npr008 = 2;
    npr008_fmt = fmt_20040;
    goto L30008;
L20040:
    if (! feas) {
	goto L20041;
    }

/*     SET LAMDA TO ZERO, COSTSC=SCOSTS, PERFORM STANDARD PHASE-2. */
    if (*kprint > 1) {
	ivout_(&c__0, idum, "(' ENTER STANDARD PHASE-2')", idg, (ftnlen)27);
    }
    xlamda = zero;
    *costsc = scosts;
    npr009 = 2;
    npr009_fmt = fmt_20044;
    goto L30009;
L20044:
L20041:
    goto L20034;
/*     CHECK IF ANY BASIC VARIABLES ARE STILL CLASSIFIED AS */
/*     INFEASIBLE.  IF ANY ARE, THEN THIS MAY NOT YET BE AN */
/*     OPTIMAL POINT.  THEREFORE SET LAMDA TO ZERO AND TRY */
/*     TO PERFORM MORE SIMPLEX STEPS. */
L20033:
    i__ = 1;
    n20046 = *mrelas;
    goto L20047;
L20046:
    ++i__;
L20047:
    if (n20046 - i__ < 0) {
	goto L20048;
    }
    if (primal[i__ + *nvars__] != zero) {
	goto L20045;
    }
    goto L20046;
L20048:
    goto L20035;
L20045:
    xlamda = zero;
    npr009 = 3;
    npr009_fmt = fmt_20050;
    goto L30009;
L20050:
L20034:

L20035:
    npr011 = 0;
    npr011_fmt = fmt_20051;
    goto L30011;
L20051:
    if (! (feas && ! unbnd)) {
	goto L20052;
    }
    *info = 1;
    goto L20053;
L20052:
    if (! (! feas && ! unbnd)) {
	goto L10001;
    }
    nerr = 1;
    xermsg_("SLATEC", "DPLPMN", "IN DSPLP, THE PROBLEM APPEARS TO BE INFEASI"
	    "BLE", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)46);
    *info = -nerr;
    goto L20053;
L10001:
    if (! (feas && unbnd)) {
	goto L10002;
    }
    nerr = 2;
    xermsg_("SLATEC", "DPLPMN", "IN DSPLP, THE PROBLEM APPEARS TO HAVE NO FI"
	    "NITE SOLUTION.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)57);
    *info = -nerr;
    goto L20053;
L10002:
    if (! (! feas && unbnd)) {
	goto L10003;
    }
    nerr = 3;
    xermsg_("SLATEC", "DPLPMN", "IN DSPLP, THE PROBLEM APPEARS TO BE INFEASI"
	    "BLE AND TO HAVE NO FINITE SOLN.", &nerr, &iopt, (ftnlen)6, (
	    ftnlen)6, (ftnlen)74);
    *info = -nerr;
L10003:
L20053:

    if (! (*info == -1 || *info == -3)) {
	goto L20055;
    }
    size = dasum_(nvars__, &primal[1], &c__1) * anorm;
    size /= dasum_(nvars__, &csc[1], &c__1);
    size += dasum_(mrelas, &primal[*nvars__ + 1], &c__1);
    i__ = 1;
    n20058 = *nvars__ + *mrelas;
    goto L20059;
L20058:
    ++i__;
L20059:
    if (n20058 - i__ < 0) {
	goto L20060;
    }
    nx0066 = ind[i__];
    if (nx0066 < 1 || nx0066 > 4) {
	goto L20066;
    }
    switch (nx0066) {
	case 1:  goto L20062;
	case 2:  goto L20063;
	case 3:  goto L20064;
	case 4:  goto L20065;
    }
L20062:
    if (! (size + (d__1 = primal[i__] - bl[i__], abs(d__1)) * factor == size))
	     {
	goto L20068;
    }
    goto L20058;
L20068:
    if (! (primal[i__] > bl[i__])) {
	goto L10004;
    }
    goto L20058;
L10004:
    ind[i__] = -4;
    goto L20067;
L20063:
    if (! (size + (d__1 = primal[i__] - bu[i__], abs(d__1)) * factor == size))
	     {
	goto L20071;
    }
    goto L20058;
L20071:
    if (! (primal[i__] < bu[i__])) {
	goto L10005;
    }
    goto L20058;
L10005:
    ind[i__] = -4;
    goto L20067;
L20064:
    if (! (size + (d__1 = primal[i__] - bl[i__], abs(d__1)) * factor == size))
	     {
	goto L20074;
    }
    goto L20058;
L20074:
    if (! (primal[i__] < bl[i__])) {
	goto L10006;
    }
    ind[i__] = -4;
    goto L20075;
L10006:
    if (! (size + (d__1 = primal[i__] - bu[i__], abs(d__1)) * factor == size))
	     {
	goto L10007;
    }
    goto L20058;
L10007:
    if (! (primal[i__] > bu[i__])) {
	goto L10008;
    }
    ind[i__] = -4;
    goto L20075;
L10008:
    goto L20058;
L20075:
    goto L20067;
L20065:
    goto L20058;
L20066:
L20067:
    goto L20058;
L20060:
L20055:

    if (! (*info == -2 || *info == -3)) {
	goto L20077;
    }
    j = 1;
    n20080 = *nvars__;
    goto L20081;
L20080:
    ++j;
L20081:
    if (n20080 - j < 0) {
	goto L20082;
    }
    if (! (ibb[j] == 0)) {
	goto L20084;
    }
    nx0091 = ind[j];
    if (nx0091 < 1 || nx0091 > 4) {
	goto L20091;
    }
    switch (nx0091) {
	case 1:  goto L20087;
	case 2:  goto L20088;
	case 3:  goto L20089;
	case 4:  goto L20090;
    }
L20087:
    bu[j] = bl[j];
    ind[j] = -3;
    goto L20092;
L20088:
    bl[j] = bu[j];
    ind[j] = -3;
    goto L20092;
L20089:
    goto L20080;
L20090:
    bl[j] = zero;
    bu[j] = zero;
    ind[j] = -3;
L20091:
L20092:
L20084:
    goto L20080;
L20082:
L20077:
/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
    if (! (*kprint >= 1)) {
	goto L20093;
    }
    npr012 = 0;
    npr012_fmt = fmt_20096;
    goto L30012;
L20096:
L20093:
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */
    goto L30001;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (COMPUTE RIGHT HAND SIDE) */
L30010:
    rhs[1] = zero;
    dcopy_(mrelas, &rhs[1], &c__0, &rhs[1], &c__1);
    j = 1;
    n20098 = *nvars__ + *mrelas;
    goto L20099;
L20098:
    ++j;
L20099:
    if (n20098 - j < 0) {
	goto L20100;
    }
    nx0106 = ind[j];
    if (nx0106 < 1 || nx0106 > 4) {
	goto L20106;
    }
    switch (nx0106) {
	case 1:  goto L20102;
	case 2:  goto L20103;
	case 3:  goto L20104;
	case 4:  goto L20105;
    }
L20102:
    scalr = -bl[j];
    goto L20107;
L20103:
    scalr = -bu[j];
    goto L20107;
L20104:
    scalr = -bl[j];
    goto L20107;
L20105:
    scalr = zero;
L20106:
L20107:
    if (! (scalr != zero)) {
	goto L20108;
    }
    if (! (j <= *nvars__)) {
	goto L20111;
    }
    i__ = 0;
L20114:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (i__ <= 0)) {
	goto L20116;
    }
    goto L20115;
L20116:
    rhs[i__] += aij * scalr;
    goto L20114;
L20115:
    goto L20112;
L20111:
    rhs[j - *nvars__] -= scalr;
L20112:
L20108:
    goto L20098;
L20100:
    j = 1;
    n20119 = *nvars__ + *mrelas;
    goto L20120;
L20119:
    ++j;
L20120:
    if (n20119 - j < 0) {
	goto L20121;
    }
    scalr = zero;
    if (ind[j] == 3 && ibb[j] % 2 == 0) {
	scalr = bu[j] - bl[j];
    }
    if (! (scalr != zero)) {
	goto L20123;
    }
    if (! (j <= *nvars__)) {
	goto L20126;
    }
    i__ = 0;
L20129:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (i__ <= 0)) {
	goto L20131;
    }
    goto L20130;
L20131:
    rhs[i__] -= aij * scalr;
    goto L20129;
L20130:
    goto L20127;
L20126:
    rhs[j - *nvars__] += scalr;
L20127:
L20123:
    goto L20119;
L20121:
    switch (npr010) {
	case 0: goto L20030;
	case 1: goto L20038;
    }
/*     PROCEDURE (PERFORM SIMPLEX STEPS) */
L30009:
    npr013 = 0;
    npr013_fmt = fmt_20134;
    goto L30013;
L20134:
    npr014 = 0;
    npr014_fmt = fmt_20135;
    goto L30014;
L20135:
    if (! (*kprint > 2)) {
	goto L20136;
    }
    dvout_(mrelas, &duals[1], "(' BASIC (INTERNAL) DUAL SOLN.')", idg, (
	    ftnlen)32);
    i__1 = *nvars__ + *mrelas;
    dvout_(&i__1, &rz[1], "(' REDUCED COSTS')", idg, (ftnlen)18);
L20136:
L20139:
    npr015 = 0;
    npr015_fmt = fmt_20141;
    goto L30015;
L20141:
    if (found) {
	goto L20142;
    }
    goto L30016;
L20145:
L20142:
    if (! found) {
	goto L20146;
    }
    if (*kprint >= 3) {
	dvout_(mrelas, &ww[1], "(' SEARCH DIRECTION')", idg, (ftnlen)21);
    }
    goto L30017;
L20149:
    if (! finite) {
	goto L20150;
    }
    goto L30018;
L20153:
    npr005 = 1;
    npr005_fmt = fmt_20154;
    goto L30005;
L20154:
    goto L20151;
L20150:
    unbnd = TRUE_;
    ibb[ibasis[ienter]] = 0;
L20151:
    goto L20147;
L20146:
    goto L20140;
L20147:
    ++itlp;
    goto L30019;
L20155:
    goto L20139;
L20140:
    switch (npr009) {
	case 0: goto L20029;
	case 1: goto L20037;
	case 2: goto L20044;
	case 3: goto L20050;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (RETRIEVE SAVED DATA FROM FILE ISAVE) */
L30002:
    lpr = *nvars__ + 4;
    al__1.aerr = 0;
    al__1.aunit = *isave;
    f_rew(&al__1);
    io___74.ciunit = *isave;
    s_rsue(&io___74);
    i__1 = lpr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_uio(&c__1, (char *)&amat[i__], (ftnlen)sizeof(doublereal));
    }
    i__2 = lpr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_uio(&c__1, (char *)&imat[i__], (ftnlen)sizeof(integer));
    }
    e_rsue();
    key = 2;
    ipage = 1;
    goto L20157;
L20156:
    if (np < 0) {
	goto L20158;
    }
L20157:
    lpr1 = lpr + 1;
    io___79.ciunit = *isave;
    s_rsue(&io___79);
    i__1 = *lmx;
    for (i__ = lpr1; i__ <= i__1; ++i__) {
	do_uio(&c__1, (char *)&amat[i__], (ftnlen)sizeof(doublereal));
    }
    i__2 = *lmx;
    for (i__ = lpr1; i__ <= i__2; ++i__) {
	do_uio(&c__1, (char *)&imat[i__], (ftnlen)sizeof(integer));
    }
    e_rsue();
    dprwpg_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    np = imat[*lmx - 1];
    ++ipage;
    goto L20156;
L20158:
    nparm = *nvars__ + *mrelas;
    io___81.ciunit = *isave;
    s_rsue(&io___81);
    i__1 = nparm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_uio(&c__1, (char *)&ibasis[i__], (ftnlen)sizeof(integer));
    }
    e_rsue();
    al__1.aerr = 0;
    al__1.aunit = *isave;
    f_rew(&al__1);
    goto L20006;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (SAVE DATA ON FILE ISAVE) */

/*     SOME PAGES MAY NOT BE WRITTEN YET. */
L30020:
    if (! (amat[*lmx] == one)) {
	goto L20159;
    }
    amat[*lmx] = zero;
    key = 2;
    ipage = (i__1 = imat[*lmx - 1], abs(i__1));
    dprwpg_(&key, &ipage, &lpg, &amat[1], &imat[1]);

/*     FORCE PAGE FILE TO BE OPENED ON RESTARTS. */
L20159:
    key = (integer) amat[4];
    amat[4] = zero;
    lpr = *nvars__ + 4;
    io___82.ciunit = *isave;
    s_wsue(&io___82);
    i__1 = lpr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_uio(&c__1, (char *)&amat[i__], (ftnlen)sizeof(doublereal));
    }
    i__2 = lpr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_uio(&c__1, (char *)&imat[i__], (ftnlen)sizeof(integer));
    }
    e_wsue();
    amat[4] = (doublereal) key;
    ipage = 1;
    key = 1;
    goto L20163;
L20162:
    if (np < 0) {
	goto L20164;
    }
L20163:
    dprwpg_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    lpr1 = lpr + 1;
    io___83.ciunit = *isave;
    s_wsue(&io___83);
    i__1 = *lmx;
    for (i__ = lpr1; i__ <= i__1; ++i__) {
	do_uio(&c__1, (char *)&amat[i__], (ftnlen)sizeof(doublereal));
    }
    i__2 = *lmx;
    for (i__ = lpr1; i__ <= i__2; ++i__) {
	do_uio(&c__1, (char *)&imat[i__], (ftnlen)sizeof(integer));
    }
    e_wsue();
    np = imat[*lmx - 1];
    ++ipage;
    goto L20162;
L20164:
    nparm = *nvars__ + *mrelas;
    io___84.ciunit = *isave;
    s_wsue(&io___84);
    i__1 = nparm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_uio(&c__1, (char *)&ibasis[i__], (ftnlen)sizeof(integer));
    }
    e_wsue();
    al__2.aerr = 0;
    al__2.aunit = *isave;
    f_end(&al__2);

/*     CLOSE FILE, IPAGEF, WHERE PAGES ARE STORED. THIS IS NEEDED SO THAT */
/*     THE PAGES MAY BE RESTORED AT A CONTINUATION OF DSPLP(). */
    goto L20317;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (DECOMPOSE BASIS MATRIX) */
/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
L30004:
    if (! (*kprint >= 2)) {
	goto L20165;
    }
    ivout_(mrelas, &ibasis[1], "(' SUBSCRIPTS OF BASIC VARIABLES DURING REDE"
	    "COMPOSITION')", idg, (ftnlen)57);
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */

/*     SET RELATIVE PIVOTING FACTOR FOR USE IN LA05 () PACKAGE. */
L20165:
    uu = .1f;
    dplpdm_(mrelas, nvars__, lmx, lbm, &nredc, info, &iopt, &ibasis[1], &imat[
	    1], &ibrc[ibrc_offset], &ipr[1], &iwr[1], &ind[1], &ibb[1], &
	    anorm, eps, &uu, &gg, &amat[1], &basmat[1], &csc[1], &wr[1], &
	    singlr, &redbas);
    if (! (*info < 0)) {
	goto L20168;
    }
    goto L30001;
L20168:
    switch (npr004) {
	case 0: goto L20013;
	case 1: goto L20204;
	case 2: goto L20242;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (CLASSIFY VARIABLES) */

/*     DEFINE THE CLASSIFICATION OF THE BASIC VARIABLES */
/*     -1 VIOLATES LOWER BOUND, 0 FEASIBLE, +1 VIOLATES UPPER BOUND. */
/*     (THIS INFO IS STORED IN PRIMAL(NVARS+1)-PRIMAL(NVARS+MRELAS)) */
/*     TRANSLATE VARIABLE TO ITS UPPER BOUND, IF .GT. UPPER BOUND */
L30007:
    primal[*nvars__ + 1] = zero;
    dcopy_(mrelas, &primal[*nvars__ + 1], &c__0, &primal[*nvars__ + 1], &c__1)
	    ;
    i__ = 1;
    n20172 = *mrelas;
    goto L20173;
L20172:
    ++i__;
L20173:
    if (n20172 - i__ < 0) {
	goto L20174;
    }
    j = ibasis[i__];
    if (! (ind[j] != 4)) {
	goto L20176;
    }
    if (! (rprim[i__] < zero)) {
	goto L20179;
    }
    primal[i__ + *nvars__] = -one;
    goto L20180;
L20179:
    if (! (ind[j] == 3)) {
	goto L10009;
    }
    upbnd = bu[j] - bl[j];
    if (j <= *nvars__) {
	upbnd /= csc[j];
    }
    if (! (rprim[i__] > upbnd)) {
	goto L20182;
    }
    rprim[i__] -= upbnd;
    if (! (j <= *nvars__)) {
	goto L20185;
    }
    k = 0;
L20188:
    dpnnzr_(&k, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (k <= 0)) {
	goto L20190;
    }
    goto L20189;
L20190:
    rhs[k] -= upbnd * aij * csc[j];
    goto L20188;
L20189:
    goto L20186;
L20185:
    rhs[j - *nvars__] += upbnd;
L20186:
    primal[i__ + *nvars__] = one;
L20182:
L10009:
L20180:
L20176:
    goto L20172;
L20174:
    switch (npr007) {
	case 0: goto L20020;
	case 1: goto L20036;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (COMPUTE ERROR IN DUAL AND PRIMAL SYSTEMS) */
L30005:
    ntries = 1;
    goto L20195;
L20194:
    ++ntries;
L20195:
    if (2 - ntries < 0) {
	goto L20196;
    }
    dplpce_(mrelas, nvars__, lmx, lbm, &itlp, itbrc, &ibasis[1], &imat[1], &
	    ibrc[ibrc_offset], &ipr[1], &iwr[1], &ind[1], &ibb[1], &erdnrm, 
	    eps, tune, &gg, &amat[1], &basmat[1], &csc[1], &wr[1], &ww[1], &
	    primal[1], &erd[1], &erp[1], &singlr, &redbas);
    if (singlr) {
	goto L20198;
    }
/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
    if (! (*kprint >= 3)) {
	goto L20201;
    }
    dvout_(mrelas, &erp[1], "(' EST. ERROR IN PRIMAL COMPS.')", idg, (ftnlen)
	    32);
    dvout_(mrelas, &erd[1], "(' EST. ERROR IN DUAL COMPS.')", idg, (ftnlen)30)
	    ;
L20201:
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */
    goto L20193;
L20198:
    if (ntries == 2) {
	goto L20197;
    }
    npr004 = 1;
    npr004_fmt = fmt_20204;
    goto L30004;
L20204:
    goto L20194;
L20196:
L20197:
    nerr = 26;
    xermsg_("SLATEC", "DPLPMN", "IN DSPLP, MOVED TO A SINGULAR POINT. THIS S"
	    "HOULD NOT HAPPEN.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)
	    60);
    *info = -nerr;
    goto L30001;
L20193:
    switch (npr005) {
	case 0: goto L20018;
	case 1: goto L20154;
	case 2: goto L20243;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (CHECK FEASIBILITY) */

/*     SEE IF NEARBY FEASIBLE POINT SATISFIES THE CONSTRAINT */
/*     EQUATIONS. */

/*     COPY RHS INTO WW(*), THEN UPDATE WW(*). */
L30008:
    dcopy_(mrelas, &rhs[1], &c__1, &ww[1], &c__1);
    j = 1;
    n20206 = *mrelas;
    goto L20207;
L20206:
    ++j;
L20207:
    if (n20206 - j < 0) {
	goto L20208;
    }
    ibas = ibasis[j];
    xval = rprim[j];

/*     ALL VARIABLES BOUNDED BELOW HAVE ZERO AS THAT BOUND. */
    if (ind[ibas] <= 3) {
	xval = max(zero,xval);
    }

/*     IF THE VARIABLE HAS AN UPPER BOUND, COMPUTE THAT BOUND. */
    if (! (ind[ibas] == 3)) {
	goto L20210;
    }
    upbnd = bu[ibas] - bl[ibas];
    if (ibas <= *nvars__) {
	upbnd /= csc[ibas];
    }
    xval = min(upbnd,xval);
L20210:

/*     SUBTRACT XVAL TIMES COLUMN VECTOR FROM RIGHT-HAND SIDE IN WW(*) */
    if (! (xval != zero)) {
	goto L20213;
    }
    if (! (ibas <= *nvars__)) {
	goto L20216;
    }
    i__ = 0;
L20219:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &ibas);
    if (! (i__ <= 0)) {
	goto L20221;
    }
    goto L20220;
L20221:
    ww[i__] -= xval * aij * csc[ibas];
    goto L20219;
L20220:
    goto L20217;
L20216:
    if (! (ind[ibas] == 2)) {
	goto L20224;
    }
    ww[ibas - *nvars__] -= xval;
    goto L20225;
L20224:
    ww[ibas - *nvars__] += xval;
L20225:
L20217:
L20213:
    goto L20206;

/*   COMPUTE NORM OF DIFFERENCE AND CHECK FOR FEASIBILITY. */
L20208:
    resnrm = dasum_(mrelas, &ww[1], &c__1);
    feas = resnrm <= *tolls * (rprnrm * anorm + rhsnrm);

/*     TRY AN ABSOLUTE ERROR TEST IF THE RELATIVE TEST FAILS. */
    if (! feas) {
	feas = resnrm <= *tolabs;
    }
    if (! feas) {
	goto L20227;
    }
    primal[*nvars__ + 1] = zero;
    dcopy_(mrelas, &primal[*nvars__ + 1], &c__0, &primal[*nvars__ + 1], &c__1)
	    ;
L20227:
    switch (npr008) {
	case 0: goto L20024;
	case 1: goto L20032;
	case 2: goto L20040;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (INITIALIZE REDUCED COSTS AND STEEPEST EDGE WEIGHTS) */
L30014:
    dpincw_(mrelas, nvars__, lmx, lbm, npp, &jstrt, &ibasis[1], &imat[1], &
	    ibrc[ibrc_offset], &ipr[1], &iwr[1], &ind[1], &ibb[1], costsc, &
	    gg, &erdnrm, &dulnrm, &amat[1], &basmat[1], &csc[1], &wr[1], &ww[
	    1], &rz[1], &rg[1], &costs[1], &colnrm[1], &duals[1], stpedg);

    switch (npr014) {
	case 0: goto L20135;
	case 1: goto L20246;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (CHECK AND RETURN WITH EXCESS ITERATIONS) */
L30019:
    if (! (itlp > *mxitlp)) {
	goto L20230;
    }
    nerr = 25;
    npr011 = 1;
    npr011_fmt = fmt_20233;
    goto L30011;
/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
L20233:
    if (! (*kprint >= 1)) {
	goto L20234;
    }
    npr012 = 1;
    npr012_fmt = fmt_20237;
    goto L30012;
L20237:
L20234:
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */
    idum[0] = 0;
    if (*savedt) {
	idum[0] = *isave;
    }
    s_wsfi(&io___100);
    do_fio(&c__1, (char *)&(*mxitlp), (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___102);
    do_fio(&c__1, (char *)&idum[0], (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 27, a__1[0] = "IN DSPLP, MAX ITERATIONS = ";
    i__3[1] = 8, a__1[1] = xern1;
    i__3[2] = 46, a__1[2] = " TAKEN.  UP-TO-DATE RESULTS SAVED ON FILE NO. ";
    i__3[3] = 8, a__1[3] = xern2;
    i__3[4] = 29, a__1[4] = ".   IF FILE NO. = 0, NO SAVE.";
    s_cat(ch__1, a__1, i__3, &c__5, (ftnlen)118);
    xermsg_("SLATEC", "DPLPMN", ch__1, &nerr, &iopt, (ftnlen)6, (ftnlen)6, (
	    ftnlen)118);
    *info = -nerr;
    goto L30001;
L20230:
    goto L20155;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (REDECOMPOSE BASIS MATRIX AND TRY AGAIN) */
L30016:
    if (redbas) {
	goto L20239;
    }
    npr004 = 2;
    npr004_fmt = fmt_20242;
    goto L30004;
L20242:
    npr005 = 2;
    npr005_fmt = fmt_20243;
    goto L30005;
L20243:
    npr006 = 3;
    npr006_fmt = fmt_20244;
    goto L30006;
L20244:
    npr013 = 1;
    npr013_fmt = fmt_20245;
    goto L30013;
L20245:
    npr014 = 1;
    npr014_fmt = fmt_20246;
    goto L30014;
L20246:

/*     ERASE NON-CYCLING MARKERS NEAR COMPLETION. */
L20239:
    i__ = *mrelas + 1;
    n20247 = *mrelas + *nvars__;
    goto L20248;
L20247:
    ++i__;
L20248:
    if (n20247 - i__ < 0) {
	goto L20249;
    }
    ibasis[i__] = (i__1 = ibasis[i__], abs(i__1));
    goto L20247;
L20249:
    npr015 = 1;
    npr015_fmt = fmt_20251;
    goto L30015;
L20251:
    goto L20145;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (COMPUTE NEW PRIMAL) */

/*     COPY RHS INTO WW(*), SOLVE SYSTEM. */
L30006:
    dcopy_(mrelas, &rhs[1], &c__1, &ww[1], &c__1);
    trans = FALSE_;
    la05bd_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], &gg, &ww[1], &trans);
    dcopy_(mrelas, &ww[1], &c__1, &rprim[1], &c__1);
    rprnrm = dasum_(mrelas, &rprim[1], &c__1);
    switch (npr006) {
	case 0: goto L20019;
	case 1: goto L20031;
	case 2: goto L20039;
	case 3: goto L20244;
	case 4: goto L20275;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (COMPUTE NEW DUALS) */

/*     SOLVE FOR DUAL VARIABLES. FIRST COPY COSTS INTO DUALS(*). */
L30013:
    i__ = 1;
    n20252 = *mrelas;
    goto L20253;
L20252:
    ++i__;
L20253:
    if (n20252 - i__ < 0) {
	goto L20254;
    }
    j = ibasis[i__];
    if (! (j <= *nvars__)) {
	goto L20256;
    }
    duals[i__] = *costsc * costs[j] * csc[j] + xlamda * primal[i__ + *nvars__]
	    ;
    goto L20257;
L20256:
    duals[i__] = xlamda * primal[i__ + *nvars__];
L20257:
    goto L20252;

L20254:
    trans = TRUE_;
    la05bd_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], &gg, &duals[1], &trans);
    dulnrm = dasum_(mrelas, &duals[1], &c__1);
    switch (npr013) {
	case 0: goto L20134;
	case 1: goto L20245;
	case 2: goto L20267;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (FIND VARIABLE TO ENTER BASIS AND GET SEARCH DIRECTION) */
L30015:
    dplpfe_(mrelas, nvars__, lmx, lbm, &ienter, &ibasis[1], &imat[1], &ibrc[
	    ibrc_offset], &ipr[1], &iwr[1], &ind[1], &ibb[1], &erdnrm, eps, &
	    gg, &dulnrm, &dirnrm, &amat[1], &basmat[1], &csc[1], &wr[1], &ww[
	    1], &bl[1], &bu[1], &rz[1], &rg[1], &colnrm[1], &duals[1], &found)
	    ;
    switch (npr015) {
	case 0: goto L20141;
	case 1: goto L20251;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (CHOOSE VARIABLE TO LEAVE BASIS) */
L30017:
    dplpfl_(mrelas, nvars__, &ienter, &ileave, &ibasis[1], &ind[1], &ibb[1], &
	    theta, &dirnrm, &rprnrm, &csc[1], &ww[1], &bl[1], &bu[1], &erp[1],
	     &rprim[1], &primal[1], &finite, &zerolv);
    goto L20149;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (MAKE MOVE AND UPDATE) */
L30018:
    dplpmu_(mrelas, nvars__, lmx, lbm, &nredc, info, &ienter, &ileave, &iopt, 
	    npp, &jstrt, &ibasis[1], &imat[1], &ibrc[ibrc_offset], &ipr[1], &
	    iwr[1], &ind[1], &ibb[1], &anorm, eps, &uu, &gg, &rprnrm, &erdnrm,
	     &dulnrm, &theta, costsc, &xlamda, &rhsnrm, &amat[1], &basmat[1], 
	    &csc[1], &wr[1], &rprim[1], &ww[1], &bu[1], &bl[1], &rhs[1], &erd[
	    1], &erp[1], &rz[1], &rg[1], &colnrm[1], &costs[1], &primal[1], &
	    duals[1], &singlr, &redbas, &zerolv, stpedg);
    if (! (*info == -26)) {
	goto L20259;
    }
    goto L30001;
/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
L20259:
    if (! (*kprint >= 2)) {
	goto L20263;
    }
    goto L30021;
L20266:
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */
L20263:
    goto L20153;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE(RESCALE AND REARRANGE VARIABLES) */

/*     RESCALE THE DUAL VARIABLES. */
L30011:
    npr013 = 2;
    npr013_fmt = fmt_20267;
    goto L30013;
L20267:
    if (! (*costsc != zero)) {
	goto L20268;
    }
    i__ = 1;
    n20271 = *mrelas;
    goto L20272;
L20271:
    ++i__;
L20272:
    if (n20271 - i__ < 0) {
	goto L20273;
    }
    duals[i__] /= *costsc;
    goto L20271;
L20273:
L20268:
    npr006 = 4;
    npr006_fmt = fmt_20275;
    goto L30006;

/*     REAPPLY COLUMN SCALING TO PRIMAL. */
L20275:
    i__ = 1;
    n20276 = *mrelas;
    goto L20277;
L20276:
    ++i__;
L20277:
    if (n20276 - i__ < 0) {
	goto L20278;
    }
    j = ibasis[i__];
    if (! (j <= *nvars__)) {
	goto L20280;
    }
    scalr = csc[j];
    if (ind[j] == 2) {
	scalr = -scalr;
    }
    rprim[i__] *= scalr;
L20280:
    goto L20276;

/*     REPLACE TRANSLATED BASIC VARIABLES INTO ARRAY PRIMAL(*) */
L20278:
    primal[1] = zero;
    i__1 = *nvars__ + *mrelas;
    dcopy_(&i__1, &primal[1], &c__0, &primal[1], &c__1);
    j = 1;
    n20283 = *nvars__ + *mrelas;
    goto L20284;
L20283:
    ++j;
L20284:
    if (n20283 - j < 0) {
	goto L20285;
    }
    ibas = (i__1 = ibasis[j], abs(i__1));
    xval = zero;
    if (j <= *mrelas) {
	xval = rprim[j];
    }
    if (ind[ibas] == 1) {
	xval += bl[ibas];
    }
    if (ind[ibas] == 2) {
	xval = bu[ibas] - xval;
    }
    if (! (ind[ibas] == 3)) {
	goto L20287;
    }
    if (ibb[ibas] % 2 == 0) {
	xval = bu[ibas] - bl[ibas] - xval;
    }
    xval += bl[ibas];
L20287:
    primal[ibas] = xval;
    goto L20283;

/*     COMPUTE DUALS FOR INDEPENDENT VARIABLES WITH BOUNDS. */
/*     OTHER ENTRIES ARE ZERO. */
L20285:
    j = 1;
    n20290 = *nvars__;
    goto L20291;
L20290:
    ++j;
L20291:
    if (n20290 - j < 0) {
	goto L20292;
    }
    rzj = zero;
    if (! ((doublereal) ibb[j] > zero && ind[j] != 4)) {
	goto L20294;
    }
    rzj = costs[j];
    i__ = 0;
L20297:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (i__ <= 0)) {
	goto L20299;
    }
    goto L20298;
L20299:
    rzj -= aij * duals[i__];
    goto L20297;
L20298:
L20294:
    duals[*mrelas + j] = rzj;
    goto L20290;
L20292:
    switch (npr011) {
	case 0: goto L20051;
	case 1: goto L20233;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* ++  CODE FOR OUTPUT=YES IS ACTIVE */
/*     PROCEDURE (PRINT PROLOGUE) */
L30003:
    idum[0] = *mrelas;
    ivout_(&c__1, idum, "('1NUM. OF DEPENDENT VARS., MRELAS')", idg, (ftnlen)
	    36);
    idum[0] = *nvars__;
    ivout_(&c__1, idum, "(' NUM. OF INDEPENDENT VARS., NVARS')", idg, (ftnlen)
	    37);
    ivout_(&c__1, idum, "(' DIMENSION OF COSTS(*)=')", idg, (ftnlen)27);
    idum[0] = *nvars__ + *mrelas;
    ivout_(&c__1, idum, "(' DIMENSIONS OF BL(*),BU(*),IND(*)'        /' PRIM"
	    "AL(*),DUALS(*) =')", idg, (ftnlen)69);
    ivout_(&c__1, idum, "(' DIMENSION OF IBASIS(*)=')", idg, (ftnlen)28);
    idum[0] = *lprg + 1;
    ivout_(&c__1, idum, "(' DIMENSION OF PRGOPT(*)=')", idg, (ftnlen)28);
    ivout_(&c__0, idum, "(' 1-NVARS=INDEPENDENT VARIABLE INDICES.'/         "
	    "            ' (NVARS+1)-(NVARS+MRELAS)=DEPENDENT VARIABLE INDICE"
	    "S.'/        ' CONSTRAINT INDICATORS ARE 1-4 AND MEAN')", idg, (
	    ftnlen)169);
    ivout_(&c__0, idum, "(' 1=VARIABLE HAS ONLY LOWER BOUND.'/              "
	    "            ' 2=VARIABLE HAS ONLY UPPER BOUND.'/                "
	    "            ' 3=VARIABLE HAS BOTH BOUNDS.'/                     "
	    "            ' 4=VARIABLE HAS NO BOUNDS, IT IS FREE.')", idg, (
	    ftnlen)232);
    dvout_(nvars__, &costs[1], "(' ARRAY OF COSTS')", idg, (ftnlen)19);
    i__1 = *nvars__ + *mrelas;
    ivout_(&i__1, &ind[1], "(' CONSTRAINT INDICATORS')", idg, (ftnlen)26);
    i__1 = *nvars__ + *mrelas;
    dvout_(&i__1, &bl[1], "(' LOWER BOUNDS FOR VARIABLES  (IGNORE UNUSED ENT"
	    "RIES.)')", idg, (ftnlen)57);
    i__1 = *nvars__ + *mrelas;
    dvout_(&i__1, &bu[1], "(' UPPER BOUNDS FOR VARIABLES  (IGNORE UNUSED ENT"
	    "RIES.)')", idg, (ftnlen)57);
    if (! (*kprint >= 2)) {
	goto L20302;
    }
    ivout_(&c__0, idum, "('0NON-BASIC INDICES THAT ARE NEGATIVE SHOW VARIABL"
	    "ES'         ' EXCHANGED AT A ZERO'/' STEP LENGTH')", idg, (ftnlen)
	    101);
    ivout_(&c__0, idum, "(' WHEN COL. NO. LEAVING=COL. NO. ENTERING, THE ENT"
	    "ERING '     'VARIABLE MOVED'/' TO ITS BOUND.  IT REMAINS NON-BAS"
	    "IC.'/     ' WHEN COL. NO. OF BASIS EXCHANGED IS NEGATIVE, THE LE"
	    "AVING'/   ' VARIABLE IS AT ITS UPPER BOUND.')", idg, (ftnlen)224);
L20302:
    goto L20011;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (PRINT SUMMARY) */
L30012:
    idum[0] = *info;
    ivout_(&c__1, idum, "(' THE OUTPUT VALUE OF INFO IS')", idg, (ftnlen)32);
    if (! (*minprb)) {
	goto L20305;
    }
    ivout_(&c__0, idum, "(' THIS IS A MINIMIZATION PROBLEM.')", idg, (ftnlen)
	    36);
    goto L20306;
L20305:
    ivout_(&c__0, idum, "(' THIS IS A MAXIMIZATION PROBLEM.')", idg, (ftnlen)
	    36);
L20306:
    if (! (*stpedg)) {
	goto L20308;
    }
    ivout_(&c__0, idum, "(' STEEPEST EDGE PRICING WAS USED.')", idg, (ftnlen)
	    36);
    goto L20309;
L20308:
    ivout_(&c__0, idum, "(' MINIMUM REDUCED COST PRICING WAS USED.')", idg, (
	    ftnlen)43);
L20309:
    rdum[0] = ddot_(nvars__, &costs[1], &c__1, &primal[1], &c__1);
    dvout_(&c__1, rdum, "(' OUTPUT VALUE OF THE OBJECTIVE FUNCTION')", idg, (
	    ftnlen)43);
    i__1 = *nvars__ + *mrelas;
    dvout_(&i__1, &primal[1], "(' THE OUTPUT INDEPENDENT AND DEPENDENT VARIA"
	    "BLES')", idg, (ftnlen)51);
    i__1 = *mrelas + *nvars__;
    dvout_(&i__1, &duals[1], "(' THE OUTPUT DUAL VARIABLES')", idg, (ftnlen)
	    30);
    i__1 = *nvars__ + *mrelas;
    ivout_(&i__1, &ibasis[1], "(' VARIABLE INDICES IN POSITIONS 1-MRELAS ARE"
	    " BASIC.')", idg, (ftnlen)54);
    idum[0] = itlp;
    ivout_(&c__1, idum, "(' NO. OF ITERATIONS')", idg, (ftnlen)22);
    idum[0] = nredc;
    ivout_(&c__1, idum, "(' NO. OF FULL REDECOMPS')", idg, (ftnlen)26);
    switch (npr012) {
	case 0: goto L20096;
	case 1: goto L20237;
    }
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (PRINT ITERATION SUMMARY) */
L30021:
    idum[0] = itlp + 1;
    ivout_(&c__1, idum, "('0ITERATION NUMBER')", idg, (ftnlen)21);
    idum[0] = ibasis[abs(ileave)];
    ivout_(&c__1, idum, "(' INDEX OF VARIABLE ENTERING THE BASIS')", idg, (
	    ftnlen)41);
    idum[0] = ileave;
    ivout_(&c__1, idum, "(' COLUMN OF THE BASIS EXCHANGED')", idg, (ftnlen)34)
	    ;
    idum[0] = ibasis[ienter];
    ivout_(&c__1, idum, "(' INDEX OF VARIABLE LEAVING THE BASIS')", idg, (
	    ftnlen)40);
    rdum[0] = theta;
    dvout_(&c__1, rdum, "(' LENGTH OF THE EXCHANGE STEP')", idg, (ftnlen)32);
    if (! (*kprint >= 3)) {
	goto L20311;
    }
    dvout_(mrelas, &rprim[1], "(' BASIC (INTERNAL) PRIMAL SOLN.')", idg, (
	    ftnlen)34);
    i__1 = *nvars__ + *mrelas;
    ivout_(&i__1, &ibasis[1], "(' VARIABLE INDICES IN POSITIONS 1-MRELAS ARE"
	    " BASIC.')", idg, (ftnlen)54);
    i__1 = *nvars__ + *mrelas;
    ivout_(&i__1, &ibb[1], "(' IBB ARRAY')", idg, (ftnlen)14);
    dvout_(mrelas, &rhs[1], "(' TRANSLATED RHS')", idg, (ftnlen)19);
    dvout_(mrelas, &duals[1], "(' BASIC (INTERNAL) DUAL SOLN.')", idg, (
	    ftnlen)32);
L20311:
    goto L20266;
/* ++  CODE FOR OUTPUT=NO IS INACTIVE */
/* ++  END */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (RETURN TO USER) */
L30001:
    if (! (*savedt)) {
	goto L20314;
    }
    goto L30020;
L20317:
L20314:
    if (imat[*lmx - 1] != -1) {
	sclosm_(ipagef);
    }

/*     THIS TEST IS THERE ONLY TO AVOID DIAGNOSTICS ON SOME FORTRAN */
/*     COMPILERS. */
    return 0;
} /* dplpmn_ */

#undef mxitlp
#undef sizeup
#undef intopt
#undef kprint
#undef usrbas
#undef stpedg
#undef cstscp
#undef costsc
#undef contin
#undef minprb
#undef savedt
#undef colscp
#undef tolabs
#undef asmall
#undef ipagef
#undef tolls
#undef isave
#undef itbrc
#undef ropt
#undef lopt
#undef tune
#undef lprg
#undef abig
#undef npp
#undef eps
#undef idg



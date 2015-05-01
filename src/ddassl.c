/* ddassl.f -- translated by f2c (version 12.02.01).
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
static integer c_n998 = -998;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__3 = 3;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__17 = 17;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c_n999 = -999;

/* DECK DDASSL */
/* Subroutine */ int ddassl_(U_fp res, integer *neq, doublereal *t, 
	doublereal *y, doublereal *yprime, doublereal *tout, integer *info, 
	doublereal *rtol, doublereal *atol, integer *idid, doublereal *rwork, 
	integer *lrw, integer *iwork, integer *liw, doublereal *rpar, integer 
	*ipar, U_fp jac)
{
    /* System generated locals */
    address a__1[4], a__2[5], a__3[6], a__4[3], a__5[2];
    integer i__1, i__2[4], i__3[5], i__4[6], i__5[3], i__6[2];
    doublereal d__1, d__2;
    char ch__1[118], ch__2[81], ch__3[128], ch__4[62], ch__5[110], ch__6[121],
	     ch__7[90], ch__8[132], ch__9[126], ch__10[85], ch__11[98], 
	    ch__12[21], ch__13[30], ch__14[61], ch__15[71], ch__16[32], 
	    ch__17[51], ch__18[78], ch__19[66], ch__20[49], ch__21[27];

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal r__;
    static integer le;
    static doublereal ho, rh, tn;
    static integer lpd, lwm, lwt;
    static logical done;
    static integer lphi;
    static doublereal hmax, hmin;
    static char xern1[8], xern2[8], xern3[16], xern4[16];
    static integer mband, lenpd;
    static doublereal atoli;
    static integer msave, itemp, leniw, nzflg, ntemp, lenrw;
    static doublereal tdist;
    static integer mxord;
    static doublereal rtoli;
    extern doublereal d1mach_(integer *);
    static doublereal tnext, tstop;
    extern /* Subroutine */ int ddaini_(doublereal *, doublereal *, 
	    doublereal *, integer *, U_fp, U_fp, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *);
    extern doublereal ddanrm_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    extern /* Subroutine */ int ddatrp_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *), ddastp_(doublereal *, doublereal *, doublereal *, 
	    integer *, U_fp, U_fp, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *), ddawts_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), xermsg_(char *, char *, char *, integer 
	    *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal uround, ypnorm;

    /* Fortran I/O blocks */
    static icilist io___10 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___34 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___35 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___36 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___37 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___39 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___40 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___41 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___42 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___43 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___44 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___45 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___46 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___47 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___48 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___49 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___50 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___51 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___52 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___53 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___54 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___56 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___57 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___58 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___59 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___60 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___61 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___62 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___63 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___64 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___65 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___66 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___67 = { 0, xern4, 0, "(1P,D15.6)", 16, 1 };
    static icilist io___68 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___69 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___70 = { 0, xern3, 0, "(1P,D15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  DDASSL */
/* ***PURPOSE  This code solves a system of differential/algebraic */
/*            equations of the form G(T,Y,YPRIME) = 0. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***CATEGORY  I1A2 */
/* ***TYPE      DOUBLE PRECISION (SDASSL-S, DDASSL-D) */
/* ***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DASSL, */
/*             DIFFERENTIAL/ALGEBRAIC, IMPLICIT DIFFERENTIAL SYSTEMS */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/*             Computing and Mathematics Research Division */
/*             Lawrence Livermore National Laboratory */
/*             L - 316, P.O. Box 808, */
/*             Livermore, CA.    94550 */
/* ***DESCRIPTION */

/* *Usage: */

/*      EXTERNAL RES, JAC */
/*      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR */
/*      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL, */
/*     *   RWORK(LRW), RPAR */

/*      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, */
/*     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC) */


/* *Arguments: */
/*  (In the following, all real arrays should be type DOUBLE PRECISION.) */

/*  RES:EXT     This is a subroutine which you provide to define the */
/*              differential/algebraic system. */

/*  NEQ:IN      This is the number of equations to be solved. */

/*  T:INOUT     This is the current value of the independent variable. */

/*  Y(*):INOUT  This array contains the solution components at T. */

/*  YPRIME(*):INOUT  This array contains the derivatives of the solution */
/*              components at T. */

/*  TOUT:IN     This is a point at which a solution is desired. */

/*  INFO(N):IN  The basic task of the code is to solve the system from T */
/*              to TOUT and return an answer at TOUT.  INFO is an integer */
/*              array which is used to communicate exactly how you want */
/*              this task to be carried out.  (See below for details.) */
/*              N must be greater than or equal to 15. */

/*  RTOL,ATOL:INOUT  These quantities represent relative and absolute */
/*              error tolerances which you provide to indicate how */
/*              accurately you wish the solution to be computed.  You */
/*              may choose them to be both scalars or else both vectors. */
/*              Caution:  In Fortran 77, a scalar is not the same as an */
/*                        array of length 1.  Some compilers may object */
/*                        to using scalars for RTOL,ATOL. */

/*  IDID:OUT    This scalar quantity is an indicator reporting what the */
/*              code did.  You must monitor this integer variable to */
/*              decide  what action to take next. */

/*  RWORK:WORK  A real work array of length LRW which provides the */
/*              code with needed storage space. */

/*  LRW:IN      The length of RWORK.  (See below for required length.) */

/*  IWORK:WORK  An integer work array of length LIW which provides the */
/*              code with needed storage space. */

/*  LIW:IN      The length of IWORK.  (See below for required length.) */

/*  RPAR,IPAR:IN  These are real and integer parameter arrays which */
/*              you can use for communication between your calling */
/*              program and the RES subroutine (and the JAC subroutine) */

/*  JAC:EXT     This is the name of a subroutine which you may choose */
/*              to provide for defining a matrix of partial derivatives */
/*              described below. */

/*  Quantities which may be altered by DDASSL are: */
/*     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL, */
/*     IDID, RWORK(*) AND IWORK(*) */

/* *Description */

/*  Subroutine DDASSL uses the backward differentiation formulas of */
/*  orders one through five to solve a system of the above form for Y and */
/*  YPRIME.  Values for Y and YPRIME at the initial time must be given as */
/*  input.  These values must be consistent, (that is, if T,Y,YPRIME are */
/*  the given initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The */
/*  subroutine solves the system from T to TOUT.  It is easy to continue */
/*  the solution to get results at additional TOUT.  This is the interval */
/*  mode of operation.  Intermediate results can also be obtained easily */
/*  by using the intermediate-output capability. */

/*  The following detailed description is divided into subsections: */
/*    1. Input required for the first call to DDASSL. */
/*    2. Output after any return from DDASSL. */
/*    3. What to do to continue the integration. */
/*    4. Error messages. */


/*  -------- INPUT -- WHAT TO DO ON THE FIRST CALL TO DDASSL ------------ */

/*  The first call of the code is defined to be the start of each new */
/*  problem. Read through the descriptions of all the following items, */
/*  provide sufficient storage space for designated arrays, set */
/*  appropriate variables for the initialization of the problem, and */
/*  give information about how you want the problem to be solved. */


/*  RES -- Provide a subroutine of the form */
/*             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR) */
/*         to define the system of differential/algebraic */
/*         equations which is to be solved. For the given values */
/*         of T,Y and YPRIME, the subroutine should */
/*         return the residual of the differential/algebraic */
/*         system */
/*             DELTA = G(T,Y,YPRIME) */
/*         (DELTA(*) is a vector of length NEQ which is */
/*         output for RES.) */

/*         Subroutine RES must not alter T,Y or YPRIME. */
/*         You must declare the name RES in an external */
/*         statement in your program that calls DDASSL. */
/*         You must dimension Y,YPRIME and DELTA in RES. */

/*         IRES is an integer flag which is always equal to */
/*         zero on input. Subroutine RES should alter IRES */
/*         only if it encounters an illegal value of Y or */
/*         a stop condition. Set IRES = -1 if an input value */
/*         is illegal, and DDASSL will try to solve the problem */
/*         without getting IRES = -1. If IRES = -2, DDASSL */
/*         will return control to the calling program */
/*         with IDID = -11. */

/*         RPAR and IPAR are real and integer parameter arrays which */
/*         you can use for communication between your calling program */
/*         and subroutine RES. They are not altered by DDASSL. If you */
/*         do not need RPAR or IPAR, ignore these parameters by treat- */
/*         ing them as dummy arguments. If you do choose to use them, */
/*         dimension them in your calling program and in RES as arrays */
/*         of appropriate length. */

/*  NEQ -- Set it to the number of differential equations. */
/*         (NEQ .GE. 1) */

/*  T -- Set it to the initial point of the integration. */
/*         T must be defined as a variable. */

/*  Y(*) -- Set this vector to the initial values of the NEQ solution */
/*         components at the initial point. You must dimension Y of */
/*         length at least NEQ in your calling program. */

/*  YPRIME(*) -- Set this vector to the initial values of the NEQ */
/*         first derivatives of the solution components at the initial */
/*         point.  You must dimension YPRIME at least NEQ in your */
/*         calling program. If you do not know initial values of some */
/*         of the solution components, see the explanation of INFO(11). */

/*  TOUT -- Set it to the first point at which a solution */
/*         is desired. You can not take TOUT = T. */
/*         integration either forward in T (TOUT .GT. T) or */
/*         backward in T (TOUT .LT. T) is permitted. */

/*         The code advances the solution from T to TOUT using */
/*         step sizes which are automatically selected so as to */
/*         achieve the desired accuracy. If you wish, the code will */
/*         return with the solution and its derivative at */
/*         intermediate steps (intermediate-output mode) so that */
/*         you can monitor them, but you still must provide TOUT in */
/*         accord with the basic aim of the code. */

/*         The first step taken by the code is a critical one */
/*         because it must reflect how fast the solution changes near */
/*         the initial point. The code automatically selects an */
/*         initial step size which is practically always suitable for */
/*         the problem. By using the fact that the code will not step */
/*         past TOUT in the first step, you could, if necessary, */
/*         restrict the length of the initial step size. */

/*         For some problems it may not be permissible to integrate */
/*         past a point TSTOP because a discontinuity occurs there */
/*         or the solution or its derivative is not defined beyond */
/*         TSTOP. When you have declared a TSTOP point (SEE INFO(4) */
/*         and RWORK(1)), you have told the code not to integrate */
/*         past TSTOP. In this case any TOUT beyond TSTOP is invalid */
/*         input. */

/*  INFO(*) -- Use the INFO array to give the code more details about */
/*         how you want your problem solved.  This array should be */
/*         dimensioned of length 15, though DDASSL uses only the first */
/*         eleven entries.  You must respond to all of the following */
/*         items, which are arranged as questions.  The simplest use */
/*         of the code corresponds to answering all questions as yes, */
/*         i.e. setting all entries of INFO to 0. */

/*       INFO(1) - This parameter enables the code to initialize */
/*              itself. You must set it to indicate the start of every */
/*              new problem. */

/*          **** Is this the first call for this problem ... */
/*                Yes - Set INFO(1) = 0 */
/*                 No - Not applicable here. */
/*                      See below for continuation calls.  **** */

/*       INFO(2) - How much accuracy you want of your solution */
/*              is specified by the error tolerances RTOL and ATOL. */
/*              The simplest use is to take them both to be scalars. */
/*              To obtain more flexibility, they can both be vectors. */
/*              The code must be told your choice. */

/*          **** Are both error tolerances RTOL, ATOL scalars ... */
/*                Yes - Set INFO(2) = 0 */
/*                      and input scalars for both RTOL and ATOL */
/*                 No - Set INFO(2) = 1 */
/*                      and input arrays for both RTOL and ATOL **** */

/*       INFO(3) - The code integrates from T in the direction */
/*              of TOUT by steps. If you wish, it will return the */
/*              computed solution and derivative at the next */
/*              intermediate step (the intermediate-output mode) or */
/*              TOUT, whichever comes first. This is a good way to */
/*              proceed if you want to see the behavior of the solution. */
/*              If you must have solutions at a great many specific */
/*              TOUT points, this code will compute them efficiently. */

/*          **** Do you want the solution only at */
/*                TOUT (and not at the next intermediate step) ... */
/*                 Yes - Set INFO(3) = 0 */
/*                  No - Set INFO(3) = 1 **** */

/*       INFO(4) - To handle solutions at a great many specific */
/*              values TOUT efficiently, this code may integrate past */
/*              TOUT and interpolate to obtain the result at TOUT. */
/*              Sometimes it is not possible to integrate beyond some */
/*              point TSTOP because the equation changes there or it is */
/*              not defined past TSTOP. Then you must tell the code */
/*              not to go past. */

/*           **** Can the integration be carried out without any */
/*                restrictions on the independent variable T ... */
/*                 Yes - Set INFO(4)=0 */
/*                  No - Set INFO(4)=1 */
/*                       and define the stopping point TSTOP by */
/*                       setting RWORK(1)=TSTOP **** */

/*       INFO(5) - To solve differential/algebraic problems it is */
/*              necessary to use a matrix of partial derivatives of the */
/*              system of differential equations. If you do not */
/*              provide a subroutine to evaluate it analytically (see */
/*              description of the item JAC in the call list), it will */
/*              be approximated by numerical differencing in this code. */
/*              although it is less trouble for you to have the code */
/*              compute partial derivatives by numerical differencing, */
/*              the solution will be more reliable if you provide the */
/*              derivatives via JAC. Sometimes numerical differencing */
/*              is cheaper than evaluating derivatives in JAC and */
/*              sometimes it is not - this depends on your problem. */

/*           **** Do you want the code to evaluate the partial */
/*                derivatives automatically by numerical differences ... */
/*                   Yes - Set INFO(5)=0 */
/*                    No - Set INFO(5)=1 */
/*                  and provide subroutine JAC for evaluating the */
/*                  matrix of partial derivatives **** */

/*       INFO(6) - DDASSL will perform much better if the matrix of */
/*              partial derivatives, DG/DY + CJ*DG/DYPRIME, */
/*              (here CJ is a scalar determined by DDASSL) */
/*              is banded and the code is told this. In this */
/*              case, the storage needed will be greatly reduced, */
/*              numerical differencing will be performed much cheaper, */
/*              and a number of important algorithms will execute much */
/*              faster. The differential equation is said to have */
/*              half-bandwidths ML (lower) and MU (upper) if equation i */
/*              involves only unknowns Y(J) with */
/*                             I-ML .LE. J .LE. I+MU */
/*              for all I=1,2,...,NEQ. Thus, ML and MU are the widths */
/*              of the lower and upper parts of the band, respectively, */
/*              with the main diagonal being excluded. If you do not */
/*              indicate that the equation has a banded matrix of partial */
/*              derivatives, the code works with a full matrix of NEQ**2 */
/*              elements (stored in the conventional way). Computations */
/*              with banded matrices cost less time and storage than with */
/*              full matrices if 2*ML+MU .LT. NEQ. If you tell the */
/*              code that the matrix of partial derivatives has a banded */
/*              structure and you want to provide subroutine JAC to */
/*              compute the partial derivatives, then you must be careful */
/*              to store the elements of the matrix in the special form */
/*              indicated in the description of JAC. */

/*          **** Do you want to solve the problem using a full */
/*               (dense) matrix (and not a special banded */
/*               structure) ... */
/*                Yes - Set INFO(6)=0 */
/*                 No - Set INFO(6)=1 */
/*                       and provide the lower (ML) and upper (MU) */
/*                       bandwidths by setting */
/*                       IWORK(1)=ML */
/*                       IWORK(2)=MU **** */


/*        INFO(7) -- You can specify a maximum (absolute value of) */
/*              stepsize, so that the code */
/*              will avoid passing over very */
/*              large regions. */

/*          ****  Do you want the code to decide */
/*                on its own maximum stepsize? */
/*                Yes - Set INFO(7)=0 */
/*                 No - Set INFO(7)=1 */
/*                      and define HMAX by setting */
/*                      RWORK(2)=HMAX **** */

/*        INFO(8) -- Differential/algebraic problems */
/*              may occasionally suffer from */
/*              severe scaling difficulties on the */
/*              first step. If you know a great deal */
/*              about the scaling of your problem, you can */
/*              help to alleviate this problem by */
/*              specifying an initial stepsize HO. */

/*          ****  Do you want the code to define */
/*                its own initial stepsize? */
/*                Yes - Set INFO(8)=0 */
/*                 No - Set INFO(8)=1 */
/*                      and define HO by setting */
/*                      RWORK(3)=HO **** */

/*        INFO(9) -- If storage is a severe problem, */
/*              you can save some locations by */
/*              restricting the maximum order MAXORD. */
/*              the default value is 5. for each */
/*              order decrease below 5, the code */
/*              requires NEQ fewer locations, however */
/*              it is likely to be slower. In any */
/*              case, you must have 1 .LE. MAXORD .LE. 5 */
/*          ****  Do you want the maximum order to */
/*                default to 5? */
/*                Yes - Set INFO(9)=0 */
/*                 No - Set INFO(9)=1 */
/*                      and define MAXORD by setting */
/*                      IWORK(3)=MAXORD **** */

/*        INFO(10) --If you know that the solutions to your equations */
/*               will always be nonnegative, it may help to set this */
/*               parameter. However, it is probably best to */
/*               try the code without using this option first, */
/*               and only to use this option if that doesn't */
/*               work very well. */
/*           ****  Do you want the code to solve the problem without */
/*                 invoking any special nonnegativity constraints? */
/*                  Yes - Set INFO(10)=0 */
/*                   No - Set INFO(10)=1 */

/*        INFO(11) --DDASSL normally requires the initial T, */
/*               Y, and YPRIME to be consistent. That is, */
/*               you must have G(T,Y,YPRIME) = 0 at the initial */
/*               time. If you do not know the initial */
/*               derivative precisely, you can let DDASSL try */
/*               to compute it. */
/*          ****   Are the initial T, Y, YPRIME consistent? */
/*                 Yes - Set INFO(11) = 0 */
/*                  No - Set INFO(11) = 1, */
/*                       and set YPRIME to an initial approximation */
/*                       to YPRIME.  (If you have no idea what */
/*                       YPRIME should be, set it to zero. Note */
/*                       that the initial Y should be such */
/*                       that there must exist a YPRIME so that */
/*                       G(T,Y,YPRIME) = 0.) */

/*  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL */
/*         error tolerances to tell the code how accurately you */
/*         want the solution to be computed.  They must be defined */
/*         as variables because the code may change them.  You */
/*         have two choices -- */
/*               Both RTOL and ATOL are scalars. (INFO(2)=0) */
/*               Both RTOL and ATOL are vectors. (INFO(2)=1) */
/*         in either case all components must be non-negative. */

/*         The tolerances are used by the code in a local error */
/*         test at each step which requires roughly that */
/*               ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL */
/*         for each vector component. */
/*         (More specifically, a root-mean-square norm is used to */
/*         measure the size of vectors, and the error test uses the */
/*         magnitude of the solution at the beginning of the step.) */

/*         The true (global) error is the difference between the */
/*         true solution of the initial value problem and the */
/*         computed approximation.  Practically all present day */
/*         codes, including this one, control the local error at */
/*         each step and do not even attempt to control the global */
/*         error directly. */
/*         Usually, but not always, the true accuracy of the */
/*         computed Y is comparable to the error tolerances. This */
/*         code will usually, but not always, deliver a more */
/*         accurate solution if you reduce the tolerances and */
/*         integrate again.  By comparing two such solutions you */
/*         can get a fairly reliable idea of the true error in the */
/*         solution at the bigger tolerances. */

/*         Setting ATOL=0. results in a pure relative error test on */
/*         that component.  Setting RTOL=0. results in a pure */
/*         absolute error test on that component.  A mixed test */
/*         with non-zero RTOL and ATOL corresponds roughly to a */
/*         relative error test when the solution component is much */
/*         bigger than ATOL and to an absolute error test when the */
/*         solution component is smaller than the threshhold ATOL. */

/*         The code will not attempt to compute a solution at an */
/*         accuracy unreasonable for the machine being used.  It will */
/*         advise you if you ask for too much accuracy and inform */
/*         you as to the maximum accuracy it believes possible. */

/*  RWORK(*) --  Dimension this real work array of length LRW in your */
/*         calling program. */

/*  LRW -- Set it to the declared length of the RWORK array. */
/*               You must have */
/*                    LRW .GE. 40+(MAXORD+4)*NEQ+NEQ**2 */
/*               for the full (dense) JACOBIAN case (when INFO(6)=0), or */
/*                    LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ */
/*               for the banded user-defined JACOBIAN case */
/*               (when INFO(5)=1 and INFO(6)=1), or */
/*                     LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ */
/*                           +2*(NEQ/(ML+MU+1)+1) */
/*               for the banded finite-difference-generated JACOBIAN case */
/*               (when INFO(5)=0 and INFO(6)=1) */

/*  IWORK(*) --  Dimension this integer work array of length LIW in */
/*         your calling program. */

/*  LIW -- Set it to the declared length of the IWORK array. */
/*               You must have LIW .GE. 20+NEQ */

/*  RPAR, IPAR -- These are parameter arrays, of real and integer */
/*         type, respectively.  You can use them for communication */
/*         between your program that calls DDASSL and the */
/*         RES subroutine (and the JAC subroutine).  They are not */
/*         altered by DDASSL.  If you do not need RPAR or IPAR, */
/*         ignore these parameters by treating them as dummy */
/*         arguments.  If you do choose to use them, dimension */
/*         them in your calling program and in RES (and in JAC) */
/*         as arrays of appropriate length. */

/*  JAC -- If you have set INFO(5)=0, you can ignore this parameter */
/*         by treating it as a dummy argument.  Otherwise, you must */
/*         provide a subroutine of the form */
/*             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR) */
/*         to define the matrix of partial derivatives */
/*             PD=DG/DY+CJ*DG/DYPRIME */
/*         CJ is a scalar which is input to JAC. */
/*         For the given values of T,Y,YPRIME, the */
/*         subroutine must evaluate the non-zero partial */
/*         derivatives for each equation and each solution */
/*         component, and store these values in the */
/*         matrix PD.  The elements of PD are set to zero */
/*         before each call to JAC so only non-zero elements */
/*         need to be defined. */

/*         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ. */
/*         You must declare the name JAC in an EXTERNAL statement in */
/*         your program that calls DDASSL.  You must dimension Y, */
/*         YPRIME and PD in JAC. */

/*         The way you must store the elements into the PD matrix */
/*         depends on the structure of the matrix which you */
/*         indicated by INFO(6). */
/*               *** INFO(6)=0 -- Full (dense) matrix *** */
/*                   Give PD a first dimension of NEQ. */
/*                   When you evaluate the (non-zero) partial derivative */
/*                   of equation I with respect to variable J, you must */
/*                   store it in PD according to */
/*                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)" */
/*               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU */
/*                   upper diagonal bands (refer to INFO(6) description */
/*                   of ML and MU) *** */
/*                   Give PD a first dimension of 2*ML+MU+1. */
/*                   when you evaluate the (non-zero) partial derivative */
/*                   of equation I with respect to variable J, you must */
/*                   store it in PD according to */
/*                   IROW = I - J + ML + MU + 1 */
/*                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)" */

/*         RPAR and IPAR are real and integer parameter arrays */
/*         which you can use for communication between your calling */
/*         program and your JACOBIAN subroutine JAC. They are not */
/*         altered by DDASSL. If you do not need RPAR or IPAR, */
/*         ignore these parameters by treating them as dummy */
/*         arguments. If you do choose to use them, dimension */
/*         them in your calling program and in JAC as arrays of */
/*         appropriate length. */


/*  OPTIONALLY REPLACEABLE NORM ROUTINE: */

/*     DDASSL uses a weighted norm DDANRM to measure the size */
/*     of vectors such as the estimated error in each step. */
/*     A FUNCTION subprogram */
/*       DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR) */
/*       DIMENSION V(NEQ),WT(NEQ) */
/*     is used to define this norm. Here, V is the vector */
/*     whose norm is to be computed, and WT is a vector of */
/*     weights.  A DDANRM routine has been included with DDASSL */
/*     which computes the weighted root-mean-square norm */
/*     given by */
/*       DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2) */
/*     this norm is suitable for most problems. In some */
/*     special cases, it may be more convenient and/or */
/*     efficient to define your own norm by writing a function */
/*     subprogram to be called instead of DDANRM. This should, */
/*     however, be attempted only after careful thought and */
/*     consideration. */


/*  -------- OUTPUT -- AFTER ANY RETURN FROM DDASSL --------------------- */

/*  The principal aim of the code is to return a computed solution at */
/*  TOUT, although it is also possible to obtain intermediate results */
/*  along the way. To find out whether the code achieved its goal */
/*  or if the integration process was interrupted before the task was */
/*  completed, you must check the IDID parameter. */


/*  T -- The solution was successfully advanced to the */
/*               output value of T. */

/*  Y(*) -- Contains the computed solution approximation at T. */

/*  YPRIME(*) -- Contains the computed derivative */
/*               approximation at T. */

/*  IDID -- Reports what the code did. */

/*                     *** Task completed *** */
/*                Reported by positive values of IDID */

/*           IDID = 1 -- A step was successfully taken in the */
/*                   intermediate-output mode. The code has not */
/*                   yet reached TOUT. */

/*           IDID = 2 -- The integration to TSTOP was successfully */
/*                   completed (T=TSTOP) by stepping exactly to TSTOP. */

/*           IDID = 3 -- The integration to TOUT was successfully */
/*                   completed (T=TOUT) by stepping past TOUT. */
/*                   Y(*) is obtained by interpolation. */
/*                   YPRIME(*) is obtained by interpolation. */

/*                    *** Task interrupted *** */
/*                Reported by negative values of IDID */

/*           IDID = -1 -- A large amount of work has been expended. */
/*                   (About 500 steps) */

/*           IDID = -2 -- The error tolerances are too stringent. */

/*           IDID = -3 -- The local error test cannot be satisfied */
/*                   because you specified a zero component in ATOL */
/*                   and the corresponding computed solution */
/*                   component is zero. Thus, a pure relative error */
/*                   test is impossible for this component. */

/*           IDID = -6 -- DDASSL had repeated error test */
/*                   failures on the last attempted step. */

/*           IDID = -7 -- The corrector could not converge. */

/*           IDID = -8 -- The matrix of partial derivatives */
/*                   is singular. */

/*           IDID = -9 -- The corrector could not converge. */
/*                   there were repeated error test failures */
/*                   in this step. */

/*           IDID =-10 -- The corrector could not converge */
/*                   because IRES was equal to minus one. */

/*           IDID =-11 -- IRES equal to -2 was encountered */
/*                   and control is being returned to the */
/*                   calling program. */

/*           IDID =-12 -- DDASSL failed to compute the initial */
/*                   YPRIME. */



/*           IDID = -13,..,-32 -- Not applicable for this code */

/*                    *** Task terminated *** */
/*                Reported by the value of IDID=-33 */

/*           IDID = -33 -- The code has encountered trouble from which */
/*                   it cannot recover. A message is printed */
/*                   explaining the trouble and control is returned */
/*                   to the calling program. For example, this occurs */
/*                   when invalid input is detected. */

/*  RTOL, ATOL -- These quantities remain unchanged except when */
/*               IDID = -2. In this case, the error tolerances have been */
/*               increased by the code to values which are estimated to */
/*               be appropriate for continuing the integration. However, */
/*               the reported solution at T was obtained using the input */
/*               values of RTOL and ATOL. */

/*  RWORK, IWORK -- Contain information which is usually of no */
/*               interest to the user but necessary for subsequent calls. */
/*               However, you may find use for */

/*               RWORK(3)--Which contains the step size H to be */
/*                       attempted on the next step. */

/*               RWORK(4)--Which contains the current value of the */
/*                       independent variable, i.e., the farthest point */
/*                       integration has reached. This will be different */
/*                       from T only when interpolation has been */
/*                       performed (IDID=3). */

/*               RWORK(7)--Which contains the stepsize used */
/*                       on the last successful step. */

/*               IWORK(7)--Which contains the order of the method to */
/*                       be attempted on the next step. */

/*               IWORK(8)--Which contains the order of the method used */
/*                       on the last step. */

/*               IWORK(11)--Which contains the number of steps taken so */
/*                        far. */

/*               IWORK(12)--Which contains the number of calls to RES */
/*                        so far. */

/*               IWORK(13)--Which contains the number of evaluations of */
/*                        the matrix of partial derivatives needed so */
/*                        far. */

/*               IWORK(14)--Which contains the total number */
/*                        of error test failures so far. */

/*               IWORK(15)--Which contains the total number */
/*                        of convergence test failures so far. */
/*                        (includes singular iteration matrix */
/*                        failures.) */


/*  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------ */
/*                    (CALLS AFTER THE FIRST) */

/*  This code is organized so that subsequent calls to continue the */
/*  integration involve little (if any) additional effort on your */
/*  part. You must monitor the IDID parameter in order to determine */
/*  what to do next. */

/*  Recalling that the principal task of the code is to integrate */
/*  from T to TOUT (the interval mode), usually all you will need */
/*  to do is specify a new TOUT upon reaching the current TOUT. */

/*  Do not alter any quantity not specifically permitted below, */
/*  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*) */
/*  or the differential equation in subroutine RES. Any such */
/*  alteration constitutes a new problem and must be treated as such, */
/*  i.e., you must start afresh. */

/*  You cannot change from vector to scalar error control or vice */
/*  versa (INFO(2)), but you can change the size of the entries of */
/*  RTOL, ATOL. Increasing a tolerance makes the equation easier */
/*  to integrate. Decreasing a tolerance will make the equation */
/*  harder to integrate and should generally be avoided. */

/*  You can switch from the intermediate-output mode to the */
/*  interval mode (INFO(3)) or vice versa at any time. */

/*  If it has been necessary to prevent the integration from going */
/*  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the */
/*  code will not integrate to any TOUT beyond the currently */
/*  specified TSTOP. Once TSTOP has been reached you must change */
/*  the value of TSTOP or set INFO(4)=0. You may change INFO(4) */
/*  or TSTOP at any time but you must supply the value of TSTOP in */
/*  RWORK(1) whenever you set INFO(4)=1. */

/*  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2) */
/*  unless you are going to restart the code. */

/*                 *** Following a completed task *** */
/*  If */
/*     IDID = 1, call the code again to continue the integration */
/*                  another step in the direction of TOUT. */

/*     IDID = 2 or 3, define a new TOUT and call the code again. */
/*                  TOUT must be different from T. You cannot change */
/*                  the direction of integration without restarting. */

/*                 *** Following an interrupted task *** */
/*               To show the code that you realize the task was */
/*               interrupted and that you want to continue, you */
/*               must take appropriate action and set INFO(1) = 1 */
/*  If */
/*    IDID = -1, The code has taken about 500 steps. */
/*                  If you want to continue, set INFO(1) = 1 and */
/*                  call the code again. An additional 500 steps */
/*                  will be allowed. */

/*    IDID = -2, The error tolerances RTOL, ATOL have been */
/*                  increased to values the code estimates appropriate */
/*                  for continuing. You may want to change them */
/*                  yourself. If you are sure you want to continue */
/*                  with relaxed error tolerances, set INFO(1)=1 and */
/*                  call the code again. */

/*    IDID = -3, A solution component is zero and you set the */
/*                  corresponding component of ATOL to zero. If you */
/*                  are sure you want to continue, you must first */
/*                  alter the error criterion to use positive values */
/*                  for those components of ATOL corresponding to zero */
/*                  solution components, then set INFO(1)=1 and call */
/*                  the code again. */

/*    IDID = -4,-5  --- Cannot occur with this code. */

/*    IDID = -6, Repeated error test failures occurred on the */
/*                  last attempted step in DDASSL. A singularity in the */
/*                  solution may be present. If you are absolutely */
/*                  certain you want to continue, you should restart */
/*                  the integration. (Provide initial values of Y and */
/*                  YPRIME which are consistent) */

/*    IDID = -7, Repeated convergence test failures occurred */
/*                  on the last attempted step in DDASSL. An inaccurate */
/*                  or ill-conditioned JACOBIAN may be the problem. If */
/*                  you are absolutely certain you want to continue, you */
/*                  should restart the integration. */

/*    IDID = -8, The matrix of partial derivatives is singular. */
/*                  Some of your equations may be redundant. */
/*                  DDASSL cannot solve the problem as stated. */
/*                  It is possible that the redundant equations */
/*                  could be removed, and then DDASSL could */
/*                  solve the problem. It is also possible */
/*                  that a solution to your problem either */
/*                  does not exist or is not unique. */

/*    IDID = -9, DDASSL had multiple convergence test */
/*                  failures, preceded by multiple error */
/*                  test failures, on the last attempted step. */
/*                  It is possible that your problem */
/*                  is ill-posed, and cannot be solved */
/*                  using this code. Or, there may be a */
/*                  discontinuity or a singularity in the */
/*                  solution. If you are absolutely certain */
/*                  you want to continue, you should restart */
/*                  the integration. */

/*    IDID =-10, DDASSL had multiple convergence test failures */
/*                  because IRES was equal to minus one. */
/*                  If you are absolutely certain you want */
/*                  to continue, you should restart the */
/*                  integration. */

/*    IDID =-11, IRES=-2 was encountered, and control is being */
/*                  returned to the calling program. */

/*    IDID =-12, DDASSL failed to compute the initial YPRIME. */
/*                  This could happen because the initial */
/*                  approximation to YPRIME was not very good, or */
/*                  if a YPRIME consistent with the initial Y */
/*                  does not exist. The problem could also be caused */
/*                  by an inaccurate or singular iteration matrix. */

/*    IDID = -13,..,-32  --- Cannot occur with this code. */


/*                 *** Following a terminated task *** */

/*  If IDID= -33, you cannot continue the solution of this problem. */
/*                  An attempt to do so will result in your */
/*                  run being terminated. */


/*  -------- ERROR MESSAGES --------------------------------------------- */

/*      The SLATEC error print routine XERMSG is called in the event of */
/*   unsuccessful completion of a task.  Most of these are treated as */
/*   "recoverable errors", which means that (unless the user has directed */
/*   otherwise) control will be returned to the calling program for */
/*   possible action after the message has been printed. */

/*   In the event of a negative value of IDID other than -33, an appro- */
/*   priate message is printed and the "error number" printed by XERMSG */
/*   is the value of IDID.  There are quite a number of illegal input */
/*   errors that can lead to a returned value IDID=-33.  The conditions */
/*   and their printed "error numbers" are as follows: */

/*   Error number       Condition */

/*        1       Some element of INFO vector is not zero or one. */
/*        2       NEQ .le. 0 */
/*        3       MAXORD not in range. */
/*        4       LRW is less than the required length for RWORK. */
/*        5       LIW is less than the required length for IWORK. */
/*        6       Some element of RTOL is .lt. 0 */
/*        7       Some element of ATOL is .lt. 0 */
/*        8       All elements of RTOL and ATOL are zero. */
/*        9       INFO(4)=1 and TSTOP is behind TOUT. */
/*       10       HMAX .lt. 0.0 */
/*       11       TOUT is behind T. */
/*       12       INFO(8)=1 and H0=0.0 */
/*       13       Some element of WT is .le. 0.0 */
/*       14       TOUT is too close to T to start integration. */
/*       15       INFO(4)=1 and TSTOP is behind T. */
/*       16       --( Not used in this version )-- */
/*       17       ML illegal.  Either .lt. 0 or .gt. NEQ */
/*       18       MU illegal.  Either .lt. 0 or .gt. NEQ */
/*       19       TOUT = T. */

/*   If DDASSL is called again without any action taken to remove the */
/*   cause of an unsuccessful return, XERMSG will be called with a fatal */
/*   error flag, which will cause unconditional termination of the */
/*   program.  There are two such fatal errors: */

/*   Error number -998:  The last step was terminated with a negative */
/*       value of IDID other than -33, and no appropriate action was */
/*       taken. */

/*   Error number -999:  The previous call was terminated because of */
/*       illegal input (IDID=-33) and there is illegal input in the */
/*       present call, as well.  (Suspect infinite loop.) */

/*  --------------------------------------------------------------------- */

/* ***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC */
/*                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637, */
/*                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982. */
/* ***ROUTINES CALLED  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS, */
/*                    XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   880387  Code changes made.  All common statements have been */
/*           replaced by a DATA statement, which defines pointers into */
/*           RWORK, and PARAMETER statements which define pointers */
/*           into IWORK.  As well the documentation has gone through */
/*           grammatical changes. */
/*   881005  The prologue has been changed to mixed case. */
/*           The subordinate routines had revision dates changed to */
/*           this date, although the documentation for these routines */
/*           is all upper case.  No code changes. */
/*   890511  Code changes made.  The DATA statement in the declaration */
/*           section of DDASSL was replaced with a PARAMETER */
/*           statement.  Also the statement S = 100.D0 was removed */
/*           from the top of the Newton iteration in DDASTP. */
/*           The subordinate routines had revision dates changed to */
/*           this date. */
/*   890517  The revision date syntax was replaced with the revision */
/*           history syntax.  Also the "DECK" comment was added to */
/*           the top of all subroutines.  These changes are consistent */
/*           with new SLATEC guidelines. */
/*           The subordinate routines had revision dates changed to */
/*           this date.  No code changes. */
/*   891013  Code changes made. */
/*           Removed all occurrences of FLOAT or DBLE.  All operations */
/*           are now performed with "mixed-mode" arithmetic. */
/*           Also, specific function names were replaced with generic */
/*           function names to be consistent with new SLATEC guidelines. */
/*           In particular: */
/*              Replaced DSQRT with SQRT everywhere. */
/*              Replaced DABS with ABS everywhere. */
/*              Replaced DMIN1 with MIN everywhere. */
/*              Replaced MIN0 with MIN everywhere. */
/*              Replaced DMAX1 with MAX everywhere. */
/*              Replaced MAX0 with MAX everywhere. */
/*              Replaced DSIGN with SIGN everywhere. */
/*           Also replaced REVISION DATE with REVISION HISTORY in all */
/*           subordinate routines. */
/*   901004  Miscellaneous changes to prologue to complete conversion */
/*           to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch) */
/*   901009  Corrected GAMS classification code and converted subsidiary */
/*           routines to 4.0 format.  No code changes.  (F.N.Fritsch) */
/*   901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens, AFWL) */
/*   901019  Code changes made. */
/*           Merged SLATEC 4.0 changes with previous changes made */
/*           by C. Ulrich.  Below is a history of the changes made by */
/*           C. Ulrich. (Changes in subsidiary routines are implied */
/*           by this history) */
/*           891228  Bug was found and repaired inside the DDASSL */
/*                   and DDAINI routines.  DDAINI was incorrectly */
/*                   returning the initial T with Y and YPRIME */
/*                   computed at T+H.  The routine now returns T+H */
/*                   rather than the initial T. */
/*                   Cosmetic changes made to DDASTP. */
/*           900904  Three modifications were made to fix a bug (inside */
/*                   DDASSL) re interpolation for continuation calls and */
/*                   cases where TN is very close to TSTOP: */

/*                   1) In testing for whether H is too large, just */
/*                      compare H to (TSTOP - TN), rather than */
/*                      (TSTOP - TN) * (1-4*UROUND), and set H to */
/*                      TSTOP - TN.  This will force DDASTP to step */
/*                      exactly to TSTOP under certain situations */
/*                      (i.e. when H returned from DDASTP would otherwise */
/*                      take TN beyond TSTOP). */

/*                   2) Inside the DDASTP loop, interpolate exactly to */
/*                      TSTOP if TN is very close to TSTOP (rather than */
/*                      interpolating to within roundoff of TSTOP). */

/*                   3) Modified IDID description for IDID = 2 to say */
/*                      that the solution is returned by stepping exactly */
/*                      to TSTOP, rather than TOUT.  (In some cases the */
/*                      solution is actually obtained by extrapolating */
/*                      over a distance near unit roundoff to TSTOP, */
/*                      but this small distance is deemed acceptable in */
/*                      these circumstances.) */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue, removed unreferenced labels, */
/*           and improved XERMSG calls.  (FNF) */
/*   901030  Added ERROR MESSAGES section and reworked other sections to */
/*           be of more uniform format.  (FNF) */
/*   910624  Fixed minor bug related to HMAX (six lines after label */
/*           525).  (LRP) */
/* ***END PROLOGUE  DDASSL */

/* **End */

/*     Declare arguments. */


/*     Declare externals. */


/*     Declare local variables. */

/*       Auxiliary variables for conversion of values to be included in */
/*       error messages. */

/*     SET POINTERS INTO IWORK */

/*     SET RELATIVE OFFSET INTO RWORK */

/*     SET POINTERS INTO RWORK */

/* ***FIRST EXECUTABLE STATEMENT  DDASSL */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwork;
    --rwork;
    --atol;
    --rtol;
    --info;
    --yprime;
    --y;

    /* Function Body */
    if (info[1] != 0) {
	goto L100;
    }

/* ----------------------------------------------------------------------- */
/*     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY. */
/*     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS. */
/* ----------------------------------------------------------------------- */

/*     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO */
/*     ARE EITHER ZERO OR ONE. */
    for (i__ = 2; i__ <= 11; ++i__) {
	if (info[i__] != 0 && info[i__] != 1) {
	    goto L701;
	}
/* L10: */
    }

    if (*neq <= 0) {
	goto L702;
    }

/*     CHECK AND COMPUTE MAXIMUM ORDER */
    mxord = 5;
    if (info[9] == 0) {
	goto L20;
    }
    mxord = iwork[3];
    if (mxord < 1 || mxord > 5) {
	goto L703;
    }
L20:
    iwork[3] = mxord;

/*     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU. */
    if (info[6] != 0) {
	goto L40;
    }
/* Computing 2nd power */
    i__1 = *neq;
    lenpd = i__1 * i__1;
    lenrw = (iwork[3] + 4) * *neq + 40 + lenpd;
    if (info[5] != 0) {
	goto L30;
    }
    iwork[4] = 2;
    goto L60;
L30:
    iwork[4] = 1;
    goto L60;
L40:
    if (iwork[1] < 0 || iwork[1] >= *neq) {
	goto L717;
    }
    if (iwork[2] < 0 || iwork[2] >= *neq) {
	goto L718;
    }
    lenpd = ((iwork[1] << 1) + iwork[2] + 1) * *neq;
    if (info[5] != 0) {
	goto L50;
    }
    iwork[4] = 5;
    mband = iwork[1] + iwork[2] + 1;
    msave = *neq / mband + 1;
    lenrw = (iwork[3] + 4) * *neq + 40 + lenpd + (msave << 1);
    goto L60;
L50:
    iwork[4] = 4;
    lenrw = (iwork[3] + 4) * *neq + 40 + lenpd;

/*     CHECK LENGTHS OF RWORK AND IWORK */
L60:
    leniw = *neq + 20;
    iwork[16] = lenpd;
    if (*lrw < lenrw) {
	goto L704;
    }
    if (*liw < leniw) {
	goto L705;
    }

/*     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T */
    if (*tout == *t) {
	goto L719;
    }

/*     CHECK HMAX */
    if (info[7] == 0) {
	goto L70;
    }
    hmax = rwork[2];
    if (hmax <= 0.) {
	goto L710;
    }
L70:

/*     INITIALIZE COUNTERS */
    iwork[11] = 0;
    iwork[12] = 0;
    iwork[13] = 0;

    iwork[10] = 0;
    *idid = 1;
    goto L200;

/* ----------------------------------------------------------------------- */
/*     THIS BLOCK IS FOR CONTINUATION CALLS */
/*     ONLY. HERE WE CHECK INFO(1), AND IF THE */
/*     LAST STEP WAS INTERRUPTED WE CHECK WHETHER */
/*     APPROPRIATE ACTION WAS TAKEN. */
/* ----------------------------------------------------------------------- */

L100:
    if (info[1] == 1) {
	goto L110;
    }
    if (info[1] != -1) {
	goto L701;
    }

/*     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED */
/*     BY AN ERROR CONDITION FROM DDASTP, AND */
/*     APPROPRIATE ACTION WAS NOT TAKEN. THIS */
/*     IS A FATAL ERROR. */
    s_wsfi(&io___10);
    do_fio(&c__1, (char *)&(*idid), (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 57, a__1[0] = "THE LAST STEP TERMINATED WITH A NEGATIVE VALUE "
	    "OF IDID = ";
    i__2[1] = 8, a__1[1] = xern1;
    i__2[2] = 39, a__1[2] = " AND NO APPROPRIATE ACTION WAS TAKEN.  ";
    i__2[3] = 14, a__1[3] = "RUN TERMINATED";
    s_cat(ch__1, a__1, i__2, &c__4, (ftnlen)118);
    xermsg_("SLATEC", "DDASSL", ch__1, &c_n998, &c__2, (ftnlen)6, (ftnlen)6, (
	    ftnlen)118);
    return 0;
L110:
    iwork[10] = iwork[11];

/* ----------------------------------------------------------------------- */
/*     THIS BLOCK IS EXECUTED ON ALL CALLS. */
/*     THE ERROR TOLERANCE PARAMETERS ARE */
/*     CHECKED, AND THE WORK ARRAY POINTERS */
/*     ARE SET. */
/* ----------------------------------------------------------------------- */

L200:
/*     CHECK RTOL,ATOL */
    nzflg = 0;
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (info[2] == 1) {
	    rtoli = rtol[i__];
	}
	if (info[2] == 1) {
	    atoli = atol[i__];
	}
	if (rtoli > 0. || atoli > 0.) {
	    nzflg = 1;
	}
	if (rtoli < 0.) {
	    goto L706;
	}
	if (atoli < 0.) {
	    goto L707;
	}
/* L210: */
    }
    if (nzflg == 0) {
	goto L708;
    }

/*     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED */
/*     IN DATA STATEMENT. */
    le = *neq + 41;
    lwt = le + *neq;
    lphi = lwt + *neq;
    lpd = lphi + (iwork[3] + 1) * *neq;
    lwm = lpd;
    ntemp = iwork[16] + 1;
    if (info[1] == 1) {
	goto L400;
    }

/* ----------------------------------------------------------------------- */
/*     THIS BLOCK IS EXECUTED ON THE INITIAL CALL */
/*     ONLY. SET THE INITIAL STEP SIZE, AND */
/*     THE ERROR WEIGHT VECTOR, AND PHI. */
/*     COMPUTE INITIAL YPRIME, IF NECESSARY. */
/* ----------------------------------------------------------------------- */

    tn = *t;
    *idid = 1;

/*     SET ERROR WEIGHT VECTOR WT */
    ddawts_(neq, &info[2], &rtol[1], &atol[1], &y[1], &rwork[lwt], &rpar[1], &
	    ipar[1]);
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[lwt + i__ - 1] <= 0.) {
	    goto L713;
	}
/* L305: */
    }

/*     COMPUTE UNIT ROUNDOFF AND HMIN */
    uround = d1mach_(&c__4);
    rwork[9] = uround;
/* Computing MAX */
    d__1 = abs(*t), d__2 = abs(*tout);
    hmin = uround * 4. * max(d__1,d__2);

/*     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH */
    tdist = (d__1 = *tout - *t, abs(d__1));
    if (tdist < hmin) {
	goto L714;
    }

/*     CHECK HO, IF THIS WAS INPUT */
    if (info[8] == 0) {
	goto L310;
    }
    ho = rwork[3];
    if ((*tout - *t) * ho < 0.) {
	goto L711;
    }
    if (ho == 0.) {
	goto L712;
    }
    goto L320;
L310:

/*     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER */
/*     DDASTP OR DDAINI, DEPENDING ON INFO(11) */
    ho = tdist * .001;
    ypnorm = ddanrm_(neq, &yprime[1], &rwork[lwt], &rpar[1], &ipar[1]);
    if (ypnorm > .5 / ho) {
	ho = .5 / ypnorm;
    }
    d__1 = *tout - *t;
    ho = d_sign(&ho, &d__1);
/*     ADJUST HO IF NECESSARY TO MEET HMAX BOUND */
L320:
    if (info[7] == 0) {
	goto L330;
    }
    rh = abs(ho) / rwork[2];
    if (rh > 1.) {
	ho /= rh;
    }
/*     COMPUTE TSTOP, IF APPLICABLE */
L330:
    if (info[4] == 0) {
	goto L340;
    }
    tstop = rwork[1];
    if ((tstop - *t) * ho < 0.) {
	goto L715;
    }
    if ((*t + ho - tstop) * ho > 0.) {
	ho = tstop - *t;
    }
    if ((tstop - *tout) * ho < 0.) {
	goto L709;
    }

/*     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, IF APPLICABLE */
L340:
    if (info[11] == 0) {
	goto L350;
    }
    ddaini_(&tn, &y[1], &yprime[1], neq, (U_fp)res, (U_fp)jac, &ho, &rwork[
	    lwt], idid, &rpar[1], &ipar[1], &rwork[lphi], &rwork[41], &rwork[
	    le], &rwork[lwm], &iwork[1], &hmin, &rwork[9], &info[10], &ntemp);
    if (*idid < 0) {
	goto L390;
    }

/*     LOAD H WITH HO.  STORE H IN RWORK(LH) */
L350:
    h__ = ho;
    rwork[3] = h__;

/*     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2) */
    itemp = lphi + *neq;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rwork[lphi + i__ - 1] = y[i__];
/* L370: */
	rwork[itemp + i__ - 1] = h__ * yprime[i__];
    }

L390:
    goto L500;

/* ------------------------------------------------------- */
/*     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS */
/*     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE */
/*     TAKING A STEP. */
/*     ADJUST H IF NECESSARY TO MEET HMAX BOUND */
/* ------------------------------------------------------- */

L400:
    uround = rwork[9];
    done = FALSE_;
    tn = rwork[4];
    h__ = rwork[3];
    if (info[7] == 0) {
	goto L410;
    }
    rh = abs(h__) / rwork[2];
    if (rh > 1.) {
	h__ /= rh;
    }
L410:
    if (*t == *tout) {
	goto L719;
    }
    if ((*t - *tout) * h__ > 0.) {
	goto L711;
    }
    if (info[4] == 1) {
	goto L430;
    }
    if (info[3] == 1) {
	goto L420;
    }
    if ((tn - *tout) * h__ < 0.) {
	goto L490;
    }
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *t = *tout;
    *idid = 3;
    done = TRUE_;
    goto L490;
L420:
    if ((tn - *t) * h__ <= 0.) {
	goto L490;
    }
    if ((tn - *tout) * h__ > 0.) {
	goto L425;
    }
    ddatrp_(&tn, &tn, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &rwork[
	    29]);
    *t = tn;
    *idid = 1;
    done = TRUE_;
    goto L490;
L425:
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *t = *tout;
    *idid = 3;
    done = TRUE_;
    goto L490;
L430:
    if (info[3] == 1) {
	goto L440;
    }
    tstop = rwork[1];
    if ((tn - tstop) * h__ > 0.) {
	goto L715;
    }
    if ((tstop - *tout) * h__ < 0.) {
	goto L709;
    }
    if ((tn - *tout) * h__ < 0.) {
	goto L450;
    }
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *t = *tout;
    *idid = 3;
    done = TRUE_;
    goto L490;
L440:
    tstop = rwork[1];
    if ((tn - tstop) * h__ > 0.) {
	goto L715;
    }
    if ((tstop - *tout) * h__ < 0.) {
	goto L709;
    }
    if ((tn - *t) * h__ <= 0.) {
	goto L450;
    }
    if ((tn - *tout) * h__ > 0.) {
	goto L445;
    }
    ddatrp_(&tn, &tn, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &rwork[
	    29]);
    *t = tn;
    *idid = 1;
    done = TRUE_;
    goto L490;
L445:
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *t = *tout;
    *idid = 3;
    done = TRUE_;
    goto L490;
L450:
/*     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP */
    if ((d__1 = tn - tstop, abs(d__1)) > uround * 100. * (abs(tn) + abs(h__)))
	     {
	goto L460;
    }
    ddatrp_(&tn, &tstop, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *idid = 2;
    *t = tstop;
    done = TRUE_;
    goto L490;
L460:
    tnext = tn + h__;
    if ((tnext - tstop) * h__ <= 0.) {
	goto L490;
    }
    h__ = tstop - tn;
    rwork[3] = h__;

L490:
    if (done) {
	goto L580;
    }

/* ------------------------------------------------------- */
/*     THE NEXT BLOCK CONTAINS THE CALL TO THE */
/*     ONE-STEP INTEGRATOR DDASTP. */
/*     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS. */
/*     CHECK FOR TOO MANY STEPS. */
/*     UPDATE WT. */
/*     CHECK FOR TOO MUCH ACCURACY REQUESTED. */
/*     COMPUTE MINIMUM STEPSIZE. */
/* ------------------------------------------------------- */

L500:
/*     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME */
    if (*idid == -12) {
	goto L527;
    }

/*     CHECK FOR TOO MANY STEPS */
    if (iwork[11] - iwork[10] < 500) {
	goto L510;
    }
    *idid = -1;
    goto L527;

/*     UPDATE WT */
L510:
    ddawts_(neq, &info[2], &rtol[1], &atol[1], &rwork[lphi], &rwork[lwt], &
	    rpar[1], &ipar[1]);
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + lwt - 1] > 0.) {
	    goto L520;
	}
	*idid = -3;
	goto L527;
L520:
	;
    }

/*     TEST FOR TOO MUCH ACCURACY REQUESTED. */
    r__ = ddanrm_(neq, &rwork[lphi], &rwork[lwt], &rpar[1], &ipar[1]) * 100. *
	     uround;
    if (r__ <= 1.) {
	goto L525;
    }
/*     MULTIPLY RTOL AND ATOL BY R AND RETURN */
    if (info[2] == 1) {
	goto L523;
    }
    rtol[1] = r__ * rtol[1];
    atol[1] = r__ * atol[1];
    *idid = -2;
    goto L527;
L523:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rtol[i__] = r__ * rtol[i__];
/* L524: */
	atol[i__] = r__ * atol[i__];
    }
    *idid = -2;
    goto L527;
L525:

/*     COMPUTE MINIMUM STEPSIZE */
/* Computing MAX */
    d__1 = abs(tn), d__2 = abs(*tout);
    hmin = uround * 4. * max(d__1,d__2);

/*     TEST H VS. HMAX */
    if (info[7] != 0) {
	rh = abs(h__) / rwork[2];
	if (rh > 1.) {
	    h__ /= rh;
	}
    }

    ddastp_(&tn, &y[1], &yprime[1], neq, (U_fp)res, (U_fp)jac, &h__, &rwork[
	    lwt], &info[1], idid, &rpar[1], &ipar[1], &rwork[lphi], &rwork[41]
	    , &rwork[le], &rwork[lwm], &iwork[1], &rwork[11], &rwork[17], &
	    rwork[23], &rwork[29], &rwork[35], &rwork[5], &rwork[6], &rwork[7]
	    , &rwork[8], &hmin, &rwork[9], &iwork[6], &iwork[5], &iwork[7], &
	    iwork[8], &iwork[9], &info[10], &ntemp);
L527:
    if (*idid < 0) {
	goto L600;
    }

/* -------------------------------------------------------- */
/*     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN */
/*     FROM DDASTP (IDID=1).  TEST FOR STOP CONDITIONS. */
/* -------------------------------------------------------- */

    if (info[4] != 0) {
	goto L540;
    }
    if (info[3] != 0) {
	goto L530;
    }
    if ((tn - *tout) * h__ < 0.) {
	goto L500;
    }
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *idid = 3;
    *t = *tout;
    goto L580;
L530:
    if ((tn - *tout) * h__ >= 0.) {
	goto L535;
    }
    *t = tn;
    *idid = 1;
    goto L580;
L535:
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *idid = 3;
    *t = *tout;
    goto L580;
L540:
    if (info[3] != 0) {
	goto L550;
    }
    if ((tn - *tout) * h__ < 0.) {
	goto L542;
    }
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *t = *tout;
    *idid = 3;
    goto L580;
L542:
    if ((d__1 = tn - tstop, abs(d__1)) <= uround * 100. * (abs(tn) + abs(h__))
	    ) {
	goto L545;
    }
    tnext = tn + h__;
    if ((tnext - tstop) * h__ <= 0.) {
	goto L500;
    }
    h__ = tstop - tn;
    goto L500;
L545:
    ddatrp_(&tn, &tstop, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *idid = 2;
    *t = tstop;
    goto L580;
L550:
    if ((tn - *tout) * h__ >= 0.) {
	goto L555;
    }
    if ((d__1 = tn - tstop, abs(d__1)) <= uround * 100. * (abs(tn) + abs(h__))
	    ) {
	goto L552;
    }
    *t = tn;
    *idid = 1;
    goto L580;
L552:
    ddatrp_(&tn, &tstop, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *idid = 2;
    *t = tstop;
    goto L580;
L555:
    ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[29]);
    *t = *tout;
    *idid = 3;
    goto L580;

/* -------------------------------------------------------- */
/*     ALL SUCCESSFUL RETURNS FROM DDASSL ARE MADE FROM */
/*     THIS BLOCK. */
/* -------------------------------------------------------- */

L580:
    rwork[4] = tn;
    rwork[3] = h__;
    return 0;

/* ----------------------------------------------------------------------- */
/*     THIS BLOCK HANDLES ALL UNSUCCESSFUL */
/*     RETURNS OTHER THAN FOR ILLEGAL INPUT. */
/* ----------------------------------------------------------------------- */

L600:
    itemp = -(*idid);
    switch (itemp) {
	case 1:  goto L610;
	case 2:  goto L620;
	case 3:  goto L630;
	case 4:  goto L690;
	case 5:  goto L690;
	case 6:  goto L640;
	case 7:  goto L650;
	case 8:  goto L660;
	case 9:  goto L670;
	case 10:  goto L675;
	case 11:  goto L680;
	case 12:  goto L685;
    }

/*     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE */
/*     REACHING TOUT */
L610:
    s_wsfi(&io___34);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 15, a__1[0] = "AT CURRENT T = ";
    i__2[1] = 16, a__1[1] = xern3;
    i__2[2] = 25, a__1[2] = " 500 STEPS TAKEN ON THIS ";
    i__2[3] = 25, a__1[3] = "CALL BEFORE REACHING TOUT";
    s_cat(ch__2, a__1, i__2, &c__4, (ftnlen)81);
    xermsg_("SLATEC", "DDASSL", ch__2, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)81);
    goto L690;

/*     TOO MUCH ACCURACY FOR MACHINE PRECISION */
L620:
    s_wsfi(&io___35);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 7, a__2[0] = "AT T = ";
    i__3[1] = 16, a__2[1] = xern3;
    i__3[2] = 33, a__2[2] = " TOO MUCH ACCURACY REQUESTED FOR ";
    i__3[3] = 54, a__2[3] = "PRECISION OF MACHINE. RTOL AND ATOL WERE INCREA"
	    "SED TO ";
    i__3[4] = 18, a__2[4] = "APPROPRIATE VALUES";
    s_cat(ch__3, a__2, i__3, &c__5, (ftnlen)128);
    xermsg_("SLATEC", "DDASSL", ch__3, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)128);
    goto L690;

/*     WT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM) */
L630:
    s_wsfi(&io___36);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 7, a__1[0] = "AT T = ";
    i__2[1] = 16, a__1[1] = xern3;
    i__2[2] = 36, a__1[2] = " SOME ELEMENT OF WT HAS BECOME .LE. ";
    i__2[3] = 3, a__1[3] = "0.0";
    s_cat(ch__4, a__1, i__2, &c__4, (ftnlen)62);
    xermsg_("SLATEC", "DDASSL", ch__4, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)62);
    goto L690;

/*     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN */
L640:
    s_wsfi(&io___37);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___39);
    do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 7, a__2[0] = "AT T = ";
    i__3[1] = 16, a__2[1] = xern3;
    i__3[2] = 18, a__2[2] = " AND STEPSIZE H = ";
    i__3[3] = 16, a__2[3] = xern4;
    i__3[4] = 53, a__2[4] = " THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H"
	    ")=HMIN";
    s_cat(ch__5, a__2, i__3, &c__5, (ftnlen)110);
    xermsg_("SLATEC", "DDASSL", ch__5, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)110);
    goto L690;

/*     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN */
L650:
    s_wsfi(&io___40);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___41);
    do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__4[0] = 7, a__3[0] = "AT T = ";
    i__4[1] = 16, a__3[1] = xern3;
    i__4[2] = 18, a__3[2] = " AND STEPSIZE H = ";
    i__4[3] = 16, a__3[3] = xern4;
    i__4[4] = 53, a__3[4] = " THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR"
	    " WITH ";
    i__4[5] = 11, a__3[5] = "ABS(H)=HMIN";
    s_cat(ch__6, a__3, i__4, &c__6, (ftnlen)121);
    xermsg_("SLATEC", "DDASSL", ch__6, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)121);
    goto L690;

/*     THE ITERATION MATRIX IS SINGULAR */
L660:
    s_wsfi(&io___42);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___43);
    do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 7, a__2[0] = "AT T = ";
    i__3[1] = 16, a__2[1] = xern3;
    i__3[2] = 18, a__2[2] = " AND STEPSIZE H = ";
    i__3[3] = 16, a__2[3] = xern4;
    i__3[4] = 33, a__2[4] = " THE ITERATION MATRIX IS SINGULAR";
    s_cat(ch__7, a__2, i__3, &c__5, (ftnlen)90);
    xermsg_("SLATEC", "DDASSL", ch__7, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)90);
    goto L690;

/*     CORRECTOR FAILURE PRECEDED BY ERROR TEST FAILURES. */
L670:
    s_wsfi(&io___44);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___45);
    do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__4[0] = 7, a__3[0] = "AT T = ";
    i__4[1] = 16, a__3[1] = xern3;
    i__4[2] = 18, a__3[2] = " AND STEPSIZE H = ";
    i__4[3] = 16, a__3[3] = xern4;
    i__4[4] = 57, a__3[4] = " THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE E"
	    "RROR TEST ";
    i__4[5] = 18, a__3[5] = "FAILED REPEATEDLY.";
    s_cat(ch__8, a__3, i__4, &c__6, (ftnlen)132);
    xermsg_("SLATEC", "DDASSL", ch__8, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)132);
    goto L690;

/*     CORRECTOR FAILURE BECAUSE IRES = -1 */
L675:
    s_wsfi(&io___46);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___47);
    do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__4[0] = 7, a__3[0] = "AT T = ";
    i__4[1] = 16, a__3[1] = xern3;
    i__4[2] = 18, a__3[2] = " AND STEPSIZE H = ";
    i__4[3] = 16, a__3[3] = xern4;
    i__4[4] = 57, a__3[4] = " THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES "
	    "WAS EQUAL ";
    i__4[5] = 12, a__3[5] = "TO MINUS ONE";
    s_cat(ch__9, a__3, i__4, &c__6, (ftnlen)126);
    xermsg_("SLATEC", "DDASSL", ch__9, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)126);
    goto L690;

/*     FAILURE BECAUSE IRES = -2 */
L680:
    s_wsfi(&io___48);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___49);
    do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 7, a__2[0] = "AT T = ";
    i__3[1] = 16, a__2[1] = xern3;
    i__3[2] = 18, a__2[2] = " AND STEPSIZE H = ";
    i__3[3] = 16, a__2[3] = xern4;
    i__3[4] = 28, a__2[4] = " IRES WAS EQUAL TO MINUS TWO";
    s_cat(ch__10, a__2, i__3, &c__5, (ftnlen)85);
    xermsg_("SLATEC", "DDASSL", ch__10, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)85);
    goto L690;

/*     FAILED TO COMPUTE INITIAL YPRIME */
L685:
    s_wsfi(&io___50);
    do_fio(&c__1, (char *)&tn, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___51);
    do_fio(&c__1, (char *)&ho, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 7, a__2[0] = "AT T = ";
    i__3[1] = 16, a__2[1] = xern3;
    i__3[2] = 18, a__2[2] = " AND STEPSIZE H = ";
    i__3[3] = 16, a__2[3] = xern4;
    i__3[4] = 41, a__2[4] = " THE INITIAL YPRIME COULD NOT BE COMPUTED";
    s_cat(ch__11, a__2, i__3, &c__5, (ftnlen)98);
    xermsg_("SLATEC", "DDASSL", ch__11, idid, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)98);
    goto L690;

L690:
    info[1] = -1;
    *t = tn;
    rwork[4] = tn;
    rwork[3] = h__;
    return 0;

/* ----------------------------------------------------------------------- */
/*     THIS BLOCK HANDLES ALL ERROR RETURNS DUE */
/*     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING */
/*     DDASTP. FIRST THE ERROR MESSAGE ROUTINE IS */
/*     CALLED. IF THIS HAPPENS TWICE IN */
/*     SUCCESSION, EXECUTION IS TERMINATED */

/* ----------------------------------------------------------------------- */
L701:
    xermsg_("SLATEC", "DDASSL", "SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR "
	    "ONE", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)46);
    goto L750;

L702:
    s_wsfi(&io___52);
    do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__5[0] = 6, a__4[0] = "NEQ = ";
    i__5[1] = 8, a__4[1] = xern1;
    i__5[2] = 7, a__4[2] = " .LE. 0";
    s_cat(ch__12, a__4, i__5, &c__3, (ftnlen)21);
    xermsg_("SLATEC", "DDASSL", ch__12, &c__2, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)21);
    goto L750;

L703:
    s_wsfi(&io___53);
    do_fio(&c__1, (char *)&mxord, (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__5[0] = 9, a__4[0] = "MAXORD = ";
    i__5[1] = 8, a__4[1] = xern1;
    i__5[2] = 13, a__4[2] = " NOT IN RANGE";
    s_cat(ch__13, a__4, i__5, &c__3, (ftnlen)30);
    xermsg_("SLATEC", "DDASSL", ch__13, &c__3, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)30);
    goto L750;

L704:
    s_wsfi(&io___54);
    do_fio(&c__1, (char *)&lenrw, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___56);
    do_fio(&c__1, (char *)&(*lrw), (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 29, a__1[0] = "RWORK LENGTH NEEDED, LENRW = ";
    i__2[1] = 8, a__1[1] = xern1;
    i__2[2] = 16, a__1[2] = ", EXCEEDS LRW = ";
    i__2[3] = 8, a__1[3] = xern2;
    s_cat(ch__14, a__1, i__2, &c__4, (ftnlen)61);
    xermsg_("SLATEC", "DDASSL", ch__14, &c__4, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)61);
    goto L750;

L705:
    s_wsfi(&io___57);
    do_fio(&c__1, (char *)&leniw, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___58);
    do_fio(&c__1, (char *)&(*liw), (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 29, a__1[0] = "IWORK LENGTH NEEDED, LENIW = ";
    i__2[1] = 8, a__1[1] = xern1;
    i__2[2] = 16, a__1[2] = ", EXCEEDS LIW = ";
    i__2[3] = 8, a__1[3] = xern2;
    s_cat(ch__14, a__1, i__2, &c__4, (ftnlen)61);
    xermsg_("SLATEC", "DDASSL", ch__14, &c__5, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)61);
    goto L750;

L706:
    xermsg_("SLATEC", "DDASSL", "SOME ELEMENT OF RTOL IS .LT. 0", &c__6, &
	    c__1, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    goto L750;

L707:
    xermsg_("SLATEC", "DDASSL", "SOME ELEMENT OF ATOL IS .LT. 0", &c__7, &
	    c__1, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    goto L750;

L708:
    xermsg_("SLATEC", "DDASSL", "ALL ELEMENTS OF RTOL AND ATOL ARE ZERO", &
	    c__8, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)38);
    goto L750;

L709:
    s_wsfi(&io___59);
    do_fio(&c__1, (char *)&tstop, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___60);
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 24, a__1[0] = "INFO(4) = 1 AND TSTOP = ";
    i__2[1] = 16, a__1[1] = xern3;
    i__2[2] = 15, a__1[2] = " BEHIND TOUT = ";
    i__2[3] = 16, a__1[3] = xern4;
    s_cat(ch__15, a__1, i__2, &c__4, (ftnlen)71);
    xermsg_("SLATEC", "DDASSL", ch__15, &c__9, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)71);
    goto L750;

L710:
    s_wsfi(&io___61);
    do_fio(&c__1, (char *)&hmax, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__5[0] = 7, a__4[0] = "HMAX = ";
    i__5[1] = 16, a__4[1] = xern3;
    i__5[2] = 9, a__4[2] = " .LT. 0.0";
    s_cat(ch__16, a__4, i__5, &c__3, (ftnlen)32);
    xermsg_("SLATEC", "DDASSL", ch__16, &c__10, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)32);
    goto L750;

L711:
    s_wsfi(&io___62);
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___63);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 7, a__1[0] = "TOUT = ";
    i__2[1] = 16, a__1[1] = xern3;
    i__2[2] = 12, a__1[2] = " BEHIND T = ";
    i__2[3] = 16, a__1[3] = xern4;
    s_cat(ch__17, a__1, i__2, &c__4, (ftnlen)51);
    xermsg_("SLATEC", "DDASSL", ch__17, &c__11, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)51);
    goto L750;

L712:
    xermsg_("SLATEC", "DDASSL", "INFO(8)=1 AND H0=0.0", &c__12, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)20);
    goto L750;

L713:
    xermsg_("SLATEC", "DDASSL", "SOME ELEMENT OF WT IS .LE. 0.0", &c__13, &
	    c__1, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    goto L750;

L714:
    s_wsfi(&io___64);
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___65);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 7, a__2[0] = "TOUT = ";
    i__3[1] = 16, a__2[1] = xern3;
    i__3[2] = 18, a__2[2] = " TOO CLOSE TO T = ";
    i__3[3] = 16, a__2[3] = xern4;
    i__3[4] = 21, a__2[4] = " TO START INTEGRATION";
    s_cat(ch__18, a__2, i__3, &c__5, (ftnlen)78);
    xermsg_("SLATEC", "DDASSL", ch__18, &c__14, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)78);
    goto L750;

L715:
    s_wsfi(&io___66);
    do_fio(&c__1, (char *)&tstop, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___67);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 22, a__1[0] = "INFO(4)=1 AND TSTOP = ";
    i__2[1] = 16, a__1[1] = xern3;
    i__2[2] = 12, a__1[2] = " BEHIND T = ";
    i__2[3] = 16, a__1[3] = xern4;
    s_cat(ch__19, a__1, i__2, &c__4, (ftnlen)66);
    xermsg_("SLATEC", "DDASSL", ch__19, &c__15, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)66);
    goto L750;

L717:
    s_wsfi(&io___68);
    do_fio(&c__1, (char *)&iwork[1], (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__5[0] = 5, a__4[0] = "ML = ";
    i__5[1] = 8, a__4[1] = xern1;
    i__5[2] = 36, a__4[2] = " ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ";
    s_cat(ch__20, a__4, i__5, &c__3, (ftnlen)49);
    xermsg_("SLATEC", "DDASSL", ch__20, &c__17, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)49);
    goto L750;

L718:
    s_wsfi(&io___69);
    do_fio(&c__1, (char *)&iwork[2], (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__5[0] = 5, a__4[0] = "MU = ";
    i__5[1] = 8, a__4[1] = xern1;
    i__5[2] = 36, a__4[2] = " ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ";
    s_cat(ch__20, a__4, i__5, &c__3, (ftnlen)49);
    xermsg_("SLATEC", "DDASSL", ch__20, &c__18, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)49);
    goto L750;

L719:
    s_wsfi(&io___70);
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__6[0] = 11, a__5[0] = "TOUT = T = ";
    i__6[1] = 16, a__5[1] = xern3;
    s_cat(ch__21, a__5, i__6, &c__2, (ftnlen)27);
    xermsg_("SLATEC", "DDASSL", ch__21, &c__19, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)27);
    goto L750;

L750:
    *idid = -33;
    if (info[1] == -1) {
	xermsg_("SLATEC", "DDASSL", "REPEATED OCCURRENCES OF ILLEGAL INPUT$$"
		"RUN TERMINATED. APPARENT INFINITE LOOP", &c_n999, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)77);
    }

    info[1] = -1;
    return 0;
/* -----------END OF SUBROUTINE DDASSL------------------------------------ */
} /* ddassl_ */


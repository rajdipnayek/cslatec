/* ddriv3.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__4 = 4;
static doublereal c_b107 = 1.;
static integer c__0 = 0;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;

/* DECK DDRIV3 */
/* Subroutine */ int ddriv3_(integer *n, doublereal *t, doublereal *y, S_fp f,
	 integer *nstate, doublereal *tout, integer *ntask, integer *nroot, 
	doublereal *eps, doublereal *ewt, integer *ierror, integer *mint, 
	integer *miter, integer *impl, integer *ml, integer *mu, integer *
	mxord, doublereal *hmax, doublereal *work, integer *lenw, integer *
	iwork, integer *leniw, U_fp jacobn, S_fp fa, integer *nde, integer *
	mxstep, D_fp g, S_fp users, integer *ierflg)
{
    /* System generated locals */
    address a__1[3], a__2[5], a__3[6], a__4[7], a__5[8], a__6[4];
    integer i__1[3], i__2[5], i__3, i__4[6], i__5[7], i__6[8], i__7, i__8[4];
    doublereal d__1, d__2, d__3;
    char ch__1[54], ch__2[51], ch__3[63], ch__4[57], ch__5[74], ch__6[53], 
	    ch__7[52], ch__8[98], ch__9[75], ch__10[156], ch__11[155], ch__12[
	    128], ch__13[145], ch__14[196], ch__15[96], ch__16[119], ch__17[
	    64];

    /* Local variables */
    static doublereal h__;
    static integer i__, j;
    static doublereal ae;
    static integer ia;
    static doublereal re;
    static char rl1[16], rl2[16];
    static doublereal big, sum;
    static integer ifac, info, npar;
    static doublereal gnow, size;
    static integer iywt;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dgbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), dgefa_(doublereal *, 
	    integer *, integer *, integer *, integer *);
    static integer iflag;
    extern /* Subroutine */ int dgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    static integer idfdy;
    extern /* Subroutine */ int dgesl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *);
    static doublereal hsign, hused, glast;
    extern /* Subroutine */ int ddntp_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    , ddstp_(doublereal *, S_fp, S_fp, doublereal *, integer *, 
	    integer *, U_fp, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, logical *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *);
    static integer ignow;
    extern /* Subroutine */ int ddzro_(doublereal *, D_fp, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal tlast;
    static integer iroot;
    extern doublereal d1mach_(integer *);
    static doublereal troot;
    static integer isave1, isave2;
    static char intgr1[8], intgr2[8];
    static integer lenchk, ndecom, matdim, liwchk, jstate, maxord, imxerr;
    static doublereal uround;
    static integer nstepl, itroot, jtroot;
    static logical convrg;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___6 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___8 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___10 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___11 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, intgr2, 0, "(I8)", 8, 1 };
    static icilist io___14 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___15 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___17 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___28 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___41 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___43 = { 0, rl2, 0, "(D16.8)", 16, 1 };
    static icilist io___44 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___45 = { 0, rl2, 0, "(D16.8)", 16, 1 };
    static icilist io___51 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___52 = { 0, rl2, 0, "(D16.8)", 16, 1 };
    static icilist io___53 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___54 = { 0, rl2, 0, "(D16.8)", 16, 1 };
    static icilist io___55 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___56 = { 0, intgr1, 0, "(I8)", 8, 1 };
    static icilist io___57 = { 0, rl2, 0, "(D16.8)", 16, 1 };
    static icilist io___66 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___67 = { 0, rl1, 0, "(D16.8)", 16, 1 };
    static icilist io___68 = { 0, rl1, 0, "(D16.8)", 16, 1 };


/* ***BEGIN PROLOGUE  DDRIV3 */
/* ***PURPOSE  The function of DDRIV3 is to solve N ordinary differential */
/*            equations of the form dY(I)/dT = F(Y(I),T), given the */
/*            initial conditions Y(I) = YI.  The program has options to */
/*            allow the solution of both stiff and non-stiff differential */
/*            equations.  Other important options are available.  DDRIV3 */
/*            uses double precision arithmetic. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***CATEGORY  I1A2, I1A1B */
/* ***TYPE      DOUBLE PRECISION (SDRIV3-S, DDRIV3-D, CDRIV3-C) */
/* ***KEYWORDS  DOUBLE PRECISION, GEAR'S METHOD, INITIAL VALUE PROBLEMS, */
/*             ODE, ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, STIFF */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*  I.  ABSTRACT  ....................................................... */

/*    The primary function of DDRIV3 is to solve N ordinary differential */
/*    equations of the form dY(I)/dT = F(Y(I),T), given the initial */
/*    conditions Y(I) = YI.  The program has options to allow the */
/*    solution of both stiff and non-stiff differential equations.  In */
/*    addition, DDRIV3 may be used to solve: */
/*      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is */
/*         a non-singular matrix depending on Y and T. */
/*      2. The hybrid differential/algebraic initial value problem, */
/*         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may */
/*         depend upon Y and T) some of whose components will be zero */
/*         corresponding to those equations which are algebraic rather */
/*         than differential. */
/*    DDRIV3 is to be called once for each output point of T. */

/*  II.  PARAMETERS  .................................................... */
/*       (REMEMBER--To run DDRIV3 correctly in double precision, ALL */
/*       non-integer arguments in the call sequence, including */
/*       arrays, MUST be declared double precision.) */

/*    The user should use parameter names in the call sequence of DDRIV3 */
/*    for those quantities whose value may be altered by DDRIV3.  The */
/*    parameters in the call sequence are: */

/*    N      = (Input) The number of dependent functions whose solution */
/*             is desired.  N must not be altered during a problem. */

/*    T      = The independent variable.  On input for the first call, T */
/*             is the initial point.  On output, T is the point at which */
/*             the solution is given. */

/*    Y      = The vector of dependent variables.  Y is used as input on */
/*             the first call, to set the initial values.  On output, Y */
/*             is the computed solution vector.  This array Y is passed */
/*             in the call sequence of the user-provided routines F, */
/*             JACOBN, FA, USERS, and G.  Thus parameters required by */
/*             those routines can be stored in this array in components */
/*             N+1 and above.  (Note: Changes by the user to the first */
/*             N components of this array will take effect only after a */
/*             restart, i.e., after setting NSTATE to 1 .) */

/*    F      = A subroutine supplied by the user.  The name must be */
/*             declared EXTERNAL in the user's calling program.  This */
/*             subroutine is of the form: */
/*                   SUBROUTINE F (N, T, Y, YDOT) */
/*                   DOUBLE PRECISION Y(*), YDOT(*) */
/*                     . */
/*                     . */
/*                   YDOT(1) = ... */
/*                     . */
/*                     . */
/*                   YDOT(N) = ... */
/*                   END (Sample) */
/*             This computes YDOT = F(Y,T), the right hand side of the */
/*             differential equations.  Here Y is a vector of length at */
/*             least N.  The actual length of Y is determined by the */
/*             user's declaration in the program which calls DDRIV3. */
/*             Thus the dimensioning of Y in F, while required by FORTRAN */
/*             convention, does not actually allocate any storage.  When */
/*             this subroutine is called, the first N components of Y are */
/*             intermediate approximations to the solution components. */
/*             The user should not alter these values.  Here YDOT is a */
/*             vector of length N.  The user should only compute YDOT(I) */
/*             for I from 1 to N.  Normally a return from F passes */
/*             control back to  DDRIV3.  However, if the user would like */
/*             to abort the calculation, i.e., return control to the */
/*             program which calls DDRIV3, he should set N to zero. */
/*             DDRIV3 will signal this by returning a value of NSTATE */
/*             equal to 6 .  Altering the value of N in F has no effect */
/*             on the value of N in the call sequence of DDRIV3. */

/*    NSTATE = An integer describing the status of integration.  The */
/*             meaning of NSTATE is as follows: */
/*               1  (Input) Means the first call to the routine.  This */
/*                  value must be set by the user.  On all subsequent */
/*                  calls the value of NSTATE should be tested by the */
/*                  user, but must not be altered.  (As a convenience to */
/*                  the user who may wish to put out the initial */
/*                  conditions, DDRIV3 can be called with NSTATE=1, and */
/*                  TOUT=T.  In this case the program will return with */
/*                  NSTATE unchanged, i.e., NSTATE=1.) */
/*               2  (Output) Means a successful integration.  If a normal */
/*                  continuation is desired (i.e., a further integration */
/*                  in the same direction), simply advance TOUT and call */
/*                  again.  All other parameters are automatically set. */
/*               3  (Output)(Unsuccessful) Means the integrator has taken */
/*                  MXSTEP steps without reaching TOUT.  The user can */
/*                  continue the integration by simply calling DDRIV3 */
/*                  again. */
/*               4  (Output)(Unsuccessful) Means too much accuracy has */
/*                  been requested.  EPS has been increased to a value */
/*                  the program estimates is appropriate.  The user can */
/*                  continue the integration by simply calling DDRIV3 */
/*                  again. */
/*               5  (Output) A root was found at a point less than TOUT. */
/*                  The user can continue the integration toward TOUT by */
/*                  simply calling DDRIV3 again. */
/*               6  (Output)(Unsuccessful) N has been set to zero in */
/*                  SUBROUTINE F. */
/*               7  (Output)(Unsuccessful) N has been set to zero in */
/*                  FUNCTION G.  See description of G below. */
/*               8  (Output)(Unsuccessful) N has been set to zero in */
/*                  SUBROUTINE JACOBN.  See description of JACOBN below. */
/*               9  (Output)(Unsuccessful) N has been set to zero in */
/*                  SUBROUTINE FA.  See description of FA below. */
/*              10  (Output)(Unsuccessful) N has been set to zero in */
/*                  SUBROUTINE USERS.  See description of USERS below. */
/*              11  (Output)(Successful) For NTASK = 2 or 3, T is beyond */
/*                  TOUT.  The solution was obtained by interpolation. */
/*                  The user can continue the integration by simply */
/*                  advancing TOUT and calling DDRIV3 again. */
/*              12  (Output)(Unsuccessful) The solution could not be */
/*                  obtained.  The value of IERFLG (see description */
/*                  below) for a "Recoverable" situation indicates the */
/*                  type of difficulty encountered: either an illegal */
/*                  value for a parameter or an inability to continue the */
/*                  solution.  For this condition the user should take */
/*                  corrective action and reset NSTATE to 1 before */
/*                  calling DDRIV3 again.  Otherwise the program will */
/*                  terminate the run. */

/*    TOUT   = (Input) The point at which the solution is desired.  The */
/*             position of TOUT relative to T on the first call */
/*             determines the direction of integration. */

/*    NTASK  = (Input) An index specifying the manner of returning the */
/*             solution, according to the following: */
/*               NTASK = 1  Means DDRIV3 will integrate past TOUT and */
/*                          interpolate the solution.  This is the most */
/*                          efficient mode. */
/*               NTASK = 2  Means DDRIV3 will return the solution after */
/*                          each internal integration step, or at TOUT, */
/*                          whichever comes first.  In the latter case, */
/*                          the program integrates exactly to TOUT. */
/*               NTASK = 3  Means DDRIV3 will adjust its internal step to */
/*                          reach TOUT exactly (useful if a singularity */
/*                          exists beyond TOUT.) */

/*    NROOT  = (Input) The number of equations whose roots are desired. */
/*             If NROOT is zero, the root search is not active.  This */
/*             option is useful for obtaining output at points which are */
/*             not known in advance, but depend upon the solution, e.g., */
/*             when some solution component takes on a specified value. */
/*             The root search is carried out using the user-written */
/*             function G (see description of G below.)  DDRIV3 attempts */
/*             to find the value of T at which one of the equations */
/*             changes sign.  DDRIV3 can find at most one root per */
/*             equation per internal integration step, and will then */
/*             return the solution either at TOUT or at a root, whichever */
/*             occurs first in the direction of integration.  The initial */
/*             point is never reported as a root.  The index of the */
/*             equation whose root is being reported is stored in the */
/*             sixth element of IWORK. */
/*             NOTE: NROOT is never altered by this program. */

/*    EPS    = On input, the requested relative accuracy in all solution */
/*             components.  EPS = 0 is allowed.  On output, the adjusted */
/*             relative accuracy if the input value was too small.  The */
/*             value of EPS should be set as large as is reasonable, */
/*             because the amount of work done by DDRIV3 increases as EPS */
/*             decreases. */

/*    EWT    = (Input) Problem zero, i.e., the smallest, nonzero, */
/*             physically meaningful value for the solution.  (Array, */
/*             possibly of length one.  See following description of */
/*             IERROR.)  Setting EWT smaller than necessary can adversely */
/*             affect the running time. */

/*    IERROR = (Input) Error control indicator.  A value of 3 is */
/*             suggested for most problems.  Other choices and detailed */
/*             explanations of EWT and IERROR are given below for those */
/*             who may need extra flexibility. */

/*             These last three input quantities EPS, EWT and IERROR */
/*             control the accuracy of the computed solution.  EWT and */
/*             IERROR are used internally to compute an array YWT.  One */
/*             step error estimates divided by YWT(I) are kept less than */
/*             EPS in root mean square norm. */
/*                 IERROR (Set by the user) = */
/*                 1  Means YWT(I) = 1. (Absolute error control) */
/*                                   EWT is ignored. */
/*                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control) */
/*                                   EWT is ignored. */
/*                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)). */
/*                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)). */
/*                    This choice is useful when the solution components */
/*                    have differing scales. */
/*                 5  Means YWT(I) = EWT(I). */
/*             If IERROR is 3, EWT need only be dimensioned one. */
/*             If IERROR is 4 or 5, the user must dimension EWT at least */
/*             N, and set its values. */

/*    MINT   = (Input) The integration method indicator. */
/*               MINT = 1  Means the Adams methods, and is used for */
/*                         non-stiff problems. */
/*               MINT = 2  Means the stiff methods of Gear (i.e., the */
/*                         backward differentiation formulas), and is */
/*                         used for stiff problems. */
/*               MINT = 3  Means the program dynamically selects the */
/*                         Adams methods when the problem is non-stiff */
/*                         and the Gear methods when the problem is */
/*                         stiff.  When using the Adams methods, the */
/*                         program uses a value of MITER=0; when using */
/*                         the Gear methods, the program uses the value */
/*                         of MITER provided by the user.  Only a value */
/*                         of IMPL = 0 and a value of MITER = 1, 2, 4, or */
/*                         5 is allowed for this option.  The user may */
/*                         not alter the value of MINT or MITER without */
/*                         restarting, i.e., setting NSTATE to 1. */

/*    MITER  = (Input) The iteration method indicator. */
/*               MITER = 0  Means functional iteration.  This value is */
/*                          suggested for non-stiff problems. */
/*               MITER = 1  Means chord method with analytic Jacobian. */
/*                          In this case, the user supplies subroutine */
/*                          JACOBN (see description below). */
/*               MITER = 2  Means chord method with Jacobian calculated */
/*                          internally by finite differences. */
/*               MITER = 3  Means chord method with corrections computed */
/*                          by the user-written routine USERS (see */
/*                          description of USERS below.)  This option */
/*                          allows all matrix algebra and storage */
/*                          decisions to be made by the user.  When using */
/*                          a value of MITER = 3, the subroutine FA is */
/*                          not required, even if IMPL is not 0.  For */
/*                          further information on using this option, see */
/*                          Section IV-E below. */
/*               MITER = 4  Means the same as MITER = 1 but the A and */
/*                          Jacobian matrices are assumed to be banded. */
/*               MITER = 5  Means the same as MITER = 2 but the A and */
/*                          Jacobian matrices are assumed to be banded. */

/*    IMPL   = (Input) The implicit method indicator. */
/*               IMPL = 0    Means solving dY(I)/dT = F(Y(I),T). */
/*               IMPL = 1    Means solving A*dY(I)/dT = F(Y(I),T), non- */
/*                           singular A (see description of FA below.) */
/*                           Only MINT = 1 or 2, and MITER = 1, 2, 3, 4, */
/*                           or 5 are allowed for this option. */
/*               IMPL = 2,3  Means solving certain systems of hybrid */
/*                           differential/algebraic equations (see */
/*                           description of FA below.)  Only MINT = 2 and */
/*                           MITER = 1, 2, 3, 4, or 5, are allowed for */
/*                           this option. */
/*               The value of IMPL must not be changed during a problem. */

/*    ML     = (Input) The lower half-bandwidth in the case of a banded */
/*             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero */
/*             A(R,C).) */

/*    MU     = (Input) The upper half-bandwidth in the case of a banded */
/*             A or Jacobian matrix.  (I.e., maximum(C-R).) */

/*    MXORD  = (Input) The maximum order desired. This is .LE. 12 for */
/*             the Adams methods and .LE. 5 for the Gear methods.  Normal */
/*             value is 12 and 5, respectively.  If MINT is 3, the */
/*             maximum order used will be MIN(MXORD, 12) when using the */
/*             Adams methods, and MIN(MXORD, 5) when using the Gear */
/*             methods.  MXORD must not be altered during a problem. */

/*    HMAX   = (Input) The maximum magnitude of the step size that will */
/*             be used for the problem.  This is useful for ensuring that */
/*             important details are not missed.  If this is not the */
/*             case, a large value, such as the interval length, is */
/*             suggested. */

/*    WORK */
/*    LENW   = (Input) */
/*             WORK is an array of LENW double precision words used */
/*             internally for temporary storage.  The user must allocate */
/*             space for this array in the calling program by a statement */
/*             such as */
/*                       DOUBLE PRECISION WORK(...) */
/*             The following table gives the required minimum value for */
/*             the length of WORK, depending on the value of IMPL and */
/*             MITER.  LENW should be set to the value used.  The */
/*             contents of WORK should not be disturbed between calls to */
/*             DDRIV3. */

/*      IMPL =   0            1               2             3 */
/*              --------------------------------------------------------- */
/* MITER =  0   (MXORD+4)*N   Not allowed     Not allowed   Not allowed */
/*              + 2*NROOT */
/*              + 250 */

/*         1,2  N*N +         2*N*N +         N*N +         N*(N + NDE) */
/*              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N */
/*              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT */
/*              + 250         + 250           + 250         + 250 */

/*          3   (MXORD+4)*N   (MXORD+4)*N     (MXORD+4)*N   (MXORD+4)*N */
/*              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT */
/*              + 250         + 250           + 250         + 250 */

/*         4,5  (2*ML+MU+1)   2*(2*ML+MU+1)   (2*ML+MU+1)   (2*ML+MU+1)* */
/*              *N +          *N +            *N +          (N+NDE) + */
/*              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N */
/*              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT */
/*              + 250         + 250           + 250         + 250 */
/*              --------------------------------------------------------- */

/*    IWORK */
/*    LENIW  = (Input) */
/*             IWORK is an integer array of length LENIW used internally */
/*             for temporary storage.  The user must allocate space for */
/*             this array in the calling program by a statement such as */
/*                       INTEGER IWORK(...) */
/*             The length of IWORK should be at least */
/*               50      if MITER is 0 or 3, or */
/*               N+50    if MITER is 1, 2, 4, or 5, or MINT is 3, */
/*             and LENIW should be set to the value used.  The contents */
/*             of IWORK should not be disturbed between calls to DDRIV3. */

/*    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4. */
/*             If this is the case, the name must be declared EXTERNAL in */
/*             the user's calling program.  Given a system of N */
/*             differential equations, it is meaningful to speak about */
/*             the partial derivative of the I-th right hand side with */
/*             respect to the J-th dependent variable.  In general there */
/*             are N*N such quantities.  Often however the equations can */
/*             be ordered so that the I-th differential equation only */
/*             involves dependent variables with index near I, e.g., I+1, */
/*             I-2.  Such a system is called banded.  If, for all I, the */
/*             I-th equation depends on at most the variables */
/*               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU) */
/*             then we call ML+MU+1 the bandwidth of the system.  In a */
/*             banded system many of the partial derivatives above are */
/*             automatically zero.  For the cases MITER = 1, 2, 4, and 5, */
/*             some of these partials are needed.  For the cases */
/*             MITER = 2 and 5 the necessary derivatives are */
/*             approximated numerically by DDRIV3, and we only ask the */
/*             user to tell DDRIV3 the value of ML and MU if the system */
/*             is banded.  For the cases MITER = 1 and 4 the user must */
/*             derive these partials algebraically and encode them in */
/*             subroutine JACOBN.  By computing these derivatives the */
/*             user can often save 20-30 per cent of the computing time. */
/*             Usually, however, the accuracy is not much affected and */
/*             most users will probably forego this option.  The optional */
/*             user-written subroutine JACOBN has the form: */
/*                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU) */
/*                   DOUBLE PRECISION Y(*), DFDY(MATDIM,*) */
/*                     . */
/*                     . */
/*                     Calculate values of DFDY */
/*                     . */
/*                     . */
/*                   END (Sample) */
/*             Here Y is a vector of length at least N.  The actual */
/*             length of Y is determined by the user's declaration in the */
/*             program which calls DDRIV3.  Thus the dimensioning of Y in */
/*             JACOBN, while required by FORTRAN convention, does not */
/*             actually allocate any storage.  When this subroutine is */
/*             called, the first N components of Y are intermediate */
/*             approximations to the solution components.  The user */
/*             should not alter these values.  If the system is not */
/*             banded (MITER=1), the partials of the I-th equation with */
/*             respect to the J-th dependent function are to be stored in */
/*             DFDY(I,J).  Thus partials of the I-th equation are stored */
/*             in the I-th row of DFDY.  If the system is banded */
/*             (MITER=4), then the partials of the I-th equation with */
/*             respect to Y(J) are to be stored in DFDY(K,J), where */
/*             K=I-J+MU+1 .  Normally a return from JACOBN passes control */
/*             back to DDRIV3.  However, if the user would like to abort */
/*             the calculation, i.e., return control to the program which */
/*             calls DDRIV3, he should set N to zero.  DDRIV3 will signal */
/*             this by returning a value of NSTATE equal to +8(-8). */
/*             Altering the value of N in JACOBN has no effect on the */
/*             value of N in the call sequence of DDRIV3. */

/*    FA     = A subroutine supplied by the user if IMPL is not zero, and */
/*             MITER is not 3.  If so, the name must be declared EXTERNAL */
/*             in the user's calling program.  This subroutine computes */
/*             the array A, where A*dY(I)/dT = F(Y(I),T). */
/*             There are three cases: */

/*               IMPL=1. */
/*               Subroutine FA is of the form: */
/*                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE) */
/*                   DOUBLE PRECISION Y(*), A(MATDIM,*) */
/*                     . */
/*                     . */
/*                     Calculate ALL values of A */
/*                     . */
/*                     . */
/*                   END (Sample) */
/*               In this case A is assumed to be a nonsingular matrix, */
/*               with the same structure as DFDY (see JACOBN description */
/*               above).  Programming considerations prevent complete */
/*               generality.  If MITER is 1 or 2, A is assumed to be full */
/*               and the user must compute and store all values of */
/*               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed */
/*               to be banded with lower and upper half bandwidth ML and */
/*               MU.  The left hand side of the I-th equation is a linear */
/*               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... , */
/*               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the */
/*               I-th equation, the coefficient of dY(J)/dT is to be */
/*               stored in A(K,J), where K=I-J+MU+1. */
/*               NOTE: The array A will be altered between calls to FA. */

/*               IMPL=2. */
/*               Subroutine FA is of the form: */
/*                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE) */
/*                   DOUBLE PRECISION Y(*), A(*) */
/*                     . */
/*                     . */
/*                     Calculate non-zero values of A(1),...,A(NDE) */
/*                     . */
/*                     . */
/*                   END (Sample) */
/*               In this case it is assumed that the system is ordered by */
/*               the user so that the differential equations appear */
/*               first, and the algebraic equations appear last.  The */
/*               algebraic equations must be written in the form: */
/*               0 = F(Y(I),T).  When using this option it is up to the */
/*               user to provide initial values for the Y(I) that satisfy */
/*               the algebraic equations as well as possible.  It is */
/*               further assumed that A is a vector of length NDE.  All */
/*               of the components of A, which may depend on T, Y(I), */
/*               etc., must be set by the user to non-zero values. */

/*               IMPL=3. */
/*               Subroutine FA is of the form: */
/*                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE) */
/*                   DOUBLE PRECISION Y(*), A(MATDIM,*) */
/*                     . */
/*                     . */
/*                     Calculate ALL values of A */
/*                     . */
/*                     . */
/*                   END (Sample) */
/*               In this case A is assumed to be a nonsingular NDE by NDE */
/*               matrix with the same structure as DFDY (see JACOBN */
/*               description above).  Programming considerations prevent */
/*               complete generality.  If MITER is 1 or 2, A is assumed */
/*               to be full and the user must compute and store all */
/*               values of A(I,J), I,J=1, ... ,NDE.  If MITER is 4 or 5, */
/*               A is assumed to be banded with lower and upper half */
/*               bandwidths ML and MU.  The left hand side of the I-th */
/*               equation is a linear combination of dY(I-ML)/dT, */
/*               dY(I-ML+1)/dT, ... , dY(I)/dT, ... , dY(I+MU-1)/dT, */
/*               dY(I+MU)/dT.  Thus in the I-th equation, the coefficient */
/*               of dY(J)/dT is to be stored in A(K,J), where K=I-J+MU+1. */
/*               It is assumed that the system is ordered by the user so */
/*               that the differential equations appear first, and the */
/*               algebraic equations appear last.  The algebraic */
/*               equations must be written in the form 0 = F(Y(I),T). */
/*               When using this option it is up to the user to provide */
/*               initial values for the Y(I) that satisfy the algebraic */
/*               equations as well as possible. */
/*               NOTE: For IMPL = 3, the array A will be altered between */
/*               calls to FA. */
/*             Here Y is a vector of length at least N.  The actual */
/*             length of Y is determined by the user's declaration in the */
/*             program which calls DDRIV3.  Thus the dimensioning of Y in */
/*             FA, while required by FORTRAN convention, does not */
/*             actually allocate any storage.  When this subroutine is */
/*             called, the first N components of Y are intermediate */
/*             approximations to the solution components.  The user */
/*             should not alter these values.  FA is always called */
/*             immediately after calling F, with the same values of T */
/*             and Y.  Normally a return from FA passes control back to */
/*             DDRIV3.  However, if the user would like to abort the */
/*             calculation, i.e., return control to the program which */
/*             calls DDRIV3, he should set N to zero.  DDRIV3 will signal */
/*             this by returning a value of NSTATE equal to +9(-9). */
/*             Altering the value of N in FA has no effect on the value */
/*             of N in the call sequence of DDRIV3. */

/*    NDE    = (Input) The number of differential equations.  This is */
/*             required only for IMPL = 2 or 3, with NDE .LT. N. */

/*    MXSTEP = (Input) The maximum number of internal steps allowed on */
/*             one call to DDRIV3. */

/*    G      = A double precision FORTRAN function supplied by the user */
/*             if NROOT is not 0.  In this case, the name must be */
/*             declared EXTERNAL in the user's calling program.  G is */
/*             repeatedly called with different values of IROOT to obtain */
/*             the value of each of the NROOT equations for which a root */
/*             is desired.  G is of the form: */
/*                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT) */
/*                   DOUBLE PRECISION Y(*) */
/*                   GO TO (10, ...), IROOT */
/*              10   G = ... */
/*                     . */
/*                     . */
/*                   END (Sample) */
/*             Here, Y is a vector of length at least N, whose first N */
/*             components are the solution components at the point T. */
/*             The user should not alter these values.  The actual length */
/*             of Y is determined by the user's declaration in the */
/*             program which calls DDRIV3.  Thus the dimensioning of Y in */
/*             G, while required by FORTRAN convention, does not actually */
/*             allocate any storage.  Normally a return from G passes */
/*             control back to  DDRIV3.  However, if the user would like */
/*             to abort the calculation, i.e., return control to the */
/*             program which calls DDRIV3, he should set N to zero. */
/*             DDRIV3 will signal this by returning a value of NSTATE */
/*             equal to +7(-7).  In this case, the index of the equation */
/*             being evaluated is stored in the sixth element of IWORK. */
/*             Altering the value of N in G has no effect on the value of */
/*             N in the call sequence of DDRIV3. */

/*    USERS  = A subroutine supplied by the user, if MITER is 3. */
/*             If this is the case, the name must be declared EXTERNAL in */
/*             the user's calling program.  The routine USERS is called */
/*             by DDRIV3 when certain linear systems must be solved.  The */
/*             user may choose any method to form, store and solve these */
/*             systems in order to obtain the solution result that is */
/*             returned to DDRIV3.  In particular, this allows sparse */
/*             matrix methods to be used.  The call sequence for this */
/*             routine is: */

/*                SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, */
/*               8                  IMPL, N, NDE, IFLAG) */
/*                DOUBLE PRECISION Y(*), YH(*), YWT(*), SAVE1(*), */
/*               8     SAVE2(*), T, H, EL */

/*             The input variable IFLAG indicates what action is to be */
/*             taken.  Subroutine USERS should perform the following */
/*             operations, depending on the value of IFLAG and IMPL. */

/*               IFLAG = 0 */
/*                 IMPL = 0.  USERS is not called. */
/*                 IMPL = 1, 2 or 3.  Solve the system A*X = SAVE2, */
/*                   returning the result in SAVE2.  The array SAVE1 can */
/*                   be used as a work array.  For IMPL = 1, there are N */
/*                   components to the system, and for IMPL = 2 or 3, */
/*                   there are NDE components to the system. */

/*               IFLAG = 1 */
/*                 IMPL = 0.  Compute, decompose and store the matrix */
/*                   (I - H*EL*J), where I is the identity matrix and J */
/*                   is the Jacobian matrix of the right hand side.  The */
/*                   array SAVE1 can be used as a work array. */
/*                 IMPL = 1, 2 or 3. Compute, decompose and store the */
/*                   matrix (A - H*EL*J).  The array SAVE1 can be used as */
/*                   a work array. */

/*               IFLAG = 2 */
/*                 IMPL = 0.   Solve the system */
/*                     (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1, */
/*                   returning the result in SAVE2. */
/*                 IMPL = 1, 2 or 3.  Solve the system */
/*                   (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1) */
/*                   returning the result in SAVE2. */
/*                 The array SAVE1 should not be altered. */
/*             If IFLAG is 0 and IMPL is 1 or 2 and the matrix A is */
/*             singular, or if IFLAG is 1 and one of the matrices */
/*             (I - H*EL*J), (A - H*EL*J) is singular, the INTEGER */
/*             variable IFLAG is to be set to -1 before RETURNing. */
/*             Normally a return from USERS passes control back to */
/*             DDRIV3.  However, if the user would like to abort the */
/*             calculation, i.e., return control to the program which */
/*             calls DDRIV3, he should set N to zero.  DDRIV3 will signal */
/*             this by returning a value of NSTATE equal to +10(-10). */
/*             Altering the value of N in USERS has no effect on the */
/*             value of N in the call sequence of DDRIV3. */

/*    IERFLG = An error flag.  The error number associated with a */
/*             diagnostic message (see Section III-A below) is the same */
/*             as the corresponding value of IERFLG.  The meaning of */
/*             IERFLG: */
/*               0  The routine completed successfully. (No message is */
/*                  issued.) */
/*               3  (Warning) The number of steps required to reach TOUT */
/*                  exceeds MXSTEP. */
/*               4  (Warning) The value of EPS is too small. */
/*              11  (Warning) For NTASK = 2 or 3, T is beyond TOUT. */
/*                  The solution was obtained by interpolation. */
/*              15  (Warning) The integration step size is below the */
/*                  roundoff level of T.  (The program issues this */
/*                  message as a warning but does not return control to */
/*                  the user.) */
/*              22  (Recoverable) N is not positive. */
/*              23  (Recoverable) MINT is less than 1 or greater than 3 . */
/*              24  (Recoverable) MITER is less than 0 or greater than */
/*                  5 . */
/*              25  (Recoverable) IMPL is less than 0 or greater than 3 . */
/*              26  (Recoverable) The value of NSTATE is less than 1 or */
/*                  greater than 12 . */
/*              27  (Recoverable) EPS is less than zero. */
/*              28  (Recoverable) MXORD is not positive. */
/*              29  (Recoverable) For MINT = 3, either MITER = 0 or 3, or */
/*                  IMPL = 0 . */
/*              30  (Recoverable) For MITER = 0, IMPL is not 0 . */
/*              31  (Recoverable) For MINT = 1, IMPL is 2 or 3 . */
/*              32  (Recoverable) Insufficient storage has been allocated */
/*                  for the WORK array. */
/*              33  (Recoverable) Insufficient storage has been allocated */
/*                  for the IWORK array. */
/*              41  (Recoverable) The integration step size has gone */
/*                  to zero. */
/*              42  (Recoverable) The integration step size has been */
/*                  reduced about 50 times without advancing the */
/*                  solution.  The problem setup may not be correct. */
/*              43  (Recoverable)  For IMPL greater than 0, the matrix A */
/*                  is singular. */
/*             999  (Fatal) The value of NSTATE is 12 . */

/*  III.  OTHER COMMUNICATION TO THE USER  .............................. */

/*    A. The solver communicates to the user through the parameters */
/*       above.  In addition it writes diagnostic messages through the */
/*       standard error handling program XERMSG.  A complete description */
/*       of XERMSG is given in "Guide to the SLATEC Common Mathematical */
/*       Library" by Kirby W. Fong et al..  At installations which do not */
/*       have this error handling package the short but serviceable */
/*       routine, XERMSG, available with this package, can be used.  That */
/*       program uses the file named OUTPUT to transmit messages. */

/*    B. The first three elements of WORK and the first five elements of */
/*       IWORK will contain the following statistical data: */
/*         AVGH     The average step size used. */
/*         HUSED    The step size last used (successfully). */
/*         AVGORD   The average order used. */
/*         IMXERR   The index of the element of the solution vector that */
/*                  contributed most to the last error test. */
/*         NQUSED   The order last used (successfully). */
/*         NSTEP    The number of steps taken since last initialization. */
/*         NFE      The number of evaluations of the right hand side. */
/*         NJE      The number of evaluations of the Jacobian matrix. */

/*  IV.  REMARKS  ....................................................... */

/*    A. Other routines used: */
/*         DDNTP, DDZRO, DDSTP, DDNTL, DDPST, DDCOR, DDCST, */
/*         DDPSC, and DDSCL; */
/*         DGEFA, DGESL, DGBFA, DGBSL, and DNRM2 (from LINPACK) */
/*         D1MACH (from the Bell Laboratories Machine Constants Package) */
/*         XERMSG (from the SLATEC Common Math Library) */
/*       The last seven routines above, not having been written by the */
/*       present authors, are not explicitly part of this package. */

/*    B. On any return from DDRIV3 all information necessary to continue */
/*       the calculation is contained in the call sequence parameters, */
/*       including the work arrays.  Thus it is possible to suspend one */
/*       problem, integrate another, and then return to the first. */

/*    C. If this package is to be used in an overlay situation, the user */
/*       must declare in the primary overlay the variables in the call */
/*       sequence to DDRIV3. */

/*    D. Changing parameters during an integration. */
/*       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may */
/*       be altered by the user between calls to DDRIV3.  For example, if */
/*       too much accuracy has been requested (the program returns with */
/*       NSTATE = 4 and an increased value of EPS) the user may wish to */
/*       increase EPS further.  In general, prudence is necessary when */
/*       making changes in parameters since such changes are not */
/*       implemented until the next integration step, which is not */
/*       necessarily the next call to DDRIV3.  This can happen if the */
/*       program has already integrated to a point which is beyond the */
/*       new point TOUT. */

/*    E. As the price for complete control of matrix algebra, the DDRIV3 */
/*       USERS option puts all responsibility for Jacobian matrix */
/*       evaluation on the user.  It is often useful to approximate */
/*       numerically all or part of the Jacobian matrix.  However this */
/*       must be done carefully.  The FORTRAN sequence below illustrates */
/*       the method we recommend.  It can be inserted directly into */
/*       subroutine USERS to approximate Jacobian elements in rows I1 */
/*       to I2 and columns J1 to J2. */
/*              DOUBLE PRECISION DFDY(N,N), EPSJ, H, R, D1MACH, */
/*             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N) */
/*              UROUND = D1MACH(4) */
/*              EPSJ = SQRT(UROUND) */
/*              DO 30 J = J1,J2 */
/*                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J))) */
/*                IF (R .EQ. 0.D0) R = YWT(J) */
/*                YJ = Y(J) */
/*                Y(J) = Y(J) + R */
/*                CALL F (N, T, Y, SAVE1) */
/*                IF (N .EQ. 0) RETURN */
/*                Y(J) = YJ */
/*                DO 20 I = I1,I2 */
/*         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R */
/*         30     CONTINUE */
/*       Many problems give rise to structured sparse Jacobians, e.g., */
/*       block banded.  It is possible to approximate them with fewer */
/*       function evaluations than the above procedure uses; see Curtis, */
/*       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13, */
/*       pp. 117-119. */

/*    F. When any of the routines JACOBN, FA, G, or USERS, is not */
/*       required, difficulties associated with unsatisfied externals can */
/*       be avoided by using the name of the routine which calculates the */
/*       right hand side of the differential equations in place of the */
/*       corresponding name in the call sequence of DDRIV3. */

/* ***REFERENCES  C. W. Gear, Numerical Initial Value Problems in */
/*                 Ordinary Differential Equations, Prentice-Hall, 1971. */
/* ***ROUTINES CALLED  D1MACH, DDNTP, DDSTP, DDZRO, DGBFA, DGBSL, DGEFA, */
/*                    DGESL, DNRM2, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  DDRIV3 */
/* ***FIRST EXECUTABLE STATEMENT  DDRIV3 */
    /* Parameter adjustments */
    --iwork;
    --work;
    --ewt;
    --y;

    /* Function Body */
    if (*nstate == 12) {
	*ierflg = 999;
	xermsg_("SLATEC", "DDRIV3", "Illegal input.  The value of NSTATE is "
		"12 .", ierflg, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
	return 0;
    } else if (*nstate < 1 || *nstate > 12) {
	s_wsfi(&io___2);
	do_fio(&c__1, (char *)&(*nstate), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 26;
/* Writing concatenation */
	i__1[0] = 44, a__1[0] = "Illegal input.  Improper value for NSTATE(= "
		;
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 2, a__1[2] = ").";
	s_cat(ch__1, a__1, i__1, &c__3, (ftnlen)54);
	xermsg_("SLATEC", "DDRIV3", ch__1, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)54);
	*nstate = 12;
	return 0;
    }
    npar = *n;
    if (*eps < 0.) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*eps), (ftnlen)sizeof(doublereal));
	e_wsfi();
	*ierflg = 27;
/* Writing concatenation */
	i__1[0] = 21, a__1[0] = "Illegal input.  EPS, ";
	i__1[1] = 16, a__1[1] = rl1;
	i__1[2] = 14, a__1[2] = ", is negative.";
	s_cat(ch__2, a__1, i__1, &c__3, (ftnlen)51);
	xermsg_("SLATEC", "DDRIV3", ch__2, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)51);
	*nstate = 12;
	return 0;
    }
    if (*n <= 0) {
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 22;
/* Writing concatenation */
	i__1[0] = 37, a__1[0] = "Illegal input.  Number of equations, ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 18, a__1[2] = ", is not positive.";
	s_cat(ch__3, a__1, i__1, &c__3, (ftnlen)63);
	xermsg_("SLATEC", "DDRIV3", ch__3, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)63);
	*nstate = 12;
	return 0;
    }
    if (*mxord <= 0) {
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*mxord), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 28;
/* Writing concatenation */
	i__1[0] = 31, a__1[0] = "Illegal input.  Maximum order, ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 18, a__1[2] = ", is not positive.";
	s_cat(ch__4, a__1, i__1, &c__3, (ftnlen)57);
	xermsg_("SLATEC", "DDRIV3", ch__4, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)57);
	*nstate = 12;
	return 0;
    }
    if (*mint < 1 || *mint > 3) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&(*mint), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 23;
/* Writing concatenation */
	i__1[0] = 64, a__1[0] = "Illegal input.  Improper value for the inte"
		"gration method flag, ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 2, a__1[2] = " .";
	s_cat(ch__5, a__1, i__1, &c__3, (ftnlen)74);
	xermsg_("SLATEC", "DDRIV3", ch__5, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)74);
	*nstate = 12;
	return 0;
    } else if (*miter < 0 || *miter > 5) {
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&(*miter), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 24;
/* Writing concatenation */
	i__1[0] = 43, a__1[0] = "Illegal input.  Improper value for MITER(= ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 2, a__1[2] = ").";
	s_cat(ch__6, a__1, i__1, &c__3, (ftnlen)53);
	xermsg_("SLATEC", "DDRIV3", ch__6, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)53);
	*nstate = 12;
	return 0;
    } else if (*impl < 0 || *impl > 3) {
	s_wsfi(&io___10);
	do_fio(&c__1, (char *)&(*impl), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 25;
/* Writing concatenation */
	i__1[0] = 42, a__1[0] = "Illegal input.  Improper value for IMPL(= ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 2, a__1[2] = ").";
	s_cat(ch__7, a__1, i__1, &c__3, (ftnlen)52);
	xermsg_("SLATEC", "DDRIV3", ch__7, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)52);
	*nstate = 12;
	return 0;
    } else if (*mint == 3 && (*miter == 0 || *miter == 3 || *impl != 0)) {
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&(*miter), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___13);
	do_fio(&c__1, (char *)&(*impl), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 29;
/* Writing concatenation */
	i__2[0] = 50, a__2[0] = "Illegal input.  For MINT = 3, the value of "
		"MITER, ";
	i__2[1] = 8, a__2[1] = intgr1;
	i__2[2] = 15, a__2[2] = ", and/or IMPL, ";
	i__2[3] = 8, a__2[3] = intgr2;
	i__2[4] = 17, a__2[4] = ", is not allowed.";
	s_cat(ch__8, a__2, i__2, &c__5, (ftnlen)98);
	xermsg_("SLATEC", "DDRIV3", ch__8, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)98);
	*nstate = 12;
	return 0;
    } else if (*impl >= 1 && *impl <= 3 && *miter == 0) {
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*impl), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 30;
/* Writing concatenation */
	i__1[0] = 50, a__1[0] = "Illegal input.  For MITER = 0, the value of"
		" IMPL, ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 17, a__1[2] = ", is not allowed.";
	s_cat(ch__9, a__1, i__1, &c__3, (ftnlen)75);
	xermsg_("SLATEC", "DDRIV3", ch__9, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)75);
	*nstate = 12;
	return 0;
    } else if ((*impl == 2 || *impl == 3) && *mint == 1) {
	s_wsfi(&io___15);
	do_fio(&c__1, (char *)&(*impl), (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 31;
/* Writing concatenation */
	i__1[0] = 49, a__1[0] = "Illegal input.  For MINT = 1, the value of "
		"IMPL, ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 17, a__1[2] = ", is not allowed.";
	s_cat(ch__5, a__1, i__1, &c__3, (ftnlen)74);
	xermsg_("SLATEC", "DDRIV3", ch__5, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)74);
	*nstate = 12;
	return 0;
    }
    if (*miter == 0 || *miter == 3) {
	liwchk = 50;
    } else if (*miter == 1 || *miter == 2 || *miter == 4 || *miter == 5) {
	liwchk = *n + 50;
    }
    if (*leniw < liwchk) {
	s_wsfi(&io___17);
	do_fio(&c__1, (char *)&liwchk, (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 33;
/* Writing concatenation */
	i__1[0] = 146, a__1[0] = "Illegal input.  Insufficient storage alloc"
		"ated for the IWORK array.  Based on the value of the input p"
		"arameters involved, the required storage is ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 2, a__1[2] = " .";
	s_cat(ch__10, a__1, i__1, &c__3, (ftnlen)156);
	xermsg_("SLATEC", "DDRIV3", ch__10, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)156);
	*nstate = 12;
	return 0;
    }
/*                                                Allocate the WORK array */
/*                                         IYH is the index of YH in WORK */
    if (*mint == 1 || *mint == 3) {
	maxord = min(*mxord,12);
    } else if (*mint == 2) {
	maxord = min(*mxord,5);
    }
    idfdy = (maxord + 1) * *n + 251;
/*                                             IDFDY is the index of DFDY */

    if (*miter == 0 || *miter == 3) {
	iywt = idfdy;
    } else if (*miter == 1 || *miter == 2) {
	iywt = idfdy + *n * *n;
    } else if (*miter == 4 || *miter == 5) {
	iywt = idfdy + ((*ml << 1) + *mu + 1) * *n;
    }
/*                                               IYWT is the index of YWT */
    isave1 = iywt + *n;
/*                                           ISAVE1 is the index of SAVE1 */
    isave2 = isave1 + *n;
/*                                           ISAVE2 is the index of SAVE2 */
    ignow = isave2 + *n;
/*                                             IGNOW is the index of GNOW */
    itroot = ignow + *nroot;
/*                                           ITROOT is the index of TROOT */
    ifac = itroot + *nroot;
/*                                               IFAC is the index of FAC */
    if (*miter == 2 || *miter == 5 || *mint == 3) {
	ia = ifac + *n;
    } else {
	ia = ifac;
    }
/*                                                   IA is the index of A */
    if (*impl == 0 || *miter == 3) {
	lenchk = ia - 1;
    } else if (*impl == 1 && (*miter == 1 || *miter == 2)) {
	lenchk = ia - 1 + *n * *n;
    } else if (*impl == 1 && (*miter == 4 || *miter == 5)) {
	lenchk = ia - 1 + ((*ml << 1) + *mu + 1) * *n;
    } else if (*impl == 2 && *miter != 3) {
	lenchk = ia - 1 + *n;
    } else if (*impl == 3 && (*miter == 1 || *miter == 2)) {
	lenchk = ia - 1 + *n * *nde;
    } else if (*impl == 3 && (*miter == 4 || *miter == 5)) {
	lenchk = ia - 1 + ((*ml << 1) + *mu + 1) * *nde;
    }
    if (*lenw < lenchk) {
	s_wsfi(&io___28);
	do_fio(&c__1, (char *)&lenchk, (ftnlen)sizeof(integer));
	e_wsfi();
	*ierflg = 32;
/* Writing concatenation */
	i__1[0] = 145, a__1[0] = "Illegal input.  Insufficient storage alloc"
		"ated for the WORK array.  Based on the value of the input pa"
		"rameters involved, the required storage is ";
	i__1[1] = 8, a__1[1] = intgr1;
	i__1[2] = 2, a__1[2] = " .";
	s_cat(ch__11, a__1, i__1, &c__3, (ftnlen)155);
	xermsg_("SLATEC", "DDRIV3", ch__11, ierflg, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)155);
	*nstate = 12;
	return 0;
    }
    if (*miter == 0 || *miter == 3) {
	matdim = 1;
    } else if (*miter == 1 || *miter == 2) {
	matdim = *n;
    } else if (*miter == 4 || *miter == 5) {
	matdim = (*ml << 1) + *mu + 1;
    }
    if (*impl == 0 || *impl == 1) {
	ndecom = *n;
    } else if (*impl == 2 || *impl == 3) {
	ndecom = *nde;
    }
    if (*nstate == 1) {
/*                                                  Initialize parameters */
	if (*mint == 1 || *mint == 3) {
	    iwork[20] = min(*mxord,12);
	} else if (*mint == 2) {
	    iwork[20] = min(*mxord,5);
	}
	iwork[19] = *mxord;
	if (*mint == 1 || *mint == 2) {
	    iwork[16] = *mint;
	    iwork[18] = *miter;
	    iwork[10] = *mint;
	    iwork[11] = *miter;
	} else if (*mint == 3) {
	    iwork[16] = 1;
	    iwork[18] = 0;
	    iwork[10] = iwork[16];
	    iwork[11] = iwork[18];
	    iwork[17] = *miter;
	}
	work[161] = *hmax;
	uround = d1mach_(&c__4);
	work[206] = uround;
	work[205] = d1mach_(&c__1);
	if (*nroot != 0) {
	    re = uround;
	    ae = work[205];
	}
	h__ = (*tout - *t) * (1. - uround * 4.);
/* Computing MIN */
	d__2 = abs(h__);
	d__1 = min(d__2,*hmax);
	h__ = d_sign(&d__1, &h__);
	work[160] = h__;
	hsign = d_sign(&c_b107, &h__);
	work[163] = hsign;
	iwork[9] = 0;
	work[1] = 0.;
	work[2] = 0.;
	work[3] = 0.;
	iwork[1] = 0;
	iwork[2] = 0;
	iwork[3] = 0;
	iwork[22] = 0;
	iwork[4] = 0;
	iwork[5] = 0;
	iwork[6] = 0;
	work[166] = *t;
	iwork[7] = 0;
	iwork[21] = 0;
/*                                                 Set initial conditions */
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L30: */
	    work[i__ + 250] = y[i__];
	}
	if (*t == *tout) {
	    return 0;
	}
	goto L180;
    } else {
	uround = work[206];
	if (*nroot != 0) {
	    re = uround;
	    ae = work[205];
	}
    }
/*                                             On a continuation, check */
/*                                             that output points have */
/*                                             been or will be overtaken. */
    if (iwork[7] == 1) {
	convrg = TRUE_;
    } else {
	convrg = FALSE_;
    }
    *t = work[166];
    h__ = work[160];
    hsign = work[163];
    if (iwork[9] == 0) {
	goto L180;
    }

/*                                   IWORK(IJROOT) flags unreported */
/*                                   roots, and is set to the value of */
/*                                   NTASK when a root was last selected. */
/*                                   It is set to zero when all roots */
/*                                   have been reported.  IWORK(INROOT) */
/*                                   contains the index and WORK(ITOUT) */
/*                                   contains the value of the root last */
/*                                   selected to be reported. */
/*                                   IWORK(INRTLD) contains the value of */
/*                                   NROOT and IWORK(INDTRT) contains */
/*                                   the value of ITROOT when the array */
/*                                   of roots was last calculated. */
    if (*nroot != 0) {
	if (iwork[8] > 0) {
/*                                      TOUT has just been reported. */
/*                                      If TROOT .LE. TOUT, report TROOT. */
	    if (*nstate != 5) {
		if (*tout * hsign >= work[167] * hsign) {
		    troot = work[167];
		    ddntp_(&h__, &c__0, n, &iwork[12], t, &troot, &work[251], 
			    &y[1]);
		    *t = troot;
		    *nstate = 5;
		    *ierflg = 0;
		    goto L580;
		}
/*                                         A root has just been reported. */
/*                                         Select the next root. */
	    } else {
		troot = *t;
		iroot = 0;
		i__3 = iwork[13];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    jtroot = i__ + iwork[14] - 1;
		    if (work[jtroot] * hsign <= troot * hsign) {

/*                                              Check for multiple roots. */

			if (work[jtroot] == work[167] && i__ > iwork[6]) {
			    iroot = i__;
			    troot = work[jtroot];
			    goto L60;
			}
			if (work[jtroot] * hsign > work[167] * hsign) {
			    iroot = i__;
			    troot = work[jtroot];
			}
		    }
/* L50: */
		}
L60:
		iwork[6] = iroot;
		work[167] = troot;
		iwork[8] = *ntask;
		if (*ntask == 1) {
		    if (iroot == 0) {
			iwork[8] = 0;
		    } else {
			if (*tout * hsign >= troot * hsign) {
			    ddntp_(&h__, &c__0, n, &iwork[12], t, &troot, &
				    work[251], &y[1]);
			    *nstate = 5;
			    *t = troot;
			    *ierflg = 0;
			    goto L580;
			}
		    }
		} else if (*ntask == 2 || *ntask == 3) {

/*                                     If there are no more roots, or the */
/*                                     user has altered TOUT to be less */
/*                                     than a root, set IJROOT to zero. */

		    if (iroot == 0 || *tout * hsign < troot * hsign) {
			iwork[8] = 0;
		    } else {
			ddntp_(&h__, &c__0, n, &iwork[12], t, &troot, &work[
				251], &y[1]);
			*nstate = 5;
			*ierflg = 0;
			*t = troot;
			goto L580;
		    }
		}
	    }
	}
    }

    if (*ntask == 1) {
	*nstate = 2;
	if (*t * hsign >= *tout * hsign) {
	    ddntp_(&h__, &c__0, n, &iwork[12], t, tout, &work[251], &y[1]);
	    *t = *tout;
	    *ierflg = 0;
	    goto L580;
	}
    } else if (*ntask == 2) {
/*                                                      Check if TOUT has */
/*                                                      been reset .LT. T */
	if (*t * hsign > *tout * hsign) {
	    s_wsfi(&io___41);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___43);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    *ierflg = 11;
/* Writing concatenation */
	    i__4[0] = 38, a__3[0] = "While integrating exactly to TOUT, T, ";
	    i__4[1] = 16, a__3[1] = rl1;
	    i__4[2] = 19, a__3[2] = ", was beyond TOUT, ";
	    i__4[3] = 16, a__3[3] = rl2;
	    i__4[4] = 25, a__3[4] = " .  Solution obtained by ";
	    i__4[5] = 14, a__3[5] = "interpolation.";
	    s_cat(ch__12, a__3, i__4, &c__6, (ftnlen)128);
	    xermsg_("SLATEC", "DDRIV3", ch__12, ierflg, &c__0, (ftnlen)6, (
		    ftnlen)6, (ftnlen)128);
	    *nstate = 11;
	    ddntp_(&h__, &c__0, n, &iwork[12], t, tout, &work[251], &y[1]);
	    *t = *tout;
	    goto L580;
	}
/*                                   Determine if TOUT has been overtaken */

/* Computing MAX */
	d__2 = abs(*t), d__3 = abs(*tout);
	if ((d__1 = *tout - *t, abs(d__1)) <= uround * 20. * max(d__2,d__3)) {
	    *t = *tout;
	    *nstate = 2;
	    *ierflg = 0;
	    goto L560;
	}
/*                                             If there are no more roots */
/*                                             to report, report T. */
	if (*nstate == 5) {
	    *nstate = 2;
	    *ierflg = 0;
	    goto L560;
	}
	*nstate = 2;
/*                                                       See if TOUT will */
/*                                                       be overtaken. */
	if ((*t + h__) * hsign > *tout * hsign) {
	    h__ = *tout - *t;
	    if ((*t + h__) * hsign > *tout * hsign) {
		h__ *= 1. - uround * 4.;
	    }
	    work[160] = h__;
	    if (h__ == 0.) {
		goto L670;
	    }
	    iwork[9] = -1;
	}
    } else if (*ntask == 3) {
	*nstate = 2;
	if (*t * hsign > *tout * hsign) {
	    s_wsfi(&io___44);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___45);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    *ierflg = 11;
/* Writing concatenation */
	    i__4[0] = 38, a__3[0] = "While integrating exactly to TOUT, T, ";
	    i__4[1] = 16, a__3[1] = rl1;
	    i__4[2] = 19, a__3[2] = ", was beyond TOUT, ";
	    i__4[3] = 16, a__3[3] = rl2;
	    i__4[4] = 25, a__3[4] = " .  Solution obtained by ";
	    i__4[5] = 14, a__3[5] = "interpolation.";
	    s_cat(ch__12, a__3, i__4, &c__6, (ftnlen)128);
	    xermsg_("SLATEC", "DDRIV3", ch__12, ierflg, &c__0, (ftnlen)6, (
		    ftnlen)6, (ftnlen)128);
	    *nstate = 11;
	    ddntp_(&h__, &c__0, n, &iwork[12], t, tout, &work[251], &y[1]);
	    *t = *tout;
	    goto L580;
	}
/* Computing MAX */
	d__2 = abs(*t), d__3 = abs(*tout);
	if ((d__1 = *tout - *t, abs(d__1)) <= uround * 20. * max(d__2,d__3)) {
	    *t = *tout;
	    *ierflg = 0;
	    goto L560;
	}
	if ((*t + h__) * hsign > *tout * hsign) {
	    h__ = *tout - *t;
	    if ((*t + h__) * hsign > *tout * hsign) {
		h__ *= 1. - uround * 4.;
	    }
	    work[160] = h__;
	    if (h__ == 0.) {
		goto L670;
	    }
	    iwork[9] = -1;
	}
    }
/*                         Implement changes in MINT, MITER, and/or HMAX. */

    if ((*mint != iwork[10] || *miter != iwork[11]) && *mint != 3 && iwork[10]
	     != 3) {
	iwork[9] = -1;
    }
    if (*hmax != work[161]) {
/* Computing MIN */
	d__2 = abs(h__);
	d__1 = min(d__2,*hmax);
	h__ = d_sign(&d__1, &h__);
	if (h__ != work[160]) {
	    iwork[9] = -1;
	    work[160] = h__;
	}
	work[161] = *hmax;
    }

L180:
    nstepl = iwork[3];
    i__3 = *n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L190: */
	y[i__] = work[i__ + 250];
    }
    if (*nroot != 0) {
	i__3 = *nroot;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    work[i__ + ignow - 1] = (*g)(&npar, t, &y[1], &i__);
	    if (npar == 0) {
		iwork[6] = i__;
		*nstate = 7;
		return 0;
	    }
/* L200: */
	}
    }
    if (*ierror == 1) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L230: */
	    work[i__ + iywt - 1] = 1.;
	}
	goto L410;
    } else if (*ierror == 5) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L250: */
	    work[i__ + iywt - 1] = ewt[i__];
	}
	goto L410;
    }
/*                                       Reset YWT array.  Looping point. */
L260:
    if (*ierror == 2) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (y[i__] == 0.) {
		goto L290;
	    }
/* L280: */
	    work[i__ + iywt - 1] = (d__1 = y[i__], abs(d__1));
	}
	goto L410;
L290:
	if (iwork[9] == 0) {
	    (*f)(&npar, t, &y[1], &work[isave2]);
	    if (npar == 0) {
		*nstate = 6;
		return 0;
	    }
	    ++iwork[4];
	    if (*miter == 3 && *impl != 0) {
		iflag = 0;
		(*users)(&y[1], &work[251], &work[iywt], &work[isave1], &work[
			isave2], t, &h__, &work[4], impl, &npar, &ndecom, &
			iflag);
		if (iflag == -1) {
		    goto L690;
		}
		if (npar == 0) {
		    *nstate = 10;
		    return 0;
		}
	    } else if (*impl == 1) {
		if (*miter == 1 || *miter == 2) {
		    (*fa)(&npar, t, &y[1], &work[ia], &matdim, ml, mu, &
			    ndecom);
		    if (npar == 0) {
			*nstate = 9;
			return 0;
		    }
		    dgefa_(&work[ia], &matdim, n, &iwork[51], &info);
		    if (info != 0) {
			goto L690;
		    }
		    dgesl_(&work[ia], &matdim, n, &iwork[51], &work[isave2], &
			    c__0);
		} else if (*miter == 4 || *miter == 5) {
		    (*fa)(&npar, t, &y[1], &work[ia + *ml], &matdim, ml, mu, &
			    ndecom);
		    if (npar == 0) {
			*nstate = 9;
			return 0;
		    }
		    dgbfa_(&work[ia], &matdim, n, ml, mu, &iwork[51], &info);
		    if (info != 0) {
			goto L690;
		    }
		    dgbsl_(&work[ia], &matdim, n, ml, mu, &iwork[51], &work[
			    isave2], &c__0);
		}
	    } else if (*impl == 2) {
		(*fa)(&npar, t, &y[1], &work[ia], &matdim, ml, mu, &ndecom);
		if (npar == 0) {
		    *nstate = 9;
		    return 0;
		}
		i__3 = ndecom;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (work[i__ + ia - 1] == 0.) {
			goto L690;
		    }
/* L340: */
		    work[i__ + isave2 - 1] /= work[i__ + ia - 1];
		}
	    } else if (*impl == 3) {
		if (*miter == 1 || *miter == 2) {
		    (*fa)(&npar, t, &y[1], &work[ia], &matdim, ml, mu, &
			    ndecom);
		    if (npar == 0) {
			*nstate = 9;
			return 0;
		    }
		    dgefa_(&work[ia], &matdim, nde, &iwork[51], &info);
		    if (info != 0) {
			goto L690;
		    }
		    dgesl_(&work[ia], &matdim, nde, &iwork[51], &work[isave2],
			     &c__0);
		} else if (*miter == 4 || *miter == 5) {
		    (*fa)(&npar, t, &y[1], &work[ia + *ml], &matdim, ml, mu, &
			    ndecom);
		    if (npar == 0) {
			*nstate = 9;
			return 0;
		    }
		    dgbfa_(&work[ia], &matdim, nde, ml, mu, &iwork[51], &info)
			    ;
		    if (info != 0) {
			goto L690;
		    }
		    dgbsl_(&work[ia], &matdim, nde, ml, mu, &iwork[51], &work[
			    isave2], &c__0);
		}
	    }
	}
	i__3 = *n;
	for (j = i__; j <= i__3; ++j) {
	    if (y[j] != 0.) {
		work[j + iywt - 1] = (d__1 = y[j], abs(d__1));
	    } else {
		if (iwork[9] == 0) {
		    work[j + iywt - 1] = (d__1 = h__ * work[j + isave2 - 1], 
			    abs(d__1));
		} else {
		    work[j + iywt - 1] = (d__1 = work[j + 251 + *n - 1], abs(
			    d__1));
		}
	    }
	    if (work[j + iywt - 1] == 0.) {
		work[j + iywt - 1] = uround;
	    }
/* L360: */
	}
    } else if (*ierror == 3) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L380: */
/* Computing MAX */
	    d__2 = ewt[1], d__3 = (d__1 = y[i__], abs(d__1));
	    work[i__ + iywt - 1] = max(d__2,d__3);
	}
    } else if (*ierror == 4) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L400: */
/* Computing MAX */
	    d__2 = ewt[i__], d__3 = (d__1 = y[i__], abs(d__1));
	    work[i__ + iywt - 1] = max(d__2,d__3);
	}
    }

L410:
    i__3 = *n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L420: */
	work[i__ + isave2 - 1] = y[i__] / work[i__ + iywt - 1];
    }
    sum = dnrm2_(n, &work[isave2], &c__1) / sqrt((doublereal) (*n));
    sum = max(1.,sum);
    if (*eps < sum * uround) {
	*eps = sum * uround * (uround * 10. + 1.);
	s_wsfi(&io___51);
	do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	e_wsfi();
	s_wsfi(&io___52);
	do_fio(&c__1, (char *)&(*eps), (ftnlen)sizeof(doublereal));
	e_wsfi();
	*ierflg = 4;
/* Writing concatenation */
	i__5[0] = 6, a__4[0] = "At T, ";
	i__5[1] = 16, a__4[1] = rl1;
	i__5[2] = 39, a__4[2] = ", the requested accuracy, EPS, was not ";
	i__5[3] = 53, a__4[3] = "obtainable with the machine precision.  EPS"
		" has been ";
	i__5[4] = 13, a__4[4] = "increased to ";
	i__5[5] = 16, a__4[5] = rl2;
	i__5[6] = 2, a__4[6] = " .";
	s_cat(ch__13, a__4, i__5, &c__7, (ftnlen)145);
	xermsg_("SLATEC", "DDRIV3", ch__13, ierflg, &c__0, (ftnlen)6, (ftnlen)
		6, (ftnlen)145);
	*nstate = 4;
	goto L560;
    }
    if (abs(h__) >= uround * abs(*t)) {
	iwork[21] = 0;
    } else if (iwork[21] == 0) {
	s_wsfi(&io___53);
	do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	e_wsfi();
	s_wsfi(&io___54);
	do_fio(&c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
	e_wsfi();
	*ierflg = 15;
/* Writing concatenation */
	i__6[0] = 6, a__5[0] = "At T, ";
	i__6[1] = 16, a__5[1] = rl1;
	i__6[2] = 17, a__5[2] = ", the step size, ";
	i__6[3] = 16, a__5[3] = rl2;
	i__6[4] = 13, a__5[4] = ", is smaller ";
	i__6[5] = 58, a__5[5] = "than the roundoff level of T.  This may occ"
		"ur if there is ";
	i__6[6] = 47, a__5[6] = "an abrupt change in the right hand side of "
		"the ";
	i__6[7] = 23, a__5[7] = "differential equations.";
	s_cat(ch__14, a__5, i__6, &c__8, (ftnlen)196);
	xermsg_("SLATEC", "DDRIV3", ch__14, ierflg, &c__0, (ftnlen)6, (ftnlen)
		6, (ftnlen)196);
	iwork[21] = 1;
    }
    if (*ntask != 2) {
	if (iwork[3] - nstepl == *mxstep) {
	    s_wsfi(&io___55);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___56);
	    do_fio(&c__1, (char *)&(*mxstep), (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___57);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    *ierflg = 3;
/* Writing concatenation */
	    i__6[0] = 6, a__5[0] = "At T, ";
	    i__6[1] = 16, a__5[1] = rl1;
	    i__6[2] = 2, a__5[2] = ", ";
	    i__6[3] = 8, a__5[3] = intgr1;
	    i__6[4] = 23, a__5[4] = " steps have been taken ";
	    i__6[5] = 23, a__5[5] = "without reaching TOUT, ";
	    i__6[6] = 16, a__5[6] = rl2;
	    i__6[7] = 2, a__5[7] = " .";
	    s_cat(ch__15, a__5, i__6, &c__8, (ftnlen)96);
	    xermsg_("SLATEC", "DDRIV3", ch__15, ierflg, &c__0, (ftnlen)6, (
		    ftnlen)6, (ftnlen)96);
	    *nstate = 3;
	    goto L560;
	}
    }

/*     CALL DDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM, */
/*    8            MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND, */
/*    8            USERS,  AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD, */
/*    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG, */
/*    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, JSTEPL, NQ, NWAIT, */
/*    8            RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, */
/*    8            MXRDSV) */

    ddstp_(eps, (S_fp)f, (S_fp)fa, &work[161], impl, ierror, (U_fp)jacobn, &
	    matdim, &iwork[20], &iwork[16], &iwork[18], ml, mu, &npar, &
	    ndecom, &work[iywt], &uround, (S_fp)users, &work[1], &work[3], &
	    work[160], &hused, &iwork[9], &iwork[10], &iwork[11], &iwork[4], &
	    iwork[5], &iwork[2], &iwork[3], &work[166], &y[1], &work[251], &
	    work[ia], &convrg, &work[idfdy], &work[4], &work[ifac], &work[162]
	    , &iwork[51], &jstate, &iwork[22], &iwork[12], &iwork[15], &work[
	    164], &work[165], &work[isave1], &work[isave2], &work[168], &work[
	    204], mint, &iwork[17], &iwork[19]);
    *t = work[166];
    h__ = work[160];
    if (convrg) {
	iwork[7] = 1;
    } else {
	iwork[7] = 0;
    }
    switch (jstate) {
	case 1:  goto L470;
	case 2:  goto L670;
	case 3:  goto L680;
	case 4:  goto L690;
	case 5:  goto L690;
	case 6:  goto L660;
	case 7:  goto L660;
	case 8:  goto L660;
	case 9:  goto L660;
	case 10:  goto L660;
    }
L470:
    iwork[9] = 1;
/*                                 Determine if a root has been overtaken */
    if (*nroot != 0) {
	iroot = 0;
	i__3 = *nroot;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    glast = work[i__ + ignow - 1];
	    gnow = (*g)(&npar, t, &y[1], &i__);
	    if (npar == 0) {
		iwork[6] = i__;
		*nstate = 7;
		return 0;
	    }
	    work[i__ + ignow - 1] = gnow;
	    if (glast * gnow > 0.) {
		work[i__ + itroot - 1] = *t + h__;
	    } else {
		if (gnow == 0.) {
		    work[i__ + itroot - 1] = *t;
		    iroot = i__;
		} else {
		    if (glast == 0.) {
			work[i__ + itroot - 1] = *t + h__;
		    } else {
			if (abs(hused) >= uround * abs(*t)) {
			    tlast = *t - hused;
			    iroot = i__;
			    troot = *t;
			    ddzro_(&ae, (D_fp)g, &h__, &npar, &iwork[12], &
				    iroot, &re, t, &work[251], &uround, &
				    troot, &tlast, &gnow, &glast, &y[1]);
			    i__7 = *n;
			    for (j = 1; j <= i__7; ++j) {
/* L480: */
				y[j] = work[j + 250];
			    }
			    if (npar == 0) {
				iwork[6] = i__;
				*nstate = 7;
				return 0;
			    }
			    work[i__ + itroot - 1] = troot;
			} else {
			    work[i__ + itroot - 1] = *t;
			    iroot = i__;
			}
		    }
		}
	    }
/* L500: */
	}
	if (iroot == 0) {
	    iwork[8] = 0;
/*                                                  Select the first root */
	} else {
	    iwork[8] = *ntask;
	    iwork[13] = *nroot;
	    iwork[14] = itroot;
	    troot = *t + h__;
	    i__3 = *nroot;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		if (work[i__ + itroot - 1] * hsign < troot * hsign) {
		    troot = work[i__ + itroot - 1];
		    iroot = i__;
		}
/* L510: */
	    }
	    iwork[6] = iroot;
	    work[167] = troot;
	    if (troot * hsign <= *tout * hsign) {
		ddntp_(&h__, &c__0, n, &iwork[12], t, &troot, &work[251], &y[
			1]);
		*nstate = 5;
		*t = troot;
		*ierflg = 0;
		goto L580;
	    }
	}
    }
/*                               Test for NTASK condition to be satisfied */
    *nstate = 2;
    if (*ntask == 1) {
	if (*t * hsign < *tout * hsign) {
	    goto L260;
	}
	ddntp_(&h__, &c__0, n, &iwork[12], t, tout, &work[251], &y[1]);
	*t = *tout;
	*ierflg = 0;
	goto L580;
/*                               TOUT is assumed to have been attained */
/*                               exactly if T is within twenty roundoff */
/*                               units of TOUT, relative to MAX(TOUT, T). */

    } else if (*ntask == 2) {
/* Computing MAX */
	d__2 = abs(*t), d__3 = abs(*tout);
	if ((d__1 = *tout - *t, abs(d__1)) <= uround * 20. * max(d__2,d__3)) {
	    *t = *tout;
	} else {
	    if ((*t + h__) * hsign > *tout * hsign) {
		h__ = *tout - *t;
		if ((*t + h__) * hsign > *tout * hsign) {
		    h__ *= 1. - uround * 4.;
		}
		work[160] = h__;
		if (h__ == 0.) {
		    goto L670;
		}
		iwork[9] = -1;
	    }
	}
    } else if (*ntask == 3) {
/* Computing MAX */
	d__2 = abs(*t), d__3 = abs(*tout);
	if ((d__1 = *tout - *t, abs(d__1)) <= uround * 20. * max(d__2,d__3)) {
	    *t = *tout;
	} else {
	    if ((*t + h__) * hsign > *tout * hsign) {
		h__ = *tout - *t;
		if ((*t + h__) * hsign > *tout * hsign) {
		    h__ *= 1. - uround * 4.;
		}
		work[160] = h__;
		if (h__ == 0.) {
		    goto L670;
		}
		iwork[9] = -1;
	    }
	    goto L260;
	}
    }
    *ierflg = 0;
/*                                      All returns are made through this */
/*                                      section.  IMXERR is determined. */
L560:
    i__3 = *n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L570: */
	y[i__] = work[i__ + 250];
    }
L580:
    if (iwork[9] == 0) {
	return 0;
    }
    big = 0.;
    imxerr = 1;
    i__3 = *n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*                                            SIZE = ABS(ERROR(I)/YWT(I)) */
	size = (d__1 = work[i__ + isave1 - 1] / work[i__ + iywt - 1], abs(
		d__1));
	if (big < size) {
	    big = size;
	    imxerr = i__;
	}
/* L590: */
    }
    iwork[1] = imxerr;
    work[2] = hused;
    return 0;

L660:
    *nstate = jstate;
    return 0;
/*                                        Fatal errors are processed here */

L670:
    s_wsfi(&io___66);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfi();
    *ierflg = 41;
/* Writing concatenation */
    i__8[0] = 6, a__6[0] = "At T, ";
    i__8[1] = 16, a__6[1] = rl1;
    i__8[2] = 38, a__6[2] = ", the attempted step size has gone to ";
    i__8[3] = 59, a__6[3] = "zero.  Often this occurs if the problem setup i"
	    "s incorrect.";
    s_cat(ch__16, a__6, i__8, &c__4, (ftnlen)119);
    xermsg_("SLATEC", "DDRIV3", ch__16, ierflg, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)119);
    *nstate = 12;
    return 0;

L680:
    s_wsfi(&io___67);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfi();
    *ierflg = 42;
/* Writing concatenation */
    i__2[0] = 6, a__2[0] = "At T, ";
    i__2[1] = 16, a__2[1] = rl1;
    i__2[2] = 42, a__2[2] = ", the step size has been reduced about 50 ";
    i__2[3] = 57, a__2[3] = "times without advancing the solution.  Often th"
	    "is occurs ";
    i__2[4] = 34, a__2[4] = "if the problem setup is incorrect.";
    s_cat(ch__11, a__2, i__2, &c__5, (ftnlen)155);
    xermsg_("SLATEC", "DDRIV3", ch__11, ierflg, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)155);
    *nstate = 12;
    return 0;

L690:
    s_wsfi(&io___68);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfi();
    *ierflg = 43;
/* Writing concatenation */
    i__1[0] = 6, a__1[0] = "At T, ";
    i__1[1] = 16, a__1[1] = rl1;
    i__1[2] = 42, a__1[2] = ", while solving A*YDOT = F, A is singular.";
    s_cat(ch__17, a__1, i__1, &c__3, (ftnlen)64);
    xermsg_("SLATEC", "DDRIV3", ch__17, ierflg, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)64);
    *nstate = 12;
    return 0;
} /* ddriv3_ */


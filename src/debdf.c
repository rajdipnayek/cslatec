/* debdf.f -- translated by f2c (version 12.02.01).
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
    real told, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iquit, init, iyh, iewt, iacor, isavf, iwm, ksteps, ibegin, itol, 
	    iinteg, itstop, ijac, iband, iowns[6], ier, jstart, kflag, l, 
	    meth, miter, maxord, n, nq, nst, nfe, nje, nqu;
} debdf1_;

#define debdf1_1 debdf1_

/* Table of constant values */

static integer c__1 = 1;
static integer c__13 = 13;
static integer c__2 = 2;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__17 = 17;
static integer c__18 = 18;

/* DECK DEBDF */
/* Subroutine */ int debdf_(U_fp f, integer *neq, real *t, real *y, real *
	tout, integer *info, real *rtol, real *atol, integer *idid, real *
	rwork, integer *lrw, integer *iwork, integer *liw, real *rpar, 
	integer *ipar, U_fp jac)
{
    /* System generated locals */
    address a__1[5], a__2[2], a__3[4];
    integer i__1[5], i__2[2], i__3[4];
    char ch__1[222], ch__2[212], ch__3[133], ch__4[155], ch__5[162], ch__6[
	    219], ch__7[196], ch__8[268], ch__9[109], ch__10[117], ch__11[90];

    /* Local variables */
    static integer ml, mu;
    extern /* Subroutine */ int lsod_(U_fp, integer *, real *, real *, real *,
	     real *, real *, integer *, real *, real *, real *, real *, real *
	    , real *, real *, integer *, U_fp, logical *, real *, real *, 
	    real *, real *, integer *);
    static integer ilrw;
    static char xern1[8], xern2[8], xern3[16];
    static integer icomi, icomr, idelsn;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer itstar, iinout;
    static logical intout;
    static integer iypout;

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___4 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___8 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___15 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___16 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___17 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DEBDF */
/* ***PURPOSE  Solve an initial value problem in ordinary differential */
/*            equations using backward differentiation formulas.  It is */
/*            intended primarily for stiff problems. */
/* ***LIBRARY   SLATEC (DEPAC) */
/* ***CATEGORY  I1A2 */
/* ***TYPE      SINGLE PRECISION (DEBDF-S, DDEBDF-D) */
/* ***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DEPAC, */
/*             INITIAL VALUE PROBLEMS, ODE, */
/*             ORDINARY DIFFERENTIAL EQUATIONS, STIFF */
/* ***AUTHOR  Shampine, L. F., (SNLA) */
/*           Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   This is the backward differentiation code in the package of */
/*   differential equation solvers DEPAC, consisting of the codes */
/*   DERKF, DEABM, and DEBDF.  Design of the package was by */
/*   L. F. Shampine and H. A. Watts.  It is documented in */
/*        SAND-79-2374 , DEPAC - Design of a User Oriented Package of ODE */
/*                              Solvers. */
/*   DEBDF is a driver for a modification of the code LSODE written by */
/*             A. C. Hindmarsh */
/*             Lawrence Livermore Laboratory */
/*             Livermore, California 94550 */

/* ********************************************************************** */
/* **             DEPAC PACKAGE OVERVIEW           ** */
/* ********************************************************************** */

/*        You have a choice of three differential equation solvers from */
/*        DEPAC.  The following brief descriptions are meant to aid you */
/*        in choosing the most appropriate code for your problem. */

/*        DERKF is a fifth order Runge-Kutta code.  It is the simplest of */
/*        the three choices, both algorithmically and in the use of the */
/*        code.  DERKF is primarily designed to solve non-stiff and mild- */
/*        ly stiff differential equations when derivative evaluations are */
/*        not expensive.  It should generally not be used to get high */
/*        accuracy results nor answers at a great many specific points. */
/*        Because DERKF has very low overhead costs, it will usually */
/*        result in the least expensive integration when solving */
/*        problems requiring a modest amount of accuracy and having */
/*        equations that are not costly to evaluate.  DERKF attempts to */
/*        discover when it is not suitable for the task posed. */

/*        DEABM is a variable order (one through twelve) Adams code. */
/*        Its complexity lies somewhere between that of DERKF and DEBDF. */
/*        DEABM is primarily designed to solve non-stiff and mildly */
/*        stiff differential equations when derivative evaluations are */
/*        expensive, high accuracy results are needed or answers at */
/*        many specific points are required.  DEABM attempts to discover */
/*        when it is not suitable for the task posed. */

/*        DEBDF is a variable order (one through five) backward */
/*        differentiation formula code.  It is the most complicated of */
/*        the three choices.  DEBDF is primarily designed to solve stiff */
/*        differential equations at crude to moderate tolerances. */
/*        If the problem is very stiff at all, DERKF and DEABM will be */
/*        quite inefficient compared to DEBDF.  However, DEBDF will be */
/*        inefficient compared to DERKF and DEABM on non-stiff problems */
/*        because it uses much more storage, has a much larger overhead, */
/*        and the low order formulas will not give high accuracies */
/*        efficiently. */

/*        The concept of stiffness cannot be described in a few words. */
/*        If you do not know the problem to be stiff, try either DERKF */
/*        or DEABM.  Both of these codes will inform you of stiffness */
/*        when the cost of solving such problems becomes important. */

/* ********************************************************************** */
/* ** ABSTRACT ** */
/* ********************************************************************** */

/*   Subroutine DEBDF uses the backward differentiation formulas of */
/*   orders one through five to integrate a system of NEQ first order */
/*   ordinary differential equations of the form */
/*                         DU/DX = F(X,U) */
/*   when the vector Y(*) of initial values for U(*) at X=T is given. The */
/*   subroutine integrates from T to TOUT.  It is easy to continue the */
/*   integration to get results at additional TOUT.  This is the interval */
/*   mode of operation.  It is also easy for the routine to return with */
/*   The solution at each intermediate step on the way to TOUT.  This is */
/*   the intermediate-output mode of operation. */

/* ********************************************************************** */
/* ** DESCRIPTION OF THE ARGUMENTS TO DEBDF (AN OVERVIEW) ** */
/* ********************************************************************** */

/*   The Parameters are: */

/*      F -- This is the name of a subroutine which you provide to */
/*             define the differential equations. */

/*      NEQ -- This is the number of (first order) differential */
/*             equations to be integrated. */

/*      T -- This is a value of the independent variable. */

/*      Y(*) -- This array contains the solution components at T. */

/*      TOUT -- This is a point at which a solution is desired. */

/*      INFO(*) -- The basic task of the code is to integrate the */
/*             differential equations from T to TOUT and return an */
/*             answer at TOUT.  INFO(*) is an INTEGER array which is used */
/*             to communicate exactly how you want this task to be */
/*             carried out. */

/*      RTOL, ATOL -- These quantities */
/*             represent relative and absolute error tolerances which you */
/*             provide to indicate how accurately you wish the solution */
/*             to be computed.  You may choose them to be both scalars */
/*             or else both vectors. */

/*      IDID -- This scalar quantity is an indicator reporting what */
/*             the code did.  You must monitor this INTEGER variable to */
/*             decide what action to take next. */

/*      RWORK(*), LRW -- RWORK(*) is a REAL work array of */
/*             length LRW which provides the code with needed storage */
/*             space. */

/*      IWORK(*), LIW -- IWORK(*) is an INTEGER work array of length LIW */
/*             which provides the code with needed storage space and an */
/*             across call flag. */

/*      RPAR, IPAR -- These are REAL and INTEGER parameter */
/*             arrays which you can use for communication between your */
/*             calling program and the F subroutine (and the JAC */
/*             subroutine). */

/*      JAC -- This is the name of a subroutine which you may choose to */
/*             provide for defining the Jacobian matrix of partial */
/*             derivatives DF/DU. */

/*  Quantities which are used as input items are */
/*             NEQ, T, Y(*), TOUT, INFO(*), */
/*             RTOL, ATOL, RWORK(1), LRW, */
/*             IWORK(1), IWORK(2), and LIW. */

/*  Quantities which may be altered by the code are */
/*             T, Y(*), INFO(1), RTOL, ATOL, */
/*             IDID, RWORK(*) and IWORK(*). */

/* ********************************************************************** */
/* * INPUT -- What To Do On The First Call To DEBDF * */
/* ********************************************************************** */

/*   The first call of the code is defined to be the start of each new */
/*   problem.  Read through the descriptions of all the following items, */
/*   provide sufficient storage space for designated arrays, set */
/*   appropriate variables for the initialization of the problem, and */
/*   give information about how you want the problem to be solved. */


/*      F -- provide a subroutine of the form */
/*                               F(X,U,UPRIME,RPAR,IPAR) */
/*             to define the system of first order differential equations */
/*             which is to be solved. For the given values of X and the */
/*             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must */
/*             evaluate the NEQ components of the system of differential */
/*             equations  DU/DX=F(X,U)  and store the derivatives in the */
/*             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for */
/*             equations I=1,...,NEQ. */

/*             Subroutine F must not alter X or U(*).  You must declare */
/*             the name F in an external statement in your program that */
/*             calls DEBDF.  You must dimension U and UPRIME in F. */

/*             RPAR and IPAR are REAL and INTEGER parameter arrays which */
/*             you can use for communication between your calling program */
/*             and subroutine F.  They are not used or altered by DEBDF. */
/*             If you do not need RPAR or IPAR, ignore these parameters */
/*             by treating them as dummy arguments.  If you do choose to */
/*             use them, dimension them in your calling program and in F */
/*             as arrays of appropriate length. */

/*      NEQ -- Set it to the number of differential equations. */
/*             (NEQ .GE. 1) */

/*      T -- Set it to the initial point of the integration. */
/*             You must use a program variable for T because the code */
/*             changes its value. */

/*      Y(*) -- Set this vector to the initial values of the NEQ solution */
/*             components at the initial point.  You must dimension Y at */
/*             least NEQ in your calling program. */

/*      TOUT -- Set it to the first point at which a solution is desired. */
/*             You can take TOUT = T, in which case the code */
/*             will evaluate the derivative of the solution at T and */
/*             return.  Integration either forward in T  (TOUT .GT. T) */
/*             or backward in T  (TOUT .LT. T)  is permitted. */

/*             The code advances the solution from T to TOUT using */
/*             step sizes which are automatically selected so as to */
/*             achieve the desired accuracy.  If you wish, the code will */
/*             return with the solution and its derivative following */
/*             each intermediate step (intermediate-output mode) so that */
/*             you can monitor them, but you still must provide TOUT in */
/*             accord with the basic aim of the code. */

/*             The first step taken by the code is a critical one */
/*             because it must reflect how fast the solution changes near */
/*             the initial point.  The code automatically selects an */
/*             initial step size which is practically always suitable for */
/*             the problem.  By using the fact that the code will not */
/*             step past TOUT in the first step, you could, if necessary, */
/*             restrict the length of the initial step size. */

/*             For some problems it may not be permissible to integrate */
/*             past a point TSTOP because a discontinuity occurs there */
/*             or the solution or its derivative is not defined beyond */
/*             TSTOP.  When you have declared a TSTOP point (see INFO(4) */
/*             and RWORK(1)), you have told the code not to integrate */
/*             past TSTOP.  In this case any TOUT beyond TSTOP is invalid */
/*             input. */

/*      INFO(*) -- Use the INFO array to give the code more details about */
/*             how you want your problem solved.  This array should be */
/*             dimensioned of length 15 to accommodate other members of */
/*             DEPAC or possible future extensions, though DEBDF uses */
/*             only the first six entries.  You must respond to all of */
/*             the following items which are arranged as questions.  The */
/*             simplest use of the code corresponds to answering all */
/*             questions as YES ,i.e. setting all entries of INFO to 0. */

/*        INFO(1) -- This parameter enables the code to initialize */
/*               itself.  You must set it to indicate the start of every */
/*               new problem. */

/*            **** Is this the first call for this problem ... */
/*                  YES -- Set INFO(1) = 0 */
/*                   NO -- Not applicable here. */
/*                         See below for continuation calls.  **** */

/*        INFO(2) -- How much accuracy you want of your solution */
/*               is specified by the error tolerances RTOL and ATOL. */
/*               The simplest use is to take them both to be scalars. */
/*               To obtain more flexibility, they can both be vectors. */
/*               The code must be told your choice. */

/*            **** Are both error tolerances RTOL, ATOL scalars ... */
/*                  YES -- Set INFO(2) = 0 */
/*                         and input scalars for both RTOL and ATOL */
/*                   NO -- Set INFO(2) = 1 */
/*                         and input arrays for both RTOL and ATOL **** */

/*        INFO(3) -- The code integrates from T in the direction */
/*               of TOUT by steps.  If you wish, it will return the */
/*               computed solution and derivative at the next */
/*               intermediate step (the intermediate-output mode) or */
/*               TOUT, whichever comes first.  This is a good way to */
/*               proceed if you want to see the behavior of the solution. */
/*               If you must have solutions at a great many specific */
/*               TOUT points, this code will compute them efficiently. */

/*            **** Do you want the solution only at */
/*                 TOUT (and NOT at the next intermediate step) ... */
/*                  YES -- Set INFO(3) = 0 */
/*                   NO -- Set INFO(3) = 1 **** */

/*        INFO(4) -- To handle solutions at a great many specific */
/*               values TOUT efficiently, this code may integrate past */
/*               TOUT and interpolate to obtain the result at TOUT. */
/*               Sometimes it is not possible to integrate beyond some */
/*               point TSTOP because the equation changes there or it is */
/*               not defined past TSTOP.  Then you must tell the code */
/*               not to go past. */

/*            **** Can the integration be carried out without any */
/*                 restrictions on the independent variable T ... */
/*                  YES -- Set INFO(4)=0 */
/*                   NO -- Set INFO(4)=1 */
/*                         and define the stopping point TSTOP by */
/*                         setting RWORK(1)=TSTOP **** */

/*        INFO(5) -- To solve stiff problems it is necessary to use the */
/*               Jacobian matrix of partial derivatives of the system */
/*               of differential equations.  If you do not provide a */
/*               subroutine to evaluate it analytically (see the */
/*               description of the item JAC in the call list), it will */
/*               be approximated by numerical differencing in this code. */
/*               Although it is less trouble for you to have the code */
/*               compute partial derivatives by numerical differencing, */
/*               the solution will be more reliable if you provide the */
/*               derivatives via JAC.  Sometimes numerical differencing */
/*               is cheaper than evaluating derivatives in JAC and */
/*               sometimes it is not - this depends on your problem. */

/*               If your problem is linear, i.e. has the form */
/*               DU/DX = F(X,U) = J(X)*U + G(X)   for some matrix J(X) */
/*               and vector G(X), the Jacobian matrix  DF/DU = J(X). */
/*               Since you must provide a subroutine to evaluate F(X,U) */
/*               analytically, it is little extra trouble to provide */
/*               subroutine JAC for evaluating J(X) analytically. */
/*               Furthermore, in such cases, numerical differencing is */
/*               much more expensive than analytic evaluation. */

/*            **** Do you want the code to evaluate the partial */
/*                 derivatives automatically by numerical differences ... */
/*                  YES -- Set INFO(5)=0 */
/*                   NO -- Set INFO(5)=1 */
/*                         and provide subroutine JAC for evaluating the */
/*                         Jacobian matrix **** */

/*        INFO(6) -- DEBDF will perform much better if the Jacobian */
/*               matrix is banded and the code is told this.  In this */
/*               case, the storage needed will be greatly reduced, */
/*               numerical differencing will be performed more cheaply, */
/*               and a number of important algorithms will execute much */
/*               faster.  The differential equation is said to have */
/*               half-bandwidths ML (lower) and MU (upper) if equation I */
/*               involves only unknowns Y(J) with */
/*                              I-ML .LE. J .LE. I+MU */
/*               for all I=1,2,...,NEQ.  Thus, ML and MU are the widths */
/*               of the lower and upper parts of the band, respectively, */
/*               with the main diagonal being excluded.  If you do not */
/*               indicate that the equation has a banded Jacobian, */
/*               the code works with a full matrix of NEQ**2 elements */
/*               (stored in the conventional way).  Computations with */
/*               banded matrices cost less time and storage than with */
/*               full matrices if  2*ML+MU .LT. NEQ.  If you tell the */
/*               code that the Jacobian matrix has a banded structure and */
/*               you want to provide subroutine JAC to compute the */
/*               partial derivatives, then you must be careful to store */
/*               the elements of the Jacobian matrix in the special form */
/*               indicated in the description of JAC. */

/*            **** Do you want to solve the problem using a full */
/*                 (dense) Jacobian matrix (and not a special banded */
/*                 structure) ... */
/*                  YES -- Set INFO(6)=0 */
/*                   NO -- Set INFO(6)=1 */
/*                         and provide the lower (ML) and upper (MU) */
/*                         bandwidths by setting */
/*                         IWORK(1)=ML */
/*                         IWORK(2)=MU **** */

/*      RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL) */
/*             error tolerances to tell the code how accurately you want */
/*             the solution to be computed.  They must be defined as */
/*             program variables because the code may change them.  You */
/*             have two choices -- */
/*                  Both RTOL and ATOL are scalars. (INFO(2)=0) */
/*                  Both RTOL and ATOL are vectors. (INFO(2)=1) */
/*             In either case all components must be non-negative. */

/*             The tolerances are used by the code in a local error test */
/*             at each step which requires roughly that */
/*                     ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL */
/*             for each vector component. */
/*             (More specifically, a root-mean-square norm is used to */
/*             measure the size of vectors, and the error test uses the */
/*             magnitude of the solution at the beginning of the step.) */

/*             The true (global) error is the difference between the true */
/*             solution of the initial value problem and the computed */
/*             approximation.  Practically all present day codes, */
/*             including this one, control the local error at each step */
/*             and do not even attempt to control the global error */
/*             directly.  Roughly speaking, they produce a solution Y(T) */
/*             which satisfies the differential equations with a */
/*             residual R(T),    DY(T)/DT = F(T,Y(T)) + R(T)   , */
/*             and, almost always, R(T) is bounded by the error */
/*             tolerances.  Usually, but not always, the true accuracy of */
/*             the computed Y is comparable to the error tolerances. This */
/*             code will usually, but not always, deliver a more accurate */
/*             solution if you reduce the tolerances and integrate again. */
/*             By comparing two such solutions you can get a fairly */
/*             reliable idea of the true error in the solution at the */
/*             bigger tolerances. */

/*             Setting ATOL=0. results in a pure relative error test on */
/*             that component.  Setting RTOL=0. results in a pure abso- */
/*             lute error test on that component.  A mixed test with non- */
/*             zero RTOL and ATOL corresponds roughly to a relative error */
/*             test when the solution component is much bigger than ATOL */
/*             and to an absolute error test when the solution component */
/*             is smaller than the threshold ATOL. */

/*             Proper selection of the absolute error control parameters */
/*             ATOL  requires you to have some idea of the scale of the */
/*             solution components.  To acquire this information may mean */
/*             that you will have to solve the problem more than once. In */
/*             the absence of scale information, you should ask for some */
/*             relative accuracy in all the components (by setting  RTOL */
/*             values non-zero) and perhaps impose extremely small */
/*             absolute error tolerances to protect against the danger of */
/*             a solution component becoming zero. */

/*             The code will not attempt to compute a solution at an */
/*             accuracy unreasonable for the machine being used.  It will */
/*             advise you if you ask for too much accuracy and inform */
/*             you as to the maximum accuracy it believes possible. */

/*      RWORK(*) -- Dimension this REAL work array of length LRW in your */
/*             calling program. */

/*      RWORK(1) -- If you have set INFO(4)=0, you can ignore this */
/*             optional input parameter.  Otherwise you must define a */
/*             stopping point TSTOP by setting   RWORK(1) = TSTOP. */
/*             (For some problems it may not be permissible to integrate */
/*             past a point TSTOP because a discontinuity occurs there */
/*             or the solution or its derivative is not defined beyond */
/*             TSTOP.) */

/*      LRW -- Set it to the declared length of the RWORK array. */
/*             You must have */
/*                  LRW .GE. 250+10*NEQ+NEQ**2 */
/*             for the full (dense) Jacobian case (when INFO(6)=0),  or */
/*                  LRW .GE. 250+10*NEQ+(2*ML+MU+1)*NEQ */
/*             for the banded Jacobian case (when INFO(6)=1). */

/*      IWORK(*) -- Dimension this INTEGER work array of length LIW in */
/*             your calling program. */

/*      IWORK(1), IWORK(2) -- If you have set INFO(6)=0, you can ignore */
/*             these optional input parameters. Otherwise you must define */
/*             the half-bandwidths ML (lower) and MU (upper) of the */
/*             Jacobian matrix by setting    IWORK(1) = ML   and */
/*             IWORK(2) = MU.  (The code will work with a full matrix */
/*             of NEQ**2 elements unless it is told that the problem has */
/*             a banded Jacobian, in which case the code will work with */
/*             a matrix containing at most  (2*ML+MU+1)*NEQ  elements.) */

/*      LIW -- Set it to the declared length of the IWORK array. */
/*             You must have LIW .GE. 56+NEQ. */

/*      RPAR, IPAR -- These are parameter arrays, of REAL and INTEGER */
/*             type, respectively.  You can use them for communication */
/*             between your program that calls DEBDF and the  F */
/*             subroutine (and the JAC subroutine).  They are not used or */
/*             altered by DEBDF.  If you do not need RPAR or IPAR, ignore */
/*             these parameters by treating them as dummy arguments.  If */
/*             you do choose to use them, dimension them in your calling */
/*             program and in F (and in JAC) as arrays of appropriate */
/*             length. */

/*      JAC -- If you have set INFO(5)=0, you can ignore this parameter */
/*             by treating it as a dummy argument. (For some compilers */
/*             you may have to write a dummy subroutine named  JAC  in */
/*             order to avoid problems associated with missing external */
/*             routine names.)  Otherwise, you must provide a subroutine */
/*             of the form */
/*                          JAC(X,U,PD,NROWPD,RPAR,IPAR) */
/*             to define the Jacobian matrix of partial derivatives DF/DU */
/*             of the system of differential equations   DU/DX = F(X,U). */
/*             For the given values of X and the vector */
/*             U(*)=(U(1),U(2),...,U(NEQ)), the subroutine must evaluate */
/*             the non-zero partial derivatives  DF(I)/DU(J)  for each */
/*             differential equation I=1,...,NEQ and each solution */
/*             component J=1,...,NEQ , and store these values in the */
/*             matrix PD.  The elements of PD are set to zero before each */
/*             call to JAC so only non-zero elements need to be defined. */

/*             Subroutine JAC must not alter X, U(*), or NROWPD.  You */
/*             must declare the name JAC in an EXTERNAL statement in your */
/*             program that calls DEBDF.  NROWPD is the row dimension of */
/*             the PD matrix and is assigned by the code.  Therefore you */
/*             must dimension PD in JAC according to */
/*                              DIMENSION PD(NROWPD,1) */
/*             You must also dimension U in JAC. */

/*             The way you must store the elements into the PD matrix */
/*             depends on the structure of the Jacobian which you */
/*             indicated by INFO(6). */
/*             *** INFO(6)=0 -- Full (Dense) Jacobian *** */
/*                 When you evaluate the (non-zero) partial derivative */
/*                 of equation I with respect to variable J, you must */
/*                 store it in PD according to */
/*                                PD(I,J) = * DF(I)/DU(J) * */
/*             *** INFO(6)=1 -- Banded Jacobian with ML Lower and MU */
/*                 Upper Diagonal Bands (refer to INFO(6) description of */
/*                 ML and MU) *** */
/*                 When you evaluate the (non-zero) partial derivative */
/*                 of equation I with respect to variable J, you must */
/*                 store it in PD according to */
/*                                IROW = I - J + ML + MU + 1 */
/*                                PD(IROW,J) = * DF(I)/DU(J) * */

/*             RPAR and IPAR are REAL and INTEGER parameter */
/*             arrays which you can use for communication between your */
/*             calling program and your Jacobian subroutine JAC.  They */
/*             are not altered by DEBDF.  If you do not need RPAR or */
/*             IPAR, ignore these parameters by treating them as dummy */
/*             arguments.  If you do choose to use them, dimension them */
/*             in your calling program and in JAC as arrays of */
/*             appropriate length. */

/* ********************************************************************** */
/* * OUTPUT -- After any return from DDEBDF * */
/* ********************************************************************** */

/*   The principal aim of the code is to return a computed solution at */
/*   TOUT, although it is also possible to obtain intermediate results */
/*   along the way.  To find out whether the code achieved its goal */
/*   or if the integration process was interrupted before the task was */
/*   completed, you must check the IDID parameter. */


/*      T -- The solution was successfully advanced to the */
/*             output value of T. */

/*      Y(*) -- Contains the computed solution approximation at T. */
/*             You may also be interested in the approximate derivative */
/*             of the solution at T.  It is contained in */
/*             RWORK(21),...,RWORK(20+NEQ). */

/*      IDID -- Reports what the code did */

/*                         *** Task Completed *** */
/*                   Reported by positive values of IDID */

/*             IDID = 1 -- A step was successfully taken in the */
/*                       intermediate-output mode.  The code has not */
/*                       yet reached TOUT. */

/*             IDID = 2 -- The integration to TOUT was successfully */
/*                       completed (T=TOUT) by stepping exactly to TOUT. */

/*             IDID = 3 -- The integration to TOUT was successfully */
/*                       completed (T=TOUT) by stepping past TOUT. */
/*                       Y(*) is obtained by interpolation. */

/*                         *** Task Interrupted *** */
/*                   Reported by negative values of IDID */

/*             IDID = -1 -- A large amount of work has been expended. */
/*                       (500 steps attempted) */

/*             IDID = -2 -- The error tolerances are too stringent. */

/*             IDID = -3 -- The local error test cannot be satisfied */
/*                       because you specified a zero component in ATOL */
/*                       and the corresponding computed solution */
/*                       component is zero.  Thus, a pure relative error */
/*                       test is impossible for this component. */

/*             IDID = -4,-5  -- Not applicable for this code but used */
/*                       by other members of DEPAC. */

/*             IDID = -6 -- DEBDF had repeated convergence test failures */
/*                       on the last attempted step. */

/*             IDID = -7 -- DEBDF had repeated error test failures on */
/*                       the last attempted step. */

/*             IDID = -8,..,-32  -- Not applicable for this code but */
/*                       used by other members of DEPAC or possible */
/*                       future extensions. */

/*                         *** Task Terminated *** */
/*                   Reported by the value of IDID=-33 */

/*             IDID = -33 -- The code has encountered trouble from which */
/*                       it cannot recover.  A message is printed */
/*                       explaining the trouble and control is returned */
/*                       to the calling program.  For example, this */
/*                       occurs when invalid input is detected. */

/*      RTOL, ATOL -- These quantities remain unchanged except when */
/*             IDID = -2.  In this case, the error tolerances have been */
/*             increased by the code to values which are estimated to be */
/*             appropriate for continuing the integration.  However, the */
/*             reported solution at T was obtained using the input values */
/*             of RTOL and ATOL. */

/*      RWORK, IWORK -- Contain information which is usually of no */
/*             interest to the user but necessary for subsequent calls. */
/*             However, you may find use for */

/*             RWORK(11)--which contains the step size H to be */
/*                        attempted on the next step. */

/*             RWORK(12)--If the tolerances have been increased by the */
/*                        code (IDID = -2) , they were multiplied by the */
/*                        value in RWORK(12). */

/*             RWORK(13)--which contains the current value of the */
/*                        independent variable, i.e. the farthest point */
/*                        integration has reached.  This will be */
/*                        different from T only when interpolation has */
/*                        been performed (IDID=3). */

/*             RWORK(20+I)--which contains the approximate derivative */
/*                        of the solution component Y(I).  In DEBDF, it */
/*                        is never obtained by calling subroutine F to */
/*                        evaluate the differential equation using T and */
/*                        Y(*), except at the initial point of */
/*                        integration. */

/* ********************************************************************** */
/* ** INPUT -- What To Do To Continue The Integration ** */
/* **             (calls after the first)             ** */
/* ********************************************************************** */

/*        This code is organized so that subsequent calls to continue the */
/*        integration involve little (if any) additional effort on your */
/*        part. You must monitor the IDID parameter in order to determine */
/*        what to do next. */

/*        Recalling that the principal task of the code is to integrate */
/*        from T to TOUT (the interval mode), usually all you will need */
/*        to do is specify a new TOUT upon reaching the current TOUT. */

/*        Do not alter any quantity not specifically permitted below, */
/*        in particular do not alter NEQ, T, Y(*), RWORK(*), IWORK(*) or */
/*        the differential equation in subroutine F.  Any such alteration */
/*        constitutes a new problem and must be treated as such, i.e. */
/*        you must start afresh. */

/*        You cannot change from vector to scalar error control or vice */
/*        versa (INFO(2)) but you can change the size of the entries of */
/*        RTOL, ATOL.  Increasing a tolerance makes the equation easier */
/*        to integrate.  Decreasing a tolerance will make the equation */
/*        harder to integrate and should generally be avoided. */

/*        You can switch from the intermediate-output mode to the */
/*        interval mode (INFO(3)) or vice versa at any time. */

/*        If it has been necessary to prevent the integration from going */
/*        past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the */
/*        code will not integrate to any TOUT beyond the currently */
/*        specified TSTOP.  Once TSTOP has been reached you must change */
/*        the value of TSTOP or set INFO(4)=0.  You may change INFO(4) */
/*        or TSTOP at any time but you must supply the value of TSTOP in */
/*        RWORK(1) whenever you set INFO(4)=1. */

/*        Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2) */
/*        unless you are going to restart the code. */

/*        The parameter INFO(1) is used by the code to indicate the */
/*        beginning of a new problem and to indicate whether integration */
/*        is to be continued.  You must input the value  INFO(1) = 0 */
/*        when starting a new problem.  You must input the value */
/*        INFO(1) = 1  if you wish to continue after an interrupted task. */
/*        Do not set  INFO(1) = 0  on a continuation call unless you */
/*        want the code to restart at the current T. */

/*                         *** Following a Completed Task *** */
/*         If */
/*             IDID = 1, call the code again to continue the integration */
/*                     another step in the direction of TOUT. */

/*             IDID = 2 or 3, define a new TOUT and call the code again. */
/*                     TOUT must be different from T.  You cannot change */
/*                     the direction of integration without restarting. */

/*                         *** Following an Interrupted Task *** */
/*                     To show the code that you realize the task was */
/*                     interrupted and that you want to continue, you */
/*                     must take appropriate action and reset INFO(1) = 1 */
/*         If */
/*             IDID = -1, the code has attempted 500 steps. */
/*                     If you want to continue, set INFO(1) = 1 and */
/*                     call the code again.  An additional 500 steps */
/*                     will be allowed. */

/*             IDID = -2, the error tolerances RTOL, ATOL have been */
/*                     increased to values the code estimates appropriate */
/*                     for continuing.  You may want to change them */
/*                     yourself.  If you are sure you want to continue */
/*                     with relaxed error tolerances, set INFO(1)=1 and */
/*                     call the code again. */

/*             IDID = -3, a solution component is zero and you set the */
/*                     corresponding component of ATOL to zero.  If you */
/*                     are sure you want to continue, you must first */
/*                     alter the error criterion to use positive values */
/*                     for those components of ATOL corresponding to zero */
/*                     solution components, then set INFO(1)=1 and call */
/*                     the code again. */

/*             IDID = -4,-5  --- cannot occur with this code but used */
/*                     by other members of DEPAC. */

/*             IDID = -6, repeated convergence test failures occurred */
/*                     on the last attempted step in DEBDF.  An inaccu- */
/*                     rate Jacobian may be the problem.  If you are */
/*                     absolutely certain you want to continue, restart */
/*                     the integration at the current T by setting */
/*                     INFO(1)=0 and call the code again. */

/*             IDID = -7, repeated error test failures occurred on the */
/*                     last attempted step in DEBDF.  A singularity in */
/*                     the solution may be present.  You should re- */
/*                     examine the problem being solved.  If you are */
/*                     absolutely certain you want to continue, restart */
/*                     the integration at the current T by setting */
/*                     INFO(1)=0 and call the code again. */

/*             IDID = -8,..,-32  --- cannot occur with this code but */
/*                     used by other members of DEPAC or possible future */
/*                     extensions. */

/*                         *** Following a Terminated Task *** */
/*         If */
/*             IDID = -33, you cannot continue the solution of this */
/*                     problem.  An attempt to do so will result in your */
/*                     run being terminated. */

/* ********************************************************************** */

/*         ***** Warning ***** */

/*     If DEBDF is to be used in an overlay situation, you must save and */
/*     restore certain items used internally by DEBDF  (values in the */
/*     common block DEBDF1).  This can be accomplished as follows. */

/*     To save the necessary values upon return from DEBDF, simply call */
/*        SVCO(RWORK(22+NEQ),IWORK(21+NEQ)). */

/*     To restore the necessary values before the next call to DEBDF, */
/*     simply call    RSCO(RWORK(22+NEQ),IWORK(21+NEQ)). */

/* ***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user */
/*                 oriented package of ODE solvers, Report SAND79-2374, */
/*                 Sandia Laboratories, 1979. */
/* ***ROUTINES CALLED  LSOD, XERMSG */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891024  Changed references from VNORM to HVNRM.  (WRB) */
/*   891024  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, change Prologue */
/*           comments to agree with DDEBDF.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DEBDF */






/*        CHECK FOR AN APPARENT INFINITE LOOP */

/* ***FIRST EXECUTABLE STATEMENT  DEBDF */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwork;
    --rwork;
    --atol;
    --rtol;
    --info;
    --y;

    /* Function Body */
    if (info[1] == 0) {
	iwork[*liw] = 0;
    }

    if (iwork[*liw] >= 5) {
	if (*t == rwork[*neq + 21]) {
	    s_wsfi(&io___2);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 82, a__1[0] = "AN APPARENT INFINITE LOOP HAS BEEN DETE"
		    "CTED.$$YOU HAVE MADE REPEATED CALLS AT T = ";
	    i__1[1] = 16, a__1[1] = xern3;
	    i__1[2] = 50, a__1[2] = " AND THE INTEGRATION HAS NOT ADVANCED. "
		    " CHECK THE ";
	    i__1[3] = 48, a__1[3] = "WAY YOU HAVE SET PARAMETERS FOR THE CAL"
		    "L TO THE ";
	    i__1[4] = 26, a__1[4] = "CODE PARTICULARLY INFO(1).";
	    s_cat(ch__1, a__1, i__1, &c__5, (ftnlen)222);
	    xermsg_("SLATEC", "DEBDF", ch__1, &c__13, &c__2, (ftnlen)6, (
		    ftnlen)5, (ftnlen)222);
	    return 0;
	}
    }

    *idid = 0;

/*        CHECK VALIDITY OF INFO PARAMETERS */

    if (info[1] != 0 && info[1] != 1) {
	s_wsfi(&io___4);
	do_fio(&c__1, (char *)&info[1], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 204, a__2[0] = "INFO(1) MUST BE SET TO 0 FOR THE  START OF"
		" A NEW PROBLEM, AND MUST BE SET TO 1 FOLLOWING AN INTERRUPTE"
		"D TASK.  YOU ARE ATTEMPTING TO CONTINUE THE INTEGRATION ILLE"
		"GALLY BY CALLING THE CODE WITH  INFO(1) = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__2, a__2, i__2, &c__2, (ftnlen)212);
	xermsg_("SLATEC", "DEBDF", ch__2, &c__3, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)212);
	*idid = -33;
    }

    if (info[2] != 0 && info[2] != 1) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&info[2], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 125, a__2[0] = "INFO(2) MUST BE 0 OR 1 INDICATING SCALAR A"
		"ND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HAVE CALLED T"
		"HE CODE WITH INFO(2) = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__3, a__2, i__2, &c__2, (ftnlen)133);
	xermsg_("SLATEC", "DEBDF", ch__3, &c__4, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)133);
	*idid = -33;
    }

    if (info[3] != 0 && info[3] != 1) {
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&info[3], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 147, a__2[0] = "INFO(3) MUST BE 0 OR 1 INDICATING THE INTE"
		"RVAL OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, RESPECTIVEL"
		"Y.  YOU HAVE CALLED THE CODE WITH  INFO(3) = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__4, a__2, i__2, &c__2, (ftnlen)155);
	xermsg_("SLATEC", "DEBDF", ch__4, &c__5, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)155);
	*idid = -33;
    }

    if (info[4] != 0 && info[4] != 1) {
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&info[4], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 154, a__2[0] = "INFO(4) MUST BE 0 OR 1 INDICATING WHETHER "
		"OR NOT THE INTEGRATION INTERVAL IS TO BE RESTRICTED BY A POI"
		"NT TSTOP.  YOU HAVE CALLED THE CODE  WITH INFO(4) = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__5, a__2, i__2, &c__2, (ftnlen)162);
	xermsg_("SLATEC", "DEBDF", ch__5, &c__14, &c__1, (ftnlen)6, (ftnlen)5,
		 (ftnlen)162);
	*idid = -33;
    }

    if (info[5] != 0 && info[5] != 1) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&info[5], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 211, a__2[0] = "INFO(5) MUST BE 0 OR 1 INDICATING WHETHER "
		"THE CODE IS TOLD TO FORM THE JACOBIAN MATRIX BY NUMERICAL DI"
		"FFERENCING OR YOU PROVIDE A SUBROUTINE TO EVALUATE IT ANALYT"
		"ICALLY.  YOU HAVE CALLED THE CODE WITH INFO(5) = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__6, a__2, i__2, &c__2, (ftnlen)219);
	xermsg_("SLATEC", "DEBDF", ch__6, &c__15, &c__1, (ftnlen)6, (ftnlen)5,
		 (ftnlen)219);
	*idid = -33;
    }

    if (info[6] != 0 && info[6] != 1) {
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&info[6], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 188, a__2[0] = "INFO(6) MUST BE 0 OR 1 INDICATING WHETHER "
		"THE CODE IS TOLD TO TREAT THE JACOBIAN AS A FULL (DENSE) MAT"
		"RIX OR AS HAVING A SPECIAL BANDED STRUCTURE.  YOU HAVE CALLE"
		"D THE CODE WITH INFO(6) = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__7, a__2, i__2, &c__2, (ftnlen)196);
	xermsg_("SLATEC", "DEBDF", ch__7, &c__16, &c__1, (ftnlen)6, (ftnlen)5,
		 (ftnlen)196);
	*idid = -33;
    }

    ilrw = *neq;
    if (info[6] != 0) {

/*        CHECK BANDWIDTH PARAMETERS */

	ml = iwork[1];
	mu = iwork[2];
	ilrw = (ml << 1) + mu + 1;

	if (ml < 0 || ml >= *neq || mu < 0 || mu >= *neq) {
	    s_wsfi(&io___13);
	    do_fio(&c__1, (char *)&ml, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&mu, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 242, a__3[0] = "YOU HAVE SET INFO(6) = 1, TELLING THE "
		    "CODE THAT THE JACOBIAN MATRIX HAS A SPECIAL BANDED STRUC"
		    "TURE.  HOWEVER, THE LOWER (UPPER) BANDWIDTHS  ML (MU) VI"
		    "OLATE THE CONSTRAINTS ML,MU .GE. 0 AND  ML,MU .LT. NEQ. "
		    " YOU HAVE CALLED THE CODE WITH ML = ";
	    i__3[1] = 8, a__3[1] = xern1;
	    i__3[2] = 10, a__3[2] = " AND MU = ";
	    i__3[3] = 8, a__3[3] = xern2;
	    s_cat(ch__8, a__3, i__3, &c__4, (ftnlen)268);
	    xermsg_("SLATEC", "DEBDF", ch__8, &c__17, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)268);
	    *idid = -33;
	}
    }

/*        CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION */

    if (*lrw < (ilrw + 10) * *neq + 250) {
	s_wsfi(&io___16);
	do_fio(&c__1, (char *)&(*lrw), (ftnlen)sizeof(integer));
	e_wsfi();
	if (info[6] == 0) {
/* Writing concatenation */
	    i__2[0] = 101, a__2[0] = "LENGTH OF ARRAY RWORK MUST BE AT LEAST"
		    " 250 + 10*NEQ + NEQ*NEQ.$$YOU HAVE CALLED THE CODE WITH "
		    " LRW = ";
	    i__2[1] = 8, a__2[1] = xern1;
	    s_cat(ch__9, a__2, i__2, &c__2, (ftnlen)109);
	    xermsg_("SLATEC", "DEBDF", ch__9, &c__1, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)109);
	} else {
/* Writing concatenation */
	    i__2[0] = 109, a__2[0] = "LENGTH OF ARRAY RWORK MUST BE AT LEAST"
		    " 250 + 10*NEQ + (2*ML+MU+1)*NEQ.$$YOU HAVE CALLED THE CO"
		    "DE WITH  LRW = ";
	    i__2[1] = 8, a__2[1] = xern1;
	    s_cat(ch__10, a__2, i__2, &c__2, (ftnlen)117);
	    xermsg_("SLATEC", "DEBDF", ch__10, &c__18, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)117);
	}
	*idid = -33;
    }

    if (*liw < *neq + 56) {
	s_wsfi(&io___17);
	do_fio(&c__1, (char *)&(*liw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 82, a__2[0] = "LENGTH OF ARRAY IWORK BE AT LEAST  56 + NEQ"
		".  YOU HAVE CALLED THE CODE WITH LIW = ";
	i__2[1] = 8, a__2[1] = xern1;
	s_cat(ch__11, a__2, i__2, &c__2, (ftnlen)90);
	xermsg_("SLATEC", "DEBDF", ch__11, &c__2, &c__1, (ftnlen)6, (ftnlen)5,
		 (ftnlen)90);
	*idid = -33;
    }

/*        COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK */
/*        ARRAY AND RESTORE COMMON BLOCK DATA */

    icomi = *neq + 21;
    iinout = icomi + 33;

    iypout = 21;
    itstar = *neq + 21;
    icomr = *neq + 22;

    if (info[1] != 0) {
	intout = iwork[iinout] != -1;
    }
/*     CALL RSCO(RWORK(ICOMR),IWORK(ICOMI)) */

    debdf1_1.iyh = icomr + 218;
    debdf1_1.iewt = debdf1_1.iyh + *neq * 6;
    debdf1_1.isavf = debdf1_1.iewt + *neq;
    debdf1_1.iacor = debdf1_1.isavf + *neq;
    debdf1_1.iwm = debdf1_1.iacor + *neq;
    idelsn = debdf1_1.iwm + 2 + ilrw * *neq;

    debdf1_1.ibegin = info[1];
    debdf1_1.itol = info[2];
    debdf1_1.iinteg = info[3];
    debdf1_1.itstop = info[4];
    debdf1_1.ijac = info[5];
    debdf1_1.iband = info[6];
    rwork[itstar] = *t;

    lsod_((U_fp)f, neq, t, &y[1], tout, &rtol[1], &atol[1], idid, &rwork[
	    iypout], &rwork[debdf1_1.iyh], &rwork[debdf1_1.iyh], &rwork[
	    debdf1_1.iewt], &rwork[debdf1_1.isavf], &rwork[debdf1_1.iacor], &
	    rwork[debdf1_1.iwm], &iwork[1], (U_fp)jac, &intout, &rwork[1], &
	    rwork[12], &rwork[idelsn], &rpar[1], &ipar[1]);

    iwork[iinout] = -1;
    if (intout) {
	iwork[iinout] = 1;
    }

    if (*idid != -2) {
	++iwork[*liw];
    }
    if (*t != rwork[itstar]) {
	iwork[*liw] = 0;
    }
/*     CALL SVCO(RWORK(ICOMR),IWORK(ICOMI)) */
    rwork[11] = debdf1_1.h__;
    rwork[13] = debdf1_1.tn;
    info[1] = debdf1_1.ibegin;

    return 0;
} /* debdf_ */


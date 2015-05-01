/* sbolsm.f -- translated by f2c (version 12.02.01).
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
static integer c__31 = 31;
static integer c__3 = 3;
static integer c__32 = 32;
static integer c__33 = 33;
static integer c__4 = 4;
static integer c__34 = 34;
static integer c__35 = 35;
static integer c__6 = 6;
static integer c__36 = 36;
static integer c__37 = 37;
static integer c__24 = 24;
static integer c__5 = 5;
static integer c__25 = 25;
static integer c__0 = 0;
static integer c__26 = 26;
static integer c__27 = 27;
static integer c__2 = 2;
static integer c__28 = 28;
static integer c__29 = 29;
static integer c__30 = 30;
static integer c__38 = 38;
static integer c__7 = 7;
static integer c__23 = 23;
static real c_b185 = 0.f;
static integer c_n4 = -4;
static integer c__22 = 22;

/* DECK SBOLSM */
/* Subroutine */ int sbolsm_(real *w, integer *mdw, integer *minput, integer *
	ncols, real *bl, real *bu, integer *ind, integer *iopt, real *x, real 
	*rnorm, integer *mode, real *rw, real *ww, real *scl, integer *ibasis,
	 integer *ibb)
{
    /* System generated locals */
    address a__1[3], a__2[4], a__3[6], a__4[5], a__5[2], a__6[7];
    integer w_dim1, w_offset, i__1[3], i__2[4], i__3, i__4[6], i__5[5], i__6[
	    2], i__7[7], i__8, i__9, i__10;
    real r__1, r__2;
    char ch__1[47], ch__2[50], ch__3[79], ch__4[53], ch__5[94], ch__6[75], 
	    ch__7[83], ch__8[92], ch__9[105], ch__10[102], ch__11[61], ch__12[
	    110], ch__13[134], ch__14[44], ch__15[76];

    /* Local variables */
    static integer i__, j;
    static real t, t1, t2, sc;
    static integer ip, jp, lp;
    static real ss, wt, cl1, cl2, cl3, fac, big;
    static integer lds;
    static real bou, beta;
    static integer jbig, jmag, ioff, jcol;
    static real wbig, wmag;
    static integer mval, iter;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real xnew;
    extern /* Subroutine */ int srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *);
    static char xern1[8], xern2[8], xern3[16], xern4[16];
    extern doublereal snrm2_(integer *, real *, integer *);
    static real alpha;
    static logical found;
    static integer nsetb, itemp, igopr, jdrop, itmax, lgopr;
    extern /* Subroutine */ int srotg_(real *, real *, real *, real *), 
	    scopy_(integer *, real *, integer *, real *, integer *), sswap_(
	    integer *, real *, integer *, real *, integer *), saxpy_(integer *
	    , real *, real *, integer *, real *, integer *), ivout_(integer *,
	     integer *, char *, integer *, ftnlen);
    static integer mrows;
    extern /* Subroutine */ int smout_(integer *, integer *, integer *, real *
	    , char *, integer *, ftnlen);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int svout_(integer *, real *, char *, integer *, 
	    ftnlen);
    static integer jdrop1, jdrop2, jlarge;
    static real colabv, colblo, wlarge, tolind;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer iprint;
    static logical constr;
    static real tolsze;

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___3 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___4 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___8 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___10 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___12 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___16 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___17 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___18 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___31 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___32 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___33 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___34 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___35 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___36 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___37 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___38 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___39 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___40 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___41 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___42 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___43 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___44 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___45 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___54 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SBOLSM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SBOCLS and SBOLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SBOLSM-S, DBOLSM-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*          Solve E*X = F (least squares sense) with bounds on */
/*            selected X values. */
/*     The user must have DIMENSION statements of the form: */

/*       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS), */
/*      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS) */
/*       INTEGER IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS) */

/*     (Here NX=number of extra locations required for options 1,...,7; */
/*     NX=0 for no options; here NI=number of extra locations possibly */
/*     required for options 1-7; NI=0 for no options; NI=14 if all the */
/*     options are simultaneously in use.) */

/*    INPUT */
/*    ----- */

/*    -------------------- */
/*    W(MDW,*),MINPUT,NCOLS */
/*    -------------------- */
/*     The array W(*,*) contains the matrix [E:F] on entry. The matrix */
/*     [E:F] has MINPUT rows and NCOLS+1 columns. This data is placed in */
/*     the array W(*,*) with E occupying the first NCOLS columns and the */
/*     right side vector F in column NCOLS+1. The row dimension, MDW, of */
/*     the array W(*,*) must satisfy the inequality MDW .ge. MINPUT. */
/*     Other values of MDW are errors. The values of MINPUT and NCOLS */
/*     must be positive. Other values are errors. */

/*    ------------------ */
/*    BL(*),BU(*),IND(*) */
/*    ------------------ */
/*     These arrays contain the information about the bounds that the */
/*     solution values are to satisfy. The value of IND(J) tells the */
/*     type of bound and BL(J) and BU(J) give the explicit values for */
/*     the respective upper and lower bounds. */

/*    1.    For IND(J)=1, require X(J) .ge. BL(J). */
/*    2.    For IND(J)=2, require X(J) .le. BU(J). */
/*    3.    For IND(J)=3, require X(J) .ge. BL(J) and */
/*                                X(J) .le. BU(J). */
/*    4.    For IND(J)=4, no bounds on X(J) are required. */
/*     The values of BL(*),BL(*) are modified by the subprogram. Values */
/*     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3 */
/*     (upper and lower bounds) the condition BL(J) .gt. BU(J) is an */
/*     error. */

/*    ------- */
/*    IOPT(*) */
/*    ------- */
/*     This is the array where the user can specify nonstandard options */
/*     for SBOLSM. Most of the time this feature can be ignored by */
/*     setting the input value IOPT(1)=99. Occasionally users may have */
/*     needs that require use of the following subprogram options. For */
/*     details about how to use the options see below: IOPT(*) CONTENTS. */

/*     Option Number   Brief Statement of Purpose */
/*     ----- ------   ----- --------- -- ------- */
/*           1         Move the IOPT(*) processing pointer. */
/*           2         Change rank determination tolerance. */
/*           3         Change blow-up factor that determines the */
/*                     size of variables being dropped from active */
/*                     status. */
/*           4         Reset the maximum number of iterations to use */
/*                     in solving the problem. */
/*           5         The data matrix is triangularized before the */
/*                     problem is solved whenever (NCOLS/MINPUT) .lt. */
/*                     FAC. Change the value of FAC. */
/*           6         Redefine the weighting matrix used for */
/*                     linear independence checking. */
/*           7         Debug output is desired. */
/*          99         No more options to change. */

/*    ---- */
/*    X(*) */
/*    ---- */
/*     This array is used to pass data associated with options 1,2,3 and */
/*     5. Ignore this input parameter if none of these options are used. */
/*     Otherwise see below: IOPT(*) CONTENTS. */

/*    ---------------- */
/*    IBASIS(*),IBB(*) */
/*    ---------------- */
/*     These arrays must be initialized by the user. The values */
/*         IBASIS(J)=J, J=1,...,NCOLS */
/*         IBB(J)   =1, J=1,...,NCOLS */
/*     are appropriate except when using nonstandard features. */

/*    ------ */
/*    SCL(*) */
/*    ------ */
/*     This is the array of scaling factors to use on the columns of the */
/*     matrix E. These values must be defined by the user. To suppress */
/*     any column scaling set SCL(J)=1.0, J=1,...,NCOLS. */

/*    OUTPUT */
/*    ------ */

/*    ---------- */
/*    X(*),RNORM */
/*    ---------- */
/*     The array X(*) contains a solution (if MODE .ge. 0 or .eq. -22) */
/*     for the constrained least squares problem. The value RNORM is the */
/*     minimum residual vector length. */

/*    ---- */
/*    MODE */
/*    ---- */
/*     The sign of mode determines whether the subprogram has completed */
/*     normally, or encountered an error condition or abnormal status. */
/*     A value of MODE .ge. 0 signifies that the subprogram has completed */
/*     normally. The value of MODE (.ge. 0) is the number of variables */
/*     in an active status: not at a bound nor at the value ZERO, for */
/*     the case of free variables. A negative value of MODE will be one */
/*     of the 18 cases -38,-37,...,-22, or -1. Values .lt. -1 correspond */
/*     to an abnormal completion of the subprogram. To understand the */
/*     abnormal completion codes see below: ERROR MESSAGES for SBOLSM */
/*     An approximate solution will be returned to the user only when */
/*     maximum iterations is reached, MODE=-22. */

/*    ----------- */
/*    RW(*),WW(*) */
/*    ----------- */
/*     These are working arrays each with NCOLS entries. The array RW(*) */
/*     contains the working (scaled, nonactive) solution values. The */
/*     array WW(*) contains the working (scaled, active) gradient vector */
/*     values. */

/*    ---------------- */
/*    IBASIS(*),IBB(*) */
/*    ---------------- */
/*     These arrays contain information about the status of the solution */
/*     when MODE .ge. 0. The indices IBASIS(K), K=1,...,MODE, show the */
/*     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are */
/*     the active variables. The value (IBB(J)-1) is the number of times */
/*     variable J was reflected from its upper bound. (Normally the user */
/*     can ignore these parameters.) */

/*    IOPT(*) CONTENTS */
/*    ------- -------- */
/*     The option array allows a user to modify internal variables in */
/*     the subprogram without recompiling the source code. A central */
/*     goal of the initial software design was to do a good job for most */
/*     people. Thus the use of options will be restricted to a select */
/*     group of users. The processing of the option array proceeds as */
/*     follows: a pointer, here called LP, is initially set to the value */
/*     1. The value is updated as the options are processed.  At the */
/*     pointer position the option number is extracted and used for */
/*     locating other information that allows for options to be changed. */
/*     The portion of the array IOPT(*) that is used for each option is */
/*     fixed; the user and the subprogram both know how many locations */
/*     are needed for each option. A great deal of error checking is */
/*     done by the subprogram on the contents of the option array. */
/*     Nevertheless it is still possible to give the subprogram optional */
/*     input that is meaningless. For example, some of the options use */
/*     the location X(NCOLS+IOFF) for passing data. The user must manage */
/*     the allocation of these locations when more than one piece of */
/*     option data is being passed to the subprogram. */

/*   1 */
/*   - */
/*     Move the processing pointer (either forward or backward) to the */
/*     location IOPT(LP+1). The processing pointer is moved to location */
/*     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over */
/*     locations 3,...,NCOLS+2 of IOPT(*), */

/*       IOPT(1)=1 */
/*       IOPT(2)=NCOLS+3 */
/*       (IOPT(I), I=3,...,NCOLS+2 are not defined here.) */
/*       IOPT(NCOLS+3)=99 */
/*       CALL SBOLSM */

/*     CAUTION: Misuse of this option can yield some very hard-to-find */
/*     bugs.  Use it with care. */

/*   2 */
/*   - */
/*     The algorithm that solves the bounded least squares problem */
/*     iteratively drops columns from the active set. This has the */
/*     effect of joining a new column vector to the QR factorization of */
/*     the rectangular matrix consisting of the partially triangularized */
/*     nonactive columns. After triangularizing this matrix a test is */
/*     made on the size of the pivot element. The column vector is */
/*     rejected as dependent if the magnitude of the pivot element is */
/*     .le. TOL* magnitude of the column in components strictly above */
/*     the pivot element. Nominally the value of this (rank) tolerance */
/*     is TOL = SQRT(R1MACH(4)). To change only the value of TOL, for */
/*     example, */

/*       X(NCOLS+1)=TOL */
/*       IOPT(1)=2 */
/*       IOPT(2)=1 */
/*       IOPT(3)=99 */
/*       CALL SBOLSM */

/*     Generally, if LP is the processing pointer for IOPT(*), */

/*       X(NCOLS+IOFF)=TOL */
/*       IOPT(LP)=2 */
/*       IOPT(LP+1)=IOFF */
/*        . */
/*       CALL SBOLSM */

/*     The required length of IOPT(*) is increased by 2 if option 2 is */
/*     used; The required length of X(*) is increased by 1. A value of */
/*     IOFF .le. 0 is an error. A value of TOL .le. R1MACH(4) gives a */
/*     warning message; it is not considered an error. */

/*   3 */
/*   - */
/*     A solution component is left active (not used) if, roughly */
/*     speaking, it seems too large. Mathematically the new component is */
/*     left active if the magnitude is .ge.((vector norm of F)/(matrix */
/*     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(R1MACH(4)). */
/*     To change only the value of BLOWUP, for example, */

/*       X(NCOLS+2)=BLOWUP */
/*       IOPT(1)=3 */
/*       IOPT(2)=2 */
/*       IOPT(3)=99 */
/*       CALL SBOLSM */

/*     Generally, if LP is the processing pointer for IOPT(*), */

/*       X(NCOLS+IOFF)=BLOWUP */
/*       IOPT(LP)=3 */
/*       IOPT(LP+1)=IOFF */
/*        . */
/*       CALL SBOLSM */

/*     The required length of IOPT(*) is increased by 2 if option 3 is */
/*     used; the required length of X(*) is increased by 1. A value of */
/*     IOFF .le. 0 is an error. A value of BLOWUP .le. 0.0 is an error. */

/*   4 */
/*   - */
/*     Normally the algorithm for solving the bounded least squares */
/*     problem requires between NCOLS/3 and NCOLS drop-add steps to */
/*     converge. (this remark is based on examining a small number of */
/*     test cases.) The amount of arithmetic for such problems is */
/*     typically about twice that required for linear least squares if */
/*     there are no bounds and if plane rotations are used in the */
/*     solution method. Convergence of the algorithm, while */
/*     mathematically certain, can be much slower than indicated. To */
/*     avoid this potential but unlikely event ITMAX drop-add steps are */
/*     permitted. Nominally ITMAX=5*(MAX(MINPUT,NCOLS)). To change the */
/*     value of ITMAX, for example, */

/*       IOPT(1)=4 */
/*       IOPT(2)=ITMAX */
/*       IOPT(3)=99 */
/*       CALL SBOLSM */

/*     Generally, if LP is the processing pointer for IOPT(*), */

/*       IOPT(LP)=4 */
/*       IOPT(LP+1)=ITMAX */
/*        . */
/*       CALL SBOLSM */

/*     The value of ITMAX must be .gt. 0. Other values are errors. Use */
/*     of this option increases the required length of IOPT(*) by 2. */

/*   5 */
/*   - */
/*     For purposes of increased efficiency the MINPUT by NCOLS+1 data */
/*     matrix [E:F] is triangularized as a first step whenever MINPUT */
/*     satisfies FAC*MINPUT .gt. NCOLS. Nominally FAC=0.75. To change the */
/*     value of FAC, */

/*       X(NCOLS+3)=FAC */
/*       IOPT(1)=5 */
/*       IOPT(2)=3 */
/*       IOPT(3)=99 */
/*       CALL SBOLSM */

/*     Generally, if LP is the processing pointer for IOPT(*), */

/*       X(NCOLS+IOFF)=FAC */
/*       IOPT(LP)=5 */
/*       IOPT(LP+1)=IOFF */
/*        . */
/*       CALL SBOLSM */

/*     The value of FAC must be nonnegative. Other values are errors. */
/*     Resetting FAC=0.0 suppresses the initial triangularization step. */
/*     Use of this option increases the required length of IOPT(*) by 2; */
/*     The required length of of X(*) is increased by 1. */

/*   6 */
/*   - */
/*     The norm used in testing the magnitudes of the pivot element */
/*     compared to the mass of the column above the pivot line can be */
/*     changed. The type of change that this option allows is to weight */
/*     the components with an index larger than MVAL by the parameter */
/*     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and */
/*     WT, where LP is the processing pointer for IOPT(*), */

/*       X(NCOLS+IOFF)=WT */
/*       IOPT(LP)=6 */
/*       IOPT(LP+1)=IOFF */
/*       IOPT(LP+2)=MVAL */

/*     Use of this option increases the required length of IOPT(*) by 3. */
/*     The length of X(*) is increased by 1. Values of MVAL must be */
/*     nonnegative and not greater than MINPUT. Other values are errors. */
/*     The value of WT must be positive. Any other value is an error. If */
/*     either error condition is present a message will be printed. */

/*   7 */
/*   - */
/*     Debug output, showing the detailed add-drop steps for the */
/*     constrained least squares problem, is desired. This option is */
/*     intended to be used to locate suspected bugs. */

/*   99 */
/*   -- */
/*     There are no more options to change. */

/*     The values for options are 1,...,7,99, and are the only ones */
/*     permitted. Other values are errors. Options -99,-1,...,-7 mean */
/*     that the repective options 99,1,...,7 are left at their default */
/*     values. An example is the option to modify the (rank) tolerance: */

/*       X(NCOLS+1)=TOL */
/*       IOPT(1)=-2 */
/*       IOPT(2)=1 */
/*       IOPT(3)=99 */

/*    Error Messages for SBOLSM */
/*    ----- -------- --- --------- */
/*    -22    MORE THAN ITMAX = ... ITERATIONS SOLVING BOUNDED LEAST */
/*           SQUARES PROBLEM. */

/*    -23    THE OPTION NUMBER = ... IS NOT DEFINED. */

/*    -24    THE OFFSET = ... BEYOND POSTION NCOLS = ... MUST BE POSITIVE */
/*           FOR OPTION NUMBER 2. */

/*    -25    THE TOLERANCE FOR RANK DETERMINATION = ... IS LESS THAN */
/*           MACHINE PRECISION = .... */

/*    -26    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE */
/*           FOR OPTION NUMBER 3. */

/*    -27    THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING VARIABLES */
/*           MUST BE POSITIVE. NOW = .... */

/*    -28    THE MAXIMUM NUMBER OF ITERATIONS = ... MUST BE POSITIVE. */

/*    -29    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE */
/*           FOR OPTION NUMBER 5. */

/*    -30    THE FACTOR (NCOLS/MINPUT) WHERE PRETRIANGULARIZING IS */
/*           PERFORMED MUST BE NONNEGATIVE. NOW = .... */

/*    -31    THE NUMBER OF ROWS = ... MUST BE POSITIVE. */

/*    -32    THE NUMBER OF COLUMNS = ... MUST BE POSTIVE. */

/*    -33    THE ROW DIMENSION OF W(,) = ... MUST BE .GE. THE NUMBER OF */
/*           ROWS = .... */

/*    -34    FOR J = ... THE CONSTRAINT INDICATOR MUST BE 1-4. */

/*    -35    FOR J = ... THE LOWER BOUND = ... IS .GT. THE UPPER BOUND = */
/*           .... */

/*    -36    THE INPUT ORDER OF COLUMNS = ... IS NOT BETWEEN 1 AND NCOLS */
/*           = .... */

/*    -37    THE BOUND POLARITY FLAG IN COMPONENT J = ... MUST BE */
/*           POSITIVE. NOW = .... */

/*    -38    THE ROW SEPARATOR TO APPLY WEIGHTING (...) MUST LIE BETWEEN */
/*           0 AND MINPUT = .... WEIGHT = ... MUST BE POSITIVE. */

/* ***SEE ALSO  SBOCLS, SBOLS */
/* ***ROUTINES CALLED  IVOUT, R1MACH, SAXPY, SCOPY, SDOT, SMOUT, SNRM2, */
/*                    SROT, SROTG, SSWAP, SVOUT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   821220  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   920422  Fixed usage of MINPUT.  (WRB) */
/*   901009  Editorial changes, code now reads from top to bottom.  (RWC) */
/* ***END PROLOGUE  SBOLSM */

/*     PURPOSE */
/*     ------- */
/*     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED */
/*     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS: */

/*     SOLVE E*X =  F  (LEAST SQUARES SENSE) */
/*     WITH BOUNDS ON SELECTED X VALUES. */

/*     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN */
/*     EDITING AT THE CARD 'C++'. */
/*     CHANGE THE SUBPROGRAM NAME TO DBOLSM AND THE STRINGS */
/*     /SAXPY/ TO /DAXPY/, /SCOPY/ TO /DCOPY/, */
/*     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/, */
/*     /SROT/ TO /DROT/, /SROTG/ TO /DROTG/, /R1MACH/ TO /D1MACH/, */
/*     /SVOUT/ TO /DVOUT/, /SMOUT/ TO /DMOUT/, */
/*     /SSWAP/ TO /DSWAP/, /E0/ TO /D0/, */
/*     /REAL            / TO /DOUBLE PRECISION/. */
/* ++ */



/* ***FIRST EXECUTABLE STATEMENT  SBOLSM */

/*     Verify that the problem dimensions are defined properly. */

    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --bl;
    --bu;
    --ind;
    --iopt;
    --x;
    --rw;
    --ww;
    --scl;
    --ibasis;
    --ibb;

    /* Function Body */
    if (*minput <= 0) {
	s_wsfi(&io___2);
	do_fio(&c__1, (char *)&(*minput), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 21, a__1[0] = "THE NUMBER OF ROWS = ";
	i__1[1] = 8, a__1[1] = xern1;
	i__1[2] = 18, a__1[2] = " MUST BE POSITIVE.";
	s_cat(ch__1, a__1, i__1, &c__3, (ftnlen)47);
	xermsg_("SLATEC", "SBOLSM", ch__1, &c__31, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)47);
	*mode = -31;
	return 0;
    }

    if (*ncols <= 0) {
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 24, a__1[0] = "THE NUMBER OF COLUMNS = ";
	i__1[1] = 8, a__1[1] = xern1;
	i__1[2] = 18, a__1[2] = " MUST BE POSITIVE.";
	s_cat(ch__2, a__1, i__1, &c__3, (ftnlen)50);
	xermsg_("SLATEC", "SBOLSM", ch__2, &c__32, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)50);
	*mode = -32;
	return 0;
    }

    if (*mdw < *minput) {
	s_wsfi(&io___4);
	do_fio(&c__1, (char *)&(*mdw), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&(*minput), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 28, a__2[0] = "THE ROW DIMENSION OF W(,) = ";
	i__2[1] = 8, a__2[1] = xern1;
	i__2[2] = 35, a__2[2] = " MUST BE .GE. THE NUMBER OF ROWS = ";
	i__2[3] = 8, a__2[3] = xern2;
	s_cat(ch__3, a__2, i__2, &c__4, (ftnlen)79);
	xermsg_("SLATEC", "SBOLSM", ch__3, &c__33, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)79);
	*mode = -33;
	return 0;
    }

/*     Verify that bound information is correct. */

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (ind[j] < 1 || ind[j] > 4) {
	    s_wsfi(&io___8);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___9);
	    do_fio(&c__1, (char *)&ind[j], (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 8, a__1[0] = "FOR J = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    i__1[2] = 37, a__1[2] = " THE CONSTRAINT INDICATOR MUST BE 1-4";
	    s_cat(ch__4, a__1, i__1, &c__3, (ftnlen)53);
	    xermsg_("SLATEC", "SBOLSM", ch__4, &c__34, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)53);
	    *mode = -34;
	    return 0;
	}
/* L10: */
    }

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (ind[j] == 3) {
	    if (bu[j] < bl[j]) {
		s_wsfi(&io___10);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___12);
		do_fio(&c__1, (char *)&bl[j], (ftnlen)sizeof(real));
		e_wsfi();
		s_wsfi(&io___14);
		do_fio(&c__1, (char *)&bu[j], (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__4[0] = 8, a__3[0] = "FOR J = ";
		i__4[1] = 8, a__3[1] = xern1;
		i__4[2] = 19, a__3[2] = " THE LOWER BOUND = ";
		i__4[3] = 16, a__3[3] = xern3;
		i__4[4] = 27, a__3[4] = " IS .GT. THE UPPER BOUND = ";
		i__4[5] = 16, a__3[5] = xern4;
		s_cat(ch__5, a__3, i__4, &c__6, (ftnlen)94);
		xermsg_("SLATEC", "SBOLSM", ch__5, &c__35, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)94);
		*mode = -35;
		return 0;
	    }
	}
/* L20: */
    }

/*     Check that permutation and polarity arrays have been set. */

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (ibasis[j] < 1 || ibasis[j] > *ncols) {
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&ibasis[j], (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__2[0] = 29, a__2[0] = "THE INPUT ORDER OF COLUMNS = ";
	    i__2[1] = 8, a__2[1] = xern1;
	    i__2[2] = 30, a__2[2] = " IS NOT BETWEEN 1 AND NCOLS = ";
	    i__2[3] = 8, a__2[3] = xern2;
	    s_cat(ch__6, a__2, i__2, &c__4, (ftnlen)75);
	    xermsg_("SLATEC", "SBOLSM", ch__6, &c__36, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)75);
	    *mode = -36;
	    return 0;
	}

	if (ibb[j] <= 0) {
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___18);
	    do_fio(&c__1, (char *)&ibb[j], (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__2[0] = 41, a__2[0] = "THE BOUND POLARITY FLAG IN COMPONENT J "
		    "= ";
	    i__2[1] = 8, a__2[1] = xern1;
	    i__2[2] = 26, a__2[2] = " MUST BE POSITIVE.$$NOW = ";
	    i__2[3] = 8, a__2[3] = xern2;
	    s_cat(ch__7, a__2, i__2, &c__4, (ftnlen)83);
	    xermsg_("SLATEC", "SBOLSM", ch__7, &c__37, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)83);
	    *mode = -37;
	    return 0;
	}
/* L30: */
    }

/*     Process the option array. */

    fac = .75f;
    tolind = sqrt(r1mach_(&c__4));
    tolsze = sqrt(r1mach_(&c__4));
    itmax = max(*minput,*ncols) * 5;
    wt = 1.f;
    mval = 0;
    iprint = 0;

/*     Changes to some parameters can occur through the option array, */
/*     IOPT(*).  Process this array looking carefully for input data */
/*     errors. */

    lp = 0;
    lds = 0;

/*     Test for no more options. */

L590:
    lp += lds;
    ip = iopt[lp + 1];
    jp = abs(ip);
    if (ip == 99) {
	goto L470;
    } else if (jp == 99) {
	lds = 1;
    } else if (jp == 1) {

/*         Move the IOPT(*) processing pointer. */

	if (ip > 0) {
	    lp = iopt[lp + 2] - 1;
	    lds = 0;
	} else {
	    lds = 2;
	}
    } else if (jp == 2) {

/*         Change tolerance for rank determination. */

	if (ip > 0) {
	    ioff = iopt[lp + 2];
	    if (ioff <= 0) {
		s_wsfi(&io___31);
		do_fio(&c__1, (char *)&ioff, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___32);
		do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 13, a__4[0] = "THE OFFSET = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 25, a__4[2] = " BEYOND POSITION NCOLS = ";
		i__5[3] = 8, a__4[3] = xern2;
		i__5[4] = 38, a__4[4] = " MUST BE POSITIVE FOR OPTION NUMBER"
			" 2.";
		s_cat(ch__8, a__4, i__5, &c__5, (ftnlen)92);
		xermsg_("SLATEC", "SBOLSM", ch__8, &c__24, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)92);
		*mode = -24;
		return 0;
	    }

	    tolind = x[*ncols + ioff];
	    if (tolind < r1mach_(&c__4)) {
		s_wsfi(&io___33);
		do_fio(&c__1, (char *)&tolind, (ftnlen)sizeof(real));
		e_wsfi();
		s_wsfi(&io___34);
		r__1 = r1mach_(&c__4);
		do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__2[0] = 39, a__2[0] = "THE TOLERANCE FOR RANK DETERMINATIO"
			"N = ";
		i__2[1] = 16, a__2[1] = xern3;
		i__2[2] = 34, a__2[2] = " IS LESS THAN MACHINE PRECISION = ";
		i__2[3] = 16, a__2[3] = xern4;
		s_cat(ch__9, a__2, i__2, &c__4, (ftnlen)105);
		xermsg_("SLATEC", "SBOLSM", ch__9, &c__25, &c__0, (ftnlen)6, (
			ftnlen)6, (ftnlen)105);
		*mode = -25;
	    }
	}
	lds = 2;
    } else if (jp == 3) {

/*         Change blowup factor for allowing variables to become */
/*         inactive. */

	if (ip > 0) {
	    ioff = iopt[lp + 2];
	    if (ioff <= 0) {
		s_wsfi(&io___35);
		do_fio(&c__1, (char *)&ioff, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___36);
		do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 13, a__4[0] = "THE OFFSET = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 25, a__4[2] = " BEYOND POSITION NCOLS = ";
		i__5[3] = 8, a__4[3] = xern2;
		i__5[4] = 38, a__4[4] = " MUST BE POSITIVE FOR OPTION NUMBER"
			" 3.";
		s_cat(ch__8, a__4, i__5, &c__5, (ftnlen)92);
		xermsg_("SLATEC", "SBOLSM", ch__8, &c__26, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)92);
		*mode = -26;
		return 0;
	    }

	    tolsze = x[*ncols + ioff];
	    if (tolsze <= 0.f) {
		s_wsfi(&io___37);
		do_fio(&c__1, (char *)&tolsze, (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__6[0] = 86, a__5[0] = "THE RECIPROCAL OF THE BLOW-UP FACTO"
			"R FOR REJECTING VARIABLES MUST BE POSITIVE.$$NOW = ";
		i__6[1] = 16, a__5[1] = xern3;
		s_cat(ch__10, a__5, i__6, &c__2, (ftnlen)102);
		xermsg_("SLATEC", "SBOLSM", ch__10, &c__27, &c__1, (ftnlen)6, 
			(ftnlen)6, (ftnlen)102);
		*mode = -27;
		return 0;
	    }
	}
	lds = 2;
    } else if (jp == 4) {

/*         Change the maximum number of iterations allowed. */

	if (ip > 0) {
	    itmax = iopt[lp + 2];
	    if (itmax <= 0) {
		s_wsfi(&io___38);
		do_fio(&c__1, (char *)&itmax, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__1[0] = 35, a__1[0] = "THE MAXIMUM NUMBER OF ITERATIONS = ";
		i__1[1] = 8, a__1[1] = xern1;
		i__1[2] = 18, a__1[2] = " MUST BE POSITIVE.";
		s_cat(ch__11, a__1, i__1, &c__3, (ftnlen)61);
		xermsg_("SLATEC", "SBOLSM", ch__11, &c__28, &c__1, (ftnlen)6, 
			(ftnlen)6, (ftnlen)61);
		*mode = -28;
		return 0;
	    }
	}
	lds = 2;
    } else if (jp == 5) {

/*         Change the factor for pretriangularizing the data matrix. */

	if (ip > 0) {
	    ioff = iopt[lp + 2];
	    if (ioff <= 0) {
		s_wsfi(&io___39);
		do_fio(&c__1, (char *)&ioff, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___40);
		do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 13, a__4[0] = "THE OFFSET = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 25, a__4[2] = " BEYOND POSITION NCOLS = ";
		i__5[3] = 8, a__4[3] = xern2;
		i__5[4] = 38, a__4[4] = " MUST BE POSITIVE FOR OPTION NUMBER"
			" 5.";
		s_cat(ch__8, a__4, i__5, &c__5, (ftnlen)92);
		xermsg_("SLATEC", "SBOLSM", ch__8, &c__29, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)92);
		*mode = -29;
		return 0;
	    }

	    fac = x[*ncols + ioff];
	    if (fac < 0.f) {
		s_wsfi(&io___41);
		do_fio(&c__1, (char *)&fac, (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__6[0] = 94, a__5[0] = "THE FACTOR (NCOLS/MINPUT) WHERE PRE"
			"-TRIANGULARIZING IS PERFORMED MUST BE NON-NEGATIVE.$"
			"$NOW = ";
		i__6[1] = 16, a__5[1] = xern3;
		s_cat(ch__12, a__5, i__6, &c__2, (ftnlen)110);
		xermsg_("SLATEC", "SBOLSM", ch__12, &c__30, &c__0, (ftnlen)6, 
			(ftnlen)6, (ftnlen)110);
		*mode = -30;
		return 0;
	    }
	}
	lds = 2;
    } else if (jp == 6) {

/*         Change the weighting factor (from 1.0) to apply to components */
/*         numbered .gt. MVAL (initially set to 1.)  This trick is needed */
/*         for applications of this subprogram to the heavily weighted */
/*         least squares problem that come from equality constraints. */

	if (ip > 0) {
	    ioff = iopt[lp + 2];
	    mval = iopt[lp + 3];
	    wt = x[*ncols + ioff];
	}

	if (mval < 0 || mval > *minput || wt <= 0.f) {
	    s_wsfi(&io___42);
	    do_fio(&c__1, (char *)&mval, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___43);
	    do_fio(&c__1, (char *)&(*minput), (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___44);
	    do_fio(&c__1, (char *)&wt, (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__7[0] = 38, a__6[0] = "THE ROW SEPARATOR TO APPLY WEIGHTING (";
	    i__7[1] = 8, a__6[1] = xern1;
	    i__7[2] = 34, a__6[2] = ") MUST LIE BETWEEN 0 AND MINPUT = ";
	    i__7[3] = 8, a__6[3] = xern2;
	    i__7[4] = 12, a__6[4] = ".$$WEIGHT = ";
	    i__7[5] = 16, a__6[5] = xern3;
	    i__7[6] = 18, a__6[6] = " MUST BE POSITIVE.";
	    s_cat(ch__13, a__6, i__7, &c__7, (ftnlen)134);
	    xermsg_("SLATEC", "SBOLSM", ch__13, &c__38, &c__0, (ftnlen)6, (
		    ftnlen)6, (ftnlen)134);
	    *mode = -38;
	    return 0;
	}
	lds = 3;
    } else if (jp == 7) {

/*         Turn on debug output. */

	if (ip > 0) {
	    iprint = 1;
	}
	lds = 2;
    } else {
	s_wsfi(&io___45);
	do_fio(&c__1, (char *)&ip, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 20, a__1[0] = "THE OPTION NUMBER = ";
	i__1[1] = 8, a__1[1] = xern1;
	i__1[2] = 16, a__1[2] = " IS NOT DEFINED.";
	s_cat(ch__14, a__1, i__1, &c__3, (ftnlen)44);
	xermsg_("SLATEC", "SBOLSM", ch__14, &c__23, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)44);
	*mode = -23;
	return 0;
    }
    goto L590;

/*     Pretriangularize rectangular arrays of certain sizes for */
/*     increased efficiency. */

L470:
    if (fac * *minput > (real) (*ncols)) {
	i__3 = *ncols + 1;
	for (j = 1; j <= i__3; ++j) {
	    i__8 = j + mval + 1;
	    for (i__ = *minput; i__ >= i__8; --i__) {
		srotg_(&w[i__ - 1 + j * w_dim1], &w[i__ + j * w_dim1], &sc, &
			ss);
		w[i__ + j * w_dim1] = 0.f;
		i__9 = *ncols - j + 1;
		srot_(&i__9, &w[i__ - 1 + (j + 1) * w_dim1], mdw, &w[i__ + (j 
			+ 1) * w_dim1], mdw, &sc, &ss);
/* L480: */
	    }
/* L490: */
	}
	mrows = *ncols + mval + 1;
    } else {
	mrows = *minput;
    }

/*     Set the X(*) array to zero so all components are defined. */

    scopy_(ncols, &c_b185, &c__0, &x[1], &c__1);

/*     The arrays IBASIS(*) and IBB(*) are initialized by the calling */
/*     program and the column scaling is defined in the calling program. */
/*     'BIG' is plus infinity on this machine. */

    big = r1mach_(&c__2);
    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (ind[j] == 1) {
	    bu[j] = big;
	} else if (ind[j] == 2) {
	    bl[j] = -big;
	} else if (ind[j] == 4) {
	    bl[j] = -big;
	    bu[j] = big;
	}
/* L550: */
    }

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (bl[j] <= 0.f && 0.f <= bu[j] && (r__1 = bu[j], dabs(r__1)) < (
		r__2 = bl[j], dabs(r__2)) || bu[j] < 0.f) {
	    t = bu[j];
	    bu[j] = -bl[j];
	    bl[j] = -t;
	    scl[j] = -scl[j];
	    i__8 = mrows;
	    for (i__ = 1; i__ <= i__8; ++i__) {
		w[i__ + j * w_dim1] = -w[i__ + j * w_dim1];
/* L560: */
	    }
	}

/*         Indices in set T(=TIGHT) are denoted by negative values */
/*         of IBASIS(*). */

	if (bl[j] >= 0.f) {
	    ibasis[j] = -ibasis[j];
	    t = -bl[j];
	    bu[j] += t;
	    saxpy_(&mrows, &t, &w[j * w_dim1 + 1], &c__1, &w[(*ncols + 1) * 
		    w_dim1 + 1], &c__1);
	}
/* L570: */
    }

    nsetb = 0;
    iter = 0;

    if (iprint > 0) {
	i__3 = *ncols + 1;
	smout_(&mrows, &i__3, mdw, &w[w_offset], "(' PRETRI. INPUT MATRIX')", 
		&c_n4, (ftnlen)25);
	svout_(ncols, &bl[1], "(' LOWER BOUNDS')", &c_n4, (ftnlen)17);
	svout_(ncols, &bu[1], "(' UPPER BOUNDS')", &c_n4, (ftnlen)17);
    }

L580:
    ++iter;
    if (iter > itmax) {
	s_wsfi(&io___54);
	do_fio(&c__1, (char *)&itmax, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 18, a__1[0] = "MORE THAN ITMAX = ";
	i__1[1] = 8, a__1[1] = xern1;
	i__1[2] = 50, a__1[2] = " ITERATIONS SOLVING BOUNDED LEAST SQUARES P"
		"ROBLEM.";
	s_cat(ch__15, a__1, i__1, &c__3, (ftnlen)76);
	xermsg_("SLATEC", "SBOLSM", ch__15, &c__22, &c__1, (ftnlen)6, (ftnlen)
		6, (ftnlen)76);
	*mode = -22;

/*        Rescale and translate variables. */

	igopr = 1;
	goto L130;
    }

/*     Find a variable to become non-active. */
/*                                                 T */
/*     Compute (negative) of gradient vector, W = E *(F-E*X). */

    scopy_(ncols, &c_b185, &c__0, &ww[1], &c__1);
    i__3 = *ncols;
    for (j = nsetb + 1; j <= i__3; ++j) {
	jcol = (i__8 = ibasis[j], abs(i__8));
	i__8 = mrows - nsetb;
/* Computing MIN */
	i__9 = nsetb + 1;
/* Computing MIN */
	i__10 = nsetb + 1;
	ww[j] = sdot_(&i__8, &w[min(i__9,mrows) + j * w_dim1], &c__1, &w[min(
		i__10,mrows) + (*ncols + 1) * w_dim1], &c__1) * (r__1 = scl[
		jcol], dabs(r__1));
/* L200: */
    }

    if (iprint > 0) {
	svout_(ncols, &ww[1], "(' GRADIENT VALUES')", &c_n4, (ftnlen)20);
	ivout_(ncols, &ibasis[1], "(' INTERNAL VARIABLE ORDER')", &c_n4, (
		ftnlen)28);
	ivout_(ncols, &ibb[1], "(' BOUND POLARITY')", &c_n4, (ftnlen)19);
    }

/*     If active set = number of total rows, quit. */

L210:
    if (nsetb == mrows) {
	found = FALSE_;
	goto L120;
    }

/*     Choose an extremal component of gradient vector for a candidate */
/*     to become non-active. */

    wlarge = -big;
    wmag = -big;
    i__3 = *ncols;
    for (j = nsetb + 1; j <= i__3; ++j) {
	t = ww[j];
	if (t == big) {
	    goto L220;
	}
	itemp = ibasis[j];
	jcol = abs(itemp);
	i__8 = mval - nsetb;
/* Computing MIN */
	i__9 = nsetb + 1;
	t1 = snrm2_(&i__8, &w[min(i__9,mrows) + j * w_dim1], &c__1);
	if (itemp < 0) {
	    if (ibb[jcol] % 2 == 0) {
		t = -t;
	    }
	    if (t < 0.f) {
		goto L220;
	    }
	    if (mval > nsetb) {
		t = t1;
	    }
	    if (t > wlarge) {
		wlarge = t;
		jlarge = j;
	    }
	} else {
	    if (mval > nsetb) {
		t = t1;
	    }
	    if (dabs(t) > wmag) {
		wmag = dabs(t);
		jmag = j;
	    }
	}
L220:
	;
    }

/*     Choose magnitude of largest component of gradient for candidate. */

    jbig = 0;
    wbig = 0.f;
    if (wlarge > 0.f) {
	jbig = jlarge;
	wbig = wlarge;
    }

    if (wmag >= wbig) {
	jbig = jmag;
	wbig = wmag;
    }

    if (jbig == 0) {
	found = FALSE_;
	if (iprint > 0) {
	    ivout_(&c__0, &i__, "(' FOUND NO VARIABLE TO ENTER')", &c_n4, (
		    ftnlen)31);
	}
	goto L120;
    }

/*     See if the incoming column is sufficiently independent.  This */
/*     test is made before an elimination is performed. */

    if (iprint > 0) {
	ivout_(&c__1, &jbig, "(' TRY TO BRING IN THIS COL.')", &c_n4, (ftnlen)
		30);
    }

    if (mval <= nsetb) {
	cl1 = snrm2_(&mval, &w[jbig * w_dim1 + 1], &c__1);
	i__3 = nsetb - mval;
/* Computing MIN */
	i__8 = mval + 1;
	cl2 = dabs(wt) * snrm2_(&i__3, &w[min(i__8,mrows) + jbig * w_dim1], &
		c__1);
	i__3 = mrows - nsetb;
/* Computing MIN */
	i__8 = nsetb + 1;
	cl3 = dabs(wt) * snrm2_(&i__3, &w[min(i__8,mrows) + jbig * w_dim1], &
		c__1);
	srotg_(&cl1, &cl2, &sc, &ss);
	colabv = dabs(cl1);
	colblo = cl3;
    } else {
	cl1 = snrm2_(&nsetb, &w[jbig * w_dim1 + 1], &c__1);
	i__3 = mval - nsetb;
/* Computing MIN */
	i__8 = nsetb + 1;
	cl2 = snrm2_(&i__3, &w[min(i__8,mrows) + jbig * w_dim1], &c__1);
	i__3 = mrows - mval;
/* Computing MIN */
	i__8 = mval + 1;
	cl3 = dabs(wt) * snrm2_(&i__3, &w[min(i__8,mrows) + jbig * w_dim1], &
		c__1);
	colabv = cl1;
	srotg_(&cl2, &cl3, &sc, &ss);
	colblo = dabs(cl2);
    }

    if (colblo <= tolind * colabv) {
	ww[jbig] = big;
	if (iprint > 0) {
	    ivout_(&c__0, &i__, "(' VARIABLE IS DEPENDENT, NOT USED.')", &
		    c_n4, (ftnlen)37);
	}
	goto L210;
    }

/*     Swap matrix columns NSETB+1 and JBIG, plus pointer information, */
/*     and gradient values. */

    ++nsetb;
    if (nsetb != jbig) {
	sswap_(&mrows, &w[nsetb * w_dim1 + 1], &c__1, &w[jbig * w_dim1 + 1], &
		c__1);
	sswap_(&c__1, &ww[nsetb], &c__1, &ww[jbig], &c__1);
	itemp = ibasis[nsetb];
	ibasis[nsetb] = ibasis[jbig];
	ibasis[jbig] = itemp;
    }

/*     Eliminate entries below the pivot line in column NSETB. */

    if (mrows > nsetb) {
	i__3 = nsetb + 1;
	for (i__ = mrows; i__ >= i__3; --i__) {
	    if (i__ == mval + 1) {
		goto L230;
	    }
	    srotg_(&w[i__ - 1 + nsetb * w_dim1], &w[i__ + nsetb * w_dim1], &
		    sc, &ss);
	    w[i__ + nsetb * w_dim1] = 0.f;
	    i__8 = *ncols - nsetb + 1;
	    srot_(&i__8, &w[i__ - 1 + (nsetb + 1) * w_dim1], mdw, &w[i__ + (
		    nsetb + 1) * w_dim1], mdw, &sc, &ss);
L230:
	    ;
	}

	if (mval >= nsetb && mval < mrows) {
	    srotg_(&w[nsetb + nsetb * w_dim1], &w[mval + 1 + nsetb * w_dim1], 
		    &sc, &ss);
	    w[mval + 1 + nsetb * w_dim1] = 0.f;
	    i__3 = *ncols - nsetb + 1;
	    srot_(&i__3, &w[nsetb + (nsetb + 1) * w_dim1], mdw, &w[mval + 1 + 
		    (nsetb + 1) * w_dim1], mdw, &sc, &ss);
	}
    }

    if (w[nsetb + nsetb * w_dim1] == 0.f) {
	ww[nsetb] = big;
	--nsetb;
	if (iprint > 0) {
	    ivout_(&c__0, &i__, "(' PIVOT IS ZERO, NOT USED.')", &c_n4, (
		    ftnlen)29);
	}
	goto L210;
    }

/*     Check that new variable is moving in the right direction. */

    itemp = ibasis[nsetb];
    jcol = abs(itemp);
    xnew = w[nsetb + (*ncols + 1) * w_dim1] / w[nsetb + nsetb * w_dim1] / (
	    r__1 = scl[jcol], dabs(r__1));
    if (itemp < 0) {

/*         IF(WW(NSETB).GE.ZERO.AND.XNEW.LE.ZERO) exit(quit) */
/*         IF(WW(NSETB).LE.ZERO.AND.XNEW.GE.ZERO) exit(quit) */

	if (ww[nsetb] >= 0.f && xnew <= 0.f || ww[nsetb] <= 0.f && xnew >= 
		0.f) {
	    goto L240;
	}
    }
    found = TRUE_;
    goto L120;

L240:
    ww[nsetb] = big;
    --nsetb;
    if (iprint > 0) {
	ivout_(&c__0, &i__, "(' VARIABLE HAS BAD DIRECTION, NOT USED.')", &
		c_n4, (ftnlen)42);
    }
    goto L210;

/*     Solve the triangular system. */

L270:
    scopy_(&nsetb, &w[(*ncols + 1) * w_dim1 + 1], &c__1, &rw[1], &c__1);
    for (j = nsetb; j >= 1; --j) {
	rw[j] /= w[j + j * w_dim1];
	jcol = (i__3 = ibasis[j], abs(i__3));
	t = rw[j];
	if (ibb[jcol] % 2 == 0) {
	    rw[j] = -rw[j];
	}
	i__3 = j - 1;
	r__1 = -t;
	saxpy_(&i__3, &r__1, &w[j * w_dim1 + 1], &c__1, &rw[1], &c__1);
	rw[j] /= (r__1 = scl[jcol], dabs(r__1));
/* L280: */
    }

    if (iprint > 0) {
	svout_(&nsetb, &rw[1], "(' SOLN. VALUES')", &c_n4, (ftnlen)17);
	ivout_(&nsetb, &ibasis[1], "(' COLS. USED')", &c_n4, (ftnlen)15);
    }

    if (lgopr == 2) {
	scopy_(&nsetb, &rw[1], &c__1, &x[1], &c__1);
	i__3 = nsetb;
	for (j = 1; j <= i__3; ++j) {
	    itemp = ibasis[j];
	    jcol = abs(itemp);
	    if (itemp < 0) {
		bou = 0.f;
	    } else {
		bou = bl[jcol];
	    }

	    if (-bou != big) {
		bou /= (r__1 = scl[jcol], dabs(r__1));
	    }
	    if (x[j] <= bou) {
		jdrop1 = j;
		goto L340;
	    }

	    bou = bu[jcol];
	    if (bou != big) {
		bou /= (r__1 = scl[jcol], dabs(r__1));
	    }
	    if (x[j] >= bou) {
		jdrop2 = j;
		goto L340;
	    }
/* L450: */
	}
	goto L340;
    }

/*     See if the unconstrained solution (obtained by solving the */
/*     triangular system) satisfies the problem bounds. */

    alpha = 2.f;
    beta = 2.f;
    x[nsetb] = 0.f;
    i__3 = nsetb;
    for (j = 1; j <= i__3; ++j) {
	itemp = ibasis[j];
	jcol = abs(itemp);
	t1 = 2.f;
	t2 = 2.f;
	if (itemp < 0) {
	    bou = 0.f;
	} else {
	    bou = bl[jcol];
	}
	if (-bou != big) {
	    bou /= (r__1 = scl[jcol], dabs(r__1));
	}
	if (rw[j] <= bou) {
	    t1 = (x[j] - bou) / (x[j] - rw[j]);
	}
	bou = bu[jcol];
	if (bou != big) {
	    bou /= (r__1 = scl[jcol], dabs(r__1));
	}
	if (rw[j] >= bou) {
	    t2 = (bou - x[j]) / (rw[j] - x[j]);
	}

/*     If not, then compute a step length so that the variables remain */
/*     feasible. */

	if (t1 < alpha) {
	    alpha = t1;
	    jdrop1 = j;
	}

	if (t2 < beta) {
	    beta = t2;
	    jdrop2 = j;
	}
/* L310: */
    }

    constr = alpha < 2.f || beta < 2.f;
    if (! constr) {

/*         Accept the candidate because it satisfies the stated bounds */
/*         on the variables. */

	scopy_(&nsetb, &rw[1], &c__1, &x[1], &c__1);
	goto L580;
    }

/*     Take a step that is as large as possible with all variables */
/*     remaining feasible. */

    i__3 = nsetb;
    for (j = 1; j <= i__3; ++j) {
	x[j] += dmin(alpha,beta) * (rw[j] - x[j]);
/* L330: */
    }

    if (alpha <= beta) {
	jdrop2 = 0;
    } else {
	jdrop1 = 0;
    }

L340:
    if (jdrop1 + jdrop2 <= 0 || nsetb <= 0) {
	goto L580;
    }
/* L350: */
    jdrop = jdrop1 + jdrop2;
    itemp = ibasis[jdrop];
    jcol = abs(itemp);
    if (jdrop2 > 0) {

/*         Variable is at an upper bound.  Subtract multiple of this */
/*         column from right hand side. */

	t = bu[jcol];
	if (itemp > 0) {
	    bu[jcol] = t - bl[jcol];
	    bl[jcol] = -t;
	    itemp = -itemp;
	    scl[jcol] = -scl[jcol];
	    i__3 = jdrop;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		w[i__ + jdrop * w_dim1] = -w[i__ + jdrop * w_dim1];
/* L360: */
	    }
	} else {
	    ++ibb[jcol];
	    if (ibb[jcol] % 2 == 0) {
		t = -t;
	    }
	}

/*     Variable is at a lower bound. */

    } else {
	if ((real) itemp < 0.f) {
	    t = 0.f;
	} else {
	    t = -bl[jcol];
	    bu[jcol] += t;
	    itemp = -itemp;
	}
    }

    saxpy_(&jdrop, &t, &w[jdrop * w_dim1 + 1], &c__1, &w[(*ncols + 1) * 
	    w_dim1 + 1], &c__1);

/*     Move certain columns left to achieve upper Hessenberg form. */

    scopy_(&jdrop, &w[jdrop * w_dim1 + 1], &c__1, &rw[1], &c__1);
    i__3 = nsetb;
    for (j = jdrop + 1; j <= i__3; ++j) {
	ibasis[j - 1] = ibasis[j];
	x[j - 1] = x[j];
	scopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[(j - 1) * w_dim1 + 1], &c__1)
		;
/* L370: */
    }

    ibasis[nsetb] = itemp;
    w[nsetb * w_dim1 + 1] = 0.f;
    i__3 = mrows - jdrop;
    scopy_(&i__3, &w[nsetb * w_dim1 + 1], &c__0, &w[jdrop + 1 + nsetb * 
	    w_dim1], &c__1);
    scopy_(&jdrop, &rw[1], &c__1, &w[nsetb * w_dim1 + 1], &c__1);

/*     Transform the matrix from upper Hessenberg form to upper */
/*     triangular form. */

    --nsetb;
    i__3 = nsetb;
    for (i__ = jdrop; i__ <= i__3; ++i__) {

/*         Look for small pivots and avoid mixing weighted and */
/*         nonweighted rows. */

	if (i__ == mval) {
	    t = 0.f;
	    i__8 = nsetb;
	    for (j = i__; j <= i__8; ++j) {
		jcol = (i__9 = ibasis[j], abs(i__9));
		t1 = (r__1 = w[i__ + j * w_dim1] * scl[jcol], dabs(r__1));
		if (t1 > t) {
		    jbig = j;
		    t = t1;
		}
/* L380: */
	    }
	    goto L400;
	}
	srotg_(&w[i__ + i__ * w_dim1], &w[i__ + 1 + i__ * w_dim1], &sc, &ss);
	w[i__ + 1 + i__ * w_dim1] = 0.f;
	i__8 = *ncols - i__ + 1;
	srot_(&i__8, &w[i__ + (i__ + 1) * w_dim1], mdw, &w[i__ + 1 + (i__ + 1)
		 * w_dim1], mdw, &sc, &ss);
/* L390: */
    }
    goto L430;

/*     The triangularization is completed by giving up the Hessenberg */
/*     form and triangularizing a rectangular matrix. */

L400:
    sswap_(&mrows, &w[i__ * w_dim1 + 1], &c__1, &w[jbig * w_dim1 + 1], &c__1);
    sswap_(&c__1, &ww[i__], &c__1, &ww[jbig], &c__1);
    sswap_(&c__1, &x[i__], &c__1, &x[jbig], &c__1);
    itemp = ibasis[i__];
    ibasis[i__] = ibasis[jbig];
    ibasis[jbig] = itemp;
    jbig = i__;
    i__3 = nsetb;
    for (j = jbig; j <= i__3; ++j) {
	i__8 = mrows;
	for (i__ = j + 1; i__ <= i__8; ++i__) {
	    srotg_(&w[j + j * w_dim1], &w[i__ + j * w_dim1], &sc, &ss);
	    w[i__ + j * w_dim1] = 0.f;
	    i__9 = *ncols - j + 1;
	    srot_(&i__9, &w[j + (j + 1) * w_dim1], mdw, &w[i__ + (j + 1) * 
		    w_dim1], mdw, &sc, &ss);
/* L410: */
	}
/* L420: */
    }

/*     See if the remaining coefficients are feasible.  They should be */
/*     because of the way MIN(ALPHA,BETA) was chosen.  Any that are not */
/*     feasible will be set to their bounds and appropriately translated. */

L430:
    jdrop1 = 0;
    jdrop2 = 0;
    lgopr = 2;
    goto L270;

/*     Find a variable to become non-active. */

L120:
    if (found) {
	lgopr = 1;
	goto L270;
    }

/*     Rescale and translate variables. */

    igopr = 2;
L130:
    scopy_(&nsetb, &x[1], &c__1, &rw[1], &c__1);
    scopy_(ncols, &c_b185, &c__0, &x[1], &c__1);
    i__3 = nsetb;
    for (j = 1; j <= i__3; ++j) {
	jcol = (i__8 = ibasis[j], abs(i__8));
	x[jcol] = rw[j] * (r__1 = scl[jcol], dabs(r__1));
/* L140: */
    }

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (ibb[j] % 2 == 0) {
	    x[j] = bu[j] - x[j];
	}
/* L150: */
    }

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	jcol = ibasis[j];
	if (jcol < 0) {
	    x[-jcol] = bl[-jcol] + x[-jcol];
	}
/* L160: */
    }

    i__3 = *ncols;
    for (j = 1; j <= i__3; ++j) {
	if (scl[j] < 0.f) {
	    x[j] = -x[j];
	}
/* L170: */
    }

    i__ = max(nsetb,mval);
    i__3 = mrows - i__;
/* Computing MIN */
    i__8 = i__ + 1;
    *rnorm = snrm2_(&i__3, &w[min(i__8,mrows) + (*ncols + 1) * w_dim1], &c__1)
	    ;

    if (igopr == 2) {
	*mode = nsetb;
    }
    return 0;
} /* sbolsm_ */


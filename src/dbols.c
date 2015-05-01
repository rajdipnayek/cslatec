/* dbols.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__8 = 8;
static integer c__7 = 7;
static integer c__9 = 9;
static integer c__6 = 6;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__17 = 17;
static integer c__10 = 10;
static integer c__0 = 0;

/* DECK DBOLS */
/* Subroutine */ int dbols_(doublereal *w, integer *mdw, integer *mrows, 
	integer *ncols, doublereal *bl, doublereal *bu, integer *ind, integer 
	*iopt, doublereal *x, doublereal *rnorm, integer *mode, doublereal *
	rw, integer *iw)
{
    /* Initialized data */

    static integer igo = 0;

    /* System generated locals */
    address a__1[3], a__2[5], a__3[8], a__4[4];
    integer w_dim1, w_offset, i__1[3], i__2, i__3[5], i__4[8], i__5[4], i__6, 
	    i__7;
    doublereal d__1;
    char ch__1[32], ch__2[55], ch__3[37], ch__4[77], ch__5[36], ch__6[82], 
	    ch__7[94], ch__8[44], ch__9[79], ch__10[71], ch__11[91], ch__12[
	    75], ch__13[64], ch__14[63];

    /* Local variables */
    static integer i__, j;
    static doublereal sc;
    static integer ip, jp, lp;
    static doublereal ss;
    static integer llb;
    static doublereal one;
    static integer lds, llx, ibig, lmdw, lndw, nerr, lenx, lliw, mnew;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer lopt;
    static doublereal zero;
    static integer llrw;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static char xern1[8], xern2[8], xern3[16], xern4[16];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), drotg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer liopt, locacc;
    static logical checkl;
    static integer iscale;
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer locdim;
    extern /* Subroutine */ int dbolsm_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *), xermsg_(char *, char *, 
	    char *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer inrows;

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___10 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___12 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___34 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___35 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___36 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___37 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___38 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___39 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___40 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___41 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___42 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___43 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___44 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___45 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___46 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___47 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___48 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___49 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___50 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___51 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___52 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___55 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___56 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DBOLS */
/* ***PURPOSE  Solve the problem */
/*                 E*X = F (in the least  squares  sense) */
/*            with bounds on selected X values. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1A2A, G2E, G2H1, G2H2 */
/* ***TYPE      DOUBLE PRECISION (SBOLS-S, DBOLS-D) */
/* ***KEYWORDS  BOUNDS, CONSTRAINTS, INEQUALITY, LEAST SQUARES, LINEAR */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION **** */

/*     The user must have dimension statements of the form: */

/*       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS), */
/*      * X(NCOLS+NX), RW(5*NCOLS) */
/*       INTEGER IND(NCOLS), IOPT(1+NI), IW(2*NCOLS) */

/*     (Here NX=number of extra locations required for option 4; NX=0 */
/*     for no options; NX=NCOLS if this option is in use. Here NI=number */
/*     of extra locations required for options 1-6; NI=0 for no */
/*     options.) */

/*   INPUT */
/*   ----- */

/*    -------------------- */
/*    W(MDW,*),MROWS,NCOLS */
/*    -------------------- */
/*     The array W(*,*) contains the matrix [E:F] on entry. The matrix */
/*     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in */
/*     the array W(*,*) with E occupying the first NCOLS columns and the */
/*     right side vector F in column NCOLS+1. The row dimension, MDW, of */
/*     the array W(*,*) must satisfy the inequality MDW .ge. MROWS. */
/*     Other values of MDW are errors. The values of MROWS and NCOLS */
/*     must be positive. Other values are errors. There is an exception */
/*     to this when using option 1 for accumulation of blocks of */
/*     equations. In that case MROWS is an OUTPUT variable ONLY, and the */
/*     matrix data for [E:F] is placed in W(*,*), one block of rows at a */
/*     time.  MROWS contains the number of rows in the matrix after */
/*     triangularizing several blocks of equations. This is an OUTPUT */
/*     parameter ONLY when option 1 is used. See IOPT(*) CONTENTS */
/*     for details about option 1. */

/*    ------------------ */
/*    BL(*),BU(*),IND(*) */
/*    ------------------ */
/*     These arrays contain the information about the bounds that the */
/*     solution values are to satisfy. The value of IND(J) tells the */
/*     type of bound and BL(J) and BU(J) give the explicit values for */
/*     the respective upper and lower bounds. */

/*    1.    For IND(J)=1, require X(J) .ge. BL(J). */
/*          (the value of BU(J) is not used.) */
/*    2.    For IND(J)=2, require X(J) .le. BU(J). */
/*          (the value of BL(J) is not used.) */
/*    3.    For IND(J)=3, require X(J) .ge. BL(J) and */
/*                                X(J) .le. BU(J). */
/*    4.    For IND(J)=4, no bounds on X(J) are required. */
/*          (the values of BL(J) and BU(J) are not used.) */

/*     Values other than 1,2,3 or 4 for IND(J) are errors. In the case */
/*     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J) */
/*     is an error. */

/*    ------- */
/*    IOPT(*) */
/*    ------- */
/*     This is the array where the user can specify nonstandard options */
/*     for DBOLSM( ). Most of the time this feature can be ignored by */
/*     setting the input value IOPT(1)=99. Occasionally users may have */
/*     needs that require use of the following subprogram options. For */
/*     details about how to use the options see below: IOPT(*) CONTENTS. */

/*     Option Number   Brief Statement of Purpose */
/*     ------ ------   ----- --------- -- ------- */
/*           1         Return to user for accumulation of blocks */
/*                     of least squares equations. */
/*           2         Check lengths of all arrays used in the */
/*                     subprogram. */
/*           3         Standard scaling of the data matrix, E. */
/*           4         User provides column scaling for matrix E. */
/*           5         Provide option array to the low-level */
/*                     subprogram DBOLSM( ). */
/*           6         Move the IOPT(*) processing pointer. */
/*          99         No more options to change. */

/*    ---- */
/*    X(*) */
/*    ---- */
/*     This array is used to pass data associated with option 4. Ignore */
/*     this parameter if this option is not used. Otherwise see below: */
/*     IOPT(*) CONTENTS. */

/*    OUTPUT */
/*    ------ */

/*    ---------- */
/*    X(*),RNORM */
/*    ---------- */
/*     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for */
/*     the constrained least squares problem. The value RNORM is the */
/*     minimum residual vector length. */

/*    ---- */
/*    MODE */
/*    ---- */
/*     The sign of MODE determines whether the subprogram has completed */
/*     normally, or encountered an error condition or abnormal status. A */
/*     value of MODE .ge. 0 signifies that the subprogram has completed */
/*     normally. The value of MODE (.GE. 0) is the number of variables */
/*     in an active status: not at a bound nor at the value ZERO, for */
/*     the case of free variables. A negative value of MODE will be one */
/*     of the cases -37,-36,...,-22, or -17,...,-2. Values .lt. -1 */
/*     correspond to an abnormal completion of the subprogram. To */
/*     understand the abnormal completion codes see below: ERROR */
/*     MESSAGES for DBOLS( ). AN approximate solution will be returned */
/*     to the user only when max. iterations is reached, MODE=-22. */
/*     Values for MODE=-37,...,-22 come from the low-level subprogram */
/*     DBOLSM(). See the section ERROR MESSAGES for DBOLSM() in the */
/*     documentation for DBOLSM(). */

/*    ----------- */
/*    RW(*),IW(*) */
/*    ----------- */
/*     These are working arrays with 5*NCOLS and 2*NCOLS entries. */
/*     (normally the user can ignore the contents of these arrays, */
/*     but they must be dimensioned properly.) */

/*    IOPT(*) CONTENTS */
/*    ------- -------- */
/*     The option array allows a user to modify internal variables in */
/*     the subprogram without recompiling the source code. A central */
/*     goal of the initial software design was to do a good job for most */
/*     people. Thus the use of options will be restricted to a select */
/*     group of users. The processing of the option array proceeds as */
/*     follows: a pointer, here called LP, is initially set to the value */
/*     1. This value is updated as each option is processed. At the */
/*     pointer position the option number is extracted and used for */
/*     locating other information that allows for options to be changed. */
/*     The portion of the array IOPT(*) that is used for each option is */
/*     fixed; the user and the subprogram both know how many locations */
/*     are needed for each option. A great deal of error checking is */
/*     done by the subprogram on the contents of the option array. */
/*     Nevertheless it is still possible to give the subprogram optional */
/*     input that is meaningless. For example option 4 uses the */
/*     locations X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing */
/*     scaling data. The user must manage the allocation of these */
/*     locations. */

/*   1 */
/*   - */
/*     This option allows the user to solve problems with a large number */
/*     of rows compared to the number of variables. The idea is that the */
/*     subprogram returns to the user (perhaps many times) and receives */
/*     new least squares equations from the calling program unit. */
/*     Eventually the user signals "that's all" and then computes the */
/*     solution with one final call to subprogram DBOLS( ). The value of */
/*     MROWS is an OUTPUT variable when this option is used. Its value */
/*     is always in the range 0 .le. MROWS .le. NCOLS+1. It is equal to */
/*     the number of rows after the triangularization of the entire set */
/*     of equations. If LP is the processing pointer for IOPT(*), the */
/*     usage for the sequential processing of blocks of equations is */

/*        IOPT(LP)=1 */
/*        Move block of equations to W(*,*) starting at */
/*        the first row of W(*,*). */
/*        IOPT(LP+3)=# of rows in the block; user defined */

/*     The user now calls DBOLS( ) in a loop. The value of IOPT(LP+1) */
/*     directs the user's action. The value of IOPT(LP+2) points to */
/*     where the subsequent rows are to be placed in W(*,*). */

/*      .<LOOP */
/*      . CALL DBOLS() */
/*      . IF(IOPT(LP+1) .EQ. 1) THEN */
/*      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED */
/*      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN */
/*      .    W(*,*) STARTING AT ROW IOPT(LP+2). */
/*      . */
/*      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN */
/*      .       IOPT(LP+1)=2 */
/*      .<------CYCLE LOOP */
/*      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN */
/*      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0 */
/*      . ELSE */
/*      . ERROR CONDITION; SHOULD NOT HAPPEN. */
/*      .<END LOOP */

/*     Use of this option adds 4 to the required length of IOPT(*). */


/*   2 */
/*   - */
/*     This option is useful for checking the lengths of all arrays used */
/*     by DBOLS() against their actual requirements for this problem. */
/*     The idea is simple: the user's program unit passes the declared */
/*     dimension information of the arrays. These values are compared */
/*     against the problem-dependent needs within the subprogram. If any */
/*     of the dimensions are too small an error message is printed and a */
/*     negative value of MODE is returned, -11 to -17. The printed error */
/*     message tells how long the dimension should be. If LP is the */
/*     processing pointer for IOPT(*), */

/*        IOPT(LP)=2 */
/*        IOPT(LP+1)=Row dimension of W(*,*) */
/*        IOPT(LP+2)=Col. dimension of W(*,*) */
/*        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*) */
/*        IOPT(LP+4)=Dimension of X(*) */
/*        IOPT(LP+5)=Dimension of RW(*) */
/*        IOPT(LP+6)=Dimension of IW(*) */
/*        IOPT(LP+7)=Dimension of IOPT(*) */
/*         . */
/*        CALL DBOLS() */

/*     Use of this option adds 8 to the required length of IOPT(*). */

/*   3 */
/*   - */
/*     This option changes the type of scaling for the data matrix E. */
/*     Nominally each nonzero column of E is scaled so that the */
/*     magnitude of its largest entry is equal to the value ONE. If LP */
/*     is the processing pointer for IOPT(*), */

/*        IOPT(LP)=3 */
/*        IOPT(LP+1)=1,2 or 3 */
/*            1= Nominal scaling as noted; */
/*            2= Each nonzero column scaled to have length ONE; */
/*            3= Identity scaling; scaling effectively suppressed. */
/*         . */
/*        CALL DBOLS() */

/*     Use of this option adds 2 to the required length of IOPT(*). */

/*   4 */
/*   - */
/*     This option allows the user to provide arbitrary (positive) */
/*     column scaling for the matrix E. If LP is the processing pointer */
/*     for IOPT(*), */

/*        IOPT(LP)=4 */
/*        IOPT(LP+1)=IOFF */
/*        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) */
/*        = Positive scale factors for cols. of E. */
/*         . */
/*        CALL DBOLS() */

/*     Use of this option adds 2 to the required length of IOPT(*) and */
/*     NCOLS to the required length of X(*). */

/*   5 */
/*   - */
/*     This option allows the user to provide an option array to the */
/*     low-level subprogram DBOLSM(). If LP is the processing pointer */
/*     for IOPT(*), */

/*        IOPT(LP)=5 */
/*        IOPT(LP+1)= Position in IOPT(*) where option array */
/*                    data for DBOLSM() begins. */
/*         . */
/*        CALL DBOLS() */

/*     Use of this option adds 2 to the required length of IOPT(*). */

/*   6 */
/*   - */
/*     Move the processing pointer (either forward or backward) to the */
/*     location IOPT(LP+1). The processing point is moved to entry */
/*     LP+2 of IOPT(*) if the option is left with -6 in IOPT(LP).  For */
/*     example to skip over locations 3,...,NCOLS+2 of IOPT(*), */

/*       IOPT(1)=6 */
/*       IOPT(2)=NCOLS+3 */
/*       (IOPT(I), I=3,...,NCOLS+2 are not defined here.) */
/*       IOPT(NCOLS+3)=99 */
/*       CALL DBOLS() */

/*     CAUTION: Misuse of this option can yield some very hard */
/*     -to-find bugs.  Use it with care. */

/*   99 */
/*   -- */
/*     There are no more options to change. */

/*     Only option numbers -99, -6,-5,...,-1, 1,2,...,6, and 99 are */
/*     permitted. Other values are errors. Options -99,-1,...,-6 mean */
/*     that the respective options 99,1,...,6 are left at their default */
/*     values. An example is the option to modify the (rank) tolerance: */

/*       IOPT(1)=-3 Option is recognized but not changed */
/*       IOPT(2)=2  Scale nonzero cols. to have length ONE */
/*       IOPT(3)=99 */

/*    ERROR MESSAGES for DBOLS() */
/*    ----- -------- --- ------- */

/* WARNING IN... */
/* DBOLS(). MDW=(I1) MUST BE POSITIVE. */
/*           IN ABOVE MESSAGE, I1=         0 */
/* ERROR NUMBER =         2 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE. */
/*           IN ABOVE MESSAGE, I1=         0 */
/* ERROR NUMBER =         3 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4. */
/*           IN ABOVE MESSAGE, I1=         1 */
/*           IN ABOVE MESSAGE, I2=         0 */
/* ERROR NUMBER =         4 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2). */
/*           IN ABOVE MESSAGE, I1=         1 */
/*           IN ABOVE MESSAGE, R1=    0. */
/*           IN ABOVE MESSAGE, R2=    ABOVE MESSAGE, I1=         0 */
/* ERROR NUMBER =         6 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). ISCALE OPTION=(I1) MUST BE 1-3. */
/*           IN ABOVE MESSAGE, I1=         0 */
/* ERROR NUMBER =         7 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED  COLUMN SCALING */
/* MUST BE POSITIVE. */
/*           IN ABOVE MESSAGE, I1=         0 */
/* ERROR NUMBER =         8 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE. */
/* COMPONENT (I1) NOW = (R1). */
/*           IN ABOVE MESSAGE, I1=        ND. .LE. MDW=(I2). */
/*           IN ABOVE MESSAGE, I1=         1 */
/*           IN ABOVE MESSAGE, I2=         0 */
/* ERROR NUMBER =        10 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS().THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE.THE NUMBER OF ROWS= */
/* (I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         1 */
/* ERROR NUMBER =        11 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+1=(I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         2 */
/* ERROR NUMBER =        12 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS().THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1) MUST BE */
/* .GE. NCOLS=(I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         1 */
/* ERROR NUMBER =        13 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). THE DIMENSION OF X()=(I1) MUST BE .GE. THE REQD. LENGTH=(I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         2 */
/* ERROR NUMBER =        14 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS(). THE DIMENSION OF RW()=(I1) MUST BE .GE. 5*NCOLS=(I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         3 */
/* ERROR NUMBER =        15 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS=(I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         2 */
/* ERROR NUMBER =        16 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* WARNING IN... */
/* DBOLS() THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD. LEN.=(I2). */
/*           IN ABOVE MESSAGE, I1=         0 */
/*           IN ABOVE MESSAGE, I2=         1 */
/* ERROR NUMBER =        17 */
/* (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

/* ***REFERENCES  R. J. Hanson, Linear least squares with bounds and */
/*                 linear constraints, Report SAND82-1517, Sandia */
/*                 Laboratories, August 1982. */
/* ***ROUTINES CALLED  DBOLSM, DCOPY, DNRM2, DROT, DROTG, IDAMAX, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   821220  DATE WRITTEN */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBOLS */

/*     SOLVE LINEAR LEAST SQUARES SYSTEM WITH BOUNDS ON */
/*     SELECTED VARIABLES. */
/*     REVISED 850329-1400 */
/*     REVISED YYMMDD-HHMM */
/*     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN */
/*     EDITING AT THE CARD 'C++'. */
/*     CHANGE THIS SUBPROGRAM NAME TO DBOLS AND THE STRINGS */
/*     /SCOPY/ TO /DCOPY/, /SBOL/ TO /DBOL/, */
/*     /SNRM2/ TO /DNRM2/, /ISAMAX/ TO /IDAMAX/, */
/*     /SROTG/ TO /DROTG/, /SROT/ TO /DROT/, /E0/ TO /D0/, */
/*     /REAL            / TO /DOUBLE PRECISION/. */
/* ++ */

/*     THIS VARIABLE SHOULD REMAIN TYPE REAL. */
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
    --iw;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DBOLS */
    nerr = 0;
    *mode = 0;
    if (igo == 0) {
/*     DO(CHECK VALIDITY OF INPUT DATA) */
/*     PROCEDURE(CHECK VALIDITY OF INPUT DATA) */

/*     SEE THAT MDW IS .GT.0. GROSS CHECK ONLY. */
	if (*mdw <= 0) {
	    s_wsfi(&io___4);
	    do_fio(&c__1, (char *)&(*mdw), (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 6, a__1[0] = "MDW = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    i__1[2] = 18, a__1[2] = " MUST BE POSITIVE.";
	    s_cat(ch__1, a__1, i__1, &c__3, (ftnlen)32);
	    xermsg_("SLATEC", "DBOLS", ch__1, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)32);
/*     DO(RETURN TO USER PROGRAM UNIT) */
	    goto L190;
	}

/*     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE. */
	if (*ncols <= 0) {
	    s_wsfi(&io___5);
	    do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 8, a__1[0] = "NCOLS = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    i__1[2] = 39, a__1[2] = " THE NO. OF VARIABLES MUST BE POSITIVE.";
	    s_cat(ch__2, a__1, i__1, &c__3, (ftnlen)55);
	    xermsg_("SLATEC", "DBOLS", ch__2, &c__3, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)55);
/*     DO(RETURN TO USER PROGRAM UNIT) */
	    goto L190;
	}

/*     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED. */
	i__2 = *ncols;
	for (j = 1; j <= i__2; ++j) {
	    if (ind[j] < 1 || ind[j] > 4) {
		s_wsfi(&io___7);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___9);
		do_fio(&c__1, (char *)&ind[j], (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = "IND(";
		i__3[1] = 8, a__2[1] = xern1;
		i__3[2] = 4, a__2[2] = ") = ";
		i__3[3] = 8, a__2[3] = xern2;
		i__3[4] = 13, a__2[4] = " MUST BE 1-4.";
		s_cat(ch__3, a__2, i__3, &c__5, (ftnlen)37);
		xermsg_("SLATEC", "DBOLS", ch__3, &c__4, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)37);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
/* L10: */
	}

/*     SEE THAT BOUNDS ARE CONSISTENT. */
	i__2 = *ncols;
	for (j = 1; j <= i__2; ++j) {
	    if (ind[j] == 3) {
		if (bl[j] > bu[j]) {
		    s_wsfi(&io___10);
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    e_wsfi();
		    s_wsfi(&io___12);
		    do_fio(&c__1, (char *)&bl[j], (ftnlen)sizeof(doublereal));
		    e_wsfi();
		    s_wsfi(&io___14);
		    do_fio(&c__1, (char *)&bu[j], (ftnlen)sizeof(doublereal));
		    e_wsfi();
/* Writing concatenation */
		    i__4[0] = 9, a__3[0] = "BOUND BL(";
		    i__4[1] = 8, a__3[1] = xern1;
		    i__4[2] = 4, a__3[2] = ") = ";
		    i__4[3] = 16, a__3[3] = xern3;
		    i__4[4] = 12, a__3[4] = " IS .GT. BU(";
		    i__4[5] = 8, a__3[5] = xern1;
		    i__4[6] = 4, a__3[6] = ") = ";
		    i__4[7] = 16, a__3[7] = xern4;
		    s_cat(ch__4, a__3, i__4, &c__8, (ftnlen)77);
		    xermsg_("SLATEC", "DBOLS", ch__4, &c__5, &c__1, (ftnlen)6,
			     (ftnlen)5, (ftnlen)77);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		    goto L190;
		}
	    }
/* L20: */
	}
/*     END PROCEDURE */
/*     DO(PROCESS OPTION ARRAY) */
/*     PROCEDURE(PROCESS OPTION ARRAY) */
	zero = 0.;
	one = 1.;
	checkl = FALSE_;
	lenx = *ncols;
	iscale = 1;
	igo = 2;
	lopt = 0;
	lp = 0;
	lds = 0;
L30:
	lp += lds;
	ip = iopt[lp + 1];
	jp = abs(ip);

/*     TEST FOR NO MORE OPTIONS. */
	if (ip == 99) {
	    if (lopt == 0) {
		lopt = lp + 1;
	    }
	    goto L50;
	} else if (jp == 99) {
	    lds = 1;
	    goto L30;
	} else if (jp == 1) {
	    if (ip > 0) {

/*     SET UP DIRECTION FLAG, ROW STACKING POINTER */
/*     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS. */
		locacc = lp + 2;

/*                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION. */
/*     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2. */
/*                  IOPT(LOCACC+1)=ROW STACKING POINTER. */
/*                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS. */
/*     USER ACTION WITH THIS OPTION.. */
/*      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*). */
/*      MUST ALSO START PROCESS WITH IOPT(LOCACC)=1.) */
/*      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST */
/*       ROW OF W(*,*).  SET IOPT(LOCACC+2)=NO. OF ROWS IN BLOCK.) */
/*              LOOP */
/*              CALL DBOLS() */

/*                  IF(IOPT(LOCACC) .EQ. 1) THEN */
/*                      STACK EQUAS., STARTING AT ROW IOPT(LOCACC+1), */
/*                       INTO W(*,*). */
/*                       SET IOPT(LOCACC+2)=NO. OF EQUAS. */
/*                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2. */
/*                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN */
/*                      (PROCESS IS OVER. EXIT LOOP.) */
/*                  ELSE */
/*                      (ERROR CONDITION. SHOULD NOT HAPPEN.) */
/*                  END IF */
/*              END LOOP */
/*              SET IOPT(LOCACC-1)=-OPTION NUMBER FOR SEQ. ACCUMULATION. */
/*              CALL DBOLS( ) */
		iopt[locacc + 1] = 1;
		igo = 1;
	    }
	    lds = 4;
	    goto L30;
	} else if (jp == 2) {
	    if (ip > 0) {

/*     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS. */
		locdim = lp + 2;

/*     LMDW.GE.MROWS */
/*     LNDW.GE.NCOLS+1 */
/*     LLB .GE.NCOLS */
/*     LLX .GE.NCOLS+EXTRA REQD. IN OPTIONS. */
/*     LLRW.GE.5*NCOLS */
/*     LLIW.GE.2*NCOLS */
/*     LIOP.GE. AMOUNT REQD. FOR IOPTION ARRAY. */
		lmdw = iopt[locdim];
		lndw = iopt[locdim + 1];
		llb = iopt[locdim + 2];
		llx = iopt[locdim + 3];
		llrw = iopt[locdim + 4];
		lliw = iopt[locdim + 5];
		liopt = iopt[locdim + 6];
		checkl = TRUE_;
	    }
	    lds = 8;
	    goto L30;

/*     OPTION TO MODIFY THE COLUMN SCALING. */
	} else if (jp == 3) {
	    if (ip > 0) {
		iscale = iopt[lp + 2];

/*     SEE THAT ISCALE IS 1 THRU 3. */
		if (iscale < 1 || iscale > 3) {
		    s_wsfi(&io___34);
		    do_fio(&c__1, (char *)&iscale, (ftnlen)sizeof(integer));
		    e_wsfi();
/* Writing concatenation */
		    i__1[0] = 16, a__1[0] = "ISCALE OPTION = ";
		    i__1[1] = 8, a__1[1] = xern1;
		    i__1[2] = 12, a__1[2] = " MUST BE 1-3";
		    s_cat(ch__5, a__1, i__1, &c__3, (ftnlen)36);
		    xermsg_("SLATEC", "DBOLS", ch__5, &c__7, &c__1, (ftnlen)6,
			     (ftnlen)5, (ftnlen)36);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		    goto L190;
		}
	    }
	    lds = 2;
/*     CYCLE FOREVER */
	    goto L30;

/*     IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE */
/*     SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)). */
	} else if (jp == 4) {
	    if (ip > 0) {
		iscale = 4;
		if (iopt[lp + 2] <= 0) {
		    s_wsfi(&io___35);
		    do_fio(&c__1, (char *)&iopt[lp + 2], (ftnlen)sizeof(
			    integer));
		    e_wsfi();
/* Writing concatenation */
		    i__5[0] = 22, a__4[0] = "OFFSET PAST X(NCOLS) (";
		    i__5[1] = 8, a__4[1] = xern1;
		    i__5[2] = 40, a__4[2] = ") FOR USER-PROVIDED COLUMN SCAL"
			    "ING MUST ";
		    i__5[3] = 12, a__4[3] = "BE POSITIVE.";
		    s_cat(ch__6, a__4, i__5, &c__4, (ftnlen)82);
		    xermsg_("SLATEC", "DBOLS", ch__6, &c__8, &c__1, (ftnlen)6,
			     (ftnlen)5, (ftnlen)82);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		    goto L190;
		}
		dcopy_(ncols, &x[*ncols + iopt[lp + 2]], &c__1, &rw[1], &c__1)
			;
		lenx += *ncols;
		i__2 = *ncols;
		for (j = 1; j <= i__2; ++j) {
		    if (rw[j] <= zero) {
			s_wsfi(&io___36);
			do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
			e_wsfi();
			s_wsfi(&io___37);
			do_fio(&c__1, (char *)&rw[j], (ftnlen)sizeof(
				doublereal));
			e_wsfi();
/* Writing concatenation */
			i__5[0] = 63, a__4[0] = "EACH PROVIDED COLUMN SCALE "
				"FACTOR MUST BE POSITIVE.$$COMPONENT ";
			i__5[1] = 8, a__4[1] = xern1;
			i__5[2] = 7, a__4[2] = " NOW = ";
			i__5[3] = 16, a__4[3] = xern3;
			s_cat(ch__7, a__4, i__5, &c__4, (ftnlen)94);
			xermsg_("SLATEC", "DBOLS", ch__7, &c__9, &c__1, (
				ftnlen)6, (ftnlen)5, (ftnlen)94);
			goto L190;
		    }
/* L40: */
		}
	    }
	    lds = 2;
/*     CYCLE FOREVER */
	    goto L30;

/*     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM(). */
	} else if (jp == 5) {
	    if (ip > 0) {
		lopt = iopt[lp + 2];
	    }
	    lds = 2;
/*     CYCLE FOREVER */
	    goto L30;

/*     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS AN */
/*     INCREMENT TO SKIP. */
	} else if (jp == 6) {
	    if (ip > 0) {
		lp = iopt[lp + 2] - 1;
		lds = 0;
	    } else {
		lds = 2;
	    }
/*     CYCLE FOREVER */
	    goto L30;

/*     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION. */
	} else {
	    s_wsfi(&io___38);
	    do_fio(&c__1, (char *)&jp, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 20, a__1[0] = "THE OPTION NUMBER = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    i__1[2] = 16, a__1[2] = " IS NOT DEFINED.";
	    s_cat(ch__8, a__1, i__1, &c__3, (ftnlen)44);
	    xermsg_("SLATEC", "DBOLS", ch__8, &c__6, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)44);
/*     DO(RETURN TO USER PROGRAM UNIT) */
	    goto L190;
	}
L50:
/*     END PROCEDURE */
	if (checkl) {
/*     DO(CHECK LENGTHS OF ARRAYS) */
/*     PROCEDURE(CHECK LENGTHS OF ARRAYS) */

/*     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE */
/*     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE. */
	    if (lmdw < *mrows) {
		s_wsfi(&io___39);
		do_fio(&c__1, (char *)&lmdw, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___40);
		do_fio(&c__1, (char *)&(*mrows), (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 28, a__4[0] = "THE ROW DIMENSION OF W(,) = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 35, a__4[2] = " MUST BE .GE. THE NUMBER OF ROWS = ";
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__9, a__4, i__5, &c__4, (ftnlen)79);
		xermsg_("SLATEC", "DBOLS", ch__9, &c__11, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)79);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
	    if (lndw < *ncols + 1) {
		s_wsfi(&io___41);
		do_fio(&c__1, (char *)&lndw, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___42);
		i__2 = *ncols + 1;
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 31, a__4[0] = "THE COLUMN DIMENSION OF W(,) = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 24, a__4[2] = " MUST BE .GE. NCOLS+1 = ";
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__10, a__4, i__5, &c__4, (ftnlen)71);
		xermsg_("SLATEC", "DBOLS", ch__10, &c__12, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)71);
		goto L190;
	    }
	    if (llb < *ncols) {
		s_wsfi(&io___43);
		do_fio(&c__1, (char *)&llb, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___44);
		do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 53, a__4[0] = "THE DIMENSIONS OF THE ARRAYS BL(), "
			"BU(), AND IND() = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 22, a__4[2] = " MUST BE .GE. NCOLS = ";
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__11, a__4, i__5, &c__4, (ftnlen)91);
		xermsg_("SLATEC", "DBOLS", ch__11, &c__13, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)91);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
	    if (llx < lenx) {
		s_wsfi(&io___45);
		do_fio(&c__1, (char *)&llx, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___46);
		do_fio(&c__1, (char *)&lenx, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 23, a__4[0] = "THE DIMENSION OF X() = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 36, a__4[2] = " MUST BE .GE. THE REQUIRED LENGTH = "
			;
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__12, a__4, i__5, &c__4, (ftnlen)75);
		xermsg_("SLATEC", "DBOLS", ch__12, &c__14, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)75);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
	    if (llrw < *ncols * 5) {
		s_wsfi(&io___47);
		do_fio(&c__1, (char *)&llrw, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___48);
		i__2 = *ncols * 5;
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 24, a__4[0] = "THE DIMENSION OF RW() = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 24, a__4[2] = " MUST BE .GE. 5*NCOLS = ";
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__13, a__4, i__5, &c__4, (ftnlen)64);
		xermsg_("SLATEC", "DBOLS", ch__13, &c__15, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)64);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
	    if (lliw < *ncols << 1) {
		s_wsfi(&io___49);
		do_fio(&c__1, (char *)&lliw, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___50);
		i__2 = *ncols << 1;
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 24, a__4[0] = "THE DIMENSION OF IW() = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 24, a__4[2] = " MUST BE .GE. 2*NCOLS = ";
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__13, a__4, i__5, &c__4, (ftnlen)64);
		xermsg_("SLATEC", "DBOLS", ch__13, &c__16, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)64);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
	    if (liopt < lp + 1) {
		s_wsfi(&io___51);
		do_fio(&c__1, (char *)&liopt, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___52);
		i__2 = lp + 1;
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 26, a__4[0] = "THE DIMENSION OF IOPT() = ";
		i__5[1] = 8, a__4[1] = xern1;
		i__5[2] = 33, a__4[2] = " MUST BE .GE. THE REQUIRED LEN = ";
		i__5[3] = 8, a__4[3] = xern2;
		s_cat(ch__12, a__4, i__5, &c__4, (ftnlen)75);
		xermsg_("SLATEC", "DBOLS", ch__12, &c__17, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)75);
/*     DO(RETURN TO USER PROGRAM UNIT) */
		goto L190;
	    }
/*     END PROCEDURE */
	}
    }
    switch (igo) {
	case 1:  goto L60;
	case 2:  goto L90;
    }
    goto L180;

/*     GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES */
/*     EQUATIONS AND DIRECTIONS TO QUIT PROCESSING. */
/*     CASE 1 */
L60:
/*     DO(ACCUMULATE LEAST SQUARES EQUATIONS) */
/*     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS) */
    *mrows = iopt[locacc + 1] - 1;
    inrows = iopt[locacc + 2];
    mnew = *mrows + inrows;
    if (mnew < 0 || mnew > *mdw) {
	s_wsfi(&io___55);
	do_fio(&c__1, (char *)&mnew, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___56);
	do_fio(&c__1, (char *)&(*mdw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__5[0] = 14, a__4[0] = "NO. OF ROWS = ";
	i__5[1] = 8, a__4[1] = xern1;
	i__5[2] = 33, a__4[2] = " MUST BE .GE. 0 .AND. .LE. MDW = ";
	i__5[3] = 8, a__4[3] = xern2;
	s_cat(ch__14, a__4, i__5, &c__4, (ftnlen)63);
	xermsg_("SLATEC", "DBOLS", ch__14, &c__10, &c__1, (ftnlen)6, (ftnlen)
		5, (ftnlen)63);
/*     DO(RETURN TO USER PROGRAM UNIT) */
	goto L190;
    }
/* Computing MIN */
    i__6 = *ncols + 1;
    i__2 = min(i__6,mnew);
    for (j = 1; j <= i__2; ++j) {
	i__6 = max(*mrows,j) + 1;
	for (i__ = mnew; i__ >= i__6; --i__) {
	    i__7 = i__ - j;
	    ibig = idamax_(&i__7, &w[j + j * w_dim1], &c__1) + j - 1;

/*     PIVOT FOR INCREASED STABILITY. */
	    drotg_(&w[ibig + j * w_dim1], &w[i__ + j * w_dim1], &sc, &ss);
	    i__7 = *ncols + 1 - j;
	    drot_(&i__7, &w[ibig + (j + 1) * w_dim1], mdw, &w[i__ + (j + 1) * 
		    w_dim1], mdw, &sc, &ss);
	    w[i__ + j * w_dim1] = zero;
/* L70: */
	}
/* L80: */
    }
/* Computing MIN */
    i__2 = *ncols + 1;
    *mrows = min(i__2,mnew);
    iopt[locacc + 1] = *mrows + 1;
    igo = iopt[locacc];
/*     END PROCEDURE */
    if (igo == 2) {
	igo = 0;
    }
    goto L180;
/*     CASE 2 */
L90:
/*     DO(INITIALIZE VARIABLES AND DATA VALUES) */
/*     PROCEDURE(INITIALIZE VARIABLES AND DATA VALUES) */
    i__2 = *ncols;
    for (j = 1; j <= i__2; ++j) {
	switch (iscale) {
	    case 1:  goto L100;
	    case 2:  goto L110;
	    case 3:  goto L120;
	    case 4:  goto L130;
	}
	goto L140;
L100:
/*     CASE 1 */

/*     THIS IS THE NOMINAL SCALING. EACH NONZERO */
/*     COL. HAS MAX. NORM EQUAL TO ONE. */
	ibig = idamax_(mrows, &w[j * w_dim1 + 1], &c__1);
	rw[j] = (d__1 = w[ibig + j * w_dim1], abs(d__1));
	if (rw[j] == zero) {
	    rw[j] = one;
	} else {
	    rw[j] = one / rw[j];
	}
	goto L140;
L110:
/*     CASE 2 */

/*     THIS CHOICE OF SCALING MAKES EACH NONZERO COLUMN */
/*     HAVE EUCLIDEAN LENGTH EQUAL TO ONE. */
	rw[j] = dnrm2_(mrows, &w[j * w_dim1 + 1], &c__1);
	if (rw[j] == zero) {
	    rw[j] = one;
	} else {
	    rw[j] = one / rw[j];
	}
	goto L140;
L120:
/*     CASE 3 */

/*     THIS CASE EFFECTIVELY SUPPRESSES SCALING BY SETTING */
/*     THE SCALING MATRIX TO THE IDENTITY MATRIX. */
	rw[1] = one;
	dcopy_(ncols, &rw[1], &c__0, &rw[1], &c__1);
	goto L160;
L130:
/*     CASE 4 */
	goto L160;
L140:
/* L150: */
	;
    }
L160:
/*     END PROCEDURE */
/*     DO(SOLVE BOUNDED LEAST SQUARES PROBLEM) */
/*     PROCEDURE(SOLVE BOUNDED LEAST SQUARES PROBLEM) */

/*     INITIALIZE IBASIS(*), J=1,NCOLS, AND IBB(*), J=1,NCOLS, */
/*     TO =J,AND =1, FOR USE IN DBOLSM( ). */
    i__2 = *ncols;
    for (j = 1; j <= i__2; ++j) {
	iw[j] = j;
	iw[j + *ncols] = 1;
	rw[*ncols * 3 + j] = bl[j];
	rw[(*ncols << 2) + j] = bu[j];
/* L170: */
    }
    dbolsm_(&w[w_offset], mdw, mrows, ncols, &rw[*ncols * 3 + 1], &rw[(*ncols 
	    << 2) + 1], &ind[1], &iopt[lopt], &x[1], rnorm, mode, &rw[*ncols 
	    + 1], &rw[(*ncols << 1) + 1], &rw[1], &iw[1], &iw[*ncols + 1]);
/*     END PROCEDURE */
    igo = 0;
L180:
    return 0;
/*     PROCEDURE(RETURN TO USER PROGRAM UNIT) */
L190:
    if (*mode >= 0) {
	*mode = -nerr;
    }
    igo = 0;
    return 0;
/*     END PROCEDURE */
} /* dbols_ */


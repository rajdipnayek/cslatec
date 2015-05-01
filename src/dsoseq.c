/* dsoseq.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DSOSEQ */
/* Subroutine */ int dsoseq_(D_fp fnc, integer *n, doublereal *s, doublereal *
	rtolx, doublereal *atolx, doublereal *tolf, integer *iflag, integer *
	mxit, integer *ncjs, integer *nsrrc, integer *nsri, integer *iprint, 
	doublereal *fmax, doublereal *c__, integer *nc, doublereal *b, 
	doublereal *p, doublereal *temp, doublereal *x, doublereal *y, 
	doublereal *fac, integer *is)
{
    /* Format strings */
    static char fmt_210[] = "(\0020RESIDUAL NORM =\002,d9.2,/1x,\002SOLUTION"
	    " ITERATE (\002,i3,\002)\002,/(1x,5d26.14))";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal f, h__;
    static integer j, k, l, m, ic, kd, jk, kj, kk;
    static doublereal fp;
    static integer kn, mm;
    static doublereal re;
    static integer it, js, ls;
    static doublereal hx, yj, fn1, fn2;
    static integer km1, np1;
    static doublereal yn1, yn2, yn3;
    static integer icr, isj, mit;
    static doublereal csv;
    static integer isv, ksv;
    static doublereal uro, yns, fdif, fact, fmin;
    static integer item;
    static doublereal pmax;
    static integer loun;
    static doublereal fmxs, test, zero;
    static integer itry;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    static doublereal xnorm, ynorm, sruro;
    extern /* Subroutine */ int dsossl_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 0, 0, fmt_210, 0 };


/* ***BEGIN PROLOGUE  DSOSEQ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSOS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SOSEQS-S, DSOSEQ-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     DSOSEQ solves a system of N simultaneous nonlinear equations. */
/*     See the comments in the interfacing routine DSOS for a more */
/*     detailed description of some of the items in the calling list. */

/* ********************************************************************** */
/*   -Input- */

/*     FNC- Function subprogram which evaluates the equations */
/*     N  -number of equations */
/*     S  -Solution vector of initial guesses */
/*     RTOLX-Relative error tolerance on solution components */
/*     ATOLX-Absolute error tolerance on solution components */
/*     TOLF-Residual error tolerance */
/*     MXIT-Maximum number of allowable iterations. */
/*     NCJS-Maximum number of consecutive iterative steps to perform */
/*          using the same triangular Jacobian matrix approximation. */
/*     NSRRC-Number of consecutive iterative steps for which the */
/*          limiting precision accuracy test must be satisfied */
/*          before the routine exits with IFLAG=4. */
/*     NSRI-Number of consecutive iterative steps for which the */
/*          diverging condition test must be satisfied before */
/*          the routine exits with IFLAG=7. */
/*     IPRINT-Internal printing parameter. You must set IPRINT=-1 if you */
/*          want the intermediate solution iterates and a residual norm */
/*          to be printed. */
/*     C   -Internal work array, dimensioned at least N*(N+1)/2. */
/*     NC  -Dimension of C array. NC  .GE.  N*(N+1)/2. */
/*     B   -Internal work array, dimensioned N. */
/*     P   -Internal work array, dimensioned N. */
/*     TEMP-Internal work array, dimensioned N. */
/*     X   -Internal work array, dimensioned N. */
/*     Y   -Internal work array, dimensioned N. */
/*     FAC -Internal work array, dimensioned N. */
/*     IS  -Internal work array, dimensioned N. */

/*   -Output- */
/*     S    -Solution vector */
/*     IFLAG-Status indicator flag */
/*     MXIT-The actual number of iterations performed */
/*     FMAX-Residual norm */
/*     C   -Upper unit triangular matrix which approximates the */
/*          forward triangularization of the full Jacobian matrix. */
/*          Stored in a vector with dimension at least N*(N+1)/2. */
/*     B   -Contains the residuals (function values) divided */
/*          by the corresponding components of the P vector */
/*     P   -Array used to store the partial derivatives. After */
/*          each iteration P(K) contains the maximal derivative */
/*          occurring in the K-th reduced equation. */
/*     TEMP-Array used to store the previous solution iterate. */
/*     X   -Solution vector. Contains the values achieved on the */
/*          last iteration loop upon exit from DSOS. */
/*     Y   -Array containing the solution increments. */
/*     FAC -Array containing factors used in computing numerical */
/*          derivatives. */
/*     IS  -Records the pivotal information (column interchanges) */

/* ********************************************************************** */
/* *** Three machine dependent parameters appear in this subroutine. */

/* *** The smallest positive magnitude, zero, is defined by the function */
/* *** routine D1MACH(1). */

/* *** URO, the computer unit roundoff value, is defined by D1MACH(3) for */
/* *** machines that round or D1MACH(4) for machines that truncate. */
/* *** URO is the smallest positive number such that 1.+URO  .GT.  1. */

/* *** The output tape unit number, LOUN, is defined by the function */
/* *** I1MACH(2). */
/* ********************************************************************** */

/* ***SEE ALSO  DSOS */
/* ***ROUTINES CALLED  D1MACH, DSOSSL, I1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DSOSEQ */



/*     BEGIN BLOCK PERMITTING ...EXITS TO 430 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 410 */
/*           BEGIN BLOCK PERMITTING ...EXITS TO 390 */
/* ***FIRST EXECUTABLE STATEMENT  DSOSEQ */
    /* Parameter adjustments */
    --is;
    --fac;
    --y;
    --x;
    --temp;
    --p;
    --b;
    --c__;
    --s;

    /* Function Body */
    uro = d1mach_(&c__4);
    loun = i1mach_(&c__2);
    zero = d1mach_(&c__1);
    re = max(*rtolx,uro);
    sruro = sqrt(uro);

    *iflag = 0;
    np1 = *n + 1;
    icr = 0;
    ic = 0;
    itry = *ncjs;
    yn1 = 0.;
    yn2 = 0.;
    yn3 = 0.;
    yns = 0.;
    mit = 0;
    fn1 = 0.;
    fn2 = 0.;
    fmxs = 0.;

/*              INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND */
/*              SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	is[k] = k;
	x[k] = s[k];
	temp[k] = x[k];
/* L10: */
    }


/*              ********************************************************* */
/*              **** BEGIN PRINCIPAL ITERATION LOOP  **** */
/*              ********************************************************* */

    i__1 = *mxit;
    for (m = 1; m <= i__1; ++m) {
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 350 */
/*                    BEGIN BLOCK PERMITTING ...EXITS TO 240 */

	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    fac[k] = sruro;
/* L20: */
	}

L30:
/*                          BEGIN BLOCK PERMITTING ...EXITS TO 180 */
	kn = 1;
	*fmax = 0.;


/*                             ******** BEGIN SUBITERATION LOOP DEFINING */
/*                             THE LINEARIZATION OF EACH ******** */
/*                             EQUATION WHICH RESULTS IN THE CONSTRUCTION */
/*                             OF AN UPPER ******** TRIANGULAR MATRIX */
/*                             APPROXIMATING THE FORWARD ******** */
/*                             TRIANGULARIZATION OF THE FULL JACOBIAN */
/*                             MATRIX */

	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
/*                                BEGIN BLOCK PERMITTING ...EXITS TO 160 */
	    km1 = k - 1;

/*                                   BACK-SOLVE A TRIANGULAR LINEAR */
/*                                   SYSTEM OBTAINING IMPROVED SOLUTION */
/*                                   VALUES FOR K-1 OF THE VARIABLES FROM */
/*                                   THE FIRST K-1 EQUATIONS. THESE */
/*                                   VARIABLES ARE THEN ELIMINATED FROM */
/*                                   THE K-TH EQUATION. */

	    if (km1 == 0) {
		goto L50;
	    }
	    dsossl_(&k, n, &km1, &y[1], &c__[1], &b[1], &kn);
	    i__3 = km1;
	    for (j = 1; j <= i__3; ++j) {
		js = is[j];
		x[js] = temp[js] + y[j];
/* L40: */
	    }
L50:


/*                                   EVALUATE THE K-TH EQUATION AND THE */
/*                                   INTERMEDIATE COMPUTATION FOR THE MAX */
/*                                   NORM OF THE RESIDUAL VECTOR. */

	    f = (*fnc)(&x[1], &k);
/* Computing MAX */
	    d__1 = *fmax, d__2 = abs(f);
	    *fmax = max(d__1,d__2);

/*                                   IF WE WISH TO PERFORM SEVERAL */
/*                                   ITERATIONS USING A FIXED */
/*                                   FACTORIZATION OF AN APPROXIMATE */
/*                                   JACOBIAN,WE NEED ONLY UPDATE THE */
/*                                   CONSTANT VECTOR. */

/*                                ...EXIT */
	    if (itry < *ncjs) {
		goto L160;
	    }


	    it = 0;

/*                                   COMPUTE PARTIAL DERIVATIVES THAT ARE */
/*                                   REQUIRED IN THE LINEARIZATION OF THE */
/*                                   K-TH REDUCED EQUATION */

	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		item = is[j];
		hx = x[item];
		h__ = fac[item] * hx;
		if (abs(h__) <= zero) {
		    h__ = fac[item];
		}
		x[item] = hx + h__;
		if (km1 == 0) {
		    goto L70;
		}
		y[j] = h__;
		dsossl_(&k, n, &j, &y[1], &c__[1], &b[1], &kn);
		i__4 = km1;
		for (l = 1; l <= i__4; ++l) {
		    ls = is[l];
		    x[ls] = temp[ls] + y[l];
/* L60: */
		}
L70:
		fp = (*fnc)(&x[1], &k);
		x[item] = hx;
		fdif = fp - f;
		if (abs(fdif) > uro * abs(f)) {
		    goto L80;
		}
		fdif = 0.;
		++it;
L80:
		p[j] = fdif / h__;
/* L90: */
	    }

	    if (it <= *n - k) {
		goto L110;
	    }

/*                                      ALL COMPUTED PARTIAL DERIVATIVES */
/*                                      OF THE K-TH EQUATION ARE */
/*                                      EFFECTIVELY ZERO.TRY LARGER */
/*                                      PERTURBATIONS OF THE INDEPENDENT */
/*                                      VARIABLES. */

	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		isj = is[j];
		fact = fac[isj] * 100.;
/*           ..............................EXIT */
		if (fact > 1e10) {
		    goto L390;
		}
		fac[isj] = fact;
/* L100: */
	    }
/*                          ............EXIT */
	    goto L180;
L110:

/*                                ...EXIT */
	    if (k == *n) {
		goto L160;
	    }

/*                                   ACHIEVE A PIVOTING EFFECT BY */
/*                                   CHOOSING THE MAXIMAL DERIVATIVE */
/*                                   ELEMENT */

	    pmax = 0.;
	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		test = (d__1 = p[j], abs(d__1));
		if (test <= pmax) {
		    goto L120;
		}
		pmax = test;
		isv = j;
L120:
/* L130: */
		;
	    }
/*           ........................EXIT */
	    if (pmax == 0.) {
		goto L390;
	    }

/*                                   SET UP THE COEFFICIENTS FOR THE K-TH */
/*                                   ROW OF THE TRIANGULAR LINEAR SYSTEM */
/*                                   AND SAVE THE PARTIAL DERIVATIVE OF */
/*                                   LARGEST MAGNITUDE */

	    pmax = p[isv];
	    kk = kn;
	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		if (j != isv) {
		    c__[kk] = -p[j] / pmax;
		}
		++kk;
/* L140: */
	    }
	    p[k] = pmax;


/*                                ...EXIT */
	    if (isv == k) {
		goto L160;
	    }

/*                                   INTERCHANGE THE TWO COLUMNS OF C */
/*                                   DETERMINED BY THE PIVOTAL STRATEGY */

	    ksv = is[k];
	    is[k] = is[isv];
	    is[isv] = ksv;

	    kd = isv - k;
	    kj = k;
	    i__3 = k;
	    for (j = 1; j <= i__3; ++j) {
		csv = c__[kj];
		jk = kj + kd;
		c__[kj] = c__[jk];
		c__[jk] = csv;
		kj = kj + *n - j;
/* L150: */
	    }
L160:

	    kn = kn + np1 - k;

/*                                STORE THE COMPONENTS FOR THE CONSTANT */
/*                                VECTOR */

	    b[k] = -f / p[k];

/* L170: */
	}
/*                       ......EXIT */
	goto L190;
L180:
	goto L30;
L190:

/*                       ******** */
/*                       ******** END OF LOOP CREATING THE TRIANGULAR */
/*                       LINEARIZATION MATRIX */
/*                       ******** */


/*                        SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW */
/*                        SOLUTION APPROXIMATION AND OBTAIN THE SOLUTION */
/*                        INCREMENT NORM. */

	--kn;
	y[*n] = b[*n];
	if (*n > 1) {
	    dsossl_(n, n, n, &y[1], &c__[1], &b[1], &kn);
	}
	xnorm = 0.;
	ynorm = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    yj = y[j];
/* Computing MAX */
	    d__1 = ynorm, d__2 = abs(yj);
	    ynorm = max(d__1,d__2);
	    js = is[j];
	    x[js] = temp[js] + yj;
/* Computing MAX */
	    d__2 = xnorm, d__3 = (d__1 = x[js], abs(d__1));
	    xnorm = max(d__2,d__3);
/* L200: */
	}


/*                       PRINT INTERMEDIATE SOLUTION ITERATES AND */
/*                       RESIDUAL NORM IF DESIRED */

	if (*iprint != -1) {
	    goto L220;
	}
	mm = m - 1;
	io___48.ciunit = loun;
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&(*fmax), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&mm, (ftnlen)sizeof(integer));
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
L220:

/*                       TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE */
/*                       AND/OR ABSOLUTE ERROR COMPARISON ON SUCCESSIVE */
/*                       APPROXIMATIONS OF EACH SOLUTION VARIABLE) */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    js = is[j];
/*                    ......EXIT */
	    if ((d__2 = y[j], abs(d__2)) > re * (d__1 = x[js], abs(d__1)) + *
		    atolx) {
		goto L240;
	    }
/* L230: */
	}
	if (*fmax <= fmxs) {
	    *iflag = 1;
	}
L240:

/*                    TEST FOR CONVERGENCE TO A SOLUTION BASED ON */
/*                    RESIDUALS */

	if (*fmax <= *tolf) {
	    *iflag += 2;
	}
/*        ............EXIT */
	if (*iflag > 0) {
	    goto L410;
	}


	if (m > 1) {
	    goto L250;
	}
	fmin = *fmax;
	goto L330;
L250:
/*                       BEGIN BLOCK PERMITTING ...EXITS TO 320 */

/*                          SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM. */

	if (*fmax >= fmin) {
	    goto L270;
	}
	mit = m + 1;
	yn1 = ynorm;
	yn2 = yns;
	fn1 = fmxs;
	fmin = *fmax;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    s[j] = x[j];
/* L260: */
	}
	ic = 0;
L270:

/*                          TEST FOR LIMITING PRECISION CONVERGENCE. VERY */
/*                          SLOWLY CONVERGENT PROBLEMS MAY ALSO BE */
/*                          DETECTED. */

	if (ynorm > sruro * xnorm) {
	    goto L290;
	}
	if (*fmax < fmxs * .2 || *fmax > fmxs * 5.) {
	    goto L290;
	}
	if (ynorm < yns * .2 || ynorm > yns * 5.) {
	    goto L290;
	}
	++icr;
	if (icr >= *nsrrc) {
	    goto L280;
	}
	ic = 0;
/*                       .........EXIT */
	goto L320;
L280:
	*iflag = 4;
	*fmax = fmin;
/*     ........................EXIT */
	goto L430;
L290:
	icr = 0;

/*                          TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME. */

	if (ynorm > yns * 2. || *fmax > fmxs * 2.) {
	    goto L300;
	}
	ic = 0;
	goto L310;
L300:
	++ic;
/*                       ......EXIT */
	if (ic < *nsri) {
	    goto L320;
	}
	*iflag = 7;
/*        .....................EXIT */
	goto L410;
L310:
L320:
L330:

/*                    CHECK TO SEE IF NEXT ITERATION CAN USE THE OLD */
/*                    JACOBIAN FACTORIZATION */

	--itry;
	if (itry == 0) {
	    goto L340;
	}
	if (ynorm * 20. > xnorm) {
	    goto L340;
	}
	if (ynorm > yns * 2.) {
	    goto L340;
	}
/*                 ......EXIT */
	if (*fmax < fmxs * 2.) {
	    goto L350;
	}
L340:
	itry = *ncjs;
L350:

/*                 SAVE THE CURRENT SOLUTION APPROXIMATION AND THE */
/*                 RESIDUAL AND SOLUTION INCREMENT NORMS FOR USE IN THE */
/*                 NEXT ITERATION. */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    temp[j] = x[j];
/* L360: */
	}
	if (m != mit) {
	    goto L370;
	}
	fn2 = *fmax;
	yn3 = ynorm;
L370:
	fmxs = *fmax;
	yns = ynorm;


/* L380: */
    }

/*              ********************************************************* */
/*              **** END OF PRINCIPAL ITERATION LOOP **** */
/*              ********************************************************* */


/*               TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED. */
    m = *mxit;
    *iflag = 5;
    if (yn1 > yn2 * 10. || yn3 > yn1 * 10.) {
	*iflag = 6;
    }
    if (fn1 > fmin * 5. || fn2 > fmin * 5.) {
	*iflag = 6;
    }
    if (*fmax > fmin * 5.) {
	*iflag = 6;
    }
/*        ......EXIT */
    goto L410;
L390:


/*           A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR. */
    *iflag = 8;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s[j] = temp[j];
/* L400: */
    }
/*     ......EXIT */
    goto L430;
L410:


    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s[j] = x[j];
/* L420: */
    }
L430:


    *mxit = m;
    return 0;
} /* dsoseq_ */


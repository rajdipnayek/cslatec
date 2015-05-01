/* soseqs.f -- translated by f2c (version 12.02.01).
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

/* DECK SOSEQS */
/* Subroutine */ int soseqs_(E_fp fnc, integer *n, real *s, real *rtolx, real 
	*atolx, real *tolf, integer *iflag, integer *mxit, integer *ncjs, 
	integer *nsrrc, integer *nsri, integer *iprint, real *fmax, real *c__,
	 integer *nc, real *b, real *p, real *temp, real *x, real *y, real *
	fac, integer *is)
{
    /* Format strings */
    static char fmt_1234[] = "(\0020RESIDUAL NORM =\002,e9.2,/1x,\002SOLUTIO"
	    "N ITERATE\002,\002 (\002,i3,\002)\002,/(1x,5e26.14))";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3;

    /* Local variables */
    static real f, h__;
    static integer j, k, l, m, ic, kd, kj;
    static real fp, re;
    static integer kk, kn, jk, mm, js, it, ls;
    static real hx, yj, fn1, fn2;
    static integer km1, np1;
    static real yn1, yn2, yn3;
    static integer icr, isj, mit;
    static real csv;
    static integer isv, ksv;
    static real uro, yns, fdif, fact, fmin;
    static integer item;
    static real pmax;
    static integer loun;
    static real fmxs, zero, test;
    static integer itry;
    extern integer i1mach_(integer *);
    static real xnorm, ynorm, sruro;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int sossol_(integer *, integer *, integer *, real 
	    *, real *, real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 0, 0, fmt_1234, 0 };


/* ***BEGIN PROLOGUE  SOSEQS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SOS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SOSEQS-S, DSOSEQ-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     SOSEQS solves a system of N simultaneous nonlinear equations. */
/*     See the comments in the interfacing routine SOS for a more */
/*     detailed description of some of the items in the calling list. */

/* ******************************************************************** */

/*   -INPUT- */
/*     FNC -Function subprogram which evaluates the equations */
/*     N   -Number of equations */
/*     S   -Solution vector of initial guesses */
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
/*     IPRINT-Internal printing parameter.  You must set IPRINT=-1 if you */
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

/*   -OUTPUT- */
/*     S   -Solution vector */
/*     IFLAG-Status indicator flag */
/*     MXIT-The actual number of iterations performed */
/*     FMAX-Residual norm */
/*     C   -Upper unit triangular matrix which approximates the */
/*          forward triangularization of the full Jacobian matrix. */
/*          stored in a vector with dimension at least N*(N+1)/2. */
/*     B   -Contains the residuals (function values) divided */
/*          by the corresponding components of the P vector */
/*     P   -Array used to store the partial derivatives. After */
/*          each iteration P(K) contains the maximal derivative */
/*          occurring in the K-th reduced equation. */
/*     TEMP-Array used to store the previous solution iterate. */
/*     X   -Solution vector. Contains the values achieved on the */
/*          last iteration loop upon exit from SOS. */
/*     Y   -Array containing the solution increments. */
/*     FAC -Array containing factors used in computing numerical */
/*          derivatives. */
/*     IS  -Records the pivotal information (column interchanges) */

/* ********************************************************************** */
/* *** Three machine dependent parameters appear in this subroutine. */

/* *** The smallest positive magnitude, zero, is defined by the function */
/* *** routine R1MACH(1). */

/* *** URO, The computer unit roundoff value, is defined by R1MACH(3) for */
/* *** machines that round or R1MACH(4) for machines that truncate. */
/* *** URO is the smallest positive number such that 1.+URO  .GT.  1. */

/* *** The output tape unit number, LOUN, is defined by the function */
/* *** I1MACH(2). */
/* ********************************************************************** */

/* ***SEE ALSO  SOS */
/* ***ROUTINES CALLED  I1MACH, R1MACH, SOSSOL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SOSEQS */



/* ***FIRST EXECUTABLE STATEMENT  SOSEQS */
    /* Parameter adjustments */
    --s;
    --c__;
    --b;
    --p;
    --temp;
    --x;
    --y;
    --fac;
    --is;

    /* Function Body */
    uro = r1mach_(&c__4);
    loun = i1mach_(&c__2);
    zero = r1mach_(&c__1);
    re = dmax(*rtolx,uro);
    sruro = sqrt(uro);

    *iflag = 0;
    np1 = *n + 1;
    icr = 0;
    ic = 0;
    itry = *ncjs;
    yn1 = 0.f;
    yn2 = 0.f;
    yn3 = 0.f;
    yns = 0.f;
    mit = 0;
    fn1 = 0.f;
    fn2 = 0.f;
    fmxs = 0.f;

/*     INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND */
/*     SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	is[k] = k;
	x[k] = s[k];
	temp[k] = x[k];
/* L10: */
    }


/*    ***************************************** */
/*    **** BEGIN PRINCIPAL ITERATION LOOP  **** */
/*    ***************************************** */

    i__1 = *mxit;
    for (m = 1; m <= i__1; ++m) {

	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    fac[k] = sruro;
/* L20: */
	}

L30:
	kn = 1;
	*fmax = 0.f;


/*    ******** BEGIN SUBITERATION LOOP DEFINING THE LINEARIZATION OF EACH */
/*    ******** EQUATION WHICH RESULTS IN THE CONSTRUCTION OF AN UPPER */
/*    ******** TRIANGULAR MATRIX APPROXIMATING THE FORWARD */
/*    ******** TRIANGULARIZATION OF THE FULL JACOBIAN MATRIX */

	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    km1 = k - 1;

/*     BACK-SOLVE A TRIANGULAR LINEAR SYSTEM OBTAINING */
/*     IMPROVED SOLUTION VALUES FOR K-1 OF THE VARIABLES */
/*     FROM THE FIRST K-1 EQUATIONS. THESE VARIABLES ARE THEN */
/*     ELIMINATED FROM THE K-TH EQUATION. */

	    if (km1 == 0) {
		goto L50;
	    }
	    sossol_(&k, n, &km1, &y[1], &c__[1], &b[1], &kn);
	    i__3 = km1;
	    for (j = 1; j <= i__3; ++j) {
		js = is[j];
		x[js] = temp[js] + y[j];
/* L40: */
	    }


/*     EVALUATE THE K-TH EQUATION AND THE INTERMEDIATE COMPUTATION */
/*     FOR THE MAX NORM OF THE RESIDUAL VECTOR. */

L50:
	    f = (*fnc)(&x[1], &k);
/* Computing MAX */
	    r__1 = *fmax, r__2 = dabs(f);
	    *fmax = dmax(r__1,r__2);

/*     IF WE WISH TO PERFORM SEVERAL ITERATIONS USING A FIXED */
/*     FACTORIZATION OF AN APPROXIMATE JACOBIAN,WE NEED ONLY */
/*     UPDATE THE CONSTANT VECTOR. */

	    if (itry < *ncjs) {
		goto L160;
	    }


	    it = 0;

/*     COMPUTE PARTIAL DERIVATIVES THAT ARE REQUIRED IN THE LINEARIZATION */
/*     OF THE K-TH REDUCED EQUATION */

	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		item = is[j];
		hx = x[item];
		h__ = fac[item] * hx;
		if (dabs(h__) <= zero) {
		    h__ = fac[item];
		}
		x[item] = hx + h__;
		if (km1 == 0) {
		    goto L70;
		}
		y[j] = h__;
		sossol_(&k, n, &j, &y[1], &c__[1], &b[1], &kn);
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
		if (dabs(fdif) > uro * dabs(f)) {
		    goto L80;
		}
		fdif = 0.f;
		++it;
L80:
		p[j] = fdif / h__;
/* L90: */
	    }

	    if (it <= *n - k) {
		goto L110;
	    }

/*     ALL COMPUTED PARTIAL DERIVATIVES OF THE K-TH EQUATION */
/*     ARE EFFECTIVELY ZERO.TRY LARGER PERTURBATIONS OF THE */
/*     INDEPENDENT VARIABLES. */

	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		isj = is[j];
		fact = fac[isj] * 100.f;
		if (fact > 1e10f) {
		    goto L340;
		}
		fac[isj] = fact;
/* L100: */
	    }
	    goto L30;

L110:
	    if (k == *n) {
		goto L160;
	    }

/*     ACHIEVE A PIVOTING EFFECT BY CHOOSING THE MAXIMAL DERIVATIVE */
/*     ELEMENT */

	    pmax = 0.f;
	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		test = (r__1 = p[j], dabs(r__1));
		if (test <= pmax) {
		    goto L120;
		}
		pmax = test;
		isv = j;
L120:
		;
	    }
	    if (pmax == 0.f) {
		goto L340;
	    }

/*     SET UP THE COEFFICIENTS FOR THE K-TH ROW OF THE TRIANGULAR */
/*     LINEAR SYSTEM AND SAVE THE PARTIAL DERIVATIVE OF */
/*     LARGEST MAGNITUDE */

	    pmax = p[isv];
	    kk = kn;
	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		if (j == isv) {
		    goto L130;
		}
		c__[kk] = -p[j] / pmax;
L130:
		++kk;
/* L140: */
	    }
	    p[k] = pmax;


	    if (isv == k) {
		goto L160;
	    }

/*     INTERCHANGE THE TWO COLUMNS OF C DETERMINED BY THE */
/*     PIVOTAL STRATEGY */

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

/*     STORE THE COMPONENTS FOR THE CONSTANT VECTOR */

	    b[k] = -f / p[k];

/* L170: */
	}

/*    ******** */
/*    ******** END OF LOOP CREATING THE TRIANGULAR LINEARIZATION MATRIX */
/*    ******** */


/*     SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW SOLUTION */
/*     APPROXIMATION AND OBTAIN THE SOLUTION INCREMENT NORM. */

	--kn;
	y[*n] = b[*n];
	if (*n > 1) {
	    sossol_(n, n, n, &y[1], &c__[1], &b[1], &kn);
	}
	xnorm = 0.f;
	ynorm = 0.f;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    yj = y[j];
/* Computing MAX */
	    r__1 = ynorm, r__2 = dabs(yj);
	    ynorm = dmax(r__1,r__2);
	    js = is[j];
	    x[js] = temp[js] + yj;
/* Computing MAX */
	    r__2 = xnorm, r__3 = (r__1 = x[js], dabs(r__1));
	    xnorm = dmax(r__2,r__3);
/* L180: */
	}


/*     PRINT INTERMEDIATE SOLUTION ITERATES AND RESIDUAL NORM IF DESIRED */

	if (*iprint != -1) {
	    goto L190;
	}
	mm = m - 1;
	io___48.ciunit = loun;
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&(*fmax), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&mm, (ftnlen)sizeof(integer));
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(real));
	}
	e_wsfe();
L190:

/*     TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE AND/OR ABSOLUTE ERROR */
/*     COMPARISON ON SUCCESSIVE APPROXIMATIONS OF EACH SOLUTION VARIABLE) */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    js = is[j];
	    if ((r__2 = y[j], dabs(r__2)) > re * (r__1 = x[js], dabs(r__1)) + 
		    *atolx) {
		goto L210;
	    }
/* L200: */
	}
	if (*fmax <= fmxs) {
	    *iflag = 1;
	}

/*     TEST FOR CONVERGENCE TO A SOLUTION BASED ON RESIDUALS */

L210:
	if (*fmax > *tolf) {
	    goto L220;
	}
	*iflag += 2;
L220:
	if (*iflag > 0) {
	    goto L360;
	}


	if (m > 1) {
	    goto L230;
	}
	fmin = *fmax;
	goto L280;

/*     SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM. */

L230:
	if (*fmax >= fmin) {
	    goto L250;
	}
	mit = m + 1;
	yn1 = ynorm;
	yn2 = yns;
	fn1 = fmxs;
	fmin = *fmax;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    s[j] = x[j];
/* L240: */
	}
	ic = 0;

/*     TEST FOR LIMITING PRECISION CONVERGENCE.  VERY SLOWLY CONVERGENT */
/*     PROBLEMS MAY ALSO BE DETECTED. */

L250:
	if (ynorm > sruro * xnorm) {
	    goto L260;
	}
	if (*fmax < fmxs * .2f || *fmax > fmxs * 5.f) {
	    goto L260;
	}
	if (ynorm < yns * .2f || ynorm > yns * 5.f) {
	    goto L260;
	}
	++icr;
	if (icr < *nsrrc) {
	    goto L270;
	}
	*iflag = 4;
	*fmax = fmin;
	goto L380;
L260:
	icr = 0;

/*     TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME. */

	if (ynorm <= yns * 2.f && *fmax <= fmxs * 2.f) {
	    goto L270;
	}
	++ic;
	if (ic < *nsri) {
	    goto L280;
	}
	*iflag = 7;
	goto L360;
L270:
	ic = 0;

/*     CHECK TO SEE IF NEXT ITERATION CAN USE THE OLD JACOBIAN */
/*     FACTORIZATION */

L280:
	--itry;
	if (itry == 0) {
	    goto L290;
	}
	if (ynorm * 20.f > xnorm) {
	    goto L290;
	}
	if (ynorm > yns * 2.f) {
	    goto L290;
	}
	if (*fmax < fmxs * 2.f) {
	    goto L300;
	}
L290:
	itry = *ncjs;

/*     SAVE THE CURRENT SOLUTION APPROXIMATION AND THE RESIDUAL AND */
/*     SOLUTION INCREMENT NORMS FOR USE IN THE NEXT ITERATION. */

L300:
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    temp[j] = x[j];
/* L310: */
	}
	if (m != mit) {
	    goto L320;
	}
	fn2 = *fmax;
	yn3 = ynorm;
L320:
	fmxs = *fmax;
	yns = ynorm;


/* L330: */
    }

/*    ***************************************** */
/*    **** END OF PRINCIPAL ITERATION LOOP **** */
/*    ***************************************** */


/*     TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED. */
    m = *mxit;
    *iflag = 5;
    if (yn1 > yn2 * 10.f || yn3 > yn1 * 10.f) {
	*iflag = 6;
    }
    if (fn1 > fmin * 5.f || fn2 > fmin * 5.f) {
	*iflag = 6;
    }
    if (*fmax > fmin * 5.f) {
	*iflag = 6;
    }
    goto L360;


/*     A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR. */
L340:
    *iflag = 8;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s[j] = temp[j];
/* L350: */
    }
    goto L380;


L360:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s[j] = x[j];
/* L370: */
    }


L380:
    *mxit = m;
    return 0;
} /* soseqs_ */


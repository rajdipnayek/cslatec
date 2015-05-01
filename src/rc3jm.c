/* rc3jm.f -- translated by f2c (version 12.02.01).
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

/* DECK RC3JM */
/* Subroutine */ int rc3jm_(real *l1, real *l2, real *l3, real *m1, real *
	m2min, real *m2max, real *thrcof, integer *ndim, integer *ier)
{
    /* Initialized data */

    static real zero = 0.f;
    static real eps = .01f;
    static real one = 1.f;
    static real two = 2.f;

    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static integer i__, n;
    static real x, y, a1, c1, c2, m2, m3, x1, x2, x3, y1, y2, y3, dv, a1s, 
	    sum1, sum2, huge__;
    static integer nfin, nlim;
    static real tiny, c1old, sign1, sign2;
    static integer index;
    static real cnorm, ratio;
    static integer lstep;
    extern doublereal r1mach_(integer *);
    static integer nfinp1, nfinp2, nfinp3, nstep2;
    static real oldfac, newfac, sumbac, srhuge, thresh;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real sumfor, sumuni, srtiny;

/* ***BEGIN PROLOGUE  RC3JM */
/* ***PURPOSE  Evaluate the 3j symbol g(M2) = (L1 L2   L3  ) */
/*                                           (M1 M2 -M1-M2) */
/*            for all allowed values of M2, the other parameters */
/*            being held fixed. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C19 */
/* ***TYPE      SINGLE PRECISION (RC3JM-S, DRC3JM-D) */
/* ***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, CLEBSCH-GORDAN COEFFICIENTS, */
/*             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS, */
/*             WIGNER COEFFICIENTS */
/* ***AUTHOR  Gordon, R. G., Harvard University */
/*           Schulten, K., Max Planck Institute */
/* ***DESCRIPTION */

/* *Usage: */

/*        REAL L1, L2, L3, M1, M2MIN, M2MAX, THRCOF(NDIM) */
/*        INTEGER NDIM, IER */

/*        CALL RC3JM (L1, L2, L3, M1, M2MIN, M2MAX, THRCOF, NDIM, IER) */

/* *Arguments: */

/*     L1 :IN      Parameter in 3j symbol. */

/*     L2 :IN      Parameter in 3j symbol. */

/*     L3 :IN      Parameter in 3j symbol. */

/*     M1 :IN      Parameter in 3j symbol. */

/*     M2MIN :OUT  Smallest allowable M2 in 3j symbol. */

/*     M2MAX :OUT  Largest allowable M2 in 3j symbol. */

/*     THRCOF :OUT Set of 3j coefficients generated by evaluating the */
/*                 3j symbol for all allowed values of M2.  THRCOF(I) */
/*                 will contain g(M2MIN+I-1), I=1,2,...,M2MAX-M2MIN+1. */

/*     NDIM :IN    Declared length of THRCOF in calling program. */

/*     IER :OUT    Error flag. */
/*                 IER=0 No errors. */
/*                 IER=1 Either L1.LT.ABS(M1) or L1+ABS(M1) non-integer. */
/*                 IER=2 ABS(L1-L2).LE.L3.LE.L1+L2 not satisfied. */
/*                 IER=3 L1+L2+L3 not an integer. */
/*                 IER=4 M2MAX-M2MIN not an integer. */
/*                 IER=5 M2MAX less than M2MIN. */
/*                 IER=6 NDIM less than M2MAX-M2MIN+1. */

/* *Description: */

/*     Although conventionally the parameters of the vector addition */
/*  coefficients satisfy certain restrictions, such as being integers */
/*  or integers plus 1/2, the restrictions imposed on input to this */
/*  subroutine are somewhat weaker. See, for example, Section 27.9 of */
/*  Abramowitz and Stegun or Appendix C of Volume II of A. Messiah. */
/*  The restrictions imposed by this subroutine are */
/*       1. L1.GE.ABS(M1) and L1+ABS(M1) must be an integer; */
/*       2. ABS(L1-L2).LE.L3.LE.L1+L2; */
/*       3. L1+L2+L3 must be an integer; */
/*       4. M2MAX-M2MIN must be an integer, where */
/*          M2MAX=MIN(L2,L3-M1) and M2MIN=MAX(-L2,-L3-M1). */
/*  If the conventional restrictions are satisfied, then these */
/*  restrictions are met. */

/*     The user should be cautious in using input parameters that do */
/*  not satisfy the conventional restrictions. For example, the */
/*  the subroutine produces values of */
/*       g(M2) = (0.75 1.50   1.75  ) */
/*               (0.25  M2  -0.25-M2) */
/*  for M2=-1.5,-0.5,0.5,1.5 but none of the symmetry properties of the */
/*  3j symbol, set forth on page 1056 of Messiah, is satisfied. */

/*     The subroutine generates g(M2MIN), g(M2MIN+1), ..., g(M2MAX) */
/*  where M2MIN and M2MAX are defined above. The sequence g(M2) is */
/*  generated by a three-term recurrence algorithm with scaling to */
/*  control overflow. Both backward and forward recurrence are used to */
/*  maintain numerical stability. The two recurrence sequences are */
/*  matched at an interior point and are normalized from the unitary */
/*  property of 3j coefficients and Wigner's phase convention. */

/*    The algorithm is suited to applications in which large quantum */
/*  numbers arise, such as in molecular dynamics. */

/* ***REFERENCES  1. Abramowitz, M., and Stegun, I. A., Eds., Handbook */
/*                  of Mathematical Functions with Formulas, Graphs */
/*                  and Mathematical Tables, NBS Applied Mathematics */
/*                  Series 55, June 1964 and subsequent printings. */
/*               2. Messiah, Albert., Quantum Mechanics, Volume II, */
/*                  North-Holland Publishing Company, 1963. */
/*               3. Schulten, Klaus and Gordon, Roy G., Exact recursive */
/*                  evaluation of 3j and 6j coefficients for quantum- */
/*                  mechanical coupling of angular momenta, J Math */
/*                  Phys, v 16, no. 10, October 1975, pp. 1961-1970. */
/*               4. Schulten, Klaus and Gordon, Roy G., Semiclassical */
/*                  approximations to 3j and 6j coefficients for */
/*                  quantum-mechanical coupling of angular momenta, */
/*                  J Math Phys, v 16, no. 10, October 1975, */
/*                  pp. 1971-1988. */
/*               5. Schulten, Klaus and Gordon, Roy G., Recursive */
/*                  evaluation of 3j and 6j coefficients, Computer */
/*                  Phys Comm, v 11, 1976, pp. 269-278. */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   880515  SLATEC prologue added by G. C. Nielson, NBS; parameters */
/*           HUGE and TINY revised to depend on R1MACH. */
/*   891229  Prologue description rewritten; other prologue sections */
/*           revised; MMATCH (location of match point for recurrences) */
/*           removed from argument list; argument IER changed to serve */
/*           only as an error flag (previously, in cases without error, */
/*           it returned the number of scalings); number of error codes */
/*           increased to provide more precise error information; */
/*           program comments revised; SLATEC error handler calls */
/*           introduced to enable printing of error messages to meet */
/*           SLATEC standards. These changes were done by D. W. Lozier, */
/*           M. A. McClain and J. M. Smith of the National Institute */
/*           of Standards and Technology, formerly NBS. */
/*   910415  Mixed type expressions eliminated; variable C1 initialized; */
/*           description of THRCOF expanded. These changes were done by */
/*           D. W. Lozier. */
/* ***END PROLOGUE  RC3JM */



    /* Parameter adjustments */
    --thrcof;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  RC3JM */
    *ier = 0;
/*  HUGE is the square root of one twentieth of the largest floating */
/*  point number, approximately. */
    huge__ = sqrt(r1mach_(&c__2) / 20.f);
    srhuge = sqrt(huge__);
    tiny = 1.f / huge__;
    srtiny = 1.f / srhuge;

/*     MMATCH = ZERO */


/*  Check error conditions 1, 2, and 3. */
    r__1 = *l1 + dabs(*m1) + eps;
    if (*l1 - dabs(*m1) + eps < zero || r_mod(&r__1, &one) >= eps + eps) {
	*ier = 1;
	xermsg_("SLATEC", "RC3JM", "L1-ABS(M1) less than zero or L1+ABS(M1) "
		"not integer.", ier, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)52);
	return 0;
    } else if (*l1 + *l2 - *l3 < -eps || *l1 - *l2 + *l3 < -eps || -(*l1) + *
	    l2 + *l3 < -eps) {
	*ier = 2;
	xermsg_("SLATEC", "RC3JM", "L1, L2, L3 do not satisfy triangular con"
		"dition.", ier, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)47);
	return 0;
    } else /* if(complicated condition) */ {
	r__1 = *l1 + *l2 + *l3 + eps;
	if (r_mod(&r__1, &one) >= eps + eps) {
	    *ier = 3;
	    xermsg_("SLATEC", "RC3JM", "L1+L2+L3 not integer.", ier, &c__1, (
		    ftnlen)6, (ftnlen)5, (ftnlen)21);
	    return 0;
	}
    }


/*  Limits for M2 */
/* Computing MAX */
    r__1 = -(*l2), r__2 = -(*l3) - *m1;
    *m2min = dmax(r__1,r__2);
/* Computing MIN */
    r__1 = *l2, r__2 = *l3 - *m1;
    *m2max = dmin(r__1,r__2);

/*  Check error condition 4. */
    r__1 = *m2max - *m2min + eps;
    if (r_mod(&r__1, &one) >= eps + eps) {
	*ier = 4;
	xermsg_("SLATEC", "RC3JM", "M2MAX-M2MIN not integer.", ier, &c__1, (
		ftnlen)6, (ftnlen)5, (ftnlen)24);
	return 0;
    }
    if (*m2min < *m2max - eps) {
	goto L20;
    }
    if (*m2min < *m2max + eps) {
	goto L10;
    }

/*  Check error condition 5. */
    *ier = 5;
    xermsg_("SLATEC", "RC3JM", "M2MIN greater than M2MAX.", ier, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;


/*  This is reached in case that M2 and M3 can take only one value. */
L10:
/*     MSCALE = 0 */
    r__2 = -one;
    i__1 = (integer) ((r__1 = *l2 - *l3 - *m1, dabs(r__1)) + eps);
    thrcof[1] = pow_ri(&r__2, &i__1) / sqrt(*l1 + *l2 + *l3 + one);
    return 0;

/*  This is reached in case that M1 and M2 take more than one value. */
L20:
/*     MSCALE = 0 */
    nfin = (integer) (*m2max - *m2min + one + eps);
    if (*ndim - nfin >= 0) {
	goto L23;
    } else {
	goto L21;
    }

/*  Check error condition 6. */
L21:
    *ier = 6;
    xermsg_("SLATEC", "RC3JM", "Dimension of result array for 3j coefficient"
	    "s too small.", ier, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)56);
    return 0;



/*  Start of forward recursion from M2 = M2MIN */

L23:
    m2 = *m2min;
    thrcof[1] = srtiny;
    newfac = 0.f;
    c1 = 0.f;
    sum1 = tiny;


    lstep = 1;
L30:
    ++lstep;
    m2 += one;
    m3 = -(*m1) - m2;


    oldfac = newfac;
    a1 = (*l2 - m2 + one) * (*l2 + m2) * (*l3 + m3 + one) * (*l3 - m3);
    newfac = sqrt(a1);


    dv = (*l1 + *l2 + *l3 + one) * (*l2 + *l3 - *l1) - (*l2 - m2 + one) * (*
	    l3 + m3 + one) - (*l2 + m2 - one) * (*l3 - m3 - one);

    if (lstep - 2 <= 0) {
	goto L32;
    } else {
	goto L31;
    }

L31:
    c1old = dabs(c1);
L32:
    c1 = -dv / newfac;

    if (lstep > 2) {
	goto L60;
    }


/*  If M2 = M2MIN + 1, the third term in the recursion equation vanishes, */
/*  hence */

    x = srtiny * c1;
    thrcof[2] = x;
    sum1 += tiny * c1 * c1;
    if (lstep == nfin) {
	goto L220;
    }
    goto L30;


L60:
    c2 = -oldfac / newfac;

/*  Recursion to the next 3j coefficient */
    x = c1 * thrcof[lstep - 1] + c2 * thrcof[lstep - 2];
    thrcof[lstep] = x;
    sumfor = sum1;
    sum1 += x * x;
    if (lstep == nfin) {
	goto L100;
    }

/*  See if last unnormalized 3j coefficient exceeds SRHUGE */

    if (dabs(x) < srhuge) {
	goto L80;
    }

/*  This is reached if last 3j coefficient larger than SRHUGE, */
/*  so that the recursion series THRCOF(1), ... , THRCOF(LSTEP) */
/*  has to be rescaled to prevent overflow */

/*     MSCALE = MSCALE + 1 */
    i__1 = lstep;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((r__1 = thrcof[i__], dabs(r__1)) < srtiny) {
	    thrcof[i__] = zero;
	}
/* L70: */
	thrcof[i__] /= srhuge;
    }
    sum1 /= huge__;
    sumfor /= huge__;
    x /= srhuge;


/*  As long as ABS(C1) is decreasing, the recursion proceeds towards */
/*  increasing 3j values and, hence, is numerically stable.  Once */
/*  an increase of ABS(C1) is detected, the recursion direction is */
/*  reversed. */

L80:
    if (c1old - dabs(c1) <= 0.f) {
	goto L100;
    } else {
	goto L30;
    }


/*  Keep three 3j coefficients around MMATCH for comparison later */
/*  with backward recursion values. */

L100:
/*     MMATCH = M2 - 1 */
    nstep2 = nfin - lstep + 3;
    x1 = x;
    x2 = thrcof[lstep - 1];
    x3 = thrcof[lstep - 2];

/*  Starting backward recursion from M2MAX taking NSTEP2 steps, so */
/*  that forwards and backwards recursion overlap at the three points */
/*  M2 = MMATCH+1, MMATCH, MMATCH-1. */

    nfinp1 = nfin + 1;
    nfinp2 = nfin + 2;
    nfinp3 = nfin + 3;
    thrcof[nfin] = srtiny;
    sum2 = tiny;



    m2 = *m2max + two;
    lstep = 1;
L110:
    ++lstep;
    m2 -= one;
    m3 = -(*m1) - m2;
    oldfac = newfac;
    a1s = (*l2 - m2 + two) * (*l2 + m2 - one) * (*l3 + m3 + two) * (*l3 - m3 
	    - one);
    newfac = sqrt(a1s);
    dv = (*l1 + *l2 + *l3 + one) * (*l2 + *l3 - *l1) - (*l2 - m2 + one) * (*
	    l3 + m3 + one) - (*l2 + m2 - one) * (*l3 - m3 - one);
    c1 = -dv / newfac;
    if (lstep > 2) {
	goto L120;
    }

/*  If M2 = M2MAX + 1 the third term in the recursion equation vanishes */

    y = srtiny * c1;
    thrcof[nfin - 1] = y;
    if (lstep == nstep2) {
	goto L200;
    }
    sumbac = sum2;
    sum2 += y * y;
    goto L110;

L120:
    c2 = -oldfac / newfac;

/*  Recursion to the next 3j coefficient */

    y = c1 * thrcof[nfinp2 - lstep] + c2 * thrcof[nfinp3 - lstep];

    if (lstep == nstep2) {
	goto L200;
    }

    thrcof[nfinp1 - lstep] = y;
    sumbac = sum2;
    sum2 += y * y;


/*  See if last 3j coefficient exceeds SRHUGE */

    if (dabs(y) < srhuge) {
	goto L110;
    }

/*  This is reached if last 3j coefficient larger than SRHUGE, */
/*  so that the recursion series THRCOF(NFIN), ... , THRCOF(NFIN-LSTEP+1) */
/*  has to be rescaled to prevent overflow. */

/*     MSCALE = MSCALE + 1 */
    i__1 = lstep;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index = nfin - i__ + 1;
	if ((r__1 = thrcof[index], dabs(r__1)) < srtiny) {
	    thrcof[index] = zero;
	}
/* L111: */
	thrcof[index] /= srhuge;
    }
    sum2 /= huge__;
    sumbac /= huge__;

    goto L110;



/*  The forward recursion 3j coefficients X1, X2, X3 are to be matched */
/*  with the corresponding backward recursion values Y1, Y2, Y3. */

L200:
    y3 = y;
    y2 = thrcof[nfinp2 - lstep];
    y1 = thrcof[nfinp3 - lstep];


/*  Determine now RATIO such that YI = RATIO * XI  (I=1,2,3) holds */
/*  with minimal error. */

    ratio = (x1 * y1 + x2 * y2 + x3 * y3) / (x1 * x1 + x2 * x2 + x3 * x3);
    nlim = nfin - nstep2 + 1;

    if (dabs(ratio) < one) {
	goto L211;
    }

    i__1 = nlim;
    for (n = 1; n <= i__1; ++n) {
/* L210: */
	thrcof[n] = ratio * thrcof[n];
    }
    sumuni = ratio * ratio * sumfor + sumbac;
    goto L230;

L211:
    ++nlim;
    ratio = one / ratio;
    i__1 = nfin;
    for (n = nlim; n <= i__1; ++n) {
/* L212: */
	thrcof[n] = ratio * thrcof[n];
    }
    sumuni = sumfor + ratio * ratio * sumbac;
    goto L230;

L220:
    sumuni = sum1;


/*  Normalize 3j coefficients */

L230:
    cnorm = one / sqrt((*l1 + *l1 + one) * sumuni);

/*  Sign convention for last 3j coefficient determines overall phase */

    sign1 = r_sign(&one, &thrcof[nfin]);
    r__2 = -one;
    i__1 = (integer) ((r__1 = *l2 - *l3 - *m1, dabs(r__1)) + eps);
    sign2 = pow_ri(&r__2, &i__1);
    if (sign1 * sign2 <= 0.f) {
	goto L235;
    } else {
	goto L236;
    }
L235:
    cnorm = -cnorm;

L236:
    if (dabs(cnorm) < one) {
	goto L250;
    }

    i__1 = nfin;
    for (n = 1; n <= i__1; ++n) {
/* L240: */
	thrcof[n] = cnorm * thrcof[n];
    }
    return 0;

L250:
    thresh = tiny / dabs(cnorm);
    i__1 = nfin;
    for (n = 1; n <= i__1; ++n) {
	if ((r__1 = thrcof[n], dabs(r__1)) < thresh) {
	    thrcof[n] = zero;
	}
/* L251: */
	thrcof[n] = cnorm * thrcof[n];
    }



    return 0;
} /* rc3jm_ */


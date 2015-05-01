/* dintp.f -- translated by f2c (version 12.02.01).
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

/* DECK DINTP */
/* Subroutine */ int dintp_(doublereal *x, doublereal *y, doublereal *xout, 
	doublereal *yout, doublereal *ypout, integer *neqn, integer *kold, 
	doublereal *phi, integer *ivc, integer *iv, integer *kgi, doublereal *
	gi, doublereal *alpha, doublereal *og, doublereal *ow, doublereal *ox,
	 doublereal *oy)
{
    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__[13], g[13], h__;
    static integer i__, j, l, m;
    static doublereal w[13], hi;
    static integer iq, jq, iw;
    static doublereal xi;
    static integer kp1, kp2;
    static doublereal gdi, alp, hmu, xiq, rmu, xim1, gdif, temp1, temp2, 
	    temp3, gamma, sigma;

/* ***BEGIN PROLOGUE  DINTP */
/* ***PURPOSE  Approximate the solution at XOUT by evaluating the */
/*            polynomial computed in DSTEPS at XOUT.  Must be used in */
/*            conjunction with DSTEPS. */
/* ***LIBRARY   SLATEC (DEPAC) */
/* ***CATEGORY  I1A1B */
/* ***TYPE      DOUBLE PRECISION (SINTRP-S, DINTP-D) */
/* ***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE, */
/*             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR, */
/*             SMOOTH INTERPOLANT */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   The methods in subroutine  DSTEPS  approximate the solution near  X */
/*   by a polynomial.  Subroutine  DINTP  approximates the solution at */
/*   XOUT  by evaluating the polynomial there.  Information defining this */
/*   polynomial is passed from  DSTEPS  so  DINTP  cannot be used alone. */

/*   Subroutine DSTEPS is completely explained and documented in the text */
/*   "Computer Solution of Ordinary Differential Equations, the Initial */
/*   Value Problem"  by L. F. Shampine and M. K. Gordon. */

/*   Input to DINTP -- */

/*   The user provides storage in the calling program for the arrays in */
/*   the call list */
/*      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN) */
/*                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10) */
/*   and defines */
/*      XOUT -- point at which solution is desired. */
/*   The remaining parameters are defined in  DSTEPS  and passed to */
/*   DINTP  from that subroutine */

/*   Output from  DINTP -- */

/*      YOUT(*) -- solution at  XOUT */
/*      YPOUT(*) -- derivative of solution at  XOUT */
/*   The remaining parameters are returned unaltered from their input */
/*   values.  Integration with  DSTEPS  may be continued. */

/* ***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP */
/*                 II, Report SAND84-0293, Sandia Laboratories, 1984. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   840201  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DINTP */



/* ***FIRST EXECUTABLE STATEMENT  DINTP */
    /* Parameter adjustments */
    --y;
    --yout;
    --ypout;
    phi_dim1 = *neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --iv;
    --gi;
    --alpha;
    --og;
    --ow;
    --oy;

    /* Function Body */
    kp1 = *kold + 1;
    kp2 = *kold + 2;

    hi = *xout - *ox;
    h__ = *x - *ox;
    xi = hi / h__;
    xim1 = xi - 1.;

/*   INITIALIZE W(*) FOR COMPUTING G(*) */

    xiq = xi;
    i__1 = kp1;
    for (iq = 1; iq <= i__1; ++iq) {
	xiq = xi * xiq;
	temp1 = (doublereal) (iq * (iq + 1));
/* L10: */
	w[iq - 1] = xiq / temp1;
    }

/*   COMPUTE THE DOUBLE INTEGRAL TERM GDI */

    if (*kold <= *kgi) {
	goto L50;
    }
    if (*ivc > 0) {
	goto L20;
    }
    gdi = 1. / temp1;
    m = 2;
    goto L30;
L20:
    iw = iv[*ivc];
    gdi = ow[iw];
    m = *kold - iw + 3;
L30:
    if (m > *kold) {
	goto L60;
    }
    i__1 = *kold;
    for (i__ = m; i__ <= i__1; ++i__) {
/* L40: */
	gdi = ow[kp2 - i__] - alpha[i__] * gdi;
    }
    goto L60;
L50:
    gdi = gi[*kold];

/*   COMPUTE G(*) AND C(*) */

L60:
    g[0] = xi;
    g[1] = xi * .5 * xi;
    c__[0] = 1.;
    c__[1] = xi;
    if (*kold < 2) {
	goto L90;
    }
    i__1 = *kold;
    for (i__ = 2; i__ <= i__1; ++i__) {
	alp = alpha[i__];
	gamma = xim1 * alp + 1.;
	l = kp2 - i__;
	i__2 = l;
	for (jq = 1; jq <= i__2; ++jq) {
/* L70: */
	    w[jq - 1] = gamma * w[jq - 1] - alp * w[jq];
	}
	g[i__] = w[0];
/* L80: */
	c__[i__] = gamma * c__[i__ - 1];
    }

/*   DEFINE INTERPOLATION PARAMETERS */

L90:
    sigma = (w[1] - xim1 * w[0]) / gdi;
    rmu = xim1 * c__[kp1 - 1] / gdi;
    hmu = rmu / h__;

/*   INTERPOLATE FOR THE SOLUTION -- YOUT */
/*   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT */

    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	yout[l] = 0.;
/* L100: */
	ypout[l] = 0.;
    }
    i__1 = *kold;
    for (j = 1; j <= i__1; ++j) {
	i__ = kp2 - j;
	gdif = og[i__] - og[i__ - 1];
	temp2 = g[i__ - 1] - g[i__ - 2] - sigma * gdif;
	temp3 = c__[i__ - 1] - c__[i__ - 2] + rmu * gdif;
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
	    yout[l] += temp2 * phi[l + i__ * phi_dim1];
/* L110: */
	    ypout[l] += temp3 * phi[l + i__ * phi_dim1];
	}
/* L120: */
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	yout[l] = (1. - sigma) * oy[l] + sigma * y[l] + h__ * (yout[l] + (g[0]
		 - sigma * og[1]) * phi[l + phi_dim1]);
/* L130: */
	ypout[l] = hmu * (oy[l] - y[l]) + (ypout[l] + (c__[0] + rmu * og[1]) *
		 phi[l + phi_dim1]);
    }

    return 0;
} /* dintp_ */


/* la05ad.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__0 = 0;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c_n4 = -4;
static integer c_n1 = -1;
static integer c_n3 = -3;
static integer c__6 = 6;
static integer c_n2 = -2;
static integer c__2 = 2;
static integer c_n7 = -7;
static integer c_n5 = -5;

/* DECK LA05AD */
/* Subroutine */ int la05ad_(doublereal *a, integer *ind, integer *nz, 
	integer *ia, integer *n, integer *ip, integer *iw, doublereal *w, 
	doublereal *g, doublereal *u)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[4], a__2[6], a__3[2];
    integer ip_dim1, ip_offset, ind_dim1, ind_offset, iw_dim1, iw_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6[4], i__7[6], i__8[2];
    doublereal d__1, d__2, d__3;
    char ch__1[67], ch__2[78], ch__3[53], ch__4[18], ch__5[50];

    /* Local variables */
    static integer i__, j, k, l, k1, k2;
    static doublereal am;
    static integer kc, nc, ii, kj;
    static doublereal au;
    static integer kl, in, kk, ir, kp, kr, jp, il, kq, ks, kn, klc, kpc, mcp, 
	    kpl;
    static doublereal eps;
    static integer ipp, krl, nzc, knp, ipv;
    static doublereal amax;
    static char xern0[8], xern1[8], xern2[8];
    extern /* Subroutine */ int mc20ad_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *), la05ed_(doublereal *,
	     integer *, integer *, integer *, integer *, integer *, logical *)
	    ;
    static integer jcost, kcost;
    extern doublereal d1mach_(integer *);
    static integer idummy;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), xsetun_(integer *);

    /* Fortran I/O blocks */
    static icilist io___40 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___42 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___44 = { 0, xern0, 0, "(I8)", 8, 1 };
    static icilist io___45 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___46 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___47 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___48 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___49 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  LA05AD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LA05AS-S, LA05AD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM */
/*     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE */
/*     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING */
/*     THE FINAL LETTER =D= IN THE NAMES USED HERE. */
/*     REVISIONS MADE BY R J HANSON, SNLA, AUGUST, 1979. */
/*     REVISED SEP. 13, 1979. */

/*     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES */
/*     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL */
/*     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN */
/*     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES */
/*     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED. */

/* IP(I,1),IP(I,2) POINT TO THE START OF ROW/COL I. */
/* IW(I,1),IW(I,2) HOLD THE NUMBER OF NON-ZEROS IN ROW/COL I. */
/* DURING THE MAIN BODY OF THIS SUBROUTINE THE VECTORS IW(.,3),IW(.,5), */
/*     IW(.,7) ARE USED TO HOLD DOUBLY LINKED LISTS OF ROWS THAT HAVE */
/*     NOT BEEN PIVOTAL AND HAVE EQUAL NUMBERS OF NON-ZEROS. */
/* IW(.,4),IW(.,6),IW(.,8) HOLD SIMILAR LISTS FOR THE COLUMNS. */
/* IW(I,3),IW(I,4) HOLD FIRST ROW/COLUMN TO HAVE I NON-ZEROS */
/*     OR ZERO IF THERE ARE NONE. */
/* IW(I,5), IW(I,6) HOLD ROW/COL NUMBER OF ROW/COL PRIOR TO ROW/COL I */
/*     IN ITS LIST, OR ZERO IF NONE. */
/* IW(I,7), IW(I,8) HOLD ROW/COL NUMBER OF ROW/COL AFTER ROW/COL I */
/*     IN ITS LIST, OR ZERO IF NONE. */
/* FOR ROWS/COLS THAT HAVE BEEN PIVOTAL IW(I,5),IW(I,6) HOLD NEGATION OF */
/*     POSITION OF ROW/COL I IN THE PIVOTAL ORDERING. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  D1MACH, LA05ED, MC20AD, XERMSG, XSETUN */
/* ***COMMON BLOCKS    LA05DD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Added D1MACH to list of DOUBLE PRECISION variables. */
/*   890605  Corrected references to XERRWV.  (WRB) */
/*           (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900402  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  LA05AD */

/* EPS IS THE RELATIVE ACCURACY OF FLOATING-POINT COMPUTATION */
    /* Parameter adjustments */
    --a;
    ind_dim1 = *ia;
    ind_offset = 1 + ind_dim1;
    ind -= ind_offset;
    iw_dim1 = *n;
    iw_offset = 1 + iw_dim1;
    iw -= iw_offset;
    ip_dim1 = *n;
    ip_offset = 1 + ip_dim1;
    ip -= ip_offset;
    --w;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  LA05AD */
    if (first) {
	eps = d1mach_(&c__4) * 2.;
    }
    first = FALSE_;

/*     SET THE OUTPUT UNIT NUMBER FOR THE ERROR PROCESSOR. */
/*     THE USAGE OF THIS ERROR PROCESSOR IS DOCUMENTED IN THE */
/*     SANDIA LABS. TECH. REPT. SAND78-1189, BY R E JONES. */
    xsetun_(&la05dd_1.lp);
    if (*u > 1.) {
	*u = 1.;
    }
    if (*u < eps) {
	*u = eps;
    }
    if (*n < 1) {
	goto L670;
    }
    *g = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = 0.f;
	for (j = 1; j <= 5; ++j) {
	    iw[i__ + j * iw_dim1] = 0;
/* L40: */
	}
/* L50: */
    }

/* FLUSH OUT SMALL ENTRIES, COUNT ELEMENTS IN ROWS AND COLUMNS */
    l = 1;
    la05dd_1.lenu = *nz;
    i__1 = *nz;
    for (idummy = 1; idummy <= i__1; ++idummy) {
	if (l > la05dd_1.lenu) {
	    goto L90;
	}
	i__2 = la05dd_1.lenu;
	for (k = l; k <= i__2; ++k) {
	    if ((d__1 = a[k], abs(d__1)) <= la05dd_1.small) {
		goto L70;
	    }
	    i__ = ind[k + ind_dim1];
	    j = ind[k + (ind_dim1 << 1)];
/* Computing MAX */
	    d__2 = (d__1 = a[k], abs(d__1));
	    *g = max(d__2,*g);
	    if (i__ < 1 || i__ > *n) {
		goto L680;
	    }
	    if (j < 1 || j > *n) {
		goto L680;
	    }
	    ++iw[i__ + iw_dim1];
	    ++iw[j + (iw_dim1 << 1)];
/* L60: */
	}
	goto L90;
L70:
	l = k;
	a[l] = a[la05dd_1.lenu];
	ind[l + ind_dim1] = ind[la05dd_1.lenu + ind_dim1];
	ind[l + (ind_dim1 << 1)] = ind[la05dd_1.lenu + (ind_dim1 << 1)];
	--la05dd_1.lenu;
/* L80: */
    }

L90:
    la05dd_1.lenl = 0;
    la05dd_1.lrow = la05dd_1.lenu;
    la05dd_1.lcol = la05dd_1.lrow;
/* MCP IS THE MAXIMUM NUMBER OF COMPRESSES PERMITTED BEFORE AN */
/*     ERROR RETURN RESULTS. */
/* Computing MAX */
    i__1 = *n / 10;
    mcp = max(i__1,20);
    la05dd_1.ncp = 0;
/* CHECK FOR NULL ROW OR COLUMN AND INITIALIZE IP(I,2) TO POINT */
/*     JUST BEYOND WHERE THE LAST COMPONENT OF COLUMN I OF A WILL */
/*     BE STORED. */
    k = 1;
    i__1 = *n;
    for (ir = 1; ir <= i__1; ++ir) {
	k += iw[ir + (iw_dim1 << 1)];
	ip[ir + (ip_dim1 << 1)] = k;
	for (l = 1; l <= 2; ++l) {
	    if (iw[ir + l * iw_dim1] <= 0) {
		goto L700;
	    }
/* L100: */
	}
/* L110: */
    }
/* REORDER BY ROWS */
/* CHECK FOR DOUBLE ENTRIES WHILE USING THE NEWLY CONSTRUCTED */
/*     ROW FILE TO CONSTRUCT THE COLUMN FILE. NOTE THAT BY PUTTING */
/*    THE ENTRIES IN BACKWARDS AND DECREASING IP(J,2) EACH TIME IT */
/*     IS USED WE AUTOMATICALLY LEAVE IT POINTING TO THE FIRST ELEMENT. */
    mc20ad_(n, &la05dd_1.lenu, &a[1], &ind[(ind_dim1 << 1) + 1], &ip[
	    ip_offset], &ind[ind_dim1 + 1], &c__0);
    kl = la05dd_1.lenu;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	ir = *n + 1 - ii;
	kp = ip[ir + ip_dim1];
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    if (iw[j + iw_dim1 * 5] == ir) {
		goto L660;
	    }
	    iw[j + iw_dim1 * 5] = ir;
	    kr = ip[j + (ip_dim1 << 1)] - 1;
	    ip[j + (ip_dim1 << 1)] = kr;
	    ind[kr + ind_dim1] = ir;
/* L120: */
	}
	kl = kp - 1;
/* L130: */
    }

/* SET UP LINKED LISTS OF ROWS AND COLS WITH EQUAL NUMBERS OF NON-ZEROS. */
    for (l = 1; l <= 2; ++l) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *nz = iw[i__ + l * iw_dim1];
	    in = iw[*nz + (l + 2) * iw_dim1];
	    iw[*nz + (l + 2) * iw_dim1] = i__;
	    iw[i__ + (l + 6) * iw_dim1] = in;
	    iw[i__ + (l + 4) * iw_dim1] = 0;
	    if (in != 0) {
		iw[in + (l + 4) * iw_dim1] = i__;
	    }
/* L140: */
	}
/* L150: */
    }


/* START OF MAIN ELIMINATION LOOP. */
    i__1 = *n;
    for (ipv = 1; ipv <= i__1; ++ipv) {
/* FIND PIVOT. JCOST IS MARKOWITZ COST OF CHEAPEST PIVOT FOUND SO FAR, */
/*     WHICH IS IN ROW IPP AND COLUMN JP. */
	jcost = *n * *n;
/* LOOP ON LENGTH OF COLUMN TO BE SEARCHED */
	i__2 = *n;
	for (*nz = 1; *nz <= i__2; ++(*nz)) {
/* Computing 2nd power */
	    i__3 = *nz - 1;
	    if (jcost <= i__3 * i__3) {
		goto L250;
	    }
	    j = iw[*nz + (iw_dim1 << 2)];
/* SEARCH COLUMNS WITH NZ NON-ZEROS. */
	    i__3 = *n;
	    for (idummy = 1; idummy <= i__3; ++idummy) {
		if (j <= 0) {
		    goto L200;
		}
		kp = ip[j + (ip_dim1 << 1)];
		kl = kp + iw[j + (iw_dim1 << 1)] - 1;
		i__4 = kl;
		for (k = kp; k <= i__4; ++k) {
		    i__ = ind[k + ind_dim1];
		    kcost = (*nz - 1) * (iw[i__ + iw_dim1] - 1);
		    if (kcost >= jcost) {
			goto L180;
		    }
		    if (*nz == 1) {
			goto L170;
		    }
/* FIND LARGEST ELEMENT IN ROW OF POTENTIAL PIVOT. */
		    amax = 0.f;
		    k1 = ip[i__ + ip_dim1];
		    k2 = iw[i__ + iw_dim1] + k1 - 1;
		    i__5 = k2;
		    for (kk = k1; kk <= i__5; ++kk) {
/* Computing MAX */
			d__2 = amax, d__3 = (d__1 = a[kk], abs(d__1));
			amax = max(d__2,d__3);
			if (ind[kk + (ind_dim1 << 1)] == j) {
			    kj = kk;
			}
/* L160: */
		    }
/* PERFORM STABILITY TEST. */
		    if ((d__1 = a[kj], abs(d__1)) < amax * *u) {
			goto L180;
		    }
L170:
		    jcost = kcost;
		    ipp = i__;
		    jp = j;
/* Computing 2nd power */
		    i__5 = *nz - 1;
		    if (jcost <= i__5 * i__5) {
			goto L250;
		    }
L180:
		    ;
		}
		j = iw[j + (iw_dim1 << 3)];
/* L190: */
	    }
/* SEARCH ROWS WITH NZ NON-ZEROS. */
L200:
	    i__ = iw[*nz + iw_dim1 * 3];
	    i__3 = *n;
	    for (idummy = 1; idummy <= i__3; ++idummy) {
		if (i__ <= 0) {
		    goto L240;
		}
		amax = 0.f;
		kp = ip[i__ + ip_dim1];
		kl = kp + iw[i__ + iw_dim1] - 1;
/* FIND LARGEST ELEMENT IN THE ROW */
		i__4 = kl;
		for (k = kp; k <= i__4; ++k) {
/* Computing MAX */
		    d__2 = (d__1 = a[k], abs(d__1));
		    amax = max(d__2,amax);
/* L210: */
		}
		au = amax * *u;
		i__4 = kl;
		for (k = kp; k <= i__4; ++k) {
/* PERFORM STABILITY TEST. */
		    if ((d__1 = a[k], abs(d__1)) < au) {
			goto L220;
		    }
		    j = ind[k + (ind_dim1 << 1)];
		    kcost = (*nz - 1) * (iw[j + (iw_dim1 << 1)] - 1);
		    if (kcost >= jcost) {
			goto L220;
		    }
		    jcost = kcost;
		    ipp = i__;
		    jp = j;
/* Computing 2nd power */
		    i__5 = *nz - 1;
		    if (jcost <= i__5 * i__5) {
			goto L250;
		    }
L220:
		    ;
		}
		i__ = iw[i__ + iw_dim1 * 7];
/* L230: */
	    }
L240:
	    ;
	}

/* PIVOT FOUND. */
/* REMOVE ROWS AND COLUMNS INVOLVED IN ELIMINATION FROM ORDERING VECTORS. */
L250:
	kp = ip[jp + (ip_dim1 << 1)];
	kl = iw[jp + (iw_dim1 << 1)] + kp - 1;
	for (l = 1; l <= 2; ++l) {
	    i__2 = kl;
	    for (k = kp; k <= i__2; ++k) {
		i__ = ind[k + l * ind_dim1];
		il = iw[i__ + (l + 4) * iw_dim1];
		in = iw[i__ + (l + 6) * iw_dim1];
		if (il == 0) {
		    goto L260;
		}
		iw[il + (l + 6) * iw_dim1] = in;
		goto L270;
L260:
		*nz = iw[i__ + l * iw_dim1];
		iw[*nz + (l + 2) * iw_dim1] = in;
L270:
		if (in > 0) {
		    iw[in + (l + 4) * iw_dim1] = il;
		}
/* L280: */
	    }
	    kp = ip[ipp + ip_dim1];
	    kl = kp + iw[ipp + iw_dim1] - 1;
/* L290: */
	}
/* STORE PIVOT */
	iw[ipp + iw_dim1 * 5] = -ipv;
	iw[jp + iw_dim1 * 6] = -ipv;
/* ELIMINATE PIVOTAL ROW FROM COLUMN FILE AND FIND PIVOT IN ROW FILE. */
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    kpc = ip[j + (ip_dim1 << 1)];
	    --iw[j + (iw_dim1 << 1)];
	    klc = kpc + iw[j + (iw_dim1 << 1)];
	    i__3 = klc;
	    for (kc = kpc; kc <= i__3; ++kc) {
		if (ipp == ind[kc + ind_dim1]) {
		    goto L310;
		}
/* L300: */
	    }
L310:
	    ind[kc + ind_dim1] = ind[klc + ind_dim1];
	    ind[klc + ind_dim1] = 0;
	    if (j == jp) {
		kr = k;
	    }
/* L320: */
	}
/* BRING PIVOT TO FRONT OF PIVOTAL ROW. */
	au = a[kr];
	a[kr] = a[kp];
	a[kp] = au;
	ind[kr + (ind_dim1 << 1)] = ind[kp + (ind_dim1 << 1)];
	ind[kp + (ind_dim1 << 1)] = jp;

/* PERFORM ELIMINATION ITSELF, LOOPING ON NON-ZEROS IN PIVOT COLUMN. */
	nzc = iw[jp + (iw_dim1 << 1)];
	if (nzc == 0) {
	    goto L550;
	}
	i__2 = nzc;
	for (nc = 1; nc <= i__2; ++nc) {
	    kc = ip[jp + (ip_dim1 << 1)] + nc - 1;
	    ir = ind[kc + ind_dim1];
/* SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED. */
	    kr = ip[ir + ip_dim1];
	    krl = kr + iw[ir + iw_dim1] - 1;
	    i__3 = krl;
	    for (knp = kr; knp <= i__3; ++knp) {
		if (jp == ind[knp + (ind_dim1 << 1)]) {
		    goto L340;
		}
/* L330: */
	    }
/* BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW. */
L340:
	    am = a[knp];
	    a[knp] = a[kr];
	    a[kr] = am;
	    ind[knp + (ind_dim1 << 1)] = ind[kr + (ind_dim1 << 1)];
	    ind[kr + (ind_dim1 << 1)] = jp;
	    am = -a[kr] / a[kp];
/* COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW. */
	    if (la05dd_1.lrow + iw[ir + iw_dim1] + iw[ipp + iw_dim1] + 
		    la05dd_1.lenl <= *ia) {
		goto L350;
	    }
	    if (la05dd_1.ncp >= mcp || la05dd_1.lenu + iw[ir + iw_dim1] + iw[
		    ipp + iw_dim1] + la05dd_1.lenl > *ia) {
		goto L710;
	    }
	    la05ed_(&a[1], &ind[(ind_dim1 << 1) + 1], &ip[ip_offset], n, &iw[
		    iw_offset], ia, &c_true);
	    kp = ip[ipp + ip_dim1];
	    kr = ip[ir + ip_dim1];
L350:
	    krl = kr + iw[ir + iw_dim1] - 1;
	    kq = kp + 1;
	    kpl = kp + iw[ipp + iw_dim1] - 1;
/* PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W. */
	    if (kq > kpl) {
		goto L370;
	    }
	    i__3 = kpl;
	    for (k = kq; k <= i__3; ++k) {
		j = ind[k + (ind_dim1 << 1)];
		w[j] = a[k];
/* L360: */
	    }
L370:
	    ip[ir + ip_dim1] = la05dd_1.lrow + 1;

/* TRANSFER MODIFIED ELEMENTS. */
	    ind[kr + (ind_dim1 << 1)] = 0;
	    ++kr;
	    if (kr > krl) {
		goto L430;
	    }
	    i__3 = krl;
	    for (ks = kr; ks <= i__3; ++ks) {
		j = ind[ks + (ind_dim1 << 1)];
		au = a[ks] + am * w[j];
		ind[ks + (ind_dim1 << 1)] = 0;
/* IF ELEMENT IS VERY SMALL REMOVE IT FROM U. */
		if (abs(au) <= la05dd_1.small) {
		    goto L380;
		}
/* Computing MAX */
		d__1 = *g, d__2 = abs(au);
		*g = max(d__1,d__2);
		++la05dd_1.lrow;
		a[la05dd_1.lrow] = au;
		ind[la05dd_1.lrow + (ind_dim1 << 1)] = j;
		goto L410;
L380:
		--la05dd_1.lenu;
/* REMOVE ELEMENT FROM COL FILE. */
		k = ip[j + (ip_dim1 << 1)];
		kl = k + iw[j + (iw_dim1 << 1)] - 1;
		iw[j + (iw_dim1 << 1)] = kl - k;
		i__4 = kl;
		for (kk = k; kk <= i__4; ++kk) {
		    if (ind[kk + ind_dim1] == ir) {
			goto L400;
		    }
/* L390: */
		}
L400:
		ind[kk + ind_dim1] = ind[kl + ind_dim1];
		ind[kl + ind_dim1] = 0;
L410:
		w[j] = 0.f;
/* L420: */
	    }

/* SCAN PIVOT ROW FOR FILLS. */
L430:
	    if (kq > kpl) {
		goto L520;
	    }
	    i__3 = kpl;
	    for (ks = kq; ks <= i__3; ++ks) {
		j = ind[ks + (ind_dim1 << 1)];
		au = am * w[j];
		if (abs(au) <= la05dd_1.small) {
		    goto L500;
		}
		++la05dd_1.lrow;
		a[la05dd_1.lrow] = au;
		ind[la05dd_1.lrow + (ind_dim1 << 1)] = j;
		++la05dd_1.lenu;

/* CREATE FILL IN COLUMN FILE. */
		*nz = iw[j + (iw_dim1 << 1)];
		k = ip[j + (ip_dim1 << 1)];
		kl = k + *nz - 1;
		if (*nz == 0) {
		    goto L460;
		}
/* IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY. */
		if (kl != la05dd_1.lcol) {
		    goto L440;
		}
		if (la05dd_1.lcol + la05dd_1.lenl >= *ia) {
		    goto L460;
		}
		++la05dd_1.lcol;
		goto L450;
L440:
		if (ind[kl + 1 + ind_dim1] != 0) {
		    goto L460;
		}
L450:
		ind[kl + 1 + ind_dim1] = ir;
		goto L490;
/* NEW ENTRY HAS TO BE CREATED. */
L460:
		if (la05dd_1.lcol + la05dd_1.lenl + *nz + 1 < *ia) {
		    goto L470;
		}
/* COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY. */
		if (la05dd_1.ncp >= mcp || la05dd_1.lenu + la05dd_1.lenl + *
			nz + 1 >= *ia) {
		    goto L710;
		}
		la05ed_(&a[1], &ind[ind_offset], &ip[(ip_dim1 << 1) + 1], n, &
			iw[(iw_dim1 << 1) + 1], ia, &c_false);
		k = ip[j + (ip_dim1 << 1)];
		kl = k + *nz - 1;
/* TRANSFER OLD ENTRY INTO NEW. */
L470:
		ip[j + (ip_dim1 << 1)] = la05dd_1.lcol + 1;
		if (kl < k) {
		    goto L485;
		}
		i__4 = kl;
		for (kk = k; kk <= i__4; ++kk) {
		    ++la05dd_1.lcol;
		    ind[la05dd_1.lcol + ind_dim1] = ind[kk + ind_dim1];
		    ind[kk + ind_dim1] = 0;
/* L480: */
		}
L485:
/* ADD NEW ELEMENT. */
		++la05dd_1.lcol;
		ind[la05dd_1.lcol + ind_dim1] = ir;
L490:
/* Computing MAX */
		d__1 = *g, d__2 = abs(au);
		*g = max(d__1,d__2);
		iw[j + (iw_dim1 << 1)] = *nz + 1;
L500:
		w[j] = 0.f;
/* L510: */
	    }
L520:
	    iw[ir + iw_dim1] = la05dd_1.lrow + 1 - ip[ir + ip_dim1];

/* STORE MULTIPLIER */
	    if (la05dd_1.lenl + la05dd_1.lcol + 1 <= *ia) {
		goto L530;
	    }
/* COMPRESS COL FILE IF NECESSARY. */
	    if (la05dd_1.ncp >= mcp) {
		goto L710;
	    }
	    la05ed_(&a[1], &ind[ind_offset], &ip[(ip_dim1 << 1) + 1], n, &iw[(
		    iw_dim1 << 1) + 1], ia, &c_false);
L530:
	    k = *ia - la05dd_1.lenl;
	    ++la05dd_1.lenl;
	    a[k] = am;
	    ind[k + ind_dim1] = ipp;
	    ind[k + (ind_dim1 << 1)] = ir;
	    --la05dd_1.lenu;
/* L540: */
	}

/* INSERT ROWS AND COLUMNS INVOLVED IN ELIMINATION IN LINKED LISTS */
/*     OF EQUAL NUMBERS OF NON-ZEROS. */
L550:
	k1 = ip[jp + (ip_dim1 << 1)];
	k2 = iw[jp + (iw_dim1 << 1)] + k1 - 1;
	iw[jp + (iw_dim1 << 1)] = 0;
	for (l = 1; l <= 2; ++l) {
	    if (k2 < k1) {
		goto L570;
	    }
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
		ir = ind[k + l * ind_dim1];
		if (l == 1) {
		    ind[k + l * ind_dim1] = 0;
		}
		*nz = iw[ir + l * iw_dim1];
		if (*nz <= 0) {
		    goto L720;
		}
		in = iw[*nz + (l + 2) * iw_dim1];
		iw[ir + (l + 6) * iw_dim1] = in;
		iw[ir + (l + 4) * iw_dim1] = 0;
		iw[*nz + (l + 2) * iw_dim1] = ir;
		if (in != 0) {
		    iw[in + (l + 4) * iw_dim1] = ir;
		}
/* L560: */
	    }
L570:
	    k1 = ip[ipp + ip_dim1] + 1;
	    k2 = iw[ipp + iw_dim1] + k1 - 2;
/* L580: */
	}
/* L590: */
    }

/* RESET COLUMN FILE TO REFER TO U AND STORE ROW/COL NUMBERS IN */
/*     PIVOTAL ORDER IN IW(.,3),IW(.,4) */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = -iw[i__ + iw_dim1 * 5];
	iw[j + iw_dim1 * 3] = i__;
	j = -iw[i__ + iw_dim1 * 6];
	iw[j + (iw_dim1 << 2)] = i__;
	iw[i__ + (iw_dim1 << 1)] = 0;
/* L600: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kp = ip[i__ + ip_dim1];
	kl = iw[i__ + iw_dim1] + kp - 1;
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    ++iw[j + (iw_dim1 << 1)];
/* L610: */
	}
/* L620: */
    }
    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += iw[i__ + (iw_dim1 << 1)];
	ip[i__ + (ip_dim1 << 1)] = k;
/* L630: */
    }
    la05dd_1.lcol = k - 1;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = iw[ii + iw_dim1 * 3];
	kp = ip[i__ + ip_dim1];
	kl = iw[i__ + iw_dim1] + kp - 1;
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    kn = ip[j + (ip_dim1 << 1)] - 1;
	    ip[j + (ip_dim1 << 1)] = kn;
	    ind[kn + ind_dim1] = i__;
/* L640: */
	}
/* L650: */
    }
    return 0;

/*     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS. */

L660:
    if (la05dd_1.lp > 0) {
	s_wsfi(&io___40);
	do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___42);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__6[0] = 40, a__1[0] = "MORE THAN ONE MATRIX ENTRY.  HERE ROW = ";
	i__6[1] = 8, a__1[1] = xern1;
	i__6[2] = 11, a__1[2] = " AND COL = ";
	i__6[3] = 8, a__1[3] = xern2;
	s_cat(ch__1, a__1, i__6, &c__4, (ftnlen)67);
	xermsg_("SLATEC", "LA05AD", ch__1, &c_n4, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)67);
    }
    *g = -4.f;
    return 0;

L670:
    if (la05dd_1.lp > 0) {
	xermsg_("SLATEC", "LA05AD", "THE ORDER OF THE SYSTEM, N, IS NOT POSI"
		"TIVE.", &c_n1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)44);
    }
    *g = -1.;
    return 0;

L680:
    if (la05dd_1.lp > 0) {
	s_wsfi(&io___44);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___45);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___46);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__7[0] = 12, a__2[0] = "ELEMENT K = ";
	i__7[1] = 8, a__2[1] = xern0;
	i__7[2] = 31, a__2[2] = " IS OUT OF BOUNDS.$$HERE ROW = ";
	i__7[3] = 8, a__2[3] = xern1;
	i__7[4] = 11, a__2[4] = " AND COL = ";
	i__7[5] = 8, a__2[5] = xern2;
	s_cat(ch__2, a__2, i__7, &c__6, (ftnlen)78);
	xermsg_("SLATEC", "LA05AD", ch__2, &c_n3, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)78);
    }
    *g = -3.f;
    return 0;

L700:
    if (la05dd_1.lp > 0) {
	s_wsfi(&io___47);
	do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__8[0] = 45, a__3[0] = "ROW OR COLUMN HAS NO ELEMENTS.  HERE INDEX "
		"= ";
	i__8[1] = 8, a__3[1] = xern1;
	s_cat(ch__3, a__3, i__8, &c__2, (ftnlen)53);
	xermsg_("SLATEC", "LA05AD", ch__3, &c_n2, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)53);
    }
    *g = -2.f;
    return 0;

L710:
    if (la05dd_1.lp > 0) {
	xermsg_("SLATEC", "LA05AD", "LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE"
		" TOO SMALL.", &c_n7, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)50);
    }
    *g = -7.f;
    return 0;

L720:
    ++ipv;
    iw[ipv + iw_dim1] = ir;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = -iw[i__ + (l + 4) * iw_dim1];
	if (ii > 0) {
	    iw[ii + iw_dim1] = i__;
	}
/* L730: */
    }

    if (la05dd_1.lp > 0) {
	s_copy(xern1, "ROWS", (ftnlen)8, (ftnlen)4);
	if (l == 2) {
	    s_copy(xern1, "COLUMNS", (ftnlen)8, (ftnlen)7);
	}
/* Writing concatenation */
	i__8[0] = 10, a__3[0] = "DEPENDANT ";
	i__8[1] = 8, a__3[1] = xern1;
	s_cat(ch__4, a__3, i__8, &c__2, (ftnlen)18);
	xermsg_("SLATEC", "LA05AD", ch__4, &c_n5, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)18);

L740:
	s_wsfi(&io___48);
	do_fio(&c__1, (char *)&iw[i__ + iw_dim1], (ftnlen)sizeof(integer));
	e_wsfi();
	s_copy(xern2, " ", (ftnlen)8, (ftnlen)1);
	if (i__ + 1 <= ipv) {
	    s_wsfi(&io___49);
	    do_fio(&c__1, (char *)&iw[i__ + 1 + iw_dim1], (ftnlen)sizeof(
		    integer));
	    e_wsfi();
	}
/* Writing concatenation */
	i__6[0] = 29, a__1[0] = "DEPENDENT VECTOR INDICES ARE ";
	i__6[1] = 8, a__1[1] = xern1;
	i__6[2] = 5, a__1[2] = " AND ";
	i__6[3] = 8, a__1[3] = xern2;
	s_cat(ch__5, a__1, i__6, &c__4, (ftnlen)50);
	xermsg_("SLATEC", "LA05AD", ch__5, &c_n5, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)50);
	i__ += 2;
	if (i__ <= ipv) {
	    goto L740;
	}
    }
    *g = -5.f;
    return 0;
} /* la05ad_ */


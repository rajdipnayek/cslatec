/* dpsixn.f -- translated by f2c (version 12.02.01).
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

/* DECK DPSIXN */
doublereal dpsixn_(integer *n)
{
    /* Initialized data */

    static doublereal c__[100] = { -.577215664901532861,.422784335098467139,
	    .922784335098467139,1.25611766843180047,1.50611766843180047,
	    1.70611766843180047,1.87278433509846714,2.01564147795561,
	    2.14064147795561,2.25175258906672111,2.35175258906672111,
	    2.44266167997581202,2.52599501330914535,2.60291809023222227,
	    2.6743466616607937,2.74101332832746037,2.80351332832746037,
	    2.86233685773922507,2.91789241329478063,2.97052399224214905,
	    3.02052399224214905,3.06814303986119667,3.11359758531574212,
	    3.15707584618530734,3.19874251285197401,3.23874251285197401,
	    3.27720405131351247,3.31424108835054951,3.34995537406483522,
	    3.38443813268552488,3.41777146601885821,3.45002953053498724,
	    3.48127953053498724,3.51158256083801755,3.5409943255438999,
	    3.56956575411532847,3.59734353189310625,3.62437055892013327,
	    3.65068634839381748,3.67632737403484313,3.70132737403484313,
	    3.72571761793728215,3.74952714174680596,3.77278295570029433,
	    3.79551022842756706,3.81773245064978928,3.83947158108457189,
	    3.86074817682925274,3.88158151016258607,3.9019896734278922,
	    3.9219896734278922,3.9415975165651471,3.96082828579591633,
	    3.97969621032421822,3.99821472884273674,4.01639654702455492,
	    4.03425368988169777,4.05179754953082058,4.06903892884116541,
	    4.08598808138353829,4.10265474805020496,4.11904819067315578,
	    4.13517722293122029,4.15105023880423617,4.16667523880423617,
	    4.18205985418885155,4.1972113693403667,4.21213674247469506,
	    4.22684262482763624,4.24133537845082464,4.25562109273653893,
	    4.26970559977879245,4.28359448866768134,4.29729311880466764,
	    4.31080663231818115,4.32413996565151449,4.33729786038835659,
	    4.35028487337536958,4.3631053861958824,4.37576361404398366,
	    4.38826361404398366,4.40060929305632934,4.41280441500754886,
	    4.42485260777863319,4.4367573696833951,4.44852207556574804,
	    4.46014998254249223,4.47164423541605544,4.48300787177969181,
	    4.49424382683587158,4.50535493794698269,4.51634394893599368,
	    4.52721351415338499,4.537966202325428,4.54860450019776842,
	    4.55913081598724211,4.56954748265390877,4.57985676100442424,
	    4.5900608426370773,4.6001618527380874 };
    static doublereal b[6] = { .0833333333333333333,-.00833333333333333333,
	    .00396825396825396825,-.00416666666666666666,
	    .00757575757575757576,-.0210927960927960928 };

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static integer k;
    static doublereal s, fn, ax, trm, rfn2, wdtol;
    extern doublereal d1mach_(integer *);

/* ***BEGIN PROLOGUE  DPSIXN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEXINT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (PSIXN-S, DPSIXN-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     This subroutine returns values of PSI(X)=derivative of log */
/*     GAMMA(X), X.GT.0.0 at integer arguments. A table look-up is */
/*     performed for N .LE. 100, and the asymptotic expansion is */
/*     evaluated for N.GT.100. */

/* ***SEE ALSO  DEXINT */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DPSIXN */


/*             DPSIXN(N), N = 1,100 */
/*             COEFFICIENTS OF ASYMPTOTIC EXPANSION */

/* ***FIRST EXECUTABLE STATEMENT  DPSIXN */
    if (*n > 100) {
	goto L10;
    }
    ret_val = c__[*n - 1];
    return ret_val;
L10:
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    wdtol = max(d__1,1e-18);
    fn = (doublereal) (*n);
    ax = 1.;
    s = -.5 / fn;
    if (abs(s) <= wdtol) {
	goto L30;
    }
    rfn2 = 1. / (fn * fn);
    for (k = 1; k <= 6; ++k) {
	ax *= rfn2;
	trm = -b[k - 1] * ax;
	if (abs(trm) < wdtol) {
	    goto L30;
	}
	s += trm;
/* L20: */
    }
L30:
    ret_val = s + log(fn);
    return ret_val;
} /* dpsixn_ */


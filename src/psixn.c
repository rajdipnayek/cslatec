/* psixn.f -- translated by f2c (version 12.02.01).
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

/* DECK PSIXN */
doublereal psixn_(integer *n)
{
    /* Initialized data */

    static real c__[100] = { -.577215664901532861f,.422784335098467139f,
	    .922784335098467139f,1.25611766843180047f,1.50611766843180047f,
	    1.70611766843180047f,1.87278433509846714f,2.01564147795561f,
	    2.14064147795561f,2.25175258906672111f,2.35175258906672111f,
	    2.44266167997581202f,2.52599501330914535f,2.60291809023222227f,
	    2.6743466616607937f,2.74101332832746037f,2.80351332832746037f,
	    2.86233685773922507f,2.91789241329478063f,2.97052399224214905f,
	    3.02052399224214905f,3.06814303986119667f,3.11359758531574212f,
	    3.15707584618530734f,3.19874251285197401f,3.23874251285197401f,
	    3.27720405131351247f,3.31424108835054951f,3.34995537406483522f,
	    3.38443813268552488f,3.41777146601885821f,3.45002953053498724f,
	    3.48127953053498724f,3.51158256083801755f,3.5409943255438999f,
	    3.56956575411532847f,3.59734353189310625f,3.62437055892013327f,
	    3.65068634839381748f,3.67632737403484313f,3.70132737403484313f,
	    3.72571761793728215f,3.74952714174680596f,3.77278295570029433f,
	    3.79551022842756706f,3.81773245064978928f,3.83947158108457189f,
	    3.86074817682925274f,3.88158151016258607f,3.9019896734278922f,
	    3.9219896734278922f,3.9415975165651471f,3.96082828579591633f,
	    3.97969621032421822f,3.99821472884273674f,4.01639654702455492f,
	    4.03425368988169777f,4.05179754953082058f,4.06903892884116541f,
	    4.08598808138353829f,4.10265474805020496f,4.11904819067315578f,
	    4.13517722293122029f,4.15105023880423617f,4.16667523880423617f,
	    4.18205985418885155f,4.1972113693403667f,4.21213674247469506f,
	    4.22684262482763624f,4.24133537845082464f,4.25562109273653893f,
	    4.26970559977879245f,4.28359448866768134f,4.29729311880466764f,
	    4.31080663231818115f,4.32413996565151449f,4.33729786038835659f,
	    4.35028487337536958f,4.3631053861958824f,4.37576361404398366f,
	    4.38826361404398366f,4.40060929305632934f,4.41280441500754886f,
	    4.42485260777863319f,4.4367573696833951f,4.44852207556574804f,
	    4.46014998254249223f,4.47164423541605544f,4.48300787177969181f,
	    4.49424382683587158f,4.50535493794698269f,4.51634394893599368f,
	    4.52721351415338499f,4.537966202325428f,4.54860450019776842f,
	    4.55913081598724211f,4.56954748265390877f,4.57985676100442424f,
	    4.5900608426370773f,4.6001618527380874f };
    static real b[6] = { .0833333333333333333f,-.00833333333333333333f,
	    .00396825396825396825f,-.00416666666666666666f,
	    .00757575757575757576f,-.0210927960927960928f };

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static integer k;
    static real s, fn, ax, trm, rfn2, wdtol;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  PSIXN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to EXINT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PSIXN-S, DPSIXN-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     This subroutine returns values of PSI(X)=derivative of log */
/*     GAMMA(X), X .GT. 0.0 at integer arguments. A table look-up is */
/*     performed for N .LE. 100, and the asymptotic expansion is */
/*     evaluated for N .GT. 100. */

/* ***SEE ALSO  EXINT */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  PSIXN */

/* ----------------------------------------------------------------------- */
/*             PSIXN(N), N = 1,100 */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*             COEFFICIENTS OF ASYMPTOTIC EXPANSION */
/* ----------------------------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  PSIXN */
    if (*n > 100) {
	goto L10;
    }
    ret_val = c__[*n - 1];
    return ret_val;
L10:
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    wdtol = dmax(r__1,1e-18f);
    fn = (real) (*n);
    ax = 1.f;
    s = -.5f / fn;
    if (dabs(s) <= wdtol) {
	goto L30;
    }
    rfn2 = 1.f / (fn * fn);
    for (k = 1; k <= 6; ++k) {
	ax *= rfn2;
	trm = -b[k - 1] * ax;
	if (dabs(trm) < wdtol) {
	    goto L30;
	}
	s += trm;
/* L20: */
    }
L30:
    ret_val = s + log(fn);
    return ret_val;
} /* psixn_ */


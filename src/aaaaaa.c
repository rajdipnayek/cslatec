/* aaaaaa.f -- translated by f2c (version 12.02.01).
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

/* DECK AAAAAA */
/* Subroutine */ int aaaaaa_(char *ver, ftnlen ver_len)
{
/* ***BEGIN PROLOGUE  AAAAAA */
/* ***PURPOSE  SLATEC Common Mathematical Library disclaimer and version. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  Z */
/* ***TYPE      ALL (AAAAAA-A) */
/* ***KEYWORDS  DISCLAIMER, DOCUMENTATION, VERSION */
/* ***AUTHOR  SLATEC Common Mathematical Library Committee */
/* ***DESCRIPTION */

/*   The SLATEC Common Mathematical Library is issued by the following */

/*           Air Force Weapons Laboratory, Albuquerque */
/*           Lawrence Livermore National Laboratory, Livermore */
/*           Los Alamos National Laboratory, Los Alamos */
/*           National Institute of Standards and Technology, Washington */
/*           National Energy Research Supercomputer Center, Livermore */
/*           Oak Ridge National Laboratory, Oak Ridge */
/*           Sandia National Laboratories, Albuquerque */
/*           Sandia National Laboratories, Livermore */

/*   All questions concerning the distribution of the library should be */
/*   directed to the NATIONAL ENERGY SOFTWARE CENTER, 9700 Cass Ave., */
/*   Argonne, Illinois  60439, and not to the authors of the subprograms. */

/*                    * * * * * Notice * * * * * */

/*   This material was prepared as an account of work sponsored by the */
/*   United States Government.  Neither the United States, nor the */
/*   Department of Energy, nor the Department of Defense, nor any of */
/*   their employees, nor any of their contractors, subcontractors, or */
/*   their employees, makes any warranty, expressed or implied, or */
/*   assumes any legal liability or responsibility for the accuracy, */
/*   completeness, or usefulness of any information, apparatus, product, */
/*   or process disclosed, or represents that its use would not infringe */
/*   upon privately owned rights. */

/* *Usage: */

/*        CHARACTER * 16 VER */

/*        CALL AAAAAA (VER) */

/* *Arguments: */

/*     VER:OUT   will contain the version number of the SLATEC CML. */

/* *Description: */

/*   This routine contains the SLATEC Common Mathematical Library */
/*   disclaimer and can be used to return the library version number. */

/* ***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro */
/*                 and Lee Walton, Guide to the SLATEC Common Mathema- */
/*                 tical Library, April 10, 1990. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800424  DATE WRITTEN */
/*   890414  REVISION DATE from Version 3.2 */
/*   890713  Routine modified to return version number.  (WRB) */
/*   900330  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   921215  Updated for Version 4.0.  (WRB) */
/*   930701  Updated for Version 4.1.  (WRB) */
/* ***END PROLOGUE  AAAAAA */
/* ***FIRST EXECUTABLE STATEMENT  AAAAAA */
    s_copy(ver, " 4.1", ver_len, (ftnlen)4);
    return 0;
} /* aaaaaa_ */


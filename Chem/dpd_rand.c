/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 15.07.2004  10:31h
  $Id: dpd_rand.c,v 1.3 2004/07/26 15:20:35 tmaeke Exp $
********************************************************************/

/*******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <string.h> 

#include "dpd_rand.h"

/*Defines***********************************************************/

#define max(a,b)  ((a)>(b) ? (a):(b))
#define dmin(a,b)  ((a)<(b) ? (a):(b))


/*Types*************************************************************/

/*Vars**************************************************************/

/*Funcs*************************************************************/


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     The pseudorandom number generator.                              */
/*     This subroutine has been taken as it is.                        */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
float r2s_0_(int n__, int iseeds)
{
    /* Initialized data */

    static int idum1s = 1720212868;
    static int idum2s = 1;
    static int iys = 1720212868;
    static int ivs[32] = { 1720212868,1392842846,1031324961,718590712,
        82237802,1816996195,1529538438,1789446856,156648835,52437849,
        1441478319,36906150,1269685686,1644535938,394503142,310212663,
        1596049480,7553450,322224693,445508654,28884682,643161691,
        407948861,479214492,2124954851,612891482,112933431,1814689225,
        53445315,1904850491,1695805043,1860990862 };

    /* System generated locals */
    int i__1;
    double ret_val, r__1;

    /* Builtin functions */
    int s_rnge();

    /* Local variables */
    static int js, ks;

    /* *************************************************************(GB 12/95) */

    /*   Portable long-period (ca. 2.3 * 10^18) random number generator of     */
    /*   L'ECUYER [1] with BAYS-DURHAM shuffle [2] and added safeguards as     */
    /*   published by PRESS et al. [3] as "ran2" in 1992. In this version      */
    /*   (called "R2S" for distinction) no argument is needed, and the         */
    /*   initialization can be done by an entry point "R2INIS(iseedS)" with    */
    /*   any positive int "iseedS" as seed. The internal state corres-         */
    /*   ponding to  iseedS = 1  is used if no initialization with "R2INIS"    */
    /*   is done before the first call of "R2S".                               */

    /*   "R2S" returns a uniform random number in the interval  ]0., 1.[       */
    /*   (NOTE: The endpoint values are excluded!)                             */
    /*   *  INITIALIZATION of "R2S":           rdummy = R2INIS(iseedS)         */
    /*   *  GENERATION of a random number:     r = R2S()                       */

    /*   NOTE:                                                                 */
    /*   *  "rdummy" is a dummy variable of type REAL.                         */
    /*   *  No variable declaractions are necessary in the calling module.     */
    /*   *  Parameter RNMX=1.-EPS should approximate the largest floating      */
    /*      point value that is less than 1. (EPS for a specific machine       */
    /*      can be determined with subroutine "MACHAR" from chapt. 20.1 of     */
    /*      ref [3]). (For "IBM RISC System/6000" workstations with "AIX XL    */
    /*      Fortran/6000", subroutine MACHAR gives  EPS=1.192092896E-07 ).     */


    /*   REFERENCES:                                                           */

    /*   [1]  P. L'ECUYER, Communications of the ACM, vol. 31 (1988) 742-774.  */
    /*   [2]  in D.E. KNUTH, "Seminumerical Algorithms" (2nd. ed.), vol. 2 of  */
    /*        "The Art of Computer Programming", Addison-Wesley, Reading, MA   */
    /*        (1981) paragraphs 3.2-3.3 .                                      */
    /*   [3]  W.H. PRESS, S.A. TEUKOLSKY, W.T. VETTERLING, and B.P. FLANNERY,  */
    /*        "Numerical Recipes in FORTRAN" (2nd ed.), Cambridge University   */
    /*        Press, Cambridge (1992), chapt. 7.1 .                            */


    /*   TEST OUTPUT (first 35 numbers for iseed=1, in row-wise sequence):     */

    /*   0.285381  0.253358  0.093469  0.608497  0.903420  0.195873  0.462954  */
    /*   0.939021  0.127216  0.415931  0.533739  0.107446  0.399671  0.650371  */
    /*   0.027072  0.076975  0.068986  0.851946  0.637346  0.573097  0.902278  */
    /*   0.887676  0.372177  0.347516  0.207896  0.569131  0.364677  0.392418  */
    /*   0.366707  0.432149  0.153942  0.626891  0.749454  0.341041  0.830393  */

    /* *********************************************************************** */
    if (n__ == 1) {
        goto L_r2inis;
    }

    /* *********************************************************************** */

    /* ---- Compute "idum1S" by  idum1S = mod( IA1*idum1S, IM1 )  without */
    /*     overflow by SCHRAGE's method (see ref. [3]):                   */
    ks = idum1s / 53668;
    idum1s = (idum1s - ks * 53668) * 40014 - ks * 12211;
    if (idum1s < 0) {
        idum1s += 2147483563;
    }
    /* ---- Compute "idum2S" likewise: */
    ks = idum2s / 52774;
    idum2s = (idum2s - ks * 52774) * 40692 - ks * 3791;
    if (idum2s < 0) {
        idum2s += 2147483399;
    }
    /* ---- "jS" will be in the range [1 (1) NTAB] : */
    js = iys / 67108862 + 1;
    /* ---- Here "idum1S" is shuffled, and "idum1S" and "idum2S" are combined */
    /*     to produce output:                                                 */
    iys = ivs[js - 1] - idum2s;
    ivs[js - 1] = idum1s;
    if (iys < 1) {
        iys += 2147483562;
    }
    /* ---- Because users don't expect endpoint values: */
    /* Computing MIN */
    r__1 = iys * (double)4.6566130573917691e-10;
    ret_val = dmin(r__1,(double).99999987999999995);
    return (float) ret_val;
    /* >>>> Initialization: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

L_r2inis:
    /* ---- Be sure to prevent a negative "iseedS" or  iseedS = 0 : */
    /* Computing MAX */
    i__1 = abs(iseeds);
    idum1s = max(i__1,1);
    idum2s = idum1s;
    /* ---- Load the shuffle table (after 8 warm-ups): */
    for (js = 40; js >= 1; --js) {
        ks = idum1s / 53668;
        idum1s = (idum1s - ks * 53668) * 40014 - ks * 12211;
        if (idum1s < 0) {
            idum1s += 2147483563;
        }
        if (js <= 32) {
            ivs[js - 1] = idum1s;
        }
        /* L10: */
    }
    iys = ivs[0];
    ret_val = (double) iys;
    return (float) ret_val;
} /* r2s_ */



/* *********************************************************************** */
float r2s_()
{
    return r2s_0_(0, 0);
}


/* *********************************************************************** */
float r2inis_(int iseeds)
{
    return r2s_0_(1, iseeds);
}

/* *********************************************************************** */
void rand_test()
{
    int i,j;
    r2inis_(1);
    for (i=0; i<5; i++) {
        for (j=0; j<7; j++) {
            printf("  %8.6f", r2s_());
        }
        printf("\n");
    }
}


/*******************************************************************/

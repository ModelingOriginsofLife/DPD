/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 15.07.2004  10:31h
  $Id: dpd_rand.h,v 1.4 2004/07/30 09:13:47 tmaeke Exp $
********************************************************************/

#ifndef _dpd_rand_h_
#define _dpd_rand_h_

/*Include***********************************************************/

/*Defines***********************************************************/

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

float r2inis_(int iseeds);
float r2s_();
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     The pseudorandom number generator.                              */
/*     This subroutine has been taken as it is.                        */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
float r2s_0_(int n__, int iseeds);
void rand_test();

#ifdef __cplusplus
}
#endif

/*Types*************************************************************/

/*Vars**************************************************************/

/*******************************************************************/
#endif  /* _dpd_rand_h_ */

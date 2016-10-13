/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 22.07.2004  11:04h
  $Id: dpd_elec.h,v 1.5 2004/08/09 12:12:15 tmaeke Exp $
********************************************************************/

#ifndef _dpd_elec_h_
#define _dpd_elec_h_

/*Include***********************************************************/

/*Defines***********************************************************/

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

void change_polar(int wie);
void change_efield(int wie); 
void draw_elecs();

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* calculate forces from e-field  (the simple version)                 */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void calc_efield();

void elec_translate();  
void elec_add_electrode(char * tag, int seqno, real_t defpol);


#ifdef __cplusplus
}
#endif

/*Types*************************************************************/

/*Vars**************************************************************/

/*******************************************************************/
#endif  /* _dpd_elec_h_ */

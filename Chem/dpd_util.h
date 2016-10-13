/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 20.07.2004  17:15h
  $Id: dpd_util.h,v 1.6 2004/08/10 15:13:05 tmaeke Exp $
********************************************************************/

#ifndef _dpd_util_h_
#define _dpd_util_h_

/*Include***********************************************************/

/*Defines***********************************************************/

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

void print_info(); 
void print_info2(); 

void set_dt(real_t t);
void change_dt(int wie); 
void change_bond(int wie);
void change_rfac(int wie); 
void savedata(int no);
int loaddata(int no);
int loadinitialdata();

real_t ipow (real_t v, int pw);
void test_pow();

void stopit (char * str);
inline real_t d__nint(real_t i);
inline void Boundaryfunc(fpoint_t * r, int kp);

#ifdef __cplusplus
}
#endif

/*Types*************************************************************/

/*Vars**************************************************************/

/*******************************************************************/
#endif  /* _dpd_util_h_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 20.07.2004  11:59h
  $Id: dpd_stat.h,v 1.3 2004/07/30 09:13:47 tmaeke Exp $
********************************************************************/

#ifndef _dpd_stat_h_
#define _dpd_stat_h_

/*Include***********************************************************/

/*Defines***********************************************************/
#define VELO_ANZ 50
#define VACF_SIZE    100

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

    extern  void init_vacf(void);
    extern  void end_vacf(int);

    //extern  int r_gyr__(); 
    //extern  int snaps_();
    //extern  int vacf_();
    extern  int vcmtp_(real_t* ,real_t* ,real_t* ,real_t* );
    extern  int velo_prof(void);

#ifdef __cplusplus
}
#endif

/*Types*************************************************************/

/*Vars**************************************************************/
    // velocity profile:    
    extern real_t velo[VELO_ANZ];
    extern int veloanz[VELO_ANZ];
    extern int drawinfoupd;

    extern fpoint_t vel[VACF_SIZE];         // For the velocity autocorrelation function
    extern real_t dcs_vacf__[VACF_SIZE];
    extern real_t dcs_vacf__[VACF_SIZE];
    extern real_t dvacf[VACF_SIZE];
    extern int    it_bound__[VACF_SIZE+1];
    extern int    it_vacf__;
    extern int    ivacf_loop__[2*VACF_SIZE+1];
    extern int    ivacf_counter__;


/*******************************************************************/
#endif  /* _dpd_stat_h_ */

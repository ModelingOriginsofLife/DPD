/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 22.07.2004  13:33h
  $Id: dpd_calc.h,v 1.2 2004/07/30 09:13:47 tmaeke Exp $
********************************************************************/

#ifndef _dpd_calc_h_
#define _dpd_calc_h_

/*Include***********************************************************/

/*Defines***********************************************************/

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Integrate velocities using a time step of "hcd" for             */
/*     conservative and dissipative forces, and a time step            */
/*     of "hr" for random forces. The routine integrates               */
/*     velocities over a time increment of dt/2, and therefore         */
/*     integrates velocities in two parts. This approach is based      */
/*     on the DPD-VV integration scheme (see the reference:            */
/*     [G. Besold, I. Vattulainen, M. Karttunen, and J.M. Polson,      */
/*     Phys. Rev. E Rapid Comm. vol. 62, R7611 (2000)]).               */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int intv_(real_t pminv_, real_t hcd_, real_t hr_);

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Integrate the positions of the particles.                       */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int intr_(real_t dt_);

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Calculate all forces: conservative, dissipative, and random     */
/*     forces.                                                         */
/*     If needed (the int variable update = .true.),                   */
/*     then update the Verlet neighbor table. Otherwise use            */
/*     the present table.                                              */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int fcdr_(int maxnab);

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Calculate dissipative forces only.                              */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int fdr_();

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Check if the Verlet neighbor table should be updated.           */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int check_();

#ifdef __cplusplus
}
#endif

/*Types*************************************************************/

/*Vars**************************************************************/

/*******************************************************************/
#endif  /* _dpd_calc_h_ */

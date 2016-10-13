/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 22.07.2004  13:27h
  $Id: dpd_init.h,v 1.3 2004/08/02 15:11:07 tmaeke Exp $
********************************************************************/

#ifndef _dpd_init_h_
#define _dpd_init_h_

/*Include***********************************************************/

/*Defines***********************************************************/

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

void set_geometry();
void set_particle();
int get_random_particle(int unassigned[], int *highest_unassigned);
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Set positions and velocities.                                   */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void set_positions_velocity();

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Shift (fold) particles back into primary simulation box.        */
/*     NOTE: Calls of FOLDB are allowed only when the Verlet list      */
/*     has to be reconstructed.                                        */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int foldb_();
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Fix some constants:                                                   */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void fix_constants();

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Generate the initial configuration                              */
/*     (positions and velocities).                                     */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int genconf();

void geom_allow_cube(fpoint_t * min, fpoint_t * max);
void geom_deny_all();
void geom_allow_all(); 

#ifdef __cplusplus
}
#endif


/*Types*************************************************************/

/*Vars**************************************************************/

/*******************************************************************/
#endif  /* _dpd_init_h_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
tmm 20.07.2004  08:11h
$Id: dpd_vars.c,v 1.24 2004/09/29 12:38:27 tmaeke Exp $
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

#include "dpd_vars.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/

// simple chemistry
real_t chem_dist;
// BUGBUG, is a fixed distance followed by a prob really the way...
// BUGBUG isn't this just a rehash of the bonding style question?
chem_reax_t chem_reax[MAX_TYP][MAX_TYP];

// consts
ipoint_t izero = {{0,0,0}};
fpoint_t fzero = {{0.0,0.0,0.0}};

// dimension
int n_dim;
int sys_size;                // System size from x,y,z = +/- sys_size/2
int sys_depth; // if set, depth of system for simulating thin systems (x = y > z)
int sys_size2;
int skin_size;
real_t rho0;                // Density

bool_t verbose;
geometry_t geometry_type;          

int display_intv;               // SDL display interval
int poly_intv;    // Interval for dumping data about polymers
int reax_intv;    // Interval for dumping reaction data
int maxstp;                     // Number of time steps after equilibration
int issstp;                     // Number of steps between snapshots
int incstp;                     // Number of steps between samples
int no_calc;                    // count for that

char ctrl_file_name[MAX_NAME] = "";

char initial_position_file[MAX_NAME] = "";  // the file to read in particle positions from
// the only stuff read in is what's spit out into the .out files
// so don't assume you can just do whatever you want in there
// and be aware of what effects whatever parameter file you run it with will have

bool_t wall_at_y0 = 0;             // simulate a wall at y=0
bool_t wall_brakes = 0;            // does the wall brakes particles? vx=0

bool_t drawinfo = 0;               // draw infos in animation
bool_t do_diss = 1;                // dissipative Term rechnen
bool_t do_cons = 1;                // conservativeb Term rechnen
bool_t do_rand = 1;                // random Term rechnen

bool_t with_brick = 0;               // a brick in the way
bool_t init_with_brick = 0;
fpoint_t brick_ul = {{-1.0, -2.0, 0.0}};
fpoint_t brick_or = {{ 1.0,  2.0, 0.0}};

bool_t pressure_to_x = 0;               // accelerate a particle 
real_t pressure_ofs = 1.0;             // velocity increment (if pressure_to_x=1)
real_t pressure_to_x_minpos = -10.1;   // between these two x limits
real_t pressure_to_x_maxpos = -10.0;   // 
real_t pressure_diff_at_y = 0.0;       // above this only pressure_ofs/2 applied (diff pressure for H-Structures)

bool_t add_vx = 0;                      // give once all particles vx +=pressure_ofs
bool_t add_vy = 0;                      // give once all particles vy +=pressure_ofs


int outsiders = 0;

bool_t bounds = 0;

bool_t efeld = 0;                       // E-field on/off
real_t polarity = 1;                   // 1 or -1
real_t emin = 1; // 0.25;              // cutoff around electrode
real_t epsilonfac = 1;                 // an extra factor to dec-/increase for all particles
fpoint_t evector[F_ASIZE];             // electric-field-vector for each unit cube


int    num_electrodes = 0;
e_seq_t e_seq = {
	0, 10, { 1,1, 0,0, -1,-1, 0,0,0,0 }
};

electr_t electrodes[MAX_ELECTRODES] = {
	{ 1, {{-5.0, 8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 0 },
	{ 1, {{0.0,  8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 1 },
	{ 1, {{0.0,  3.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 2 },
	{ 1, {{0.0, -3.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 3 },
	
	{ 1, {{0.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 4 },
	{ 1, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 5 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 6 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 7 },
	
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 8 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0, 9 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0,10 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0,11 },
	
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0,12 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0,13 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0,14 },
	{ 0, {{5.0, -8.0,  0.0}}, {{1.0, 1.0, 1.0}},  1.0,  0,15 },
};

int    num_particle;            // Total particle number
particle_t  liste[F_SIZE];    

int npoint[F_SIZE];             // For the Verlet neighbor tables
int nvsize;                     // For the Verlet neighbor tables
int maxverlet = F_POINTS;
int list_verl[F_POINTS];        // For the Verlet neighbor tables
int update;                     // int variable for the Verlet tables

char allowed[F_ASIZE];          // list of allowed unit cubes for particles
int chainarray[F_ASIZE];        // chain linking update to Verlet neighbors

int nb3Dx[27] = { 0, 1,-1, 0, 1,-1, 0, 1,-1,  0, 1,-1, 0, 1,-1, 0, 1,-1,  0, 1,-1, 0, 1,-1, 0, 1,-1 };
int nb3Dy[27] = { 0, 0, 0, 1, 1, 1,-1,-1,-1,  0, 0, 0, 1, 1, 1,-1,-1,-1,  0, 0, 0, 1, 1, 1,-1,-1,-1 };
int nb3Dz[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1, 1, -1,-1,-1,-1,-1,-1,-1,-1,-1 };

real_t dtrth; 
real_t dth; 
real_t vol;

real_t vcm;               // Center-of-mass velocity
real_t pm;                // Particle mass
real_t pminv;             // Inverse particle mass

fpoint_t dl;               // Size of the simulation box in 2D
fpoint_t dlinv;            // Inverse system size 
ipoint_t dli;              // same as int
ipoint_t dldiff;

real_t dt;                // Time step for the molecular dynamics

int iseed;        // For the random number generator

real_t rfac;      // Scaling factor of the random force
real_t tfac;      // Scaling factor of temperature
real_t skin;      // Skin radius for the Verlet tables
real_t skinsq;    // Skin radius squared
real_t temp;      // Temperature
real_t tmax;      // Time of the simulation after equilibration

param_cons_t cons_parm[MAX_MOLMOL];  // Conservative force strengths && dynamic bonding/breaking
real_t beta_min;
real_t beta_min_force;
real_t d_bond__;                    // Harmonic force constant for the polymerechain
real_t spring_dist;                 // minimal energy length

int dyn_loops = FALSE;
int max_dyn_polylen;
// bond forming and breaking probabilities are handled in cons_parm

// for stiff bonds
real_t bond_angle_strength = 0;

// bookkeeping
reaction_record *reaction_list_root = NULL;
poly_pop_rec *poly_pop_list_root = NULL;

real_t gamma_;                      // Dissipative force strength

real_t sigma;                       // Random force strength

particle_dist_t particle[MOLMOL] = {
	{ "Poly ",  0, 0, { 240, 240,  240 }},
	{ "Water", 60, 0, { 0x00,0x00, 250 }},
	{ "Oil  ",  0, 0, { 250, 0x00,0x00 }},
	{ "Poly2",  0, 0, { 0x00,250, 0x00 }},
	{ "Poly3",  0, 0, { 0x00,250,  250 }}
};

polymere_t polymere[MAX_POLYMERS] = {  // polymere distibution and definition and chains
	{ 20.0, 0, 2, { TYP_OIL,  TYP_POLY2,  -1,        -1,        -1, -1, -1, -1 }},
	{ 0.0, 0, 3, { TYP_POLY,  TYP_OIL,    TYP_POLY2, -1,        -1, -1, -1, -1 }},
	{ 0.0, 0, 3, { TYP_POLY,  TYP_OIL,    TYP_POLY3, -1,        -1, -1, -1, -1 }},
	{ 0.0, 0, 3, { TYP_POLY2, TYP_POLY,   TYP_POLY2, -1,        -1, -1, -1, -1 }},
	{ 0.0, 0, 3, { TYP_POLY,  TYP_POLY2,  TYP_POLY,  -1,        -1, -1, -1, -1 }},
	{ 0.0, 0, 4, { TYP_POLY,  TYP_POLY2,  TYP_POLY,  TYP_POLY2, -1, -1, -1, -1 }},
	{ 0.0, 0, 4, { TYP_POLY,  TYP_POLY2,  TYP_POLY2, TYP_POLY,  -1, -1, -1, -1 }},
	{ 0.0, 0, 7, { TYP_POLY2, TYP_POLY,   TYP_POLY,  TYP_POLY2, TYP_WATER, TYP_OIL, TYP_WATER, -1 }}
}; 

#ifdef WITH_THREAD
SDL_mutex *hndev_mutex=NULL;              
SDL_sem   *sem_waitcalc=NULL, 
*sem_waitend =NULL,
*sem_waitdisp=NULL;  
#endif


/*Funcs*************************************************************/

/*******************************************************************/

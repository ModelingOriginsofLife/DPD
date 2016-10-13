/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
tmm 20.07.2004  08:11h
$Id: dpd_vars.h,v 1.20 2004/09/29 12:38:27 tmaeke Exp $
********************************************************************/

#ifndef _dpd_vars_h_
#define _dpd_vars_h_

//#ifndef MAC_NATIVE
//#undef __APPLE__  /* Added by Uwe to hinder removal of _main entry-point */
//#endif

#define OLD_EFIELD    // old/simple calculations

/*Include***********************************************************/
#ifdef WITH_THREAD
#include "SDL.h"
#include "SDL_thread.h"
#include "SDL_mutex.h"
#include "SDL_etc.h"
#endif


/*Defines***********************************************************/

/* Table of constant values */

#define TRUE 1
#define FALSE 0

#define SYS_SIZEMAX  128       // max. simulation space
#define F_SIZE       1000000 //(SYS_SIZEMAX*SYS_SIZEMAX*SYS_SIZEMAX)   // max. number of particles
#define F_POINTS     200000000 //40000000  // max. number of entries in Verlet-Table

// an entry in units cubed (max (SYS_SIZEMAX+2)^3)
#define F_ASIZE      ((SYS_SIZEMAX+2)*(SYS_SIZEMAX+2)*(SYS_SIZEMAX+2))  

// the different particle types
#define TYP_POLY  2
#define TYP_WATER 1
#define TYP_OIL   2
#define TYP_POLY2 3
#define TYP_POLY3 4

#define NO_TYPE  -1
#define MIN_TYP   0
#define MAX_TYP   8

#define MOLMOL 8    // defines size of alpha-array, number of different particles

#define ALPHA(nfi,nfj)   (cons_parm[(nfi) + (nfj)*MOLMOL].alpha)
#define BETA(nfi,nfj)    (cons_parm[(nfi) + (nfj)*MOLMOL].beta)
#define BETA6(nfi,nfj)   (cons_parm[(nfi) + (nfj)*MOLMOL].beta6)
#define BETA12(nfi,nfj)  (cons_parm[(nfi) + (nfj)*MOLMOL].beta12)
#define BOND_FORM(nfi,nfj)  (cons_parm[(nfi) + (nfj)*MOLMOL].bond_form)
#define BOND_BREAK(nfi,nfj) (cons_parm[(nfi) + (nfj)*MOLMOL].bond_break)

#define MAX_POLYMERS    8      // number of different chains
#define MAX_POLYLEN    8      // maximum polymerechainlength
#define MAX_REACTION 100	// size of char array which holds individual reactions
			// current theoretical max is 4 + 4 + 8 + 7, 2*MAX_POLYLEN + 7

#define MAX_MOLMOL (MOLMOL*MOLMOL)

#define MAX_ELECTRODES 16
#define MAX_SEQUENCE 10

#define MAX_NAME  200       // length of filename

/*******************************************************************/
#define Boundaryf(rij,kk) (dl.k[kk] * d__nint((rij.k[kk]) * dlinv.k[kk]))

#define Boundary(rij,kk) (dl.k[kk] * ((((rij.k[kk]) * dlinv.k[kk])<0) ? \
									  (int)(((rij.k[kk]) * dlinv.k[kk])-0.5): \
									  (int)(((rij.k[kk]) * dlinv.k[kk])+0.5)))

#define NumCube(r) ((int)(r.k[0])*dldiff.k[0] + (int)(r.k[1])*dldiff.k[1] + (int)(r.k[2])*dldiff.k[2] + F_ASIZE/2)
#define ALLOWED(r) (allowed[NumCube(r)])

#define max(a,b)  ((a)>(b) ? (a):(b))

#define DEBUGOUT(str)
//#define DEBUGOUT(str)  fprintf(stderr,"dpd " str "\n")

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
typedef float  real_t;  // try float or double
typedef int    int_t;
typedef int  index_t;
typedef short int  bool_t;
typedef int color_t[3] ;

typedef struct {
	real_t   k[3];    
} fpoint_t;

typedef struct {
	int_t k[3];
} ipoint_t;

typedef enum { geo_all, geo_hole, geo_h, geo_bell, geo_extern } geometry_t;

typedef struct {
	fpoint_t  r;         // Positions of particles
	fpoint_t  v;         // Velocities of particles (x coordinate)
	fpoint_t  fr;        // Random forces
	fpoint_t  fc;        // Conservative forces
	fpoint_t  fd;        // Dissipative forces
	fpoint_t  r0;        // Saved positions
	short int type;      // Defines the particle: Water, oil, monomere ...
	short int subtype;   // reserved for reaction
	bool_t    drawn;     // for drawing
	bool_t    inside;    // is this particle inside simulating space, false after being outflowed prior to being inflowed
	index_t   left;      // bounded neighbor left side (index)
	index_t   right;     // bounded neighbor right side (index)
	index_t   nextindex; // next link in chain, started at chainstack
	index_t   next_outside; // link to the chain of outsiders
} particle_t;

typedef struct {
	real_t alpha;
	real_t beta;
	int beta12;
	int beta6;
	real_t bond_form;
	real_t bond_break;
	// int bond_chem; // the type that the particle will change to when bonded with something
} param_cons_t;

typedef struct {
	int prod1; // what type the one of them becomes
	int prod2; // the type of the other
	// for catalysis just say A + B --> A + C
	real_t prob; // the probability that the reaction takes place
	//bool_t react;  // whether the two particles react at all
} chem_reax_t;

typedef struct {
	float percent;   // percent of all particles
	int   num;       // absolute number of polymeres of this type
	int   len;       // length of this polymere
	int   part[MAX_POLYLEN];   // list of particletypes
} polymere_t;

typedef struct {
	char  *name;         // a name for that type
	float percent;
	int   num;
	color_t color;       // the color on the screen
	real_t epsilon;     // electrical attraction/repulsion : particle specific
} particle_dist_t;

typedef struct {
	bool_t   active;     // 
	fpoint_t position;   // center of electrode
	fpoint_t size;       // position +/- size/2
	real_t   polarity;   // -1 .. 0 .. 1
	int      dutycycle;  // 0 .. 100%
	int      seq_ofs;    // ofs in seq.list
} electr_t;

typedef struct {
	int sequence;                      
	int length;
	int list[MAX_SEQUENCE];
} e_seq_t; 

/* structure for bonded particle bookkeeping */
/* there's a list with the particle in number form, a count, and a link to the next */
// all polymers are recorded in the highest value form, i.e. a 233 particle gets recorded as 332 because 332>233
// BUGBUG note that you can't use particle type 0 anymore because it breaks things (think polymer 0-1-0 as an int)
typedef struct poly_pop_rec {
	int poly;
	uint count;
	struct poly_pop_rec *next;
} poly_pop_rec;

// reaction are also listed greatest particle first, for consistency
typedef struct reaction_record {
	char reaction[MAX_NAME];
	uint count;
	struct reaction_record *next;
} reaction_record;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

extern real_t chem_dist;
extern chem_reax_t chem_reax[MAX_TYP][MAX_TYP];

// extern poly_pop_branch *poly_pop_tree_root;
extern reaction_record *reaction_list_root;
extern poly_pop_rec *poly_pop_list_root;

// consts
extern ipoint_t izero;
extern fpoint_t fzero;

// dimension
extern int n_dim;
extern int sys_size;
extern int sys_depth;
extern int sys_size2;
extern int skin_size;

// siehe genconf
extern bool_t verbose;
extern geometry_t geometry_type;         

extern int display_intv;          // SDL display interval
extern int poly_intv;    // Interval for dumping data about polymers
extern int reax_intv;    // Interval for dumping reaction data
extern int maxstp;                // Number of time steps after equilibration
extern int issstp;                // Number of steps between snapshots
extern int incstp;                // Number of steps between samples
extern int no_calc;

extern char ctrl_file_name[MAX_NAME];

extern char initial_position_file[MAX_NAME];

// some vars to turn on/off features
extern bool_t wall_at_y0;             // simulate a wall at y=0
extern bool_t wall_brakes;            // does the wall brakes particles? vx=0

extern bool_t drawinfo;               // draw infos in animation
extern bool_t do_diss;                // dissipative Term rechnen
extern bool_t do_cons;                // conservativeb Term rechnen
extern bool_t do_rand;                // random Term rechnen

extern bool_t with_brick;               // a brick in the way
extern bool_t init_with_brick;          // 
extern fpoint_t brick_ul;
extern fpoint_t brick_or;

extern bool_t pressure_to_x;                 // accelerate a particle 
extern real_t pressure_ofs;                 // velocity increment (if pressure_to_x=1)
extern real_t pressure_to_x_minpos;         // between these two x limits
extern real_t pressure_to_x_maxpos;         // 
extern real_t pressure_diff_at_y;           // above this only pressure_ofs/2 applied (diff pressure for H-Structures)
extern bool_t add_vx;                        // give once all particles vx +=pressure_ofs
extern bool_t add_vy;                        // give once all particles vy +=pressure_ofs

extern int outsiders;

extern bool_t bounds;

extern bool_t efeld;                       // E-field on/off
extern real_t polarity;                   // 1 or -1
extern real_t emin;                       // cutoff around electrode
extern real_t epsilonfac;                 // an extra factor to dec-/increase for all particles
extern electr_t electrodes[MAX_ELECTRODES];
extern int    num_electrodes;
extern e_seq_t e_seq;
extern fpoint_t evector[F_ASIZE];      // electric-field-vector for each unit cube

extern int    num_particle;       // Total particle number
extern particle_t  liste[F_SIZE];    

extern int npoint[F_SIZE];             // For the Verlet neighbor tables
extern int nvsize;                     // For the Verlet neighbor tables
extern int maxverlet;
extern int list_verl[F_POINTS];        // For the Verlet neighbor tables
extern int update;                     // int variable for the Verlet tables

extern char allowed[F_ASIZE];

// chain linking update to Verlet neighbors
extern int chainarray[F_ASIZE];        

extern int nb3Dx[27];
extern int nb3Dy[27];
extern int nb3Dz[27];

extern real_t dtrth; 
extern real_t dth; 
extern real_t vol;

extern real_t vcm;               // Center-of-mass velocity
extern real_t pm;                // Particle mass
extern real_t pminv;             // Inverse particle mass

extern fpoint_t dl;               // Size of the simulation box in 2D
extern fpoint_t dlinv;            // Inverse system size 
extern ipoint_t dli;              // same as int
extern ipoint_t dldiff;

extern real_t dt;                // Time step for the molecular dynamics

extern int iseed;                 // For the random number generator

extern real_t rfac;              // Scaling factor of the random force
extern real_t tfac;              // Scaling factor of temperature
extern real_t skin;              // Skin radius for the Verlet tables
extern real_t skinsq;            // Skin radius squared
extern real_t temp;              // Temperature
extern real_t tmax;              // Time of the simulation after equilibration

extern real_t rho0;              // Density

extern param_cons_t cons_parm[MAX_MOLMOL];  // Conservative force strengths
extern real_t beta_min;
extern real_t beta_min_force;
extern real_t d_bond__;                    // Harmonic force constant for the polymerechain
extern real_t spring_dist;                 // minimal energy length

extern int dyn_loops;
extern int max_dyn_polylen;

// extern real_t bond_angle; // pi/2
extern real_t bond_angle_strength;

extern real_t gamma_;                      // Dissipative force strength

extern real_t sigma;                       // Random force strength

extern particle_dist_t particle[MOLMOL];    // particle distribution
extern polymere_t polymere[MAX_POLYMERS];   // polymere distibution and defiition

#ifdef WITH_THREAD
extern SDL_mutex *hndev_mutex;              
extern SDL_sem   *sem_waitcalc,
*sem_waitend,
*sem_waitdisp;  
#endif


/*Funcs*************************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/

/*******************************************************************/
#endif  /* _dpd_vars_h_ */

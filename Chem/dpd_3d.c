/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/*                                                                     */
/* based on dpd_2d.f from:                                             */
/*  http://www.lce.hut.fi/research/polymer/softsimu2002/software.shtml */
/*  softsimu2002@lce.hut.fi                                            */
/*                                                                     */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/************************************************************************
$Id: dpd_3d.c,v 1.47 2004/09/29 12:38:27 tmaeke Exp $
************************************************************************/

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/* Notes: see changelog.txt                                            */
/*                                                                     */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <string.h> 
#include <time.h>

#include "SDL.h"

#include "dpd_vars.h"
#include "dpd_util.h"
#include "dpd_rand.h"
#include "dpd_file.h"
#include "dpd_stat.h"
#include "dpd_disp.h"
#include "dpd_elec.h"            
#include "dpd_init.h"
#include "dpd_calc.h"
#include "dpd_geom.h"
#include "dpd_dynamic.h"

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  the calculating thread                                                 */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int calc_all(void *t)
{
#ifdef WITH_THREAD
	if((sem_waitcalc==NULL)||(sem_waitdisp==NULL)) return(-1);
#endif
	
    update = 1;        /* Make sure that the Verlet table will be constructed */
	/*
	 fcdr_(maxverlet);  // Calculate all forces for the initial configuration.
#ifdef OLD_EFIELD
	 if (efeld) { calc_efield(); }      
#else
	 if (efeld) { e_field(); }
#endif  
	 */
    printf("Type 'h' for help\n");
	
    for (no_calc = 0; no_calc < maxstp; ++no_calc) {
		// first statistics and the like
        if (no_calc % incstp == 0) {
            velo_prof();
            /* Take samples of the center-of-mass velocity and the temperature of the system */
            vcmtp_(&pm, &pminv, &vcm, &temp);
            printf("  Step:%6d  Temp:%7.4f\n", no_calc, temp);
            
            /* Calculate the radius of gyration squar and the end-to-end distance squared */
            /* r_gyr__(&dr_gyr__, &dr_ee__); */
            
            /* Update the velocity autocorrelation functions that are being calculated */
            /* vacf_(max_vacf, &it_vacf__, &ivacf_counter__, it_bound__, ivacf_loop__, dvacf); */
        }
		
        if (issstp != 0 && no_calc % issstp == 0) savedata(no_calc);  /* save all positions into files */
		if (poly_intv != 0 && no_calc % poly_intv == 0) output_census();
		if (reax_intv != 0 && no_calc % reax_intv == 0) output_reactions();         
        snapshot();  /* save all positions into files */
		
		// now comes the heavy lifting
        intv_(pminv, dth, dtrth);    /* Integrate velocities over dt/2 */
        intr_(dt);                   /* Integrate particle positions over dt */
        
        check_();                    /* Check if we should update the  Verlet neighbor table */
		
		// if enabled this slows things down, but alleviates weird problems
		// wherein things would go nuts around the stat_intv update with some values
		// update = TRUE;
        if (update) {
			// printf("Updating Verlet list at timestep: %d\n", no_calc);
            /* Particles have wandered outside */
            /* the simulation box, so let us shift them back inside */
            foldb_();
        }
        
        fcdr_(maxverlet);            /* Calculate all forces: (conservative, dissipative, and random forces) */
		
#ifdef OLD_EFIELD
		if (efeld) { calc_efield(); }      
#else
		if (efeld) { e_field(); }
#endif  
		
        intv_(pminv, dth, dtrth);    /* Integrate velocities over dt/2 */
        fdr_();                      /* Calculate dissipative forces */
		
#ifdef WITH_THREAD
		if (no_calc % display_intv == 0) {
			//printf(" calc will post\n");
			SDL_SemPost(sem_waitcalc);  // +1  
										//printf(" calc will wait\n");
			SDL_SemWait(sem_waitdisp);  // -1
										// here was: display();
										//printf(" calc is back\n");
		}
		
		int res;
		res = SDL_SemTryWait(sem_waitend);
		if (res == 0) { 
			//printf(" calculation ends\n");
			return 0; 
		}
#else
		if (no_calc % display_intv == 0) {
			DrawIt();
		}
		if (handleevent()) {
			break;
		}
#endif
    }
    end_vacf(incstp);
#ifdef WITH_THREAD
	SDL_SemPost(sem_waitend);  // calculation ended
#endif
    return (0);
}




/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     This code is for a single homopolymer in a solvent in 2D            */
/*     by dissipative particle dynamics (DPD).                             */
/*     We do the following (main points):                                  */
/*     - generate the initial configuration                                */
/*     - equilibrate the system                                            */
/*     - do the dynamics in equilibrium                                    */
/*       -- calculate the temperature of the system                        */
/*       -- for the chain, calculate the radius of gyration squared        */
/*          as well as the end-to-end distance squared                     */
/*       -- determine the velocity autocorrelation function which          */
/*          is related to the diffusion coefficient of the chain           */
/*     - close the simulation                                              */
/*     and some other things                                               */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ 
int main(int argc, char*argv[])
{
    sys_size = 10;
	sys_depth = 1;
    n_dim = 2;  //3
    rho0 = 10;  //1       /* total density of all particles [3.0]*/
	
    d_bond__ = 10.;       /* harmonic bonding btw. nn beads in a ch */
    spring_dist = 0.01;
	
    max_dyn_polylen = 8;
	
    for (int ai=MIN_TYP; ai<MAX_TYP; ai++) {
		for (int bi=MIN_TYP; bi<MAX_TYP; bi++) {
			BOND_FORM(ai,bi) = -1; // bonds don't form/break by default
			BOND_BREAK(ai,bi) = -1;
			chem_reax[ai][bi].prob = 0; // no chemistry by default
			chem_reax[ai][bi].prob = 0;
		}
    }
	
    chem_dist = 0.5;
	
    pm = 1.;              /* particle mass (identical for all particles) */
    sigma = 3.;           /* dpd parameter */
    set_dt(0.01);         /* init integration time */
    rfac = sqrt(3.);      /* Scaling of random forces */
    iseed = 0;            /* seed for a random number generator */
    
    init_with_brick = 0;
	
    verbose = 0;
    geometry_type = geo_all;          
    geom_allow_all();
    
    display_intv = 1;     // display interval
    poly_intv = 0; // interval for dumping data about polymers
    maxstp = 200;    /* number of time steps after equilibration */
    incstp = 250;         /* how often do we calculate quantities */
    issstp = 0;         /* how often do we take snapshots */
	
    iseed = iseed + time(0);  
	
	//		int i,j;
	for (int i=0; i<MOLMOL; i++) {
		for (int j=0; j<MOLMOL; j++) {
			BETA(i,j) = 1.0;
			BETA12(i,j) = 0;
			BETA6(i,j) = -1;
			ALPHA(i,j) = (i==j) ? 3.0 : 10.0;
		}
	}
	ALPHA(1,3) = 3.0; // just to have some micelle formation
	ALPHA(3,1) = 3.0;	
	
    beta_min =  0.0000000001;
    beta_min_force = 0.1;
    
    particle[TYP_WATER].epsilon = 0.0;
    particle[TYP_OIL].epsilon   = 0.0;
    particle[TYP_POLY].epsilon  = 1.0;   
    particle[TYP_POLY2].epsilon = 1.0;
    epsilonfac = 1;
	
    memset (evector, 0, sizeof(evector));
	
    /* --------------------------------------------------------------- */
    //helpinteract();
    //test_pow();  exit(0);
    //rand_test(); exit(0);
    init_coords();
    //fix_constants();
    //genconf();         /* Generate the initial configuration */
	
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Commandline parameters                                                */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    process_cmdline(argc, argv);    
    
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Fix some constants:                                                   */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Initialize particles                                                  */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    genconf();         /* Generate the initial configuration (init.c)*/
	if (!!strcmp("",initial_position_file)) {
		loadinitialdata();
		do_census();
		print_census();
	}
    kill_bookkeeping_files();
    print_info();
	
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     SDL init setup                                                        */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
	init_sdl(0);
#ifdef WITH_THREAD
	if ( SDL_InitSubSystem(SDL_INIT_EVENTTHREAD) < 0 ) {
		fprintf(stderr, "Couldn't initialize SDL event threading: %s\n",SDL_GetError());
		exit(1);
	}
	hndev_mutex=SDL_CreateMutex();
	sem_waitdisp=SDL_CreateSemaphore(0);
	sem_waitcalc=SDL_CreateSemaphore(0);
	sem_waitend=SDL_CreateSemaphore(0);
#endif
    display();
    
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Time to start the actual DPD simulation                               */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
	
#ifdef WITH_THREAD
	if(sem_waitcalc && sem_waitdisp && sem_waitend) {
		SDL_Thread *calc_thread;
		int thread_state;
		
		calc_thread = SDL_CreateThread(calc_all,(void *)(long)1);
		while (1) {
			int res; 
			res = SDL_SemTryWait(sem_waitcalc);  // -1  calculation ready?
			if (res == 0) {
				SDL_mutexP(hndev_mutex); 
				display();
				// DrawIt();
				SDL_mutexV(hndev_mutex);
				SDL_SemPost(sem_waitdisp);  // +1  
			}
			
			res = SDL_SemTryWait(sem_waitend);  // -1  calculation ended?
			if (res == 0) {
				SDL_SemPost(sem_waitdisp);  // +1  
				break;
			}
			
			SDL_mutexP(hndev_mutex); /* Need to lock mutex before calling handleevent */
			res = handleevent();
			SDL_mutexV(hndev_mutex); /* exit critical section */
			if (res != 0) {
				SDL_SemPost(sem_waitend);
				SDL_SemPost(sem_waitdisp);
				break;
			}
			RestCPU();
		}
		
		SDL_WaitThread(calc_thread, &thread_state);
		printf ("   main ends\n");
		
		if (sem_waitdisp!=NULL)    SDL_DestroySemaphore(sem_waitdisp);
		if (sem_waitcalc!=NULL)    SDL_DestroySemaphore(sem_waitcalc);
		if (sem_waitend !=NULL)    SDL_DestroySemaphore(sem_waitend);
	}
	if(hndev_mutex!=NULL)  SDL_DestroyMutex(hndev_mutex);
#else
	calc_all(0);
#endif     
	
    if (strcmp(snapsfn,"dpd.bmp")) {
        snapshot();
    }
	if (poly_intv!=0) output_census();
	if (reax_intv!=0) output_reactions();         
	
    done_coords();
    end_sdl();
    return 0;
} /* main */



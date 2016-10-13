/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
tmm 22.07.2004  13:27h
$Id: dpd_init.c,v 1.9 2004/09/29 09:53:30 tmaeke Exp $
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
#include "dpd_util.h"
#include "dpd_rand.h"
#include "dpd_disp.h"
#include "dpd_stat.h"
#include "dpd_init.h"
#include "dpd_dynamic.h"

/*Defines***********************************************************/
#define DEF_SKIN 2

/*Types*************************************************************/

/*Vars**************************************************************/

/*Funcs*************************************************************/

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void save_frame(FILE* fi)    // for the 3D-Viewer
{            
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, -dl.k[1]/2, -dl.k[2]/2, -dl.k[0]/2, -dl.k[1]/2, +dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", +dl.k[0]/2, -dl.k[1]/2, -dl.k[2]/2, +dl.k[0]/2, -dl.k[1]/2, +dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", +dl.k[0]/2, +dl.k[1]/2, -dl.k[2]/2, +dl.k[0]/2, +dl.k[1]/2, +dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, +dl.k[1]/2, -dl.k[2]/2, -dl.k[0]/2, +dl.k[1]/2, +dl.k[2]/2);
	
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, -dl.k[1]/2, -dl.k[2]/2, -dl.k[0]/2, +dl.k[1]/2, -dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", +dl.k[0]/2, -dl.k[1]/2, -dl.k[2]/2, +dl.k[0]/2, +dl.k[1]/2, -dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", +dl.k[0]/2, -dl.k[1]/2, +dl.k[2]/2, +dl.k[0]/2, +dl.k[1]/2, +dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, -dl.k[1]/2, +dl.k[2]/2, -dl.k[0]/2, +dl.k[1]/2, +dl.k[2]/2);
	
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, -dl.k[1]/2, -dl.k[2]/2, +dl.k[0]/2, -dl.k[1]/2, -dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, -dl.k[1]/2, +dl.k[2]/2, +dl.k[0]/2, -dl.k[1]/2, +dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, +dl.k[1]/2, +dl.k[2]/2, +dl.k[0]/2, +dl.k[1]/2, +dl.k[2]/2);
    fprintf(fi, "%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", -dl.k[0]/2, +dl.k[1]/2, -dl.k[2]/2, +dl.k[0]/2, +dl.k[1]/2, -dl.k[2]/2);
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void geom_allow_cube(
					 fpoint_t * min,
					 fpoint_t * max)
{
    ipoint_t p;
    const real_t d=-0.00;
    for (p.k[2] = -dli.k[2]/2-1; p.k[2] <= dli.k[2]/2+1; p.k[2]++) {
        for (p.k[1] = -dli.k[1]/2-1; p.k[1] <= dli.k[1]/2+1; p.k[1]++) { 
            for (p.k[0] = -dli.k[0]/2-1; p.k[0] <= dli.k[0]/2+1; p.k[0]++) {
                if (  (p.k[2]-d <= max->k[2])
					  && (p.k[2]+d >= min->k[2]) 
					  && (p.k[1]-d <= max->k[1])
					  && (p.k[1]+d >= min->k[1])
					  && (p.k[0]-d <= max->k[0])
					  && (p.k[0]+d >= min->k[0]) ) {
					
                    ALLOWED(p) = 1;
                }
            }
        }
    }
}    


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void geom_deny_all()
{
    int ip;
    fprintf(stderr,"geom_deny_all()\n");
    for (ip = 0; ip < F_ASIZE; ++ip) {
        allowed[ip] = 0;
    }
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void geom_allow_all()
{
    int ip;
    for (ip = 0; ip < F_ASIZE; ++ip) {
        allowed[ip] = 1;
    }
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int get_allowed_vol()
{
    int ip,s=0;
    for (ip = 0; ip < F_ASIZE; ++ip) {
        if (allowed[ip]) s++;
    }
    //if (n_dim==2) { s = s / (SYS_SIZEMAX+2); }
    
    if (geometry_type == geo_all || geometry_type == geo_bell) {
        s = (int)(dl.k[0] * dl.k[1] * dl.k[2]);            /* Volume (area) of the system */
    }
    printf("Volume is %d\n",s);
    return s;
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void set_geometry()
{    
    FILE * fi;
    if ((fi = fopen("vertices.out", "w")) == NULL) {
        stopit("Fatal error: cannot open vertices.out");
        exit(1);
    }
    
    if (geometry_type == geo_all) {         // all allowed
        printf("      geometry 0\n");
        geom_allow_all();
        save_frame(fi);
    }
    else if (geometry_type == geo_hole) {  
        printf("      geometry 1\n");  // a hole
        ipoint_t p;
        geom_allow_all();
        save_frame(fi);
        for (p.k[0]=-2; p.k[0]<2; p.k[0]++) {
            for (p.k[1]=-2; p.k[1]<1; p.k[1]++) {
                for (p.k[2] = -dli.k[2]/2-1; p.k[2] <= dli.k[2]/2+1; p.k[2]++) {
                    ALLOWED(p) = 0;
                }
            }
        }
    }
    else if (geometry_type == geo_h) {
        printf("      geometry 2\n");  // a H-structure
        ipoint_t p;
        geom_deny_all();
        save_frame(fi);
        for (p.k[2]=-2; p.k[2]<=2; p.k[2]++) {
            for (p.k[1]=-8; p.k[1]<=-4; p.k[1]++) { 
                for (p.k[0] = -dli.k[0]/2-1; p.k[0] <= dli.k[0]/2+1; p.k[0]++) {
                    ALLOWED(p) = 1;
                }
            }
        }
        for (p.k[2]=-2; p.k[2]<=2; p.k[2]++) {
            for (p.k[1]=4; p.k[1]<=8; p.k[1]++) {  
                for (p.k[0] = -dli.k[0]/2-1; p.k[0] <= dli.k[0]/2+1; p.k[0]++) {
                    ALLOWED(p) = 1;
                }
            }
        }
        for (p.k[2]=-2; p.k[2]<=2; p.k[2]++) {
            for (p.k[1]=-3; p.k[1]<=4; p.k[1]++) { 
                for (p.k[0] = -2; p.k[0] <= 2; p.k[0]++) {
                    ALLOWED(p) = 1;
                }
            }
        }
        // Electrodes
		//        electrodes[0].position.k[0] = 0.0; electrodes[0].position.k[1] = 2.0; electrodes[0].position.k[2] = -0.0;
		//        electrodes[1].position.k[0] = 0.0; electrodes[1].position.k[1] =-2.0; electrodes[1].position.k[2] = -0.0;
		//        electrodes[2].position.k[0] = 0.0; electrodes[2].position.k[1] = 6.0; electrodes[2].position.k[2] = -0.0;
		//        electrodes[3].position.k[0] = 0.0; electrodes[3].position.k[1] =-6.0; electrodes[3].position.k[2] = -0.0;
        emin = 0.5; 
        // pressure
        pressure_to_x_minpos = -dl.k[0]/2+1;           // accelerate particles between these two x-limits 
        pressure_to_x_maxpos = -dl.k[0]/2+1+0.1;       //
        pressure_diff_at_y = 0.0;                  // above this only pressure_ofs/2 applied (diff pressure for H-Structures)
    }
    else if (geometry_type == geo_bell) {         // bell shape
        printf("      geometry 3\n");
        geom_allow_all();
        save_frame(fi);
    }
    
    fclose(fi); 
}



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void set_particle()
{
    int ip, kp, ii;
    
    /* ----------------------------------------------------------- */
    /*  calc num_particle                                          */
    /* ----------------------------------------------------------- */
    real_t psum = 0;
    for (kp=0; kp < MOLMOL; kp++) {
        psum += particle[kp].percent;
    }
    for (kp=0; kp < MAX_POLYMERS; kp++) {
        psum += polymere[kp].percent * polymere[kp].len;
    }
    int sum = 0;
    for (kp=0; kp < MOLMOL; kp++) {
        particle[kp].num = (int)(num_particle * particle[kp].percent / psum +0.5);
        sum += particle[kp].num;
		// printf("settting %d to %d\n", kp, particle[kp].num);
		if (particle[kp].num != 0) printf("   num of %-7s = %d\n", particle[kp].name, particle[kp].num);
    }
    for (kp=0; kp < MAX_POLYMERS; kp++) {
        polymere[kp].num = (int)(num_particle * polymere[kp].percent * polymere[kp].len / psum +0.5);
        sum += polymere[kp].num;
		int poly=0;
		for (int i=0; i < polymere[kp].len; i++) {
			poly = poly*10 + polymere[kp].part[i];
		}
		if (reverse_int(poly) > poly) {
			poly=reverse_int(poly);
		}
		if (polymere[kp].num != 0) {
			printf("   num of chain %d = %-4d  ", kp, polymere[kp].num);
			int lp;
			for (lp=0; lp< MAX_POLYLEN; lp++) {
				if (polymere[kp].part[lp] != NO_TYPE) printf(" %s", particle[polymere[kp].part[lp]].name);
				else break;
			}
			printf("\n");
		}
    }
    num_particle = sum;
    printf("   num particles  = %d\n", num_particle);
	
    /* ----------------------------------------------------------- */
    /*  initialize particle types                                  */
    /* ----------------------------------------------------------- */
    // clear all
    for (ip = 0; ip < num_particle; ++ip) {
        liste[ip].left = -1;
        liste[ip].right = -1;
        liste[ip].drawn = 0;
        liste[ip].inside = 0;
        liste[ip].type = TYP_WATER; // in the beginning, there was water
        liste[ip].subtype = 0;
    }
	
	// new mechanism assigns types somewhat randomly
	int unassigned[num_particle];
    for (ip = 0; ip < num_particle; ++ip) {
		unassigned[ip] = ip; // at the outset unassigned[n] == n
	}
	int highest_unassigned = num_particle;
	
    for (kp = MIN_TYP; kp<MAX_TYP; kp++) {
        for (ii = 0; ii < particle[kp].num; ++ii) {
			ip = get_random_particle(unassigned, &highest_unassigned);
            liste[ip].type = kp;
            liste[ip].inside = 1;  // this is used
        }
    }                                              
    
    /* ----------------------------------------------------------- */
    /*  initialize chained particles                               */
    /* ----------------------------------------------------------- */
    for (ii=0; ii< MAX_POLYMERS; ii++) {
        if (polymere[ii].num > 0) {  
            int jj, cc=0;
			int last_ip=-1;
            for (jj=0; jj< polymere[ii].num; jj++) {
                if (polymere[ii].part[cc] != NO_TYPE) {
					ip = get_random_particle(unassigned, &highest_unassigned);
                    liste[ip].type = polymere[ii].part[cc];
                    liste[ip].inside = 1;  // this is used
										   //printf(" chainparticle: %d := %s ", ip, particle[polymere[ii].part[cc]].name );
                    if (cc > 0) { // verketten // verkettaboutdit
                        liste[ip].left = last_ip;
                        liste[last_ip].right = ip;
                        //printf(" l=%d", ip-1);
                    }
					last_ip = ip;
                    cc++;
                    if ((polymere[ii].part[cc] == NO_TYPE) ){// || (cc>=MAX_POLYLEN)) {
                        cc = 0;
						last_ip = -1; // shouldn't be necessary, but never hurts
                        //printf("\n");
                    }
                }
            }
            //printf("eof %d\n", ii);
        }
    }
}

int get_random_particle(int unassigned[], int *highest_unassigned) {
//	int i;
//	printf("Particle list:\n");
//	for (i=0; i<*highest_unassigned; i++) {
//		printf("%d ",unassigned[i]);
//	}
//	printf("\n");
	int index, particle;
    index = int (r2s_() * (float)(*highest_unassigned)); // generate random index to list of unassigned
	particle = unassigned[index]; // use the particle there (initially same as index)
//	printf("Index is: %d\tParticle is thus: %d\n", index, particle);

	int tmp;
	// swap(unassigned[index], unassigned[highest_unassigned - 1]); // move that one to the end of the list
	tmp = unassigned[index];
	unassigned[index] = unassigned[*highest_unassigned - 1];
	unassigned[*highest_unassigned - 1] = tmp;
	
	(*highest_unassigned)--; // only look at the beginning of the list
	return particle;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Shift (fold) particles back into primary simulation box.        */
/*     NOTE: Calls of FOLDB are allowed only when the Verlet list      */
/*     has to be reconstructed.                                        */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int foldb_()
{
    int ip, kp;
    // printf("Foldback\n");
	
    for (ip = 1; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].r.k[kp] -= Boundary(liste[ip].r,kp);
        }
    }
    return 0;
} /* foldb_ */






/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Set positions and velocities.                                   */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void set_positions_velocity()
{
    /* Local variables */
    real_t tscal;
    int ip, kp;
    real_t vv, temper;
	
    /* ************************************************************ */    
    /* Generate random initial positions */
    printf("  Init: random positions\n");

	// note that most particles are shuffled around below as well
	for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].r.k[kp]  = (r2s_() - 0.5)*dl.k[kp];
        }                                           
    }
	
    if (geometry_type == geo_bell) {
        // ruedi's bell for systemology 
        // special boundaries
        real_t H, B, C;
		
        H = dl.k[2]/10*5;  
        B = 0.01 * pow(M_PI,4.0/3) / pow(H,4.0/3);
        C = (4*dl.k[1])*(4*dl.k[1]);
		
        for (ip = 0; ip < num_particle; ++ip) {
            int ok = 0;
            do {
                if (liste[ip].left >= 0) {
                    liste[ip].r = liste[liste[ip].left].r;  // same as chained partner
                    liste[ip].r.k[0] = liste[liste[ip].left].r.k[0] + (r2s_()-0.5)*0.0;
                    liste[ip].r.k[1] = liste[liste[ip].left].r.k[1] + (r2s_()-0.5)*0.0;
                    liste[ip].r.k[2] = liste[liste[ip].left].r.k[2] + (r2s_()-0.5)*0.0;
                }
                else { 
                    liste[ip].r.k[0] = (r2s_() - 0.5) * dl.k[0];
                    liste[ip].r.k[1] = (r2s_() - 0.5) * dl.k[1];
                    liste[ip].r.k[2] = (r2s_()) * dl.k[2]/2;
                }
                if (liste[ip].r.k[2] >= 0.0) {
                    real_t c2 = liste[ip].r.k[0]*liste[ip].r.k[0] + liste[ip].r.k[1]*liste[ip].r.k[1];
                    if (c2 < C) {
                        real_t dr = H/(1+B*c2*c2);
                        if (liste[ip].r.k[2] < dr) {
                            ok = 1;
                        }
                    }
                }
            } while (!ok);
        }
    }
    else {
        for (ip = 0; ip < num_particle; ++ip) {
			
            /* for the particles */
            bool_t inside_brick = 0;
            do {
                if (liste[ip].left >= 0) {   // chained ones
					int current = traverse_left(ip);
					while (liste[current].right >= 0) {
						real_t d_scale, r_dist;
						fpoint_t dr_tmp;
						int right = liste[current].right;
						
						r_dist = 0;
						for (kp = 0; kp < n_dim; kp++) {
							dr_tmp.k[kp] = r2s_() - 0.5;
							r_dist += dr_tmp.k[kp] * dr_tmp.k[kp];
							//liste[ip2].r.k[kp] = (r2s_() - 0.5) * dl.k[kp];
						}
						
						r_dist = sqrt(r_dist);
						d_scale = 0.5 / r_dist;
						
						for (kp = 0; kp < n_dim; kp++) {
							liste[right].r.k[kp]  = liste[current].r.k[kp] + dr_tmp.k[kp] * d_scale;
						}
						current = liste[current].right;
					}
				}
				else if (liste[ip].right >= 0) {
					// just ignore it and it should all work out okay
					// really though, particle pairs are set together based on
					// leftmost so you wouldn't want to reset the end of the 
					// chain by letting it fall through
				}
                else {
                    for (kp = 0; kp < n_dim; kp++) {
                        liste[ip].r.k[kp] = (r2s_() - 0.5) * dl.k[kp];
                    }
                }
                if (init_with_brick) {
                    inside_brick =    (liste[ip].r.k[0]>brick_ul.k[0])
					&&(liste[ip].r.k[0]<brick_or.k[0])
					&&(liste[ip].r.k[1]<brick_or.k[1])
					&&(liste[ip].r.k[1]>brick_ul.k[1]);
                }
            } while (! ALLOWED(liste[ip].r) ||  inside_brick);
            //if (ip % 1000 == 0) { printf("  i = %d\n", ip); }
        }
    }
	
	
    /* ************************************************************ */    
    /* Make sure that all particles are */
    /* inside the simulation box */
    foldb_();
    
	
    /* ************************************************************ */    
    printf("  Init: random velocities\n");
	
    /* Then generate random initial */
    /* velocities for the particles */
    fpoint_t vcm = {{ 0,0,0,}};
    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            vv = r2s_() - 0.5;
            liste[ip].v.k[kp] = vv;
            vcm.k[kp] += pm * vv;
        }
    }
    
    
    /* ************************************************************ */    
    /* Satisfy the condition that the */
    /* center-of-mass velocity vanishes */
	
    real_t dnum_particle_inv = 1. / ((real_t)num_particle);  /* Inverse particle number */
	
    for (kp = 0; kp < n_dim; kp++) {
        vcm.k[kp] = vcm.k[kp] * dnum_particle_inv * pminv; 
    }
	
    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].v.k[kp] = liste[ip].v.k[kp] - vcm.k[kp];
        }
    }
    
	
    /* ************************************************************ */    
    temper = 0.;   /* Finally rescale velocities such that */
    for (ip = 0; ip < num_particle; ++ip) {
        /* the kinetic temperature is equal to one */
        for (kp = 0; kp < n_dim; kp++) {
            temper += pm * (liste[ip].v.k[kp] * liste[ip].v.k[kp]);
        }
    }
	
    temper = tfac * temper;
    tscal = 1. / sqrt(temper);
    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].v.k[kp] = liste[ip].v.k[kp] * tscal;
        }
    }
    printf("  Init: done!\n");
} /* genconf */


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Fix some constants:                                                   */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void fix_constants() {
	bool_t stop=FALSE;
	// a little bit obfuscated because I wanted to keep the specific error messages but not add a bunch of braces
    if (rho0 < 0.0001 && (stop=TRUE)) printf("Error: density too low: %f.  Must be >0.0001\n", rho0); //rho0 = 0.0001;
    if (n_dim!=2 && n_dim!=3 && (stop=TRUE)) printf("Error: bad number of dimensions: %d.  Options are 2 or 3\n", n_dim); //n_dim = 2;
    if (sys_size >= SYS_SIZEMAX && (stop=TRUE)) printf("Error: system too large: %d.  Must be <%d\n", sys_size, SYS_SIZEMAX); //sys_size = SYS_SIZEMAX;
    if (sys_size <= 3 && (stop=TRUE)) printf("Error: system too small: %d.  Minimum size is 4\n", sys_size); //sys_size = 4;
    if (sys_depth < 1 && (stop=TRUE)) printf("Error: system too thin: %d.  Minimum depth is 1\n", sys_depth);
    if (sigma <= 0.0 && (stop=TRUE)) printf("Error: sigma must not be <= 0: %f\n", sigma); //sigma = 3.0;
    if (dt <= 0.0 && (stop=TRUE)) printf("Error: dt must not be <= 0.0: %f\n", dt); //dt = 0.001;
    if (rfac <= 0.0 && (stop=TRUE)) printf("Error: rfac too small: %f.  Must be >0\n", rfac); //rfac = sqrt(3);
    if (beta_min < 0.0 && (stop=TRUE)) printf("Error: beta_min too small: %f. Minimum is 0\n", beta_min); //beta_min = 0.1;
	if (stop) stopit("Error: parameter out of range");
	
    if (init_with_brick) with_brick = 1;    
    
	sys_size = sys_size / DEF_SKIN * DEF_SKIN;
    sys_size2 = sys_size/2;
    // fill_nb3D();  // fill the neighbour offset table
    dl.k[0] = sys_size;                    /* system size in x-direction */
    dl.k[1] = sys_size; // MAX_Y/MAX_X;    /* system size in y-direction */
    dl.k[2] = (n_dim == 2)? 1.0 : (sys_depth == 1) ? dl.k[1] : sys_depth;
	
    dlinv.k[0] = 1. / dl.k[0];                /* Inverse system size */
    dlinv.k[1] = 1. / dl.k[1];
    dlinv.k[2] = 1. / dl.k[2];
    
    dli.k[0] = (int)dl.k[0];
    dli.k[1] = (int)dl.k[1];
    dli.k[2] = (int)dl.k[2];
	
    dldiff.k[0] = dli.k[0]*dli.k[1];          // Offset z
    dldiff.k[1] = dli.k[1];                   // Offset y
    dldiff.k[2] = 1;                          // Offset x
	
    if (F_ASIZE < (dli.k[0]+2)*(dli.k[1]+2)*(dli.k[2]+2) ) {
        stopit("Fatal Error: F_ASIZE to small!");
    }             
    printf("    Sizeof(float) = %lu\n", (unsigned long) sizeof(real_t));
	
    set_dt(dt);
    
    vol = get_allowed_vol();
    //vol = dl.k[0] * dl.k[1] * dl.k[2];            /* Volume (area) of the system */
    
    num_particle = (int) (rho0 * vol);
    
    if (num_particle > F_SIZE) {
        printf("Num_particle=%d vol=%f\n", num_particle,vol);
        stopit("Fatal Error: num_particle > F_SIZE");
    }
    
    pminv = 1. / pm;                 /* Inverse particle mass */
    
    gamma_ = sigma * .5 * sigma;     /* The strength of the dissipative force */
	/* (using the fluctuation-dissipation theorem) */
	
    skin = DEF_SKIN;                 /* Skin radius for the Verlet neighbor table */
    skinsq = skin * skin;
    skin_size = (int) (sys_size/skin);
    
    tfac = 1. / (((real_t)num_particle) * 2. - 2.);     /* Factor for temperature control */
	
    sdl_set_size();
}    



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Generate the initial configuration                              */
/*     (positions and velocities).                                     */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int genconf()
{
	
    printf("  Init: \n");
    fix_constants();
    init_vacf();
    
    r2inis_(iseed);/* Initialize the random number generator */
		set_particle();
		
		printf("  Init: geometry\n");
		set_geometry();
		
		set_positions_velocity();
		do_census();
		print_census();

		return 0;
}    




/*******************************************************************/

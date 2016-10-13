/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
tmm 22.07.2004  13:33h
$Id: dpd_calc.c,v 1.11 2004/09/29 09:53:30 tmaeke Exp $
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
#include "dpd_elec.h"
#include "dpd_flow.h"
#include "dpd_init.h"
#include "dpd_calc.h"
#include "dpd_dynamic.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/

/*Funcs*************************************************************/

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
int intv_(
		  real_t pminv_, 
		  real_t hcd_,     // dt/2
		  real_t hr_)      // sqrt(dt)/2
{
    int ip, kp;
	
    for (ip = 0; ip < num_particle; ++ip) {
        if (! liste[ip].inside) continue;
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].v.k[kp] = liste[ip].v.k[kp] 
			+ pminv_ * ((liste[ip].fc.k[kp]*do_cons
						 + liste[ip].fd.k[kp]*do_diss) * hcd_ 
						+ liste[ip].fr.k[kp]*do_rand * hr_);
        }
		
        if (pressure_to_x) {
            if (liste[ip].r.k[0]< pressure_to_x_maxpos && liste[ip].r.k[0]> pressure_to_x_minpos) {
                if (liste[ip].r.k[1] > pressure_diff_at_y) {
                    liste[ip].v.k[0] += pressure_ofs/2;
                }
                else {
                    liste[ip].v.k[0] += pressure_ofs;  // displace to right
                }
            } 
        }
        if (add_vx) {
            liste[ip].v.k[0] += pressure_ofs;
        }
        if (add_vy) {
            liste[ip].v.k[1] += pressure_ofs;
        }
    }
    add_vx = 0;
    add_vy = 0;
    return 0;
} /* intv_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Integrate the positions of the particles.                       */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int intr_(
		  real_t dt_)  // dt
{
    int ip,kp;
    
    static int first = 1;
    static real_t H, B, C;
    
    if (first) {   //   for the bell-shape
        first = 0;
        H = dl.k[2]/10*5;  // Z
        B = 0.01 * pow(M_PI,4.0/3) / pow(H,4.0/3);
        // B = 0.03*10/dl.k[1];  // Y X 
        C = (4*dl.k[1])*(4*dl.k[1]);
    }
	
    for (ip = 0; ip < num_particle; ++ip) {
        if (! liste[ip].inside) continue;
        // new position
        fpoint_t ri = liste[ip].r;  // save old position
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].r.k[kp] = ri.k[kp] + liste[ip].v.k[kp] * dt_;
        }
        
        if (geometry_type == geo_bell) {  //  for the bell-shape
            int ok = 0;
            if (liste[ip].r.k[2] >= 0.0) {
                real_t c2 = liste[ip].r.k[0]*liste[ip].r.k[0] + liste[ip].r.k[1]*liste[ip].r.k[1];
                if (c2 < C) {
                    real_t dr = H/(1+B*c2*c2);
                    if (liste[ip].r.k[2] < dr) {
                        ok = 1;
                    }
                }
            }
            if (!ok) {
                liste[ip].r = ri;
                for (kp = 0; kp < n_dim; kp++) {
                    liste[ip].v.k[kp] = - liste[ip].v.k[kp];
                }
            }
        }
        else {
            int newcube = NumCube(liste[ip].r);
            if (! allowed[newcube]) {  // undo
                for (kp = 0; kp < n_dim; kp++) {
                    if ((int)ri.k[kp] != (int)liste[ip].r.k[kp]) {
                        liste[ip].v.k[kp] = 0; // -liste[ip].v.k[kp];
                    }
                    liste[ip].r.k[kp] = ri.k[kp];    
                }
            }
        }
        
        if (with_brick) {
            if (  ((liste[ip].r.k[0]>brick_ul.k[0])) 
				  &&((liste[ip].r.k[0]<brick_or.k[0]))
				  &&((liste[ip].r.k[1]<brick_or.k[1]))
				  &&((liste[ip].r.k[1]>brick_ul.k[1])) ) {  // inside
                
                // and comes from
                if (ri.k[0] < brick_ul.k[0]) {  // left
                    if (ri.k[1] > brick_or.k[1]) { // upper left
                        liste[ip].r.k[0] = brick_ul.k[0] - (liste[ip].r.k[0] - brick_ul.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else if (ri.k[1] < brick_ul.k[1]) { // lower left
                        liste[ip].r.k[0] = brick_ul.k[0] - (liste[ip].r.k[0] - brick_ul.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else { // left
                        liste[ip].r.k[0] = brick_ul.k[0] - (liste[ip].r.k[0] - brick_ul.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        if (wall_brakes) liste[ip].v.k[1] = 0;
                    }
                }
                else if (ri.k[0] > brick_or.k[0]) {  // right
                    if (ri.k[1] > brick_or.k[1]) { // upper right
                        liste[ip].r.k[0] = brick_or.k[0] + (brick_or.k[0] - liste[ip].r.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else if (ri.k[1] < brick_ul.k[1]) { // lower right
                        liste[ip].r.k[0] = brick_or.k[0] + (brick_or.k[0] - liste[ip].r.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else  { // right
                        liste[ip].r.k[0] = brick_or.k[0] + (brick_or.k[0] - liste[ip].r.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        if (wall_brakes) liste[ip].v.k[1] = 0;
                    }
                }
                else if (ri.k[1] > brick_or.k[1]) {  // top
                    if (ri.k[0] > brick_or.k[0]) { // top right
                        liste[ip].r.k[1] = brick_or.k[1] + (brick_or.k[1] - liste[ip].r.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }
                    else if (ri.k[0] < brick_ul.k[0]) { // top left
                        liste[ip].r.k[1] = brick_or.k[1] + (brick_or.k[1] - liste[ip].r.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }                        
                    else { // top
                        liste[ip].r.k[1] = brick_or.k[1] + (brick_or.k[1] - liste[ip].r.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        if (wall_brakes) liste[ip].v.k[0] = 0;
                    }
                }
                else if (ri.k[1] < brick_ul.k[1]) {  // bottom
                    if (ri.k[0] > brick_or.k[0]) { // bottom right
                        liste[ip].r.k[1] = brick_ul.k[1] - (liste[ip].r.k[1] - brick_ul.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }
                    else if (ri.k[0] < brick_ul.k[0]) { // bottom left
                        liste[ip].r.k[1] = brick_ul.k[1] - (liste[ip].r.k[1] - brick_ul.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }
                    else { // bottom left
                        liste[ip].r.k[1] = brick_ul.k[1] - (liste[ip].r.k[1] - brick_ul.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        if (wall_brakes) liste[ip].v.k[0] = 0;
                    }
                }
            }
        }
        if (wall_at_y0) {
            const real_t t=0.55;  // thickness of the wall
            if ((liste[ip].r.k[1]<t)&&(liste[ip].r.k[1]>-t)) { // inside wall
                if (ri.k[1]>t) {  // before outside top
                    liste[ip].r.k[1] = t+t-liste[ip].r.k[1];
                    liste[ip].v.k[1] = -liste[ip].v.k[1];
                    if (wall_brakes) liste[ip].v.k[0]=0.0;
                }
                else if (ri.k[1]< -t) { // before outside bottom
                    liste[ip].r.k[1] = -t-t-liste[ip].r.k[1];
                    liste[ip].v.k[1] = -liste[ip].v.k[1];
                    if (wall_brakes) liste[ip].v.k[0]=0.0;
                }
            }
        }
        if (inout_flow) {
            if (test_outflow(&liste[ip].r)) {
                flow_out(ip);
            }
        }
    }
    if (inout_flow) {
        flow_in_all();
    }
    return 0;
} /* intr_ */


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cell stack and chain handling                                       */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void fill_cell_chains(int num_particle)
{
    int i,ci,ch, kp;
	
    for (i=0;i<F_ASIZE;i++) chainarray[i]=-1;
    
    for (i=0;i<num_particle;i++) {
        // calculate cell coords
		
        fpoint_t rr = liste[i].r;
        ipoint_t ii = izero;
        
        for (kp = 0; kp < n_dim; kp++) {
			rr.k[kp] -= Boundary(rr,kp);
			//ii.k[kp] = (int)((rr.k[kp] + sys_size2+skin)/skin) % skin_size;  // +skin
			ii.k[kp] = (int)((int)(rr.k[kp] + sys_size2+skin)/skin) % skin_size;  // +skin
        } 
        
        ci=(ii.k[2]*skin_size+ii.k[1])*skin_size+ii.k[0];
		
        // error control
        if (ci>=F_ASIZE)                 {
            printf("Error GT init (particle %d) in indexing ix %d iy %d iz %d   ",i,ii.k[0],ii.k[1],ii.k[2]);
            printf("x %f y %f z %f  is=%d\n",liste[i].r.k[0],liste[i].r.k[1],liste[i].r.k[2],sys_size);
            ci=0;
        }
        else if (ci<0)                  {
            printf("Error LT init (particle %d) in indexing ix %d iy %d iz %d   ",i,ii.k[0],ii.k[1],ii.k[2]);
            printf("x %f y %f z %f  is=%d\n",liste[i].v.k[0],liste[i].v.k[1],liste[i].v.k[2],sys_size);
            ci=0;
        }
		
        // get cell from global stack and put on top of local chain
        ch=chainarray[ci];
        chainarray[ci]=i;
        liste[i].nextindex = ch;
    }
}


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*   The two particle interaction                                     */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int calc_ij(int ip, int jp)
{
    real_t  fcfac, fdfac, frfac;
    int kp;
    real_t omega;
    real_t rijinv,rrij, rijsq;
    int typej, typei;
    fpoint_t eij, rij;
    
    /*--------------------------*/
    /* Vector from "jp" to "ip" */
    /*--------------------------*/
    rijsq = 0;                
    for (kp = 0; kp < n_dim; kp++) {
        rij.k[kp] = liste[ip].r.k[kp] - liste[jp].r.k[kp];
        rij.k[kp] -= Boundary(rij,kp);
        rijsq += rij.k[kp]*rij.k[kp];
    }
	
    if (rijsq < 1.0) {
        /* The two particles interact */
		
        rrij = sqrt(rijsq);   /* Distance between the particles */
        if (rrij < 0.0000001) {
            rijinv = 1000000; //printf(".");
        }
        else {
            rijinv = 1.0 / rrij;
        }
        
        omega = 1.0 - rrij;
        
        typei = liste[ip].type;  /* Type of particle "ip"  */
        typej = liste[jp].type;  /* Type of particle "jp"  */

		/* ------------------- */
		/* Do simple chemistry */
		/* ------------------- */
		if (rrij < chem_dist && chem_reax[typei][typej].prob > 0 && chem_reax[typei][typej].prob > r2s_()) {
			// BUGBUG problem of getting some forces as one type and the rest as another
			// BUGBUG problem of census
			decrease_poly_pop(get_polymer(ip));
			bool_t separate = !joined(ip,jp);
			if (separate) {
				// check first to make sure they're not the same polymer
				decrease_poly_pop(get_polymer(jp));
				separate = TRUE;
			}
			liste[ip].type = chem_reax[typei][typej].prod1;
			liste[jp].type = chem_reax[typei][typej].prod2;
			increase_poly_pop(get_polymer(ip));
			if (separate) {
				increase_poly_pop(get_polymer(jp));				
			}
		}

		/* --------------------------------- */
		/* see if particles dynamically bond */
		/* --------------------------------- */
#ifdef PROBABALISTIC_BONDING
		if (BOND_FORM(typei,typej) > 0 && BOND_FORM(typei,typej) > r2s_() && particles_should_dynamically_bond(ip, jp))
#else
		if (BOND_FORM(typei,typej) > 0 && BOND_FORM(typei,typej) >= rrij && particles_should_dynamically_bond(ip, jp))
#endif
		{
			int bond_error = dynamically_bond_particles(ip, jp);
			if (bond_error != 0) { printf("error bonding\n"); }
		}
			
		/*-----------------*/
		/* calc fc-factor  */
		/*-----------------*/
/*
#               /   1          beta    \
#  Fc = alpha * | --------- - -------- |
#               \ rij^exp12   rij^exp6 /
#
# if exp12==0 and exp6==-1 then
#  Fc = alpha * (1 - beta * rij) = -alpha * beta * rij + alpha
*/

#ifdef SIMPLE_CONS
			fcfac = ALPHA(typei,typej) * (1.0 - BETA(typei,typej) * rrij);
/*
 Suffers from cutoff problem:
 
   |
   |       DPD cutoff radius (1 dpd-unit)
 F |\    :
 O | \   :
 R |--\--:----------
 C |   \ :
 E |    \: particles feel significant attractive
   |     : force immediately at cutoff distance
   |
      DISTANCE ->
 
 */
			//fcfac = omega  * ALPHA(typei,typej);
#else
		fcfac = ALPHA(typei,typej) * (1.0 * ipow(rrij, -BETA12(typei,typej)) 
									  - BETA(typei,typej) * ipow(rrij, -BETA6(typei,typej)));
#endif
		
		/*-----------------*/
		/* calc fr-factor  */
		/*-----------------*/
		frfac =  omega * rfac * (r2s_() * 2.0 - 1.0) * sigma;
		
		/*-----------------*/
		/* calc fd-factor  */
		/*-----------------*/
		fdfac = 0;
		for (kp = 0; kp < n_dim; kp++) {
			eij.k[kp] = rij.k[kp] * rijinv;
			fdfac += eij.k[kp] * (liste[ip].v.k[kp] - liste[jp].v.k[kp]); /*   * Relative velocity */
		}
		fdfac =  -(gamma_) *  omega * omega * fdfac;
		
		
		/*-----------------*/
		/* Update forces   */
		/*-----------------*/
		for (kp = 0; kp < n_dim; kp++) {
			/* positive for particle i */
			/* negative for particle j */
			real_t f;
			
			/* Calculate the pairwise conservative */
			f = fcfac * eij.k[kp];
			liste[ip].fc.k[kp] += f; //part_i->fc.k[kp] += f; 
			liste[jp].fc.k[kp] -= f;
			
			/* Calculate the pairwise dissipative force */
			f = fdfac * eij.k[kp];
			liste[ip].fd.k[kp] += f; //part_i->fd.k[kp] += f;
			liste[jp].fd.k[kp] -= f; 
			
			/* Calculate the pairwise random force */
			f = frfac * eij.k[kp];
			liste[ip].fr.k[kp] += f; //part_i->fr.k[kp] += f;
			liste[jp].fr.k[kp] -= f; 
		}
	}
	return (rijsq < skinsq);
}



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Calculate all forces: conservative, dissipative, and random     */
/*     forces.                                                         */
/*     If needed (the int variable update = .true.),                   */
/*     then update the Verlet neighbor table. Otherwise use            */
/*     the present table.                                              */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int fcdr_(int maxnab) {
    static int jpbeg;
    static int jpnab, jpend;
    int nlist, ip, jp, kip;
    int chaincell, nb, nbc;               // indices for cell chaining
	
    for (ip = 0; ip < num_particle; ++ip) {
        /* Initialize all forces */
        liste[ip].fc = fzero;
        liste[ip].fd = fzero;
        liste[ip].fr = fzero;
    }
	
    if (update) {  /* update the Verlet list  */
        fill_cell_chains(num_particle);
        
        for (ip = 0; ip < num_particle; ++ip) {
            /* Ok, let's do it. */
            /* Save the particle positions. */
            liste[ip].r0 = liste[ip].r;   // we use this in check_()
        }
        
        nlist = 0;
        for (ip = 0; ip < num_particle-1; ++ip) {
            if (liste[ip].inside) {
                /* Update the Verlet list and calculate all forces */
                npoint[ip] = nlist;
                
                fpoint_t rr = liste[ip].r;
                ipoint_t ii = izero;
				
                for (kip = 0; kip < n_dim; kip++) {
					rr.k[kip] -= Boundary(rr,kip);
					ii.k[kip] = (int)((int)(rr.k[kip] + sys_size2+skin)/skin);
                } 
				
                for(nb=0; nb < (n_dim==2 ? 9 : 27); nb++) {
                    int iix = (ii.k[0] + nb3Dx[nb]) % skin_size;
                    int iiy = (ii.k[1] + nb3Dy[nb]) % skin_size;
                    int iiz = (n_dim==2)? 0: (ii.k[2] + nb3Dz[nb]) % skin_size;
					
                    nbc=(iiz*skin_size+iiy)*skin_size+iix;
                    
                    // error control
                    if (nbc>=F_ASIZE) {  
						printf("Error GT in indexing  nb %d nbc %d  ip=%d\n",nb,nbc,ip);
						break; 
                    }
                    else if (nbc<0) {  
						printf("Error LT in indexing  nb %d nbc %d  ip=%d\n",nb,nbc,ip);
						break;
                    }
                    chaincell=chainarray[nbc]; 
                    
                    while (chaincell!=-1) {
                        jp=chaincell; 
                        chaincell=liste[chaincell].nextindex;
						
                        if ((jp>ip) && liste[jp].inside) {                                                                                                                                                                                                    
                            /* Go over all pairs within the skin-radius/Verlet list */
							
                            if (calc_ij(ip,jp)) {
                                /* The particles are neighbors */
                                list_verl[nlist] = jp;   /* put into verlet-list */
                                nlist++;
                                if (nlist == maxnab) {
                                    stopit("verlet-list too small");
                                }
                            }
                        } 
                    }  // while
                }  // for(nb...
            }  // if (liste[ip].inside)
        }  // for (ip...
        
        npoint[num_particle - 1] = nlist + 1;
        nvsize = nlist + 1;
	} else {
		/* There is no need to update the Verlet */
		/* so let us use the current list and calc all forces */
		for (ip = 0; ip < num_particle-1; ++ip) {
			if (liste[ip].inside) {
				jpbeg = npoint[ip];
				jpend = npoint[ip+1] - 1;
				
				if (jpbeg <= jpend) {
					for (jpnab = jpbeg; jpnab <= jpend; ++jpnab) {
						jp = list_verl[jpnab];
						
						if (liste[jp].inside) {
							calc_ij(ip,jp);
						}
					}
				}
			}
		}
	}

for (ip = 0; ip < num_particle; ++ip) {
	
	/***********************************************************************************/
	/* Calculate conservative forces due to harmonic springs between adjacent monomers */
	/***********************************************************************************/
	
	if (! liste[ip].inside) continue;
	jp = liste[ip].right;
	if (jp >= 0) { // particle is bonded to something
		fpoint_t rij;
	    float_t distrf = 0.0;
		// float_t spring_dist = 0.01;
		
		for (kip = 0; kip < n_dim; kip++) {
			rij.k[kip] = liste[ip].r.k[kip] - liste[jp].r.k[kip];
			rij.k[kip] -= Boundary(rij,kip);
		}
		for (kip = 0; kip < n_dim; kip++) {
			distrf += rij.k[kip] * rij.k[kip];
		}
		distrf = sqrt(distrf);
		
		if (distrf<0.0000001) { distrf=0.0000001; }   // div/0
		
		int typei = liste[ip].type;  /* Type of particle "ip"  */
        int typej = liste[jp].type;  /* Type of particle "jp"  */
		
		// totally unbreakable bonds are possible by setting break_prob <= 0
		// not currently true, but could be by rearranging parens in line below
#ifdef PROBABALISTIC_BONDING
		if ((BOND_BREAK(typei,typej) > 0 && BOND_BREAK(typei,typej) > r2s_()) || distrf > 1)
#else
			if (BOND_BREAK(typei,typej) > 0 && BOND_BREAK(typei,typej) <= distrf)
#endif
			{
				int unbond_err = dynamically_unbond_particles(ip, jp);
				if (unbond_err > 0) {
					printf("Error unbonding particles %d and %d\n", ip, jp);
				}
				continue;
			}
				
				// stiff bonds
				
				int kp = liste[jp].right;
		if (kp >= 0) { // there is a ip - jp - ip trimer (or possibly longer)
			fpoint_t rjk;
			float_t distij = distrf;
			float_t distjk = 0.0;
			
			for (kip = 0; kip < n_dim; kip++) {
				rjk.k[kip] = liste[jp].r.k[kip] - liste[kp].r.k[kip];
				rjk.k[kip] -= Boundary(rjk,kip);
			}
			for (kip = 0; kip < n_dim; kip++) {
				distjk += rjk.k[kip] * rjk.k[kip];
			}
			distjk = sqrt(distjk);
			
			if (distjk<0.0000001) { distjk=0.0000001; }   // div/0
			
			// thanks to Anders Erikson for the code
			/*
			 i is basically ij
			 i+1 means jk
			 dx[i] => rjk.k[kp]
			 HookeanSpringConstant => d_bond__
			 BondLength => spring_dist
			 dist[i] => distij
			 dist[i+1] => distjk
			 distInv[i] => 1/distij
			 distInv[i+1] => 1/distjk
			 */
			double a = bond_angle_strength * (1/distij) * (1/distjk);
			double b = a * (rij.k[0]*rjk.k[0] + rij.k[1]*rjk.k[1] + rij.k[2]*rjk.k[2]);
			
			double b0 = b * (1/distij) * (1/distij);
			double f0x =  b0 * rij.k[0] - a * rjk.k[0];
			double f0y =  b0 * rij.k[1] - a * rjk.k[1];
			double f0z =  b0 * rij.k[2] - a * rjk.k[2];
			liste[ip].fc.k[0] -= f0x; liste[jp].fc.k[0] += f0x;
			liste[ip].fc.k[1] -= f0y; liste[jp].fc.k[1] += f0y;
			liste[ip].fc.k[2] -= f0z; liste[jp].fc.k[2] += f0z;
			
			double b2 = b * (1/distjk) * (1/distjk);
			double f2x =  a * rij.k[0] - b2 * rjk.k[0];
			double f2y =  a * rij.k[1] - b2 * rjk.k[1];
			double f2z =  a * rij.k[2] - b2 * rjk.k[2];
			liste[kp].fc.k[0] -= f2x; liste[jp].fc.k[0] += f2x;
			liste[kp].fc.k[1] -= f2y; liste[jp].fc.k[1] += f2y;
			liste[kp].fc.k[2] -= f2z; liste[jp].fc.k[2] += f2z;
		}
		
		for (kip = 0; kip < n_dim; kip++) {
			liste[ip].fc.k[kip]  -= d_bond__ * (distrf - spring_dist) * (rij.k[kip]) * (1/distrf);
			liste[jp].fc.k[kip]  += d_bond__ * (distrf - spring_dist) * (rij.k[kip]) * (1/distrf);
			
			//liste[ip].fc.k[kp]  -= d_bond__ * rij.k[kp];
			//liste[jp].fc.k[kp]  += d_bond__ * rij.k[kp];
		}
	}
}
return 0;
} /* fcdr_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Calculate dissipative forces only.                              */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int fdr_()
{
    /* Local variables */
    fpoint_t  ri,rij,vi,fdi,eij,vij;
    real_t  fdfac;
    real_t omega;
    int jpbeg;
    int jpnab, jpend;
    int ip, jp, kp;
    real_t rijsq,rijinv, rrij;
	
    for (ip = 0; ip < num_particle; ++ip) {
        if (! liste[ip].inside) continue;
        liste[ip].fd = fzero;
    }
    
    for (ip = 0; ip < num_particle-1; ++ip) {
        if (! liste[ip].inside) continue;
        jpbeg = npoint[ip];
        jpend = npoint[ip+1] - 1;
		
        if (jpbeg <= jpend) {
            ri  = liste[ip].r;
            vi  = liste[ip].v;
            fdi = liste[ip].fd;
            
            for (jpnab = jpbeg; jpnab <= jpend; ++jpnab) {
                jp = list_verl[jpnab];
                rijsq = 0;
                for (kp = 0; kp < n_dim; kp++) {
                    rij.k[kp] = ri.k[kp] - liste[jp].r.k[kp];
                    rij.k[kp] -= Boundary(rij,kp);
                    rijsq += rij.k[kp]*rij.k[kp]; 
                }
                if (rijsq < 1.) {
                    rrij = sqrt(rijsq);
                    if (rrij < 0.0000001) {
                        rijinv = 1000000;
                    }
                    else {
                        rijinv = 1. / rrij;
                    }
                    omega = 1. - rrij;
					
                    fdfac = 0;
                    for (kp = 0; kp < n_dim; kp++) {
                        eij.k[kp] = rij.k[kp] * rijinv;
                        vij.k[kp] = vi.k[kp] - liste[jp].v.k[kp];
                        fdfac += eij.k[kp] * vij.k[kp];
                    }
                    fdfac *= omega * omega * (-gamma_);
					
                    for (kp = 0; kp < n_dim; kp++) {
                        real_t fdij = fdfac * eij.k[kp];
                        fdi.k[kp]          += fdij;
                        liste[jp].fd.k[kp] -= fdij;
                    }
                }
            }
            liste[ip].fd = fdi;
        }
    }
#if 0
    for (ip = 0; ip < num_particle; ++ip) {
        if (! liste[ip].inside) continue;
        /* Multiply forces with prefactors */
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].fd.k[kp] *= -(gamma_);
        }
    }
#endif
    return 0;
} /* fdr_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Check if the Verlet neighbor table should be updated.           */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int check_()
{
    int ip,kp;
    real_t dispmx = 0;
	
    for (ip = 0; ip < num_particle; ++ip) {
        if (! liste[ip].inside) continue;
        /* Computing MAX */
        for (kp = 0; kp < n_dim; kp++) {
            dispmx = max(fabs(liste[ip].r.k[kp] - liste[ip].r0.k[kp]),dispmx);
        }
    }
    dispmx = dispmx * sqrt(3.0) * 2.0;
    update = dispmx > skin - 1.;
    return 0;
} /* check_ */



/*******************************************************************/

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 20.07.2004  11:59h
  $Id: dpd_stat.c,v 1.6 2004/07/30 09:13:47 tmaeke Exp $
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
#include "dpd_stat.h"


/*Defines***********************************************************/
#define VELO_ANZ 50
#define VACF_SIZE    100

/*Types*************************************************************/

/*Vars**************************************************************/
    // velocity profile:    
    real_t velo[VELO_ANZ];
    int veloanz[VELO_ANZ];
    int drawinfoupd = 0;

    fpoint_t vel[VACF_SIZE];         // For the velocity autocorrelation function
    real_t dcs_vacf__[VACF_SIZE];
    real_t dvacf[VACF_SIZE];
    int    it_bound__[VACF_SIZE+1];
    int    it_vacf__;
    int    ivacf_loop__[2*VACF_SIZE+1];
    int    ivacf_counter__;

    int len_chain__ = 0; // dummy

/*Funcs*************************************************************/



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void init_vacf(void)
{    
    int it;
    for (it = 0; it < VACF_SIZE; ++it) {      /* Here we initialize a few variables */
        dvacf[it] = 0.;                       /* to examine the decay of the */
    }                                         /* velocity correlation function of */
                                              /* the polymer. This is related to */
    for (it = 1; it <= VACF_SIZE; ++it) {     /* the tracer diffusion coefficient */
        ivacf_loop__[it + VACF_SIZE] = it;    /* of the polymer motion. */
    }
    
    for (it = -VACF_SIZE; it <= 0; ++it) {
        ivacf_loop__[it + VACF_SIZE] = it + VACF_SIZE;
    }
    it_vacf__ = 0;
    ivacf_counter__ = 0;
    for (it = 1; it <= VACF_SIZE; ++it) {
        it_bound__[it - 1] = it;
    }
    it_bound__[VACF_SIZE] = 1;
}    


void end_vacf( int incstp)
{
    int it;
    /* -------------------------------------- */
    if (ivacf_counter__ >= VACF_SIZE) {
        /* -------------------------------------- */
        /* The purpose here is to calculate the */
        /* velocity autocorrelation function and */
        /* then write the final averaged result */
        for (it = 0; it < VACF_SIZE; ++it) {
            dvacf[it] = dvacf[it] / (real_t) (ivacf_counter__ - (VACF_SIZE-1));
        }

        dcs_vacf__[0] = dvacf[0] /6.0 * dt * (real_t) incstp;
    
        for (it = 1; it < VACF_SIZE; ++it) {
            dcs_vacf__[it] = dcs_vacf__[it - 1] + dvacf[it] /3.0 * dt * (real_t) incstp;
        }
        for (it = 0; it < VACF_SIZE; ++it) {
            /* Write the results */ /*io*/
        }
    }
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Calculate the center-of-mass velocity of the system */
/*     (to guarantee that momentum conservation is satisfied), */
/*     and the instantaneous kinetic temperature of the system */
/*     based on the velocoties of the particles. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int velo_prof(        
    void)
{
    /* Local variables */
    int ip,i,kp;
    fpoint_t vi;
    
    for (i=0; i<VELO_ANZ; i++) { velo[i] = 0; veloanz[i] = 0; }
    
    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            vi.k[kp] = liste[ip].v.k[kp];
        }
        
        i = (int) ((liste[ip].r.k[1]+dl.k[1]/2.0)/dl.k[1]*VELO_ANZ);
        if (i>=0 && i<VELO_ANZ) {
            velo[i] += vi.k[0];
            veloanz[i]++;
        }
    }
    if (drawinfo) {
        for (i=0; i<VELO_ANZ; i++) { 
            if (veloanz[i]>0) printf(" %2d:(%3d) = %7.4f\n",i,veloanz[i],velo[i]/veloanz[i]); 
        }
    }
    drawinfoupd = 1;
    return 0;
} /* velo_prof */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Calculate the center-of-mass velocity of the system */
/*     (to guarantee that momentum conservation is satisfied), */
/*     and the instantaneous kinetic temperature of the system */
/*     based on the velocoties of the particles. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int vcmtp_(        
    real_t *pm, 
    real_t *pminv, 
    real_t *vcm, 
    real_t *temp)
{
    /* Local variables */
    static int ip,kp;
    fpoint_t vi, vcml = {{0,0,0}};
    real_t tmp;

    *temp = 0.0;
    for (ip = 1; ip <= num_particle; ++ip) {
        tmp = 0.0;
        for (kp=0; kp<n_dim; kp++) {
           vi.k[kp] = liste[ip - 1].v.k[kp];
           vcml.k[kp]  += *pm * vi.k[kp];
           tmp += vi.k[kp]*vi.k[kp];
        }
        *temp += *pm * tmp;
    }
    tmp = 0.0;
    for (kp=0; kp<n_dim; kp++) {
        tmp += vcml.k[kp]*vcml.k[kp];
    }
    *vcm = (1/((real_t)num_particle)) * *pminv * sqrt(tmp);
    *temp = tfac * *temp;
    return 0;
} /* vcmtp_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Calculate */
/*     - the radius of gyration squared */
/*     - the end-to-end distance squared */
/*     of the polymer chain. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
 int r_gyr__(
    int n_solu, 
    real_t *dr_gyr__, 
    real_t *dr_ee__)
{
    /* Local variables */
    int lenc, j,kp;
    fpoint_t dr_g_ = {{0,0,0}}, 
             dr_ee_, 
             dr_cm_ = {{0,0,0}}, 
             dr_tm_, 
             disp_, 
             dr_dpl_[len_chain__];

    for (kp=0; kp<n_dim; kp++) {
        dr_dpl_[0].k[kp] = 0.;
    }
    
    /* Determine the positions of the */
    /* monomers in a chain with respect */
    /* to the first monomer */
    for (j = 2; j <= len_chain__; ++j) {
        for (kp=0; kp < n_dim; kp++) {
            disp_.k[kp] = liste[n_solu + j - 1].r.k[kp] - liste[n_solu + j - 2].r.k[kp];
            disp_.k[kp] -= Boundary(disp_,kp);
            dr_dpl_[j - 1].k[kp] = dr_dpl_[j - 2].k[kp] + disp_.k[kp];
        }
    }
    
    /* Then calculate the radius of gyration */
    /* over the whole chain */
    for (j = 1; j <= len_chain__; ++j) {
        for (kp=0; kp<n_dim; kp++) {
            dr_cm_.k[kp] += dr_dpl_[j - 1].k[kp];
        }
    }
    for (kp=0; kp<n_dim; kp++) {
        dr_cm_.k[kp] /= len_chain__;
    }
    for (j = 1; j <= len_chain__; ++j) {
        for (kp=0; kp<n_dim; kp++) {
            dr_tm_.k[kp] = dr_dpl_[j - 1].k[kp] - dr_cm_.k[kp];
            dr_tm_.k[kp] -= Boundary(dr_tm_,kp);
            dr_g_.k[kp] += dr_tm_.k[kp] * dr_tm_.k[kp];
        }
    }
    for (kp=0; kp<n_dim; kp++) {
        dr_g_.k[kp] /= len_chain__;
        *dr_gyr__ += dr_g_.k[kp];
    }
    
    /* Then calculate the squared */
    /* Radius of gyration squared (whole chai */
    /* end-to-end distance */
    lenc = len_chain__;
    for (kp=0; kp<n_dim; kp++) {
        dr_ee_.k[kp] = dr_dpl_[lenc - 1].k[kp];
        dr_ee_.k[kp] *= dr_ee_.k[kp];
        *dr_ee__ += dr_ee_.k[kp];
    }
    /* Squared end-to-end distance */
    return 0;
} /* r_gyr__ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Take a snapshot of the current configuration of */
/*     the polymer chain. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int snaps_(
    int n_solu, 
    int nconf, 
    int ifile)
{
    /* Local variables */
    int j, ip,kp;
    fpoint_t disp_, dr_dpl_[len_chain__];

    /* The starting position is the position */
    /* of one of the end monomers */
    ip = n_solu + 1;
    for (kp=0; kp<n_dim; kp++) {
        dr_dpl_[0].k[kp] = liste[ip - 1].r.k[kp];
    }

    for (j = 1; j <= len_chain__ - 1; ++j) {
        for (kp=0; kp<n_dim; kp++) {
            disp_.k[kp] = liste[ip + j - 1].r.k[kp] - liste[ip + j - 2].r.k[kp];
            disp_.k[kp] -= Boundary(disp_,kp);

            /* The positions of other monomers are */
            /* written with respect to the starting */
            /* position. Thus it is possible that */
            /* some of the monomers reside outside */
            /* the simulation box. */
            
            dr_dpl_[j].k[kp] = dr_dpl_[j - 1].k[kp] + disp_.k[kp];
        }
    }
    /*     Here we write the (x,y)-coordinates of the monomers */
    /*     of the polymer chain. The files will be called "fort.abc" */
    /*     where "abc" is a number ranging from Ifile to (Ifile + 49). */
    /*     The parameter Ifile is defined to have a value of 50, thus */
    /*     the maximum number of snapshots is also 50. */

    /*     The if-loop below is not needed on most computer architectures. */
    /*     However, on the Linux-workstations used during the summer */
    /*     school, the Fortran77 compilers (f77 and g77) limit the */
    /*     maximum unit number to be 99 (fort.99 is ok, while fort.100 */
    /*     leads to an error message and core dumped). */
    if (nconf < 50) {
        /*
        i__1 = *len_chain__;
        for (j = 1; j <= i__1; ++j) {
              do_fio(&c__1, (char *)&drx_dpl__[j])
              do_fio(&c__1, (char *)&dry_dpl__[j])
              do_fio(&c__1, (char *)&drz_dpl__[j])
        }
        f_clos(&cl__1);
        */
    }
    return 0;
} /* snaps_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Update the velocity autocorrelation function for the */
/*     study of the tracer diffusion coefficient to characterize */
/*     the motion of the polymer chain along the 2D plane. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int vacf_(
    int n_solu, 
    int maxt_vacf__, 
    int *it_vacf__, 
    int *ivacf_counter__, 
    int *it_bound__, 
    int *ivacf_loop__, 
    real_t *dvacf)
{
    /* System generated locals */
    int ivacf_loop_offset;

    /* Local variables */
    static int ihist, ip, kp, it_diff__;

    /* Parameter adjustments */
    ivacf_loop_offset = -(maxt_vacf__);

    /* Rotate the label of the register */
    *it_vacf__ = it_bound__[*it_vacf__];
    
    for (kp=0; kp<n_dim; kp++) {
        vel[*it_vacf__ - 1].k[kp] = 0.;
    }
    
    /* Calculate the velocity of the */
    /* center-of-mass of the polymer chain */
    for (ip = n_solu + 1; ip <= num_particle; ++ip) {
        for (kp=0; kp<n_dim; kp++) {
            vel[*it_vacf__ - 1].k[kp] = vel[*it_vacf__ - 1].k[kp] + liste[ip - 1].v.k[kp];
        }
    }
    
    for (kp=0; kp<n_dim; kp++) {
        vel[*it_vacf__ - 1].k[kp] = vel[*it_vacf__ - 1].k[kp] / len_chain__;
    }
    ++(*ivacf_counter__);
    
    if (*ivacf_counter__ >= maxt_vacf__) {
    
        for (it_diff__ = 0; it_diff__ <= maxt_vacf__ - 1; ++it_diff__) {
            
            ihist = ivacf_loop__[*it_vacf__ - it_diff__ - ivacf_loop_offset];
            
            /* Update the velocity autocorrelation fu */
            for (kp=0; kp<n_dim; kp++) {
                dvacf[it_diff__] += 
                    vel[*it_vacf__ - 1].k[kp] * vel[ihist - 1].k[kp];
            }
        }
    }
    return 0;
} /* vacf_ */


/*******************************************************************/

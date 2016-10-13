/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 22.07.2004  11:04h
  $Id: dpd_elec.c,v 1.9 2004/09/07 07:35:26 tmaeke Exp $
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
#include "dpd_geom.h"
#include "dpd_elec.h"
#include "dpd_disp.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/

/*Funcs*************************************************************/


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void update_electrodes() 
{
#ifdef OLD_ELEC
    int el;
    for (el=0; el < MAX_ELECTRODES; el++) {  
        if (electrodes[el].seq_ofs >= 0) {
            electrodes[el].polarity = e_seq.list[(electrodes[el].seq_ofs + e_seq.sequence) % e_seq.length];
        }
    }
#else
    elem_p p = geom_electrodes;
    while (p) {
        if (p->typ==electrode && p->load) {
            electrode_p elec = (electrode_p)p->load;
            if (elec->active==1  && elec->seq_ofs >= 0) {
                elec->polarity = e_seq.list[(elec->seq_ofs + e_seq.sequence) % e_seq.length];
            }
            int ci = 6;
            if (polarity * elec->polarity < 0) { ci=8; }
            else if (polarity * elec->polarity > 0) { ci=7; }
            polygon_set_color(elec->poly, ci);
        }
        p = p->next;
        if (p == geom_electrodes) break;
    }
#endif    
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void change_polar(int wie) 
{
    if (wie) {
        polarity = 1.0;
        e_seq.sequence = (e_seq.sequence+1) % e_seq.length;
    }
    else {
        polarity = -1.0;
        e_seq.sequence = (e_seq.sequence-1+e_seq.length) % e_seq.length;
    }
    update_electrodes();
    printf("   rfac = %7.4f,%d\n", rfac, e_seq.sequence);
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void change_efield(int wie) 
{
    if (wie) {
        epsilonfac = epsilonfac*2;
    }
    else {
        epsilonfac = epsilonfac/2;
    }
    printf("   efac = %7.4f\n", epsilonfac);
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void draw_elecs()
{        
#ifdef OLD_ELEC
    int jj;
    for (jj=0; jj < MAX_ELECTRODES; jj++) {
        if (electrodes[jj].active) {
            int c1 = 6;
            if (polarity * electrodes[jj].polarity < 0) { c1=8; }
            else if (polarity * electrodes[jj].polarity > 0) { c1=7; }
            drawRect_ofs(screen, electrodes[jj].position.k[0]-electrodes[jj].size.k[0], 
                                electrodes[jj].position.k[1]-electrodes[jj].size.k[1], 
                                electrodes[jj].position.k[0]+electrodes[jj].size.k[0], 
                                electrodes[jj].position.k[1]+electrodes[jj].size.k[1],
                                Color[c1], 150);
        }
    }
#else
    elem_list_draw(geom_electrodes, polygon);  
#endif
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* calculate forces from e-field                                       */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void e_field()   // the complex version
{
    /*  TODO:
          1. actualize "evector" for each (inside/allowed) unitcube
          2. calc force on charged particles within unitcubes
    */
}



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* calculate forces from e-field                                       */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void calc_efield()   // the simple version
{        
#ifdef OLD_ELEC
    int el, jp, kp;
    for (el=0; el < MAX_ELECTRODES; el++) {  
        if (electrodes[el].active==1  && (electrodes[el].polarity != 0.0)) {
            real_t pol = polarity * electrodes[el].polarity * epsilonfac;

            for (jp = 0; jp < num_particle; ++jp) {
            
                int typ = liste[jp].type;
                if (particle[typ].epsilon != 0.0) {

                    real_t rijsq = 0;      // calc distance from i to electrode
                    fpoint_t rij;
                    for (kp = 0; kp < n_dim; kp++) {
                        rij.k[kp] = liste[jp].r.k[kp] - electrodes[el].position.k[kp];
                        rij.k[kp] -= Boundary(rij,kp);
                        rijsq += rij.k[kp]*rij.k[kp];
                    }
                    real_t rrij = sqrt(rijsq);
                    if (rrij > emin) {
                    
                        real_t rijinv = 1. / (rrij*rrij);  // calc resulting force
                        real_t omega = pol * particle[typ].epsilon * rijinv;
                      
                        for (kp = 0; kp < n_dim; kp++) {  
                            real_t e = rij.k[kp]; // * rijinv;
                            liste[jp].fc.k[kp] += e * omega;
                        }
                    }
                }
            }
        }
    }
#else
    int jp, kp;
    
    elem_p p = geom_electrodes;
    while (p) {
        if (p->typ==electrode && p->load) {
            electrode_p elec = (electrode_p)p->load;
            if (elec->active==1  && (elec->polarity != 0.0)) {
                real_t pol = polarity * elec->polarity * epsilonfac;

                for (jp = 0; jp < num_particle; ++jp) {

                    real_t eps = particle[liste[jp].type].epsilon;
                    if (eps != 0.0) {

                        real_t rijsq = 0;      // calc distance from i to electrode
                        fpoint_t rij;
                        for (kp = 0; kp < n_dim; kp++) {
                            rij.k[kp] = liste[jp].r.k[kp] - elec->position.k[kp];
                            rij.k[kp] -= Boundary(rij,kp);
                            rijsq += rij.k[kp]*rij.k[kp];
                        }
                        real_t rrij = sqrt(rijsq);
                        if (rrij > emin) {

                            real_t rijinv = 1. / (rrij*rrij);  // calc resulting force
                            real_t omega = pol * eps * rijinv;

                            for (kp = 0; kp < n_dim; kp++) {  
                                real_t e = rij.k[kp]; // * rijinv;
                                liste[jp].fc.k[kp] += e * omega;
                            }
                        }
                    }
                }
            }
        }
        p = p->next;
        if (p == geom_electrodes) break;
    }
#endif    
}


/* ****************************************************** */
/* ****************************************************** */
void elec_translate()
{
    /* convert from elem-list to electrodes-array
    */
    num_electrodes = 0;
    elem_p p = geom_electrodes;
    int n;
    for (n=0; n<MAX_ELECTRODES;n++) {
        electrodes[n].active = 0;
    }
    while (p) {
        if (p->load && p->typ==polygon) {
            fpoint_t min = {{ 1E10, 1E10, 1E10}};
            fpoint_t max = {{-1E10,-1E10,-1E10}};
            polygon_min_max((poly_p)p->load, &min, &max);
            electrodes[num_electrodes].active = 1;
            electrodes[num_electrodes].polarity = 1.0;
            electrodes[num_electrodes].dutycycle = 100;
            electrodes[num_electrodes].seq_ofs = num_electrodes;

            int kp;
            for (kp=0; kp<3; kp++) {
                electrodes[num_electrodes].position.k[kp] = (max.k[kp]+min.k[kp])/2;
                electrodes[num_electrodes].size.k[kp] = (max.k[kp]-min.k[kp])/2;
            }
        }
        p = p->next;
        if (p == geom_electrodes) break;
        num_electrodes++;
        if (num_electrodes==MAX_ELECTRODES) break;
    }
}


/* ****************************************************** */
/* ****************************************************** */
void elec_add_electrode(
    char * tag,
    int seqno,
    real_t defpol)
{
#ifdef OLD_ELEC
    elem_p p = elem_list_findname(geom_electrodes,tag);
    if (p && p->typ==polygon && (num_electrodes<MAX_ELECTRODES)) {
        if (p->load) {
            fpoint_t min = {{ 1E10, 1E10, 1E10}};
            fpoint_t max = {{-1E10,-1E10,-1E10}};
            polygon_min_max((poly_p)p->load, &min, &max);
            electrodes[num_electrodes].active = 1;
            electrodes[num_electrodes].polarity = defpol;
            electrodes[num_electrodes].dutycycle = 100;
            electrodes[num_electrodes].seq_ofs = seqno;

            int kp;
            for (kp=0; kp<3; kp++) {
                electrodes[num_electrodes].position.k[kp] = (max.k[kp]+min.k[kp])/2;
                electrodes[num_electrodes].size.k[kp] = (max.k[kp]-min.k[kp])/2;
            }
            fprintf(stderr,"new elec: %d  seq=%d  pol=%3.0f\n", num_electrodes, seqno, defpol); 
            num_electrodes++;
        }
    }
#else
    elem_p p = elem_list_findname_qual(geom_electrodes,tag,polygon);
    if (p && p->typ==polygon) {
        if (p->load) {
            elem_p e = elec_new(tag, geom_electrodes, seqno,defpol,100,1);
            elec_add_polygon(e,(poly_p)p->load);
            fprintf(stderr,"elec_new: %d  seq=%d  pol=%3.0f\n", num_electrodes, seqno, defpol); 
            num_electrodes++;
        }
    }
#endif    
}    



/*******************************************************************/

/*******************************************************************
  tmm 03.08.2004  10:13h
  $Id: dpd_flow.c,v 1.3 2004/08/24 11:33:39 tmaeke Exp ${fn},v 1.10 2003/08/11 09:56:01 tmaeke Exp tmaeke $
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
#include "dpd_geom.h"
#include "dpd_flow.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/
    
    bool_t   inout_flow = 0;             // enables all
    index_t  outside_chain = -1;         // start of linked list of outsiders 
    int      inflows = 0;                // actual no of ins
    int      outflows = 0;               // actual no of outs
    ioflow_t inflow_list[MAX_IO_FLOW];   
    ioflow_t outflow_list[MAX_IO_FLOW];
    
    

/*Funcs*************************************************************/


/*******************************************************************/
/*******************************************************************/
void add_inflow(
    char * tag, 
    int i, 
    real_t rate)
{
    if (inflows < MAX_IO_FLOW) {
        elem_p p = elem_list_findname(geom_sensors,tag);
        if (p && p->typ==polygon && p->load) {
            fpoint_t min = {{ 1E10, 1E10, 1E10}};
            fpoint_t max = {{-1E10,-1E10,-1E10}};
            polygon_min_max((poly_p)p->load, &min, &max);
            inflow_list[inflows].mi = min;
            inflow_list[inflows].ma = max;
            int k;
            for (k=0; k<n_dim; k++) {
                inflow_list[inflows].mid.k[k] = (max.k[k]+min.k[k])/2.0; 
            }
            inflow_list[inflows].rate = (int) rate;
            inflows++;
            printf("Added Inflow %d %s\n", inflows,tag);
         }
    }
}    


/*******************************************************************/
/*******************************************************************/
void add_outflow(
    char * tag, 
    int i, 
    real_t rate)
{
    if (outflows < MAX_IO_FLOW) {
        elem_p p = elem_list_findname(geom_sensors,tag);
        if (p && p->typ==polygon && p->load) {
            fpoint_t min = {{ 1E10, 1E10, 1E10}};
            fpoint_t max = {{-1E10,-1E10,-1E10}};
            polygon_min_max((poly_p)p->load, &min, &max);
            outflow_list[outflows].mi = min;
            outflow_list[outflows].ma = max;
            outflow_list[outflows].rate = (int) rate;
            outflows++;
            printf("Added Outflow %d %s\n", outflows,tag);
         }
    }
}    

/*******************************************************************/
/*******************************************************************/
bool_t test_outflow(
    fpoint_t * where)
{
    int f;
    for (f=0; f<outflows; f++) {
        int kp,in=1;
        for (kp=0; kp<2 /*Z always   n_dim*/; kp++) {
            if (where->k[kp] > outflow_list[f].ma.k[kp]) { in=0; break; }
            if (where->k[kp] < outflow_list[f].mi.k[kp]) { in=0; break; }
        }
        // this point is inside the outflow volume
        if (in) return 1;
    }
    return 0;
}

/*******************************************************************/
void flow_out_one(
    index_t ip)
{
    if (!liste[ip].inside || (ip<0)) return;
    
    liste[ip].inside = 0;
    liste[ip].next_outside = outside_chain;
    outside_chain = ip;
    outsiders ++;
    
    // and now the rest of a polymer-chain
    flow_out_one(liste[ip].right);
}    

/*******************************************************************/
void flow_out(index_t ip)
{
	// BUGBUG this won't be loop safe
    while (liste[ip].left>=0) { ip = liste[ip].left; }  // go to left end of polymer
    flow_out_one(ip);
    //ow_out_one(liste[ip].right, -1);
}    


/*******************************************************************/
/*******************************************************************/
int flow_in(
    int inflow)
{
    if (outsiders > 0 && (outside_chain >= 0) && (inflow<inflows)) {
        int ip, c=0;
        do {  // push all part of a polymer chain inside
            ip = outside_chain;
            outside_chain = liste[ip].next_outside;
            outsiders--;
            liste[ip].inside = 1;
            liste[ip].next_outside = -1;
            liste[ip].r = inflow_list[inflow].mid;
            c++;
        } while ( (outsiders>0) || (liste[ip].left>=0)  );
        return c;
    }
    return 0;
}

/*******************************************************************/
/*******************************************************************/
bool_t flow_in_all()
{
    int f;
    //printf(" in rate ");
    for (f=0; f<inflows; f++) {
        int rate=0;
        while (rate < inflow_list[f].rate) {
            int r = flow_in(f);
            //if (r==0) return 0; // no more left
            if (r==0) break;
            rate += r;
        }
        //printf(" (%d) = %d ", f, rate); 
    }
    //printf("\n");
    return 1;
}


/*******************************************************************/


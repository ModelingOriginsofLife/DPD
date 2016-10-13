/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
tmm 20.07.2004  08:10h
$Id: dpd_file.c,v 1.23 2004/09/29 12:38:27 tmaeke Exp $
********************************************************************/

/*******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <string.h> 

#include "dpd_vars.h"
#include "dpd_disp.h"
#include "dpd_util.h"
#include "dpd_init.h"
#include "dpd_file.h"
#include "dpd_geom.h"
#include "dpd_elec.h"
#include "dpd_flow.h"
#include "dpd_dynamic.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/
const char delimiters[] = " ,;:!()\t\n\r";
static int err;

char * errors[] = {
	"no Error",             //  0
	"missing '='",          // -1
	"wrong int",            // -2 
	"out of range",         // -3
	"unknown token",        // -4
	"wrong float",          // -5
	"wrong float-triple",   // -6
	"missing float-triple", // -7
	"unknown polytyp",      // -8
};

/*Funcs*************************************************************/


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  note: strtok is not reentrant !!!!                                 */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
bool_t testequal() 
{
    char * tok = strtok (NULL, delimiters); 
    if (tok && !strcmp(tok,"=")) {
        err = 0;
        return 1;
    }
    err = -1;  // missing "="
    return 0;
}

bool_t testplus() 
{
    char * tok = strtok (NULL, delimiters); 
    if (tok && !strcmp(tok,"+")) {
        err = 0;
        return 1;
    }
    err = -1;  // missing "+"
    return 0;
}

bool_t testarrow() 
{
    char * tok = strtok (NULL, delimiters); 
    if (tok && !strcmp(tok,"->")) {
        err = 0;
        return 1;
    }
    err = -1;  // missing "->"
    return 0;
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*
bool_t testarrow() 
{
    char * tok = strtok (NULL, delimiters); 
    if (tok && !strcmp(tok,"->")) {
        err = 0;
        return 1;
    }
    err = -1;  // missing "->"
    return 0;
}
*/
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
bool_t get1flt(real_t *r)
{
    real_t flo;
    char * tok = strtok (NULL, delimiters); 
    if (tok && sscanf(tok, "%f", &flo)==1) {
        // printf("Float=%f\n",flo);
        *r = flo;
        err = 0;
        return 1;
    }
    err = -5;  // float error
    return 0;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
bool_t get3flt(real_t *x, real_t *y, real_t *z)
{
    real_t flo1;
    char * tok = strtok (NULL, delimiters); 
    err = -6;
    if (tok && sscanf(tok, "%f)", &flo1)==1) {
        //printf("tok=<%s>\n",tok);
        // printf("Float=%f\n",flo1);
        *x=flo1; 
        tok = strtok (NULL, delimiters); 
        if (tok && sscanf(tok, "%f)", &flo1)==1) {
            //printf("tok=<%s>\n",tok);
            *y=flo1; 
            tok = strtok (NULL, delimiters); 
            if (tok && sscanf(tok, "%f)", &flo1)==1) {
                //printf("tok=<%s>\n",tok);
                *z=flo1;
                err = 0;
            }
        }
        return 1;
    }
    return 0;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
bool_t get1int(int_t *r)
{
    int_t i;
    char * tok = strtok (NULL, delimiters); 
    if (tok && sscanf(tok, "%d", &i)==1) {
        // printf("Int=%d\n",i);
        *r = i;
        err = 0;
        return 1;
    }
    err = -2;  // int error
    return 0;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
bool_t get1str(char *r)
{
    char * tok = strtok (NULL, delimiters); 
    strncpy (r,tok,MAX_NAME);
    err = 0;
    return 1;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int interpret_line(char * ln)
{    
    int i=0;
    err=0;
    // int mode=1; // 0=std 1=uwe
    char line[MAX_FILE_LN];
    strcpy(line,ln);
    
    for (i=0; i<(int)strlen(line); i++) {
        if (line[i] == '#') { line[i] = '\0'; break; }
    }
    while (strlen(line)>0  && (
							   (line[strlen(line)-1]=='\r') ||
							   (line[strlen(line)-1]=='\n') ||
							   (line[strlen(line)-1]==' ')
							   )) {
        line[strlen(line)-1] = '\0';
    }
    
    if (strlen(line) > 0) {
        char *token;
		real_t flo=0,flo2=0,flo3=0;
		int ai=0,aj=0,e1=0,e2=0,ntag=0;
		
        //printf("Exam:<%s>\n",line);
        token = strtok (line, delimiters);
        // printf("token=%s\n", token);
        
        if (token != NULL && strlen(token)>0) {
			if (!strcmp(token,"sigma"))				 { testequal() && get1flt(&sigma); }
            else if (!strcmp(token,"randfac"))       { testequal() && get1flt(&rfac); }
            else if (!strcmp(token,"rho"))           { testequal() && get1flt(&rho0); }
            else if (!strcmp(token,"dt"))            { testequal() && get1flt(&dt); }
            else if (!strcmp(token,"displayintv"))   { testequal() && get1int(&display_intv); }
            else if (!strcmp(token,"seed"))          { testequal() && get1int(&iseed); }
            else if (!strcmp(token,"size"))          { testequal() && get1int(&sys_size); }
            else if (!strcmp(token,"depth"))          { testequal() && get1int(&sys_depth); }
            else if (!strcmp(token,"dim"))           { testequal() && get1int(&n_dim); }
			// geometry
			// slice_disp
            else if (!strcmp(token,"max_step"))      { testequal() && get1int(&maxstp); }
            else if (!strcmp(token,"stat_intv"))     { testequal() && get1int(&incstp); }
            else if (!strcmp(token,"save_intv"))     { testequal() && get1int(&issstp); }
			else if (!strcmp(token,"poly_intv"))     { testequal() && get1int(&poly_intv); }
			else if (!strcmp(token,"reax_intv"))     { testequal() && get1int(&reax_intv); }
            else if (!strcmp(token,"snapshotfn"))     { testequal() && get1str(snapsfn); }            
			// snaps_intv
			// brick
			// pressx
			// colors
            else if (!strcmp(token,"beta_min"))      { testequal() && get1flt(&beta_min); }
            else if (!strcmp(token,"beta_min_force")){ testequal() && get1flt(&beta_min_force); }
            else if (!strcmp(token,"bond"))          { testequal() && get1flt(&d_bond__); }
            else if (!strcmp(token,"spring_dist"))    { testequal() && get1flt(&spring_dist); }
			else if (!strcmp(token,"bond_angle_strength")) { testequal() && get1flt(&bond_angle_strength); }  // for stiff bonds
			/* else if (!strcmp(token,"bond_angle")) { testequal() && get1flt(&bond_angle); } not implemented yet */
			else if (!strcmp(token,"dyn_loops"))   { testequal() && get1int(&dyn_loops); }
			else if (!strcmp(token,"max_dynamic_polymer_length")) { testequal() && get1int(&max_dyn_polylen); }
			else if (!strcmp(token,"initial_position_file")) { testequal() && get1str(initial_position_file); }
			// else if (!strcmp(token,"chem_dist")) { testequal() && get1flt(&chem_dist); }
            else if (!strcmp(token,"serverhost") || !strcmp(token,"serverport")) {
				fprintf(stderr, "server settings no longer supported\n");
			}
			//-------------------------------------------------------------------------------------------
/*			else if (!strcmp(token,"bond_chem")) {
				int tmp1=-1, tmp2=-1;
				get1int(&tmp1) && testarrow() && get1int(&tmp2);
*/				

			//-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"geometry"))      { 
                if (testequal() && get1int(&e1)) {
                    geometry_type= (geometry_t)e1; 
                }
                else err = -3;
            }   
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"slice_disp"))         {
                if (testequal() && get1flt(&min_slice) 
					&& get1flt(&delta_slice)
					&& get1int(&e1)) {
                    slice_mode=e1;
                }
                else {
                    slice_mode=0;
                }
                delta_slice = fabs(delta_slice);
                if (delta_slice < 0.001) delta_slice=0.001;
                max_slice = min_slice+delta_slice;
            }			
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"snaps_intv")) { 
                if (testequal()) {
                    if (!get1int(&snaps_stp)) snaps_end = 0;
                    if (!get1int(&snaps_start)) snaps_start = 0;
                    if (!get1int(&snaps_end)) snaps_end = 0;
                }
            }
			//-------------------------------------------------------------------------------------------			
            else if (!strcmp(token,"brick"))         {
                if (testequal()   && get1flt(&brick_ul.k[0])   
					&& get1flt(&brick_ul.k[1])
					&& get1flt(&brick_or.k[0])     
					&& get1flt(&brick_or.k[1])    
					&& get1int(&e1)) {
                    with_brick=e1;
                    if (get1int(&e2)) {
                        init_with_brick = e2;
                    }
                } 
                else err = -3;
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"pressx"))         {
                if (testequal() && get1flt(&pressure_to_x_minpos) && get1flt(&pressure_to_x_maxpos)
					&& get1flt(&pressure_ofs)  && get1flt(&pressure_diff_at_y) && get1int(&e1)) {
                    pressure_to_x=e1;
                }
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"color")) {
                if (get1int(&ai) && testequal() && get1int(&aj) && get1int(&e1) && get1int(&e2) ) {
                    if (ai<MAX_TYP && ai>=MIN_TYP) {
                        particle[ai].color[0] = (int)aj;
                        particle[ai].color[1] = (int)e1;
                        particle[ai].color[2] = (int)e2;
                    }
                    else err = -3;
                }
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"epsilon")) {
                if (get1int(&ai) && testequal() && get1flt(&flo)) {
                    if (ai<MAX_TYP && ai>=MIN_TYP) {
                        particle[ai].epsilon = flo;
                    }
                    else err = -3;
                }
            }            
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"beta")) {
                if (get1int(&ai) && get1int(&aj) && testequal() &&
					get1flt(&flo) && get1flt(&flo2) && get1int(&e1) && get1int(&e2) ) {
                    if (ai<MAX_TYP && ai>=MIN_TYP && aj<MAX_TYP && aj>=MIN_TYP) {
                        ALPHA(ai,aj) = flo;
                        ALPHA(aj,ai) = flo;
                        BETA(ai,aj) = flo2;
                        BETA(aj,ai) = flo2;
                        BETA12(ai,aj) = e1;
                        BETA12(aj,ai) = e1;
                        BETA6(ai,aj) = e2;
                        BETA6(aj,ai) = e2;
                    }
                    else err = -3;
                }
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"num_particle")) {
                if (get1int(&ai) && testequal() && get1flt(&flo)) {
                    if (ai<MAX_TYP && ai>=MIN_TYP) {
                        particle[ai].percent = flo;
                    }
                    else err = -3;
                }
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"num_chain")) {
                if (get1int(&ai) && testequal() && get1flt(&flo) ) {
                    if (ai<MAX_POLYMERS && ai>=0) {
                        int j;
                        polymere[ai].percent = flo;
                        polymere[ai].len = 0;
						
                        for (j=0; j<MAX_POLYLEN; j++) polymere[ai].part[j]=NO_TYPE;
                        //printf("Line <%s>\nn=%d i=%d\n", line, n,i);
                        j=0;
                        while (get1int(&aj) && j<MAX_POLYLEN) {
                            polymere[ai].part[j] = aj;
                            if (polymere[ai].part[j]<=NO_TYPE) {
                                polymere[ai].part[j]=NO_TYPE;
                            }
                            else {
                                polymere[ai].len ++;
                            }
                            j++;
                        } 
                        err = 0;
                    } 
                    else err = -3; // range error
                }
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"dyn_bond")) {
                int type1=0, type2=0;
				real_t bond_form=-1, bond_break=-1;
                if (get1int(&type1) && get1int(&type2) && testequal() && get1flt(&bond_form) && get1flt(&bond_break)) {
					if (type1<MAX_TYP && type2<MAX_TYP && type1>=MIN_TYP && type2>=MIN_TYP) {
						BOND_FORM(type1,type2) = bond_form;
						BOND_FORM(type2,type1) = bond_form;
						BOND_BREAK(type1,type2) = bond_break;
						BOND_BREAK(type2,type1) = bond_break;
					}
					else {
						err = -3;
					}
				}
			}
			//-------------------------------------------------------------------------------------------
			else if (!strcmp(token,"simple_chem")) {
				int react1=-1, react2=-1, prod1=-1, prod2=-1;
				real_t chance=0;
				
				if (get1int(&react1) && testplus() && get1int(&react2) && testarrow() && get1int(&prod1) && testplus() && get1int(&prod2) && get1flt(&chance)) {
					if (react1<MAX_TYP && react2<MAX_TYP && prod1<MAX_TYP && prod2<MAX_TYP && react1>=MIN_TYP && react2>=MIN_TYP && prod1>=MIN_TYP && prod2>=MIN_TYP) {
						chem_reax[react1][react2].prod1 = prod1;
						chem_reax[react1][react2].prod2 = prod2;
						chem_reax[react1][react2].prob = chance;
						chem_reax[react2][react1].prod1 = prod2; // note that products are switched
						chem_reax[react2][react1].prod2 = prod1; // same particles should become same either way
						chem_reax[react2][react1].prob = chance;
						// add some check for a + a -> a + a != reaction?
						// no, maybe used for genetic drift and user error anyway
					}
					else {
						err = -3;
					}
				}
			}
			//-------------------------------------------------------------------------------------------
			else if (!strcmp(token,"info"))  { 
                elem_list_print(geom_sensors,4,stdout);
                elem_list_print(geom_electrodes,4,stdout);
                elem_list_print(geom_channels,4,stdout);
            }            
            else if (!strcmp(token,"coord"))  { 
                geometry_type = geo_extern;
                poly_translate(sys_size);
                //elec_translate();
            }            
            else if (!strcmp(token,"e")) {  // E "Elek_1a"     0(sequenceno)      1(polatity)  
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ai) && get1flt(&flo)) {
                    elec_add_electrode(strg, ai, flo);
                }
            }                
            else if (!strcmp(token,"i")) {  // I "segm_31"     1     1.0(rate)
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ai) && get1flt(&flo)) {
                    add_inflow(strg, ai, flo);
                    inout_flow = 1;
                }
            }
            else if (!strcmp(token,"o")) {  // I "segm_31"     1     1.0(rate)
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ai) && get1flt(&flo)) {
                    add_outflow(strg, ai, flo);
                    inout_flow = 1;
                }
            }                
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"p")) {  
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ai)) {
                    int j=0;
                    //printf(" P= %s  %d coords\n", strg,ai);
                    elem_p poly = NULL;
                    
                    if (strg[1] == 'E') { // Elec
                        poly = poly_malloc(strg,geom_electrodes);
                        poly_setcolor(poly,3); 
                    }
                    else if (strg[1] == 's') { //segm
                        poly = poly_malloc(strg,geom_sensors);
                        poly_setcolor(poly,2); 
                    }
                    else if (strg[1] == 'C') { // Chan
                        poly = poly_malloc(strg,geom_channels);
                        poly_setcolor(poly,5); 
                    }
                    else { err = -8; }  // unknown polygon
                    
                    if (err == 0) {
                        real_t flo3,flo2,flo;
                        while ( get3flt(&flo,&flo2,&flo3) && (ai>0)) {
                            //printf(" %d %d: %f,%f,%f\n", j, ai, flo, flo2, flo3 );
                            poly_add_point(poly, (fpoint_t){{flo,flo2,flo3}} );
                            j++;
                            ai--;
                        } 
                    }
                    err = 0;
                    if (ai != 0) err = -7; // missing triplet
                }
            }
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"E")) {  // E Name  ntag tag  (coords if tag=geom)  
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ntag)) {
                    while (ntag>0) {
                        char sfkt[MAX_NAME];
                        ntag--;
                        if (get1str(sfkt) && !strcmp(sfkt,"geom")) {
                            if (get1int(&ai) && (ai>0)) {   // coords:  num  x y z...
                                elem_p poly = NULL;
                                int j = 0;
                                
                                poly = poly_malloc(strg,geom_electrodes);
                                poly_setcolor(poly,3);
                                elec_add_electrode(strg, j, 1);
                                
                                if (err == 0) {
                                    while ( get1flt(&flo) && get1flt(&flo2) && get1flt(&flo3) && (ai>0)) {
                                        //printf(" %d %d: %f,%f,%f\n", j, ai, flo, flo2, flo3 );
                                        poly_add_point(poly, (fpoint_t){{flo,flo2,flo3}} );
                                        j++;
                                        ai--;
                                    } 
                                }
                                err = 0;
                                if (ai != 0) err = -7; // missing triplet
                            } 
                        }
                        else if (!strcmp(sfkt,"xxxx")) {  // another tag
                        }
                    } 
                }
            }                
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"S")) {  // S "Elek_1a"   koords 
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ntag)) {
                    while (ntag>0) {
                        char sfkt[MAX_NAME];
                        ntag--;
                        if (get1str(sfkt) && !strcmp(sfkt,"geom")) {
                            if (get1int(&ai) && (ai>0)) {   // coords:  num  x y z...
								
                                elem_p poly = NULL;
                                int j = 0;
                                
                                poly = poly_malloc(strg,geom_sensors);
                                poly_setcolor(poly,2);
                                
                                if (err == 0) {
                                    while ( get1flt(&flo) && get1flt(&flo2) && get1flt(&flo3) && (ai>0)) {
                                        //printf(" %d %d: %f,%f,%f\n", j, ai, flo, flo2, flo3 );
                                        poly_add_point(poly, (fpoint_t){{flo,flo2,flo3}} );
                                        j++;
                                        ai--;
                                    } 
                                }
                                err = 0;
                                if (ai != 0) err = -7; // missing triplet
                            }
                        }
                        else if (!strcmp(sfkt,"xxxx")) {  // another tag
                        }
                    } 
                }
            }                
            //-------------------------------------------------------------------------------------------
            else if (!strcmp(token,"P")) {  
                char strg[MAX_NAME];
                if (get1str(strg) && get1int(&ntag)) {
                    while (ntag>0) {
                        char sfkt[MAX_NAME];
                        ntag--;
                        if (get1str(sfkt) && !strcmp(sfkt,"geom")) {
                            if (get1int(&ai) && (ai>0)) {   // coords:  num  x y z...
								
                                int j=0;
                                elem_p poly = NULL;
                                poly = elem_list_findname(geom_channels,strg);
                                if (! poly) {
                                    poly = poly_malloc(strg,geom_channels);
                                    poly_setcolor(poly,5); 
                                }
								
                                if (err == 0) {
                                    while ( get1flt(&flo) && get1flt(&flo2) && get1flt(&flo3) && (ai>0)) {
                                        //printf(" %d %d: %f,%f,%f\n", j, ai, flo, flo2, flo3 );
                                        poly_add_point(poly, (fpoint_t){{flo,flo2,flo3}} );
                                        j++;
                                        ai--;
                                    } 
                                }
                                err = 0;
                                if (ai != 0) err = -7; // missing triplet
                            }
                        }
                        else if (!strcmp(sfkt,"xxxx")) {  // another tag
                        }
                    } 
                }
            }
            else err = -4;  // unknown token
        }
    }        
    return err;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void read_ctrl_file(char * fname)
{
    FILE * fi;
    if ((fi = fopen(fname,"r")) != NULL) {
        printf("Read %s\n", fname);
		strncpy(ctrl_file_name, fname, sizeof(ctrl_file_name)-1);
        int lcnt = 0;
        int err;
        while (!feof(fi)) {
            char line[MAX_FILE_LN];
            fgets(line, sizeof(line), fi);
            lcnt++;
            err = interpret_line(line);
            if (err < 0) {
                fprintf(stderr, "  Error in line %d: %d %s\n", lcnt, err, errors[abs(err)]);
				stopit("Error in ctrl file");
            }
        }
        fclose(fi);  
        //genconf();
		
    }
    else {
        fprintf(stderr, "  Error 11: cannot read ctrlfile %s\n", fname);
		stopit("Error 11: cannot read ctrlfile");
    }
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void template_info(FILE * fi) 
{    
    fprintf(fi, 
            "# control file for dpd3\n"
            "# ---------------------\n"
            "#  # = comment line\n"
            "#  you can create this file with:\n"
            "#     dpd3 -tfilename\n"
            "#  or if you already have a control file just type in:\n"
            "#     dpd3 oldfile [other options] -tnewfile\n"
            "#  so your old file is up to date with all new features\n"
            "#  (but your comments are lost, and this may not work)\n"
            "#  or push 'v' while running dpd3 to create a template from actual parameters\n"
            "#\n");
	
    fprintf(fi, 
            "sigma = %-12.3f          # factor for dissipative and random force\n"
            "randfac = %-12.5f        # factor for random force\n"
            "rho = %-12.3f            # particle density\n"
            "dt = %-12.4f             # integration interval\n"
            "#\n"
            "displayintv = %-12d    # display interval\n"
            "seed = %-12d           # random seed (chosen at random if not specified)\n"
            "#\n"
            "size = %-12d           # width (and height and possibly depth) of simulated system\n"
            "dim = %-12d            # dimension [2|3] of system\n"
			"depth = %-12d          # (if != 1 and 3-d) depth of non-cubic \"thin\" system\n"
            "#\n"
			,
			sigma,
            rfac,
            rho0,
            dt,
            display_intv,
            iseed,
            sys_size,
            n_dim,
			sys_depth
			);
	
    fprintf(fi, 
            "slice_disp = %-6.3f %-6.3f  %d # sliced display for 3d viewing: min delta active [0|1]\n"
            "#\n"
            "max_step = %-12d       # simulation steps\n"
            "stat_intv = %-12d      # interval for displaying statistics\n"
            "save_intv = %-12d      # interval for saving position & bond data (*.out)\n"
			"poly_intv = %-12d      # interval for census (.pout) file\n"
			"reax_intv = %-12d      # interval for reaction (.reax) file\n"
            "#\n"
            "snapshotfn = %-s       # filename for picture-snapshots (e.g. film/dpd_%%06d.bmp)\n"
            "snaps_intv = %d %d %d  #  interval for automatic snapshots if >0, start,end\n"
            "#\n"
			,
            min_slice, delta_slice, slice_mode,
            maxstp,
            incstp,
            issstp,
			poly_intv,
			reax_intv,
            snapsfn,
            snaps_stp, snaps_start, snaps_end
            );
	
	fprintf(fi,"# beta[i,j] = conservative interaction between different particle types \n");
    fprintf(fi,"#\n");
    fprintf(fi,"#               /   1          beta    \\\n");
    fprintf(fi,"#  Fc = alpha * | --------- - -------- |\n");
    fprintf(fi,"#               \\ rij^exp12   rij^exp6 /\n");
    fprintf(fi,"#\n");
	fprintf(fi,"# if exp12==0 and exp6==-1 or program was compiled with SIMPLE_CONS set in Makefile then\n");
	fprintf(fi,"#  Fc = alpha * (1 - beta * rij) = -alpha * beta * rij + alpha\n");  
    fprintf(fi,"#\n");
    fprintf(fi,"#    i j   alpha   beta    exp12 exp6\n");
	int i,j;
    for (i=MIN_TYP+1; i<MAX_TYP; i++) { // BUGBUG incredibly gross hack to prevent type 0 showing up
        for (j=i; j<MAX_TYP; j++) {
            fprintf(fi,"beta %d %d = %7.3f %7.3f %-4d  %-4d  # %s -- %s\n", 
					i,j, ALPHA(i,j), BETA(i,j), BETA12(i,j), BETA6(i,j),
					particle[i].name, particle[j].name );
        }
    }
    fprintf(fi,"beta_min = %-12.3f                # minimal rij for Fc calculation\n", beta_min);
    fprintf(fi,"beta_min_force = %-12.3f          # Fc if rij<beta_min\n", beta_min_force);
    fprintf(fi,"#\n");
	
    fprintf(fi,"# num_particle[i] = ratio of free particles of type i\n");
    for (i=MIN_TYP; i<MAX_TYP; i++) {
        fprintf(fi,"num_particle %d = %5.2f  # %s\n", 
				i, particle[i].percent, particle[i].name );
    }
    fprintf(fi,"#\n");
	
    fprintf(fi,"# num_chains[i] = ratio of polymers of type i\n");
    fprintf(fi,"#   maximun polymer length %d\n", MAX_POLYLEN);
    fprintf(fi,"#   maximun number of different polymers %d\n", MAX_POLYMERS);
    for (i=0; i<MAX_POLYMERS; i++) {
        fprintf(fi,"num_chain %d = %-5.2f  ", 
				i, polymere[i].percent );
        int j;
        for (j=0; j<MAX_POLYLEN; j++) {
            if (polymere[i].part[j] != NO_TYPE) {
                fprintf(fi," %2d", polymere[i].part[j]);
            }
            else {
                break;  // found 1. NO_TYPE
            }
        }
        fprintf(fi,"\n");
    }
    fprintf(fi,"#\n");
	
	fprintf(fi,"# Bond parameters\n");
	fprintf(fi,	
			"bond = %-12.3f             # bonding force between chained monomers\n"
			"spring_dist = %-12.3f      # bond minimal energy length\n"
			"bond_angle_strength = %-7.3f   # stiff bond force\n"
			"dyn_loops = %d                   # whether loop polymers can form or not [1|0]\n"
			"# much of the code cannot currently handle looping polymers, so this should be left off\n"
			"max_dynamic_polymer_length = %-2d # polymers will stop bonding at this length\n"
			,
			d_bond__, spring_dist, bond_angle_strength, dyn_loops, max_dyn_polylen);
    fprintf(fi,"#\n");
    fprintf(fi,"# Dynamic bond parameters:\n");
	fprintf(fi,"# default is no dynamic bonding\n");
	fprintf(fi,"# values <= 0 mean bonds won't form/break\n");
    fprintf(fi,"#        i j   form  break probabilities\n");
    for (i=MIN_TYP; i<MAX_TYP; i++) {
        for (j=i; j<MAX_TYP; j++) {
			if (BOND_FORM(i,j) > 0 || BOND_BREAK(i,j) > 0) {
				fprintf(fi,"dyn_bond %d %d = %6.f %6.f\n",
						i, j, BOND_FORM(i,j), BOND_BREAK(i,j));
			}
        }
    }
	fprintf(fi,"#\n");
	fprintf(fi,"# Simple chemistry:\n"
			"# For every timestep that a type-a is near a type-b particle, there is an X chance that\n"
			"# a will change to type c and b to type d.  Any or all of a-d may be the same, but all\n"
			"# reactions must have two reactants and two products.  See dyn_bond above for a-b <-> a + b\n"
			"# type reactions.\n"
			"#           a + b -> c + d\tX\n"
			);
    for (i=MIN_TYP; i<MAX_TYP; i++) {
        for (j=i; j<MAX_TYP; j++) {
			if (chem_reax[i][j].prob > 0) {
				fprintf(fi,"simple_chem %d + %d -> %d %d\t%4f\n",
						i, j, chem_reax[i][j].prod1, chem_reax[i][j].prod2, chem_reax[i][j].prob);
			}
        }
    }	
	fprintf(fi,
            "#\n"
			"# impenetrable brick\n"
			"# from x0, y0 (low, left) to x1,y1 (up, right) and active-flag, at init\n"
            "brick = %-6.3f %-6.3f %-6.3f %-6.3f  %d %d\n"
            "#\n"
			"# between xmin, xmax with velocity, diffy and active-flag\n"
            "pressx = %-6.3f %-6.3f %-6.3f %-6.3f  %d\n"
            "#\n"
			,
            brick_ul.k[0],  brick_ul.k[1],  brick_or.k[0],  brick_or.k[1], with_brick, init_with_brick,
            pressure_to_x_minpos, pressure_to_x_maxpos, pressure_ofs, pressure_diff_at_y, pressure_to_x
			);
	fprintf(fi,"# color[i] = RGB-color of particle of type i\n");
    for (i=MIN_TYP; i<MAX_TYP; i++) {
        fprintf(fi,"color %d = %-3d  %-3d  %-3d   # %s\n", 
				i, particle[i].color[0], particle[i].color[1], particle[i].color[2], particle[i].name );
    }
    fprintf(fi,"#\n");
	fprintf(fi,"# epsilon[i] = strength of e-field on particle of type i\n");
    for (i=MIN_TYP; i<MAX_TYP; i++) {
        fprintf(fi,"epsilon %d = %-12.3f  # %s\n", 
				i, particle[i].epsilon, particle[i].name );
    }
    fprintf(fi,"#\n"
			"geometry = %-12d       # one of the predefined geometries\n"
            "                              #    0        normal (all)\n"
            "                              #    1        a rectangular hole from -z to z\n"
            "                              #    2        a H-structure\n"
            "                              #    3        bell-shape\n"
            "#\n",
            geometry_type
			);
	fprintf(fi,"# There are a number of advanced options for controlling the shape of the space\n"
			"# including electrodes, blocks, etc.\n"
			"# Interested parties are advised to examine other .dat files which include these\n"
			"# and the source code itself, beginning with dpd_file.c and dpd_vars.*\n");
	fprintf(fi,"#\n");
}               


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void save_template(char * fname)
{
    FILE * fi;
    if ((fi = fopen(fname,"w")) != NULL) {
        template_info(fi);
		
        fclose(fi);
        fprintf(stderr, "  Wrote %s\n", fname);
    }
    else {
        fprintf(stderr, "  Error: cannot create templatefile %s\n", fname);
    }
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int get_i_num(char *str)
{
    int a;
    int b;
	
    b = (sscanf(str, "%i", &a) == 0);
    if (b != 0)
		a = 0;
    return a;
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
float get_f_num(char *str)
{
    float a;
    int b;
	
    b = (sscanf(str, "%f", &a) == 0);
    if (b != 0)
		a = 0;
    return a;
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void helpinvoke()
{
    printf("\n"
		   "Program call:\n"
		   "   dpd3 [options] [ctrlfiles] [options]\n"
		   "\n"
		   "   option:\n"
		   "       -gx        x=geometry\n"
		   "       -tfn       fn=write template ctrl file\n"
		   "       -ffn       fn=ctrl file\n"
		   "\n"
		   "       -dx     --disp=x     x=displayinterval [1]\n"
		   "       -h,?    --help       this help\n"
		   "       -v                   verbose\n"
		   "               --dim=x      x=dimension 2 or 3 [3]\n"
		   "               --dt=x.x     x.x=dt, integration time step [0.01]\n"
		   "       -bx.x   --bond=x.x   x.x=bond force between chained monomers [10.0]\n"
		   "       -rx.x   --rho=x.x    x.x=rho,   density of particles per unit cube [4.0]\n"
		   "               --sigma=x.x  x.x=sigma, factor for random and dissipative force [3.0]\n"
		   "               --rfac=x.x   x.x=rfac,  factor for random force [sqrt(3.0)]\n"
		   "       -sx     --size=x     x=system size [16]\n"
		   "               --seed=x     x=random seed [time]\n"
		   "\n"
		   "\n"
		   "   ESC: Program end\n"
		   "\n");
    exit(0);    
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void process_cmdline(int argc, char *argv[])
{    
    int i;
    for (i = 1; i <= argc - 1; i++) {
		
        if (argv[i][0] == '-') {
            if (!strncmp(argv[i],"--disp=",7) && strlen(argv[i])>7 ) {
                display_intv = get_i_num(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--seed=",7) && strlen(argv[i])>7 ) {
                iseed = get_i_num(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--help",6) && strlen(argv[i])==6 ) {
				helpinvoke();
            }
            else if (!strncmp(argv[i],"--size=",7) && strlen(argv[i])>7 ) {
                sys_size = get_i_num(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--dt=",5) && strlen(argv[i])>5 ) {
                dt = get_f_num(&argv[i][5]);
                if (verbose) printf("cmdline: dt = %f\n",dt);
            }
            else if (!strncmp(argv[i],"--rfac=",7) && strlen(argv[i])>7 ) {
                rfac = get_f_num(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--dim=",6) && strlen(argv[i])>6 ) {
                n_dim = get_i_num(&argv[i][6]);
            }
            else if (!strncmp(argv[i],"--rho=",6) && strlen(argv[i])>6 ) {
                rho0 = get_f_num(&argv[i][6]);
            }
            else if (!strncmp(argv[i],"--sigma=",8) && strlen(argv[i])>8 ) {
                sigma = get_f_num(&argv[i][8]);
            }
            else if (!strncmp(argv[i],"--bond=",7) && strlen(argv[i])>7 ) {
                d_bond__ = get_f_num(&argv[i][7]);
            }
            else if ((argv[i][2]==0)||(argv[i][2]!='-')) {
                switch (argv[i][1]) {
					
                    case 'g':
                        geometry_type = (geometry_t)get_i_num(&argv[i][2]);
                        break;
						
                    case 'h':
                        helpinvoke();
                        break;
						
					case '?':
						helpinvoke();
						break;
						
                    case 'd':
                        display_intv = get_i_num(&argv[i][2]);
                        break;
						
                    case 'r':
                        rho0 = get_f_num(&argv[i][2]);
                        break;
						
                    case 'b':
                        d_bond__ = get_f_num(&argv[i][2]);
                        break;
						
                    case 's':
                        sys_size = get_i_num(&argv[i][2]);
                        break;
                        
                    case 't':   
                        save_template(&argv[i][2]);
                        exit(2);
                        break;
                        
                    case 'f':
                        read_ctrl_file(&argv[i][2]);
                        break;
                        
                    case 'v':
                        verbose = 1;
                        break;
                        
                    default:
                        printf("Unknown -commandline parameter %s\n", argv[i]);
						helpinvoke();
                        stopit("Fatal error");
                }    
            }
            else {
                printf("Unknown --commandline parameter %s\n", argv[i]);
				helpinvoke();
                stopit("Fatal error");
            }
        }
        else {
            read_ctrl_file(&argv[i][0]);
        }
    }    
}

void output_census() {
	FILE* fi;
	char fname[MAX_NAME];
	strncpy(fname, "", sizeof(fname));
	strcat(fname, ctrl_file_name);
	strcat(fname,".pout");
	if ((fi = fopen(fname,"a")) != NULL) {
		fprintf(fi, "Step = %d\n", no_calc);
		poly_pop_rec *current = poly_pop_list_root;
		while (current != NULL) {
			if (current->count > 0) {
				fprintf(fi,"%d = %d\n", current->poly, current->count);
			}
			current = current->next;
		}
		fprintf(fi, "\n");
	} else {
		printf("Error opening bookkeeping file\n");
	}
	fclose(fi);
}

void output_reactions() {
	FILE* fi;
	char fname[MAX_NAME];
	strncpy(fname, "", sizeof(fname));
	strcat(fname, ctrl_file_name);
	strcat(fname,".reax");
	if ((fi = fopen(fname,"a")) != NULL) {
		fprintf(fi, "Step = %d\n", no_calc);
		reaction_record *current = reaction_list_root;
		while (current != NULL) {
			fprintf(fi, "%s\t%d\n", current->reaction, current->count);
			 //BUGBUG also need to list which particles were involved for later catalyst analysis
			current = current->next;
		}
		fprintf(fi,"\n");
	} else {
		printf("Error opening reaction file\n");
	}
	fclose(fi);
	reset_reaction_list();
}

void kill_bookkeeping_files() {
	char fname[MAX_NAME];
	
	// kill the census file
	strncpy(fname, ctrl_file_name, MAX_NAME);
	strcat(fname,".pout");
	if (unlink(fname) != 0) {
		if (errno == ENOENT) {
			// file doesn't exist, no worries
		}
		else {
			stopit("Error removing census file");
		}
	}
	
	// now kill the reaction file
	strncpy(fname, ctrl_file_name, MAX_NAME);
	strcat(fname,".reax");
	if (unlink(fname) != 0) {
		if (errno == ENOENT) {
			// file doesn't exist, no worries
		}
		else {
			stopit("Error removing reaction file");
		}
	}
}

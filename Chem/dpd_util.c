/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
tmm 20.07.2004  17:15h
$Id: dpd_util.c,v 1.7 2004/08/10 15:13:05 tmaeke Exp $
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
#include "dpd_dynamic.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/

/*Funcs*************************************************************/

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void print_info2()
{        
    printf("Wall(w):%d  Brake(b):%d  Pressure(p):%d (x:%d y:%d) Info(i):%d\n"
           "Brick(m):%d  E-Field(e):%d  Polarity(n):%1.0f  Power(q):%7.4f\n"
           "Diss(d):%d  Cons(c):%d  Rand(r):%d\n", 
		   wall_at_y0, wall_brakes, pressure_to_x, add_vx, add_vy, drawinfo,
		   with_brick, efeld, polarity, epsilonfac,
		   do_diss, do_cons, do_rand);
	print_census();
	// print_reactions();
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void print_info()
{
    printf("   Size: (%4.1f,%4.1f,%4.1f)  (%d,%d,%d)  Rho=%4.2f  %dD\n", 
		   dl.k[0], dl.k[1], dl.k[2], 
		   dli.k[0], dli.k[1], dli.k[2], 
		   rho0, n_dim);
    printf("   num_particle=%d  volume=%f  maxstep=%d\n", num_particle, vol, maxstp);
    printf("   rand-seed = %d\n", iseed);
    printf("   dt = %7.4f\n", dt);
    printf("   bond = %7.4f\n", d_bond__);
    printf("   rfac = %7.4f\n", rfac);
    printf("   outside = %d\n", outsiders);
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void set_dt(real_t t)
{        
    dt = t;                           /* time step for the integrator [0.01] */
    tmax = dt * (real_t) maxstp;     /* time of the production run */
    dth = dt * .5;                    /* time steps for the integrator */
    dtrth = sqrt(dt) * .5;
    //printf("   dt = %7.4f\n", dt);
}    

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void change_dt(int wie) 
{
    if (wie) {
        set_dt(2*dt);
    }
    else {
        set_dt(dt/2);
    }
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void change_bond(int wie) 
{
    if (wie) {
        d_bond__ *= 2;
    }
    else {
        d_bond__ /= 2;
    }
    // printf("   bond = %7.4f\n", d_bond__);
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void change_rfac(int wie) 
{
    if (wie) {
        rfac = rfac*2;
    }
    else {
        rfac = rfac/2;
    }
    // printf("   rfac = %7.4f\n", rfac);
}

/*******************************************************************/
void savedata(int no)
{
    FILE * fi;
	char *fn = (char *) malloc(MAX_NAME);
	snprintf(fn, MAX_NAME, "data_%d.out", no);
	if ((fi = fopen(fn, "w")) != NULL) {
		fprintf(fi, "%d\n", num_particle);
		int n;
		for (n=0; n<num_particle; n++) {			
			fprintf(fi, "%7.4f %7.4f %7.4f  %d %d  %d %d  %d\n", 
					liste[n].r.k[0], liste[n].r.k[1], liste[n].r.k[2], 
					n, liste[n].type, 
					liste[n].left, liste[n].right, 
					liste[n].inside);
		}
		fclose(fi); 
	} 
}


/*******************************************************************/
int loaddata(int no)
{
    static int loading = 0;
    FILE * fi;
    int i = 0;
    char buf[1000];
    float x,y,z;
    int num, col, left, right, inside;
	//    int fnum;
    if (! loading) {
        loading = 1;
        char fn[30] = "Error";
		//        for (fnum = 0; fnum < MAX_TYP; fnum++) {
		//            sprintf(fn,"data%d_%d.out", fnum, no);
		sprintf(fn,"data_%d.out", no);
		if ((fi = fopen(fn, "r")) != NULL) {
			fgets(buf, sizeof(buf), fi);
			sscanf(buf, "%d", &i);
			num_particle = i;
			
			int n;
			for (n=0; n<i; n++) {
				fgets(buf, sizeof(buf), fi);
				sscanf(buf, "%f %f %f %d %d %d %d %d", &x, &y, &z, &num, &col, &left, &right, &inside);
				liste[num].r.k[0] = x;
				liste[num].r.k[1] = y;
				liste[num].r.k[2] = z;
				liste[num].type = col;
				liste[num].left = left;
				liste[num].right = right;
				liste[num].inside = inside;
				liste[num].v.k[0] = 0;
				liste[num].v.k[1] = 0;
				liste[num].v.k[2] = 0;
				liste[num].fd = liste[num].v;  // 0.0.0
				liste[num].fc = liste[num].v;
				liste[num].fr = liste[num].v;
				// hier was tun // ich bin eine jelly donut
				if (feof(fi)) break;
			}
			fclose(fi); 
		}
		//      }
        printf("Loaded files: %s, %d particles\n", fn, num_particle);
        update = 1;
        loading = 0;
    }
    return i; 
}


/*******************************************************************/
int loadinitialdata()
{
    static int loading = 0;
    FILE * fi;
    int i = 0;
    char buf[1000];
    float x, y, z;
    int num, col, left, right, inside;
    if (! loading) {
        loading = 1;
		if ((fi = fopen(initial_position_file, "r")) != NULL) {
			fgets(buf, sizeof(buf), fi);
			sscanf(buf, "%d", &i);
			num_particle = i;
			
			int n;
			for (n=0; n<i; n++) {
				fgets(buf, sizeof(buf), fi);
				sscanf(buf, "%f %f %f %d %d %d %d %d", &x, &y, &z, &num, &col, &left, &right, &inside);
				liste[num].r.k[0] = x;
				liste[num].r.k[1] = y;
				liste[num].r.k[2] = z;
				liste[num].type = col;
				liste[num].left = left;
				liste[num].right = right;
				liste[num].inside = inside;
				liste[num].v.k[0] = 0;
				liste[num].v.k[1] = 0;
				liste[num].v.k[2] = 0;
				liste[num].fd = liste[num].v;  // 0.0.0
				liste[num].fc = liste[num].v;
				liste[num].fr = liste[num].v;
				// hier was tun // ich bin eine jelly donut
				if (feof(fi)) break;
			}
			fclose(fi); 
		}
		//      }
        printf("Loaded files: %s, %d particles\n", initial_position_file, num_particle);
        update = 1;
        loading = 0;
    }
    return i; 
}




/*******************************************************************/
/*
 void fill_nb3D() {
	 int c;
	 int size_x = 1;
	 int size_y = sys_size;
	 int size_z = sys_size*sys_size;
	 for (c=0; c<27; c++) {
		 nb3Dcell[c] = size_x * nb3Dx[c] + size_y * nb3Dy[c] + size_z * nb3Dz[c];
	 }
 }
 */


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
real_t ipow (real_t v, int pw)
{
    real_t x,y;
    switch (pw) {
		case 0:  return 1.0;   break;
		case 1:  return v;     break;
		case 2:  return v*v;   break;
		case 3:  return v*v*v; break;
		case 4:  x=v*v;   return x*x;   break;
		case 5:  x=v*v;   return x*x*v; break;
		case 6:  x=v*v;   return x*x*x; break;
		case 7:  x=v*v*v; return x*x*v; break;
		case 8:  x=v*v;   y=x*x; return y*y;     break;
		case 9:  x=v*v*v;        return x*x*x;   break;
		case 10: x=v*v;   y=x*x; return y*y*x;   break;
		case 11: x=v*v;   y=x*x; return y*y*x*v; break;
		case 12: x=v*v;   y=x*x; return y*y*y;   break;
		case -1:  return 1/v;     break;
		case -2:  return 1/(v*v);   break;
		case -3:  return 1/(v*v*v); break;
		case -4:  x=v*v;   return 1/(x*x);   break;
		case -5:  x=v*v;   return 1/(x*x*v); break;
		case -6:  x=v*v;   return 1/(x*x*x); break;
		case -7:  x=v*v*v; return 1/(x*x*v); break;
		case -8:  x=v*v;   y=x*x; return 1/(y*y);     break;
		case -9:  x=v*v*v;        return 1/(x*x*x);   break;
		case -10: x=v*v;   y=x*x; return 1/(y*y*x);   break;
		case -11: x=v*v;   y=x*x; return 1/(y*y*x*v); break;
		case -12: x=v*v;   y=x*x; return 1/(y*y*y);   break;
			
		default: printf ("Error exponent %f %d\n", v,pw);
    }
    return 1.0;
}


void test_pow()
{
    int i;
    for (i=-12; i<=13; i++) {
        printf(" i=%3d  pow(i)=%f\n", i, ipow(2.0, i) );
    }
}


/*******************************************************************/
void stopit (char * str)
{
    fprintf(stderr, str);
    exit(1);
}    


/*******************************************************************/
inline real_t d__nint(real_t i)
{
    if (i<0) return (int)(i-0.5);
	return (int)(i+0.5); 
}

/*******************************************************************/
inline void Boundaryfunc(
						 fpoint_t * r, 
						 int kp)
{
    if (bounds) return;
    
    r->k[kp] -= 
		dl.k[kp] * (     (r->k[kp]*dlinv.k[kp] <0.0) ? 
						 (int)((r->k[kp]*dlinv.k[kp])-0.5) :
						 (int)((r->k[kp]*dlinv.k[kp])+0.5)  );
}



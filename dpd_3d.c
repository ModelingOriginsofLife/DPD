/*
      dpd

      $Id: dpd_3d.c,v 1.19 2004/07/19 15:24:52 tmaeke Exp tmaeke $
      compile:
      gcc dpd_3d.c SDL_etc.c SDL_prim.c -I. -lSDL -I/usr/include/SDL  -O5 -lm -o dpd3
For mac:
gcc -v dpd_3d.c dpd_rand.c SDL_etc.c SDL_prim.c -I. -lSDL -lSDLmain -framework Cocoa -lpng -L/usr/local/lib -I/usr/local/include/SDL  -O5 -lm -o dpd3

      
      
      dpd3 -p2  --dim=2  --size=12   -i5  --rho=10 --disp=1 -c12  -g0 -m2 --sigma=3  --dt=0.16  --rfac=0.0541
      dpd3 -p3  --dim=2  --size=24   -i5  --rho=10 --disp=1 -c14  -g0 -m2 --sigma=3 --dt=0.04
      dpd3 -p2  --dim=2  --size=16   -i5  --rho=20 --disp=1 -c1000  -g0 -m2 --sigma=3 --dt=0.02 --bond=60
*/

#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <SDL/SDL.h>
#include "SDL_etc.h"
#include "SDL_prim.h"
#include "dpd_rand.h"
            
/* Table of constant values */

#define SYS_SIZEMAX  64

#define SYS_SIZE2   (SYS_SIZEMAX+2)
#define SYS_SIZESQ  (SYS_SIZE2*SYS_SIZE2)


#define INIT_DT      0.01    // Integration time

#define F_SIZE       100000
#define F_POINTS     10000000  // max. number of entries in Verlet-Table
                     // an entry in unit cube (max (SYS_SIZEMAX+2)^3)
#define F_ASIZE      ((SYS_SIZEMAX+2)*(SYS_SIZEMAX+2)*(SYS_SIZEMAX+2))  

#define STEPS        10000000  // simulation time


// graphics
#define MAX_X  256     // pixel a pixperdot size
#define MAX_Y  256
#define PIXPERDOT 3
#define SDLRAND   0

#define MAXCOLOR 256

#define MOLMOL 4    // defines size of alpha-array
#define ALPHA(nfi,nfj)  (alpha[(nfi) + (nfj)*MOLMOL])
#define MAX_MOLMOL (MOLMOL*MOLMOL)


#define MIN_TYP   0
#define TYP_POLY  0
#define TYP_WATER 1
#define TYP_OIL   2
#define TYP_POLY2 3
#define MAX_TYP   4


#define TYPOFS    10    // in color table
#define TYPSHADES 40    // number of shades per color per molecule=FClType
#define ZSHADES (TYPSHADES-16)
#define RGBSHADEFAC 6
#define TYP2COL(typ,shade) ((((typ)*TYPSHADES+(shade)) + TYPOFS)% MAXCOLOR)    

#define Boundary(rij,kk) (dl.k[kk] * d__nint((rij.k[kk]) * dlinv.k[kk]))

#define NumCube(r) ((int)(r.k[0])*dldiff.k[0] + (int)(r.k[1])*dldiff.k[1] + (int)(r.k[2])*dldiff.k[2] + F_ASIZE/2)
#define ALLOWED(r) (allowed[NumCube(r)])
 
#define max(a,b)  ((a)>(b) ? (a):(b))

#define DEBUGOUT(str)
//#define DEBUGOUT(str)  fprintf(stderr,"dpd " str "\n")

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
//    typedef float  float_t;  // try float or double
    typedef int    int_t;
    
    typedef struct {
        float_t   k[3];    
    } fpoint_t;

    typedef struct {
        int_t k[3];
    } ipoint_t;

    typedef struct {
        short int type;      // Defines the particle: Water, oil, monomere ...
        short int drawn;     // for drawing
        int       left;      // bounded neighbor left side (index)
        int       right;     // bounded neighbor right side (index)
        int       nextindex; // next link in stack
        fpoint_t  r;         // Positions of particles
        fpoint_t  v;         // Velocities of particles (x coordinate)
        fpoint_t  r0;        // Saved positions
        fpoint_t  fr;        // Random forces
        fpoint_t  fc;        // Conservative forces
        fpoint_t  fd;        // Dissipative forces
    } particle_t;

    typedef struct {                
        int ip;                 // index of particle
        int nextindex;          // next link in chain
    } chainlink_t;
                                                        


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    // dimension
    int n_dim     = 3;
    int sys_size  = 16;
    int sys_size2 = 8;
    int skin_size;

    // for SDL display
    static int ofsx = 0;
    static int ofsy = 0;
    static int ena_disp[10] = {1,1,1,1,1,1,1,1,1,1};  // switch on/off displayed elements

    // siehe genconf
    static int mere_types = 2;             // initialization pattern
    static int particle_types = 1;         // different kinds of particles
    static int interaction_types = 1;      // configuration of alpha
    static int geometry_type = 0;          // 
    
    static int snapsintv = 50;             // snapshot interval
    static int display_intv = 1;           // SDL display interval
    static int chain_len = 600;            // number of chained particles

    
    // some vars to turn on/off features
    static int xchg_yz = 0;                // Anzeige: y und z vertauschen
    
    static int wall_at_y0 = 0;              // simulate a wall at y=0
    static int wall_brakes = 0;             // does the wall brakes particles? vx=0
    
    static int drawinfo = 0;               // draw infos in animation
    static int do_diss = 1;                // dissipative Term rechnen
    static int do_cons = 1;                // conservativeb Term rechnen
    static int do_rand = 1;                // random Term rechnen
    
    static int with_brick = 0;               // a brick in the way
    float_t brick_xl = -7.0;
    float_t brick_xr = -6.0;
    float_t brick_yu = -2.0;
    float_t brick_yo =  2.0;
    
    int pressure_to_x = 0;                    // accelerate a particle 
    float_t pressure_ofs = 1.0;                 // velocity increment (if pressure_to_x=1)
    float_t pressure_to_x_minpos = -10.1;      // between these two x limits
    float_t pressure_to_x_maxpos = -10.0;      // 
    float_t pressure_diff_at_y = 0.0;           // above this only pressure_ofs/2 applied (diff pressure for H-Structures)
    int add_vx = 0;                        // give once all particles vx +=pressure_ofs
    int add_vy = 0;                        // give once all particles vy +=pressure_ofs

    static int efeld = 0;                  // E-field on/off
    fpoint_t  pol1 = { {0.0, 0.0, -1.0} };
    fpoint_t  pol2 = { {0.0, 0.0, 1.0} };
    float_t polarity = 1;                   // 1 or -1
    float_t emin = 1; // 0.25;              // cutoff around electrode
    float_t epsilon[MOLMOL];                // electrical attraction/repulsion : particle specific
    float_t epsilonfac = 1;                 // an extra factor to dec-/increase for all particles
    
    /* Local variables */
    
    particle_t  liste[F_SIZE];    

    static int npoint[F_SIZE];     // For the Verlet neighbor tables
    static int nvsize;             // For the Verlet neighbor tables
    static int nvintv;             // For the Verlet neighbor tables
    static int maxverlet = F_POINTS;
    static int list_verl[F_POINTS];        // For the Verlet neighbor tables
    
    static short int allowed[F_ASIZE];
    
    // chain linking update to Verlet neighbors
    int chainarray[F_ASIZE];        

    static int nb3Dx[27] = { 0, 1,-1, 0, 1,-1, 0, 1,-1,  0, 1,-1, 0, 1,-1, 0, 1,-1,  0, 1,-1, 0, 1,-1, 0, 1,-1 };
    static int nb3Dy[27] = { 0, 0, 0, 1, 1, 1,-1,-1,-1,  0, 0, 0, 1, 1, 1,-1,-1,-1,  0, 0, 0, 1, 1, 1,-1,-1,-1 };
    static int nb3Dz[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1, 1, -1,-1,-1,-1,-1,-1,-1,-1,-1 };


    
    static float_t dtrth; 
    static float_t dth; 
    static float_t vol;
    
    static int it;
	static int nconf;
    
    static float_t vcm;               // Center-of-mass velocity
    static float_t pm;                // Particle mass
    static float_t pminv;             // Inverse particle mass
    
    static fpoint_t dl;               // Size of the simulation box in 2D
    static fpoint_t dlinv;            // Inverse system size 
    static ipoint_t dli;              // same as int
    static ipoint_t dldiff;
    
    static float_t d_bond__;          // Harmonic force constant for the chain
    static float_t dt;                // Time step for the molecular dynamics
    static int    num_particle;       // Total particle number
    static float_t dnum_particle;     // Total particle number
    static float_t dnum_particle_inv; // Inverse particle number
    
    
    static int update;       // int variable for the Verlet tables
    static int isteps;       // Number of time steps for equilibration
    static int maxstp;       // Number of time steps after equilibration
    static int issstp;       // Number of steps between snapshots
    static int incstp;       // Number of steps between samples
    
    static int iseed;        // For the random number generator

    
    static float_t rfac;      // Scaling factor of the random force
    static float_t tfac;      // Scaling factor of temperature
    static float_t skin;      // Skin radius for the Verlet tables
    static float_t skinsq;    // Skin radius squared
    static float_t temp;      // Temperature
    static float_t tmax;      // Time of the simulation after equilibration

    static int n_water__;     // index to water particles
    static int n_oil; 
    static int n_solu;
    static int n_poly2;
    static int len_chain__; 
    static int len_mol2;
    
    static float_t rho;                 // Density
    static float_t rho0;
    static float_t alpha[MAX_MOLMOL];   // Conservative force strengths
    static float_t gamma_;              // Dissipative force strength
    static float_t sigma;               // Random force strength
    static float_t omegaofs = 0;

    // velocity profile:    
#define VELO_ANZ 50
    float_t velo[VELO_ANZ];
    int veloanz[VELO_ANZ];
    int drawinfoupd = 0;

#define VACF_SIZE    100
    //static float_t dr_gyr__;              // Radius of gyration
    //static float_t dr_ee__;               // End-to-end distance of the chain
    static fpoint_t vel[VACF_SIZE];         // For the velocity autocorrelation function
    //static int    max_vacf = VACF_SIZE;
    static float_t dcs_vacf__[VACF_SIZE];
    static float_t dcs_vacf__[VACF_SIZE];
    static float_t dvacf[VACF_SIZE];
    static int    it_bound__[VACF_SIZE+1];
    static int    it_vacf__;
    static int    ivacf_loop__[2*VACF_SIZE+1];
    static int    ivacf_counter__;



// Prototypes
    extern  int genconf();

    extern  int fcdr_(); 
    extern  int fdr_();
    extern  int intv_(float_t, float_t, float_t);
    extern  int intr_(float_t);
    extern  int check_();
    extern  int foldb_();
    extern  void fill_cell_chains(int);

    extern  void init_vacf(void);
    extern  int r_gyr__(); 
    extern  int snaps_();
    extern  int vacf_();
    extern  int vcmtp_();
    extern  int velo_prof(void);
    

/*******************************************************************/
void stopit (char * str)
{
    printf(str);
    exit(1);
}    
    

/*******************************************************************/
float_t d__nint(float_t i)
{
    if (i<0) return (int)(i-0.5);
             return (int)(i+0.5); 
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
float_t boounds (float_t i, float_t d) 
{
    return (d * d__nint((i) / d));
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void testbnd() 
{
    int i;
    float_t j = -5;
    for (i=0; i <=40; i++) {
        float_t erg = boounds(j,2);
        printf("  i:%2d  j:%7.3f  j+0.5:%7.3f  (int)j+0.5:%3d     max:2.0  bnd:%7.3f   j-bnd:%7.3f\n", 
                  i,     j,       j+0.5,       (int)(j+0.5),               erg,        j-erg);
        j = j + 0.2;
    }
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void set_dt(float_t t)
{        
    dt = t;                           /* time step for the integrator [0.01] */
    tmax = dt * (float_t) maxstp;     /* time of the production run */
    dth = dt * .5;                    /* time steps for the integrator */
    dtrth = sqrt(dt) * .5;
    printf("   dt = %7.4f\n", dt);
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
void change_rfac(int wie) 
{
    if (wie) {
        rfac = rfac*2;
    }
    else {
        rfac = rfac/2;
    }
    printf("   rfac = %7.4f\n", rfac);
}

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void change_polar(int wie) 
{
    if (wie) {
        polarity = 1.0;
    }
    else {
        polarity = -1.0;
    }
    printf("   rfac = %7.4f\n", rfac);
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


/*******************************************************************/
void savedata(int no)
{
    FILE * fi;
    int i;
    for (i = MIN_TYP; i < MAX_TYP; i++) {
        char fn[20];
        sprintf(fn,"data%d_%d.out", i, no);
        if ((fi = fopen(fn, "w")) != NULL) {
            fprintf(fi, "%d\n", num_particle);
            int n;
            for (n=0; n<num_particle; n++) {
                if (liste[n].type == i) {
                    fprintf(fi, "%7.4f %7.4f %7.4f %d %d %d\n", 
                            liste[n].r.k[0], liste[n].r.k[1], liste[n].r.k[2], n, liste[n].type, liste[n].left);
                }
            }
            fclose(fi); 
        } 
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
    int num,col,link;
    int fnum;
    if (! loading) {
        loading = 1;
        char fn[30] = "Error";
        for (fnum = 0; fnum < MAX_TYP; fnum++) {
            sprintf(fn,"data%d_%d.out", fnum,no);
            if ((fi = fopen(fn, "r")) != NULL) {
                fgets(buf, sizeof(buf), fi);
                sscanf(buf, "%d", &i);
                num_particle = i;

                int n;
                for (n=0; n<i; n++) {
                    fgets(buf, sizeof(buf), fi);
                    sscanf(buf, "%f %f %f %d %d %d", &x,&y,&z,&num,&col,&link);
                    liste[num].r.k[0] = x;
                    liste[num].r.k[1] = y;
                    liste[num].r.k[2] = z;
                    liste[num].left = link;
                    liste[num].type = col;
                    liste[num].v.k[0] = 0;
                    liste[num].v.k[1] = 0;
                    liste[num].v.k[2] = 0;
                    liste[num].fd = liste[num].v;  // 0.0.0
                    liste[num].fc = liste[num].v;
                    liste[num].fr = liste[num].v;
                    // hier was tun
                    if (feof(fi)) break;
                }
                fclose(fi); 
            }
        }
        printf("Loaded files: %s, %d particles\n", fn, num_particle);
        update = 1;
        loading = 0;
    }
    return i; 
}


/*******************************************************************/
void helpinteract()
{
    printf("\n"
        "Program call:\n"
        "   dpdc [option]\n"
        "\n"
        "   option:\n"
        "       -h       help for call program\n"
        "\n"
        "\n"
        "   Mouse left:  center at mouseposition\n"
        "   Mouse right: take a snapshot in 'dpd.bmp'\n"
        "\n"
        "   Upper case letters: turn off or halve a parameter\n"
        "   Lowe case letters: turn on or double a parameter\n" 
        "\n"
        "   W:  Wall at Y=0\n"
        "   M:  Rectangle in the flow\n"
        "   B:  Wall/Rectangle: brake tangent velocity\n"
        "   I:  Velocity profile in X-direction\n"
        "   P:  Create pressure to X at X=0  (vx += const)\n"
        "   E:  E-field\n"
        "   N:    Change polarity\n"
        "   Q:    Change field strength\n"
        "\n"
        "   C:  Conservative term\n"
        "   D:  Dissipative term\n" 
        "   R:  Random term\n"
        "\n"
        "   X:  All particles vx += %4.1f\n"
        "   Y:  All particles vy += %4.1f\n"
        "\n"
        "   F:  Random factor %7.4f\n"
        "   T:  dT  %7.4f\n"
        "\n"
        "   z:  Exchange y-z for display\n"
        "\n"
        "   S:   Save-Data\n"
        "   L:   Load-Data\n"
        "   H:   this help\n"
        "   ESC: End\n"
        "\n",  pressure_ofs, pressure_ofs, rfac, dt);
}


/*******************************************************************/
void handleevent() 
{        

    SDL_PollEvent(&event);
    if (Quit) exit(1);

    if (KeyDown) {
        int is = 1;
        //printf("Mod=%04X  Sym=%04X\n", MOD, SYM);
        if (MOD & KMOD_SHIFT) { 
           is = 0;
        }
        switch (SYM) {
            case SDLK_w : wall_at_y0 = is; break;
            case SDLK_p : pressure_to_x = is; break;
            case SDLK_b : wall_brakes = is; break;
            case SDLK_i : drawinfo = is; break;
            case SDLK_m : with_brick = is; break;
            case SDLK_d : do_diss = is; break;
            case SDLK_c : do_cons = is; break;
            case SDLK_r : do_rand = is; break;
            case SDLK_t : change_dt(is); break;
            case SDLK_f : change_rfac(is); break;
            case SDLK_h : helpinteract(); break;
            case SDLK_l : loaddata(1+is); break;
            case SDLK_x : add_vx = 1; break;
            case SDLK_y : add_vy = 1; break;
            case SDLK_e : efeld = is; break;
            case SDLK_q : change_efield(is); break;
            case SDLK_z : xchg_yz = 1-xchg_yz; break;
            case SDLK_n : change_polar(is); break;
            case SDLK_s : savedata(1+is); break;
            case SDLK_0 : ena_disp[0] = 1-ena_disp[0]; break;
            case SDLK_1 : ena_disp[1] = 1-ena_disp[1]; break;
            case SDLK_2 : ena_disp[2] = 1-ena_disp[2]; break;
            case SDLK_3 : ena_disp[3] = 1-ena_disp[3]; break;
            case SDLK_4 : ena_disp[4] = 1-ena_disp[4]; break;
            case SDLK_5 : ena_disp[5] = 1-ena_disp[5]; break;
            case SDLK_6 : ena_disp[6] = 1-ena_disp[6]; break;
            case SDLK_7 : ena_disp[7] = 1-ena_disp[7]; break;
            case SDLK_8 : ena_disp[8] = 1-ena_disp[8]; break;
            case SDLK_9 : ena_disp[9] = 1-ena_disp[9]; break;
            default: break;      
        }
        KeyDown = 0;
        //event.key.keysym.sym = 0;
        printf("Wall(w):%d  Brake(b):%d  Pressure(p):%d (x:%d y:%d) Info(i):%d\n"
               "Brick(m):%d  E-Field(e):%d  Polarity(n):%1.0f  Power(q):%7.4f\n"
               "Diss(d):%d  Cons(c):%d  Rand(r):%d\n"
               "dt(t):%7.4f  randfac(f):%7.4f  y-z:%d\n", 
                wall_at_y0, wall_brakes, pressure_to_x, add_vx, add_vy, drawinfo,
                with_brick, efeld, polarity, epsilonfac,
                do_diss, do_cons, do_rand, 
                dt, rfac, xchg_yz);
    }
    
    if (Mouse.B){
        switch (Mouse.B){
            case 1:       
                ofsx = (3*MAX_X*PIXPERDOT/2 + ofsx - Mouse.X) % (MAX_X*PIXPERDOT);
                ofsy = (3*MAX_Y*PIXPERDOT/2 + ofsy - Mouse.Y) % (MAX_Y*PIXPERDOT);
                break;
            case 2:     
                SDL_SaveBMP(screen,"dpd.bmp");   
                printf("wrote dpd.bmp\n");
                break;
            case 3:
                break;
       }
       Mouse.B=0;
       printf("Mouse at: %d %d   ofs=%d %d\n",Mouse.X,Mouse.Y, ofsx, ofsy); 
    }
}    

/*******************************************************************/
/*******************************************************************/
    static Uint32 Color[MAXCOLOR]; // color table

    static int basecolors[MAX_TYP][3] = {
        { 240, 240,  240 },
        { 0x00,0x00, 250 },
        { 250, 0x00,0x00 },
        { 0x00,250, 0x00 }
    };

    int runde;


/*******************************************************************/
void DrawIt()
{
    #define toX(val) (int)(PIXPERDOT*(MAX_X/(dl.k[0]+SDLRAND)*  (val) +MAX_X/2))
    #define toY(val) (int)(PIXPERDOT*(MAX_Y/(dl.k[1]+SDLRAND)*(-(val))+MAX_Y/2))
    #define toZ(val) (int)(ZSHADES/dl.k[2]*(val)+TYPSHADES/2)    // 0..100% 
    
    #define otoX(val) ((toX(val)+ofsx)%(MAX_X*PIXPERDOT))
    #define otoY(val) ((toY(val)+ofsy)%(MAX_Y*PIXPERDOT))
    #define otoZ(val) (toZ(val))
    
    static int drawn = 0;
    static int first = 0;
    static int dx;
    static int dy;
    static int x0;
    static int y0;
    
    if (!first) {
        first = 1;
        dx = toX(2.0)-toX(1.0);
        dy = toY(2.0)-toY(1.0);
        x0 = 100;
        y0 = 100;
    }
    
    if (++runde % display_intv == 0) { 

        ClearSurfaceRGB(screen,0,0,0);
        SDL_LockSurface(screen);

        // draw square of size 1.0x1.0
        SDL_drawRect_Alpha(screen, x0, y0, x0-dx, y0+dy, Color[5], 228);

        SDL_drawLine_Alpha(screen, otoX(0), otoY(-dl.k[1]/4), otoX(0), otoY(dl.k[1]/4), Color[5], 228);
        SDL_drawLine_Alpha(screen, otoX(-dl.k[0]/4),otoY(0),  otoX(dl.k[0]/4), otoY(0), Color[5], 228);
        
        
        if (with_brick) {
            SDL_drawRect_Alpha(screen, otoX(brick_xl), otoY(brick_yu), otoX(brick_xr), otoY(brick_yo), Color[6], 150);
        }
        
        if (efeld) {
            int c1 = 7;
            int c2 = 8;
            if (polarity < 0) { c1=8; c2=7; }
            SDL_drawRect_Alpha(screen, otoX(pol1.k[0]-emin), otoY(pol1.k[1]-emin), 
                                       otoX(pol1.k[0]+emin), otoY(pol1.k[1]+emin), Color[c1], 150);
            SDL_drawRect_Alpha(screen, otoX(pol2.k[0]-emin), otoY(pol2.k[1]-emin), 
                                       otoX(pol2.k[0]+emin), otoY(pol2.k[1]+emin), Color[c2], 150);
        }
                
        int n;
        for (n=0; n<num_particle; n++) {
            if (liste[n].drawn == drawn) {
                int col, SX,SY;

                if (ena_disp[liste[n].type]) {
                    if (xchg_yz) {
                        col = Color[TYP2COL(liste[n].type, otoZ(liste[n].r.k[1]))];
                        SX = otoX(liste[n].r.k[0]);
                        SY = otoY(liste[n].r.k[2]);

                        if (liste[n].left >= 0) {
                            int lastxx = otoX(liste[liste[n].left].r.k[0]);
                            int lastyy = otoY(liste[liste[n].left].r.k[2]);
                            if ((abs(SX-lastxx)<MAX_X) && (abs(SY-lastyy)<MAX_Y) ) {
                                SDL_drawLine_Alpha(screen,SX,SY,lastxx,lastyy, col,125);
                            }
                        }
                    }
                    else {                
                        col = Color[TYP2COL(liste[n].type, otoZ(liste[n].r.k[2]))];
                        SX = otoX(liste[n].r.k[0]);
                        SY = otoY(liste[n].r.k[1]);

                        if (liste[n].left >= 0) {
                            int lastxx = otoX(liste[liste[n].left].r.k[0]);
                            int lastyy = otoY(liste[liste[n].left].r.k[1]);
                            if ((abs(SX-lastxx)<MAX_X) && (abs(SY-lastyy)<MAX_Y) ) {
                                SDL_drawLine_Alpha(screen,SX,SY,lastxx,lastyy, col,125);
                            }
                        }
                    }
                                
                    SDL_putPixel( screen, SX+0, SY+0, col); // 1x1
                    
                    if (PIXPERDOT > 1) {  // 2x2
                        SDL_putPixel( screen, SX+1, SY+0, col);
                        SDL_putPixel( screen, SX+0, SY+1, col);
                        SDL_putPixel( screen, SX+1, SY+1, col);
                    }
                    if (PIXPERDOT > 2) {  // 3x3
                        SDL_putPixel( screen, SX+2, SY+0, col);
                        SDL_putPixel( screen, SX+2, SY+1, col);
                        SDL_putPixel( screen, SX+2, SY+2, col);
                        SDL_putPixel( screen, SX+1, SY+2, col);
                        SDL_putPixel( screen, SX+0, SY+2, col);
                    }
                    if (PIXPERDOT > 3) {  // 4x4
                        SDL_putPixel( screen, SX+0, SY+3, col);
                        SDL_putPixel( screen, SX+1, SY+3, col);
                        SDL_putPixel( screen, SX+2, SY+3, col);
                        SDL_putPixel( screen, SX+3, SY+0, col);
                        SDL_putPixel( screen, SX+3, SY+1, col);
                        SDL_putPixel( screen, SX+3, SY+2, col);
                        SDL_putPixel( screen, SX+3, SY+3, col);
                    }
                }
                liste[n].drawn = 1-drawn;
            }
        }
        drawn = 1-drawn;
        
        if (drawinfo) {
            if (1) {  // drawinfoupd) {
                drawinfoupd = 0;
                
                int i;
                float lastx=0,lasty=0;
                for (i=0; i<VELO_ANZ; i++) { 
                    float x = 0;
                    float y = 0;
                    if (veloanz[i]>0) x = velo[i]/veloanz[i]*dl.k[0]/8; 
                    y = -dl.k[1]/2 + dl.k[1]*i/(VELO_ANZ-1);
                    if (i>0) {
                        SDL_drawLine_Alpha(screen,otoX(lastx),otoY(lasty),otoX(x),otoY(y),Color[6],150);
                    }
                    lastx=x; lasty=y;
                }
            }
        }
        SDL_UnlockSurface(screen);
        SDL_UpdateRect(screen,0,0,0,0);

        char Str[255];
        sprintf(Str,"DPD %d", runde);
        SDL_WM_SetCaption(Str,Str);
    }

    handleevent();
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int maxrgb(
    int col)
{
    return (col>255 ? 255 : col);
}    

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void createcolors()
{
    int typ;
    for (typ = 0; typ < MAXCOLOR; typ++) {
        Color[typ] = SDL_MapRGB(screen->format,255,255,255);
    }
    for (typ = MIN_TYP;  typ < MAX_TYP; typ++) {
        int shade;
        for (shade=0; shade<TYPSHADES; shade++) {
             int r = maxrgb(basecolors[typ][0]+ shade*RGBSHADEFAC);
             int g = maxrgb(basecolors[typ][1]+ shade*RGBSHADEFAC);
             int b = maxrgb(basecolors[typ][2]+ shade*RGBSHADEFAC);
             Color[TYP2COL(typ,shade)] = SDL_MapRGB(screen->format, r,g,b);
             /*
             printf("  Color %3d %1d,%2d= %02X %02X %02X   %3d,%3d,%3d\n", 
                     TYP2COL(typ,shade), typ,shade,r,g,b,r,g,b);
                     */
        }
    }
    Color[0] = SDL_MapRGB(screen->format,0,0,0);           // black
    Color[1] = SDL_MapRGB(screen->format,0xff,0xff,0xff);  // white   Monomer
    Color[2] = SDL_MapRGB(screen->format,0x30,0x30,0xff);  // blue    Water
    Color[3] = SDL_MapRGB(screen->format,0xff,0x20,0x20);  // red     Oil
    Color[4] = SDL_MapRGB(screen->format,0x20,0xff,0x20);  // green      
    Color[5] = SDL_MapRGB(screen->format,0x40,0x40,0x40);  // gray
    Color[6] = SDL_MapRGB(screen->format,0xff,0xff,0x00);  // yellow
                                   
    Color[7] = SDL_MapRGB(screen->format,0x60,0x60,0xff);  // light blue  
    Color[8] = SDL_MapRGB(screen->format,0xff,0x60,0x60);  // light red
}



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int GetNum(char *str)
{
    int a;
    int b;
        
    b = (sscanf(str, "%i", &a) == 0);
    if (b != 0)
    a = 0;
    return a;
}
                        
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
float_t GetFNum(char *str)
{
    float_t a;
    int b;
        
    b = (sscanf(str, "%f", &a) == 0);
    if (b != 0)
    a = 0;
    return a;
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
  
                        
                        
/*******************************************************************/
void helpinvoke()
{
    printf("\n"
        "Program call:\n"
        "   dpdc [option]\n"
        "\n"
        "   option:\n"
        "       -mx        x=polymertypes\n"
        "       -px        x=particletypes\n"
        "       -gx        x=geometry\n"
        "       -cx        x=no of monomers (and polymers)\n"
        "       -dx     --disp=x     x=displayinterval [1]\n"
        "       -h      --help       this help\n"
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

    

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     This code is for a single homopolymer in a solvent in 2D */
/*     by dissipative particle dynamics (DPD). */
/*     We do the following (main points): */
/*     - generate the initial configuration */
/*     - equilibrate the system */
/*     - do the dynamics in equilibrium */
/*       -- calculate the temperature of the system */
/*       -- for the chain, calculate the radius of gyration squared */
/*          as well as the end-to-end distance squared */
/*       -- determine the velocity autocorrelation function which */
/*          is related to the diffusion coefficient of the chain */
/*     - close the simulation */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ 
int main(
    int argc,
    char *argv[])
{
    sys_size = 16;
    n_dim = 3;
    
    rho0 = 4;             /* total density of all particles [3.0]*/
    d_bond__ = 10.;       /* harmonic bonding btw. nn beads in a ch */
    pm = 1.;              /* particle mass (identical for all parti */
    sigma = 3.;           /* dpd parameter */
    set_dt(0.01);         /* init integration time */
    rfac = sqrt(3.);                           /* Scaling of random forces */


    iseed = 619099;       /* seed for a random number generator */
    isteps = 50;          /* number of steps to equilibrate the sys */
    maxstp = STEPS;       /* number of time steps after equilibrati */
    incstp = 250;         /* how often do we calculate quantities */
    issstp = 500;         /* how often do we take snapshots */

    iseed = iseed + time(0);  

    //testbnd();


    /* --------------------------------------------------------------- */
    if (maxstp / incstp >= 100) {
    /*io*/
    }
    
    helpinteract();

    DEBUGOUT("7");
            
    /* --------------------------------------------------------------- */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Commandline parameters                                                */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    int i;
    for (i = 1; i <= argc - 1; i++) {
    
        if (argv[i][0] == '-') {
            if (!strncmp(argv[i],"--disp=",7) && strlen(argv[i])>7 ) {
                display_intv = ((GetNum(&argv[i][7])-1) % 256)+1;
            }
            else if (!strncmp(argv[i],"--seed=",7) && strlen(argv[i])>7 ) {
                iseed = GetNum(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--help",6) && strlen(argv[i])==6 ) {
                 helpinvoke();
            }
            else if (!strncmp(argv[i],"--size=",7) && strlen(argv[i])>7 ) {
                sys_size = GetFNum(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--dt=",5) && strlen(argv[i])>5 ) {
                dt = GetFNum(&argv[i][5]);
            }
            else if (!strncmp(argv[i],"--rfac=",7) && strlen(argv[i])>7 ) {
                rfac = GetFNum(&argv[i][7]);
            }
            else if (!strncmp(argv[i],"--dim=",6) && strlen(argv[i])>6 ) {
                n_dim = GetNum(&argv[i][6]);
            }
            else if (!strncmp(argv[i],"--rho=",6) && strlen(argv[i])>6 ) {
                rho0 = GetFNum(&argv[i][6]);
            }
            else if (!strncmp(argv[i],"--sigma=",8) && strlen(argv[i])>8 ) {
                sigma = GetFNum(&argv[i][8]);
            }
            else if (!strncmp(argv[i],"--bond=",7) && strlen(argv[i])>7 ) {
                d_bond__ = GetFNum(&argv[i][7]);
            }
            else if ((argv[i][2]>='0' && argv[i][2]<='9') || (argv[i][2]==0)) {
                switch (argv[i][1]) {
            
                    case 'g':
                        geometry_type = ((GetNum(&argv[i][2])) % 4);
                        break;
            
                    case 'm':
                        mere_types = GetNum(&argv[i][2]);
                        break;
            
                    case 'p':
                        particle_types = ((GetNum(&argv[i][2])-1) % 5)+1;
                        break;
                        
                    case 'h':
                        helpinvoke();
                        break;
                 
                    case 'i':
                        interaction_types = ((GetNum(&argv[i][2])-1) % 5)+1;
                        break;

                    case 'd':
                        display_intv = ((GetNum(&argv[i][2])-1) % 256)+1;
                        break;

                    case 'c':
                        chain_len = GetNum(&argv[i][2]);
                        break;

                    case 'r':
                        rho0 = GetFNum(&argv[i][2]);
                        break;

                    case 'b':
                        d_bond__ = GetFNum(&argv[i][2]);
                        break;

                    case 's':
                        sys_size = GetFNum(&argv[i][2]);
                        break;
                        
                    default:
                        printf("Unknown commandline parameter %s\n", argv[i]);
                        stopit("Fatal error");
                }    
            }
            else {
                printf("Unknown commandline parameter %s\n", argv[i]);
                stopit("Fatal error");
            }
        }
    }    
    if (rho0 < 0.01) rho0 = 0.01; 
    if (n_dim!=2 && n_dim!=3) n_dim = 2;
    if (sys_size >= SYS_SIZEMAX) sys_size = SYS_SIZEMAX;
    if (sys_size <= 2) sys_size = 4;
    if (sigma <= 0.0) sigma = 3.0;
    if (dt <= 0) dt = 0.01;
    if (rfac <= 0) rfac = sqrt(3);
    
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Fix some constants:                                                   */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    sys_size2 = sys_size/2;
    // fill_nb3D();  // fill the neighbour offset table

    dl.k[0] = sys_size;                    /* system size in x-direction */
    dl.k[1] = sys_size; //*MAX_Y/MAX_X;    /* system size in y-direction */
    dl.k[2] = (n_dim == 2)? 1.0 : dl.k[1];

    dlinv.k[0] = 1. / dl.k[0];                /* Inverse system size */
    dlinv.k[1] = 1. / dl.k[1];
    dlinv.k[2] = 1. / dl.k[2];
    
    dli.k[0] = (int)dl.k[0];
    dli.k[1] = (int)dl.k[1];
    dli.k[2] = (int)dl.k[2];

    dldiff.k[0] = dli.k[0]*dli.k[1];          // Offset z
    dldiff.k[1] = dli.k[1];                   // Offset y
    dldiff.k[2] = 1;                          // Offset x

    printf("Size: %4.1f %4.1f %4.1f  %d %d %d = %d %d  %dD\n", 
                  dl.k[0], dl.k[1], dl.k[2], dli.k[0], dli.k[1], dli.k[2], 
                  dli.k[0]*dli.k[1]*dli.k[2], (dli.k[0]+2)*(dli.k[1]+2)*(dli.k[2]+2),
                  n_dim);
    if (F_ASIZE < (dli.k[0]+2)*(dli.k[1]+2)*(dli.k[2]+2) ) {
        stopit("Fatal Error: F_ASIZE to small!");
    }             
    printf("Sizeof(float) = %ld\n", sizeof(float_t));

    set_dt(dt);
    
    vol = dl.k[0] * dl.k[1] * dl.k[2];            /* Volume (area) of the system */

    num_particle = (int) (rho0 * vol);
    
    ALPHA(TYP_WATER, TYP_WATER)=  3.0;
    ALPHA(TYP_OIL,   TYP_OIL)  =  3.0;
    ALPHA(TYP_POLY,  TYP_POLY) =  3.0; // 0.9;
    ALPHA(TYP_POLY2, TYP_POLY2)=  3.0; //  0.9;
    
    ALPHA(TYP_WATER, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_WATER)= 10.0;
    ALPHA(TYP_WATER, TYP_POLY) = ALPHA(TYP_POLY,  TYP_WATER)= 16.0;
    ALPHA(TYP_POLY,  TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY) = 1.0; //0.2;
    ALPHA(TYP_POLY,  TYP_POLY2)= ALPHA(TYP_POLY2, TYP_POLY) = 3.0; //1.0;
    ALPHA(TYP_POLY2, TYP_WATER)= ALPHA(TYP_WATER, TYP_POLY2)= 1.0; //0.2;
    ALPHA(TYP_POLY2, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY2)= 16.0;
   
    epsilon[TYP_WATER] = 0.0;
    epsilon[TYP_OIL]   = 0.0;
    epsilon[TYP_POLY]  = 5.0;    // 50
    epsilon[TYP_POLY2] = 5.0;    //-50
    epsilonfac = 1;
    
    printf("   num_particle=%d  chn_particle=%d\n", num_particle, len_chain__);
    
    if (num_particle > F_SIZE) {
        stopit("Fatal Error: num_particle > F_SIZE");
    }
    
    dnum_particle = (float_t) num_particle;
    dnum_particle_inv = 1. / dnum_particle;           /* Inverse particle number */
    rho = dnum_particle / vol;                        /* Total density of the system */
    
    pminv = 1. / pm;                 /* Inverse particle mass */
    
    gamma_ = sigma * .5 * sigma;     /* The strength of the dissipative force */
                                     /* (using the fluctuation-dissipation the */
    skin = 2.0;                      /* Skin radius for the Verlet neighbor table */
    skinsq = skin * skin;
    skin_size = (int) (sys_size/skin);
    
    tfac = 1. / (dnum_particle * 2. - 2.);     /* Factor for temperature control */
    
    DEBUGOUT("7a");
    init_vacf();
    
    r2inis_(iseed);/* Initialize the random number generator */
    DEBUGOUT("9");
    it = 0;
    
    
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     SDL init setup                                                        */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    FULLSCREEN=0;
    WIDTH=MAX_X*PIXPERDOT;
    HEIGHT=MAX_Y*PIXPERDOT;
    BPP=16;
    SOUND=0;
    runde = 0;

    SDL_INIT_ALL();
    // set window caption
    SDL_WM_SetCaption("dpd","dpd");
    // clean screen
    ClearSurfaceRGB(screen,0,0,0);
    SDL_UpdateRect(screen,0,0,0,0);
    
    createcolors();

    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*     Time to start the actual DPD simulation                               */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /* Set time t = 0 */
    genconf();
        /* Generate the initial configuration */
    DrawIt();
    DEBUGOUT("A");
        
    /* Make sure that the Verlet table */
    update = 1;/* will be constructed */

    /* Calculate all forces for the initial configuration. */
    fcdr_(maxverlet);
        
    nvintv = 0;
    DEBUGOUT("D");
    DrawIt();
    
    /* Initialize counter for configuration io */
    nconf = 0; 

    /* Take a snapshot of the polymer chain */
    snaps_(nconf, snapsintv);

    for (it = 0; it < maxstp; ++it) {

        /* Integrate velocities over dt/2 */
        intv_(pminv, dth, dtrth); 
        
        /* Integrate particle positions over dt */
        intr_(dt);
            
        /* Check if we should update the */
        /* Verlet neighbor table */
        check_();

        if (update) {
            /* Particles have wandered outside */
            /* the simulation box, so let us shift them back inside */
            foldb_();
        }
        
        /* Calculate all forces */
        /* (conservative, dissipative, and random forces) */
        fcdr_(maxverlet);
            
        /* Integrate velocities over dt/2 */
        intv_(pminv, dth, dtrth);
            
        /* Calculate dissipative forces */
        fdr_();
            
        if (it % incstp == 0) {
            velo_prof();
            /* Take samples of the center-of-mass */
            /* velocity and the temperature of the system */
            vcmtp_(&pm, &pminv, &vcm, &temp);
            /*io */
            printf("  Step:%6d  Temp:%7.4f\n", it,temp);
            
            /* Calculate the radius of gyration squar */
            /* and the end-to-end distance squared */
            /* r_gyr__(&dr_gyr__, &dr_ee__); */
            
            /* Update the velocity autocorrelation functions that are being calculated */
            /* vacf_(max_vacf, 
                     &it_vacf__, &ivacf_counter__, it_bound__, ivacf_loop__, dvacf); */
        }

        if (it % issstp == 0) {
            /* Take a snapshot of the polymer chain to a file */
            ++nconf;
            // snaps_(nconf, snapsintv);
            savedata(0);
        }
        DrawIt();
    }
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

    /* -------------------------------------- */
    if (ivacf_counter__ >= VACF_SIZE) {
        /* -------------------------------------- */
        /* The purpose here is to calculate the */
        /* velocity autocorrelation function and */
        /* then write the final averaged result */
        for (it = 0; it < VACF_SIZE; ++it) {
            dvacf[it] = dvacf[it] / (float_t) (ivacf_counter__ - (VACF_SIZE-1));
        }

        dcs_vacf__[0] = dvacf[0] /6.0 * dt * (float_t) incstp;
    
        for (it = 1; it < VACF_SIZE; ++it) {
            dcs_vacf__[it] = dcs_vacf__[it - 1] + dvacf[it] /3.0 * dt * (float_t) incstp;
        }
        for (it = 0; it < VACF_SIZE; ++it) {
            /* Write the results */ /*io*/
        }
    }
    /*io close*/

    while (!Quit) { 
        handleevent();
        // do not let the event cycle consume the whole processor
        RestCPU();  
    }

    // SDL_SaveBMP(screen,"dpd.bmp");   
    SDL_Quit();
    
    stopit("Groovy!");/* The DPD simulation is over. */
    return 0;
} /* MAIN__ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void save_frame(FILE* fi)
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
/*     SUBROUTINE: */
/*     Generate the initial configuration */
/*     (positions and velocities). */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int genconf()
{
    /* Local variables */
    float_t tscal;
    int ip, kp;
    float_t vv, temper;
    FILE * fi;

    
    /*
          e.g.   dpd3 -p5 -i2 -m2 -g0 -b40 -c1200 --dim=2 --size=15 --rho=10 --disp=40
      
    -p    particle_types  
            1        1/1 water
            2        1/1 water mere
            3        1/2 water 1/2 oil mere
            4        2/3 water 1/3 oil mere
            5        9/10 water 1/10 oil mere
            
    -m    mere_types = mere
            1        Polymere (len_chain)
            2        Dimere (len_chain/2)
            3        Monomere (len_chain)  
            4        Trimere (len_chain)   Poly-Poly2-Poly
            5        Trimere (len_chain)   Poly-Poly2-Poly2
            6        Trimere (len_chain)   Poly2-Poly-Poly2
            7        Trimere (len_chain)   Poly2-Poly-Poly
            8        Tetramere              Poly-Poly2-Poly2-Poly
            
    -i    interaction_types
            1        standard
            2        oil-cluster mit surfactant
            
    -g    geometry_type
            0        normal (all)
            1        a rectangular hole from -z to z
            2        a H-structure
            3        bell-shape
    */

    printf("Init: particle types\n");
    if (interaction_types == 2) {
        ALPHA(TYP_WATER, TYP_WATER)=  3.0;
        ALPHA(TYP_OIL,   TYP_OIL)  =  0.2;
        
        ALPHA(TYP_POLY,  TYP_POLY) =  4.0;
        ALPHA(TYP_POLY2, TYP_POLY2)=  4.5;
        
        ALPHA(TYP_POLY,  TYP_POLY2)= ALPHA(TYP_POLY2, TYP_POLY) = 7.0;
        ALPHA(TYP_WATER, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_WATER)= 16.0;

        ALPHA(TYP_WATER, TYP_POLY) = ALPHA(TYP_POLY,  TYP_WATER)= 16.0;
        ALPHA(TYP_WATER, TYP_POLY2)= ALPHA(TYP_POLY2, TYP_WATER)= 0.2;
        
        ALPHA(TYP_POLY,  TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY) = 0.5;
        ALPHA(TYP_POLY2, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY2)= 16.0;
    }
    else if (interaction_types == 3) {
        ALPHA(TYP_WATER, TYP_WATER)=  3.0;
        ALPHA(TYP_OIL,   TYP_OIL)  =  0.2;
        
        ALPHA(TYP_POLY,  TYP_POLY) =  0.1; // 4.0;
        ALPHA(TYP_POLY2, TYP_POLY2)=  0.2; // 4.5;
        
        ALPHA(TYP_POLY,  TYP_POLY2)= ALPHA(TYP_POLY2, TYP_POLY) = 7.0;
        ALPHA(TYP_WATER, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_WATER)= 16.0;

        ALPHA(TYP_WATER, TYP_POLY) = ALPHA(TYP_POLY,  TYP_WATER)= 16.0;
        ALPHA(TYP_WATER, TYP_POLY2)= ALPHA(TYP_POLY2, TYP_WATER)= 0.2; // 2.0;
        
        ALPHA(TYP_POLY,  TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY) = 0.5;
        ALPHA(TYP_POLY2, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY2)= 16.0;
    }
    else if (interaction_types == 4) {  
        ALPHA(TYP_WATER, TYP_WATER)=  5.0;
        ALPHA(TYP_OIL,   TYP_OIL)  =  3.0;
        
        ALPHA(TYP_POLY,  TYP_POLY) =  2.5; 
        ALPHA(TYP_POLY2, TYP_POLY2)=  1.3; 
        
        ALPHA(TYP_POLY,  TYP_POLY2)= ALPHA(TYP_POLY2, TYP_POLY) = 3.0;
        ALPHA(TYP_WATER, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_WATER)= 16.0;

        ALPHA(TYP_WATER, TYP_POLY) = ALPHA(TYP_POLY,  TYP_WATER)= 16.0;
        ALPHA(TYP_WATER, TYP_POLY2)= ALPHA(TYP_POLY2, TYP_WATER)= 4.0;
        
        ALPHA(TYP_POLY,  TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY) = 2.1;
        ALPHA(TYP_POLY2, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY2)= 16.0;
        
        omegaofs = 0.0;
    }
    else if (interaction_types == 5) {  
        ALPHA(TYP_WATER, TYP_WATER)=  5; //0.1;
        ALPHA(TYP_OIL,   TYP_OIL)  =  5; //3;
        
        ALPHA(TYP_POLY,  TYP_POLY) =  20; //0.2; 
        ALPHA(TYP_POLY2, TYP_POLY2)=  5; //0.1; 
        
        ALPHA(TYP_POLY,  TYP_POLY2)= ALPHA(TYP_POLY2, TYP_POLY) = 10.0;
        ALPHA(TYP_WATER, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_WATER)= 16.0;

        ALPHA(TYP_WATER, TYP_POLY) = ALPHA(TYP_POLY,  TYP_WATER)= 20.0;  // white
        ALPHA(TYP_WATER, TYP_POLY2)= ALPHA(TYP_POLY2, TYP_WATER)= 1.0;   // green
        
        ALPHA(TYP_POLY,  TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY) = 1.1;
        ALPHA(TYP_POLY2, TYP_OIL)  = ALPHA(TYP_OIL,   TYP_POLY2)= 16.0;
        
        omegaofs = 0.0;
    }
    
    if (len_chain__<num_particle) len_chain__ = num_particle;
        
    if (particle_types == 1) {   // only water
        len_chain__ = 0;/* size of the chain(s) */
        len_mol2 = 0; 
        n_water__ = num_particle - len_chain__ - len_mol2;     /* LastNumber of water particle */
        n_poly2   = num_particle - len_chain__ - len_mol2;
        n_oil     = num_particle - len_chain__;
        n_solu =  n_oil;
    }
    else if (particle_types == 2) { // water and polymere
        len_chain__ = chain_len;/* size of the chain(s) */
        len_mol2 = 0; 
        n_water__ = num_particle - len_chain__ - len_mol2;     /* LastNumber of water particle */
        n_poly2   = num_particle - len_chain__ - len_mol2;
        n_oil     = num_particle - len_chain__;
        n_solu =  n_oil;
    }
    else if (particle_types == 3) {  // water oil and polymere
        len_chain__ = chain_len;/* size of the chain(s) */
        len_mol2 = (int)((num_particle - len_chain__)/2); 
        n_water__ = num_particle - len_chain__ - len_mol2;     
        n_poly2   = num_particle - len_chain__ - len_mol2;
        n_oil     = num_particle - len_chain__;
        n_solu =  n_oil;
    }
    else if (particle_types == 4) { 
        len_chain__ = chain_len;/* size of the chain(s) */
        len_mol2 = (int)((num_particle - len_chain__)/3); 
        n_water__ = num_particle - len_chain__ - len_mol2 - len_mol2;    
        n_poly2   = num_particle - len_chain__ - len_mol2;
        n_oil     = num_particle - len_chain__;
        n_solu =  n_oil;
    }
    else if (particle_types == 5) { 
        len_chain__ = chain_len;/* size of the chain(s) */
        len_mol2 = (int)((num_particle - len_chain__)/10); 
        n_water__ = num_particle - len_chain__ - len_mol2;    
        n_poly2   = n_water__;
        n_oil     = num_particle - len_chain__;
        n_solu =  n_oil;
    }
    else {
        stopit("Wrong no of mols");
    }
    
    for (ip = 0; ip < num_particle; ++ip) {
        liste[ip].left = -1;
        liste[ip].right = -1;
        liste[ip].drawn = 0;
        liste[ip].type = TYP_WATER;
    }
    
    
    for (ip = 0; ip < n_water__; ++ip) {           /* Fix the particle types. */
        liste[ip].type = TYP_WATER;               /* Type '1' corresponds to solvent (water */
    }                                               /* in a polymer chain. */
    for (ip = n_water__; ip < n_poly2; ++ip) {           
        liste[ip].type = TYP_WATER; // TYP_POLY2;               
    }  
    for (ip = n_poly2; ip < n_oil; ++ip) {           
        liste[ip].type = TYP_OIL;               
    }  
    for (ip = n_oil; ip < num_particle; ++ip) {
        if (ip>n_oil) {
            if (mere_types <= 1) {   
                liste[ip].type = TYP_POLY;
                liste[ip].left = ip-1;   // Polymere aus (POLY)
                liste[ip-1].right = ip;
            }
            else if (mere_types == 2) {   // Dimere  (POLY--POLY2)
                if (ip % 2 == 0) {
                    liste[ip].type = TYP_POLY;
                }
                else  if (ip % 2 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
            }
            else if (mere_types == 3) {   // Monomere  (POLY+, POLY2-)
                if (ip % 2 == 0) {
                    liste[ip].type = TYP_POLY;
                }
                else  if (ip % 2 == 1) {
                    //liste[ip].left = ip-1;
                    //liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
                epsilon[TYP_POLY2] = -epsilon[TYP_POLY]; 
            }
            else if (mere_types == 4) {   // Trimere  (POLY--POLY2-POLY)
                if (ip % 3 == 0) {
                    liste[ip].type = TYP_POLY;
                }
                else if (ip % 3 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
                else  if (ip % 3 == 2) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
            }
            else if (mere_types == 5) {   // Trimere  (POLY--POLY2-POLY2)
                if (ip % 3 == 0) {
                    liste[ip].type = TYP_POLY;
                }
                else if (ip % 3 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
                else  if (ip % 3 == 2) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
            }
            else if (mere_types == 6) {   // Trimere  (POLY2--POLY-POLY2)
                if (ip % 3 == 0) {
                    liste[ip].type = TYP_POLY2;
                }
                else if (ip % 3 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
                else  if (ip % 3 == 2) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
            }
            else if (mere_types == 7) {   // Trimere  (POLY2--POLY-POLY)
                if (ip % 3 == 0) {
                    liste[ip].type = TYP_POLY2;
                }
                else if (ip % 3 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
                else  if (ip % 3 == 2) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
            }
            else if (mere_types == 8) {   // Tetramere  (POLY-POLY2-POLY2-POLY)
                if (ip % 4 == 0) {
                    liste[ip].type = TYP_POLY;
                }
                else if (ip % 4 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
                else  if (ip % 4 == 2) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
                else  if (ip % 4 == 3) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
            }
            else if (mere_types >= 9) {   // Tetramere  (POLY2-POLY-POLY-POLY2)
                if (ip % 4 == 0) {
                    liste[ip].type = TYP_POLY2;
                }
                else if (ip % 4 == 1) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
                else  if (ip % 4 == 2) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY;
                }
                else  if (ip % 4 == 3) {
                    liste[ip].left = ip-1;
                    liste[ip-1].right = ip;
                    liste[ip].type = TYP_POLY2;
                }
            }
        }    
    }
    
    if ((fi = fopen("vertices.out", "w")) == NULL) {
        stopit("Fatal error: cannot open vertices.out");
        exit(1);
    }
    
    printf("Init: geometry\n");
    if (geometry_type == 0) {         // all allowed
        printf("      geometry 0\n");
        for (ip = 0; ip < F_ASIZE; ++ip) {
            allowed[ip] = 1;
        }
        save_frame(fi);
    }
    else if (geometry_type == 1) {  
        printf("      geometry 1\n");  // a hole
        ipoint_t p;
        for (ip = 0; ip < F_ASIZE; ++ip) {
            allowed[ip] = 1;
        }
        save_frame(fi);
        for (p.k[0]=-2; p.k[0]<2; p.k[0]++) {
            for (p.k[1]=-2; p.k[1]<1; p.k[1]++) {
                for (p.k[2] = -dli.k[2]/2-1; p.k[2] <= dli.k[2]/2+1; p.k[2]++) {
                    ALLOWED(p) = 0;
                }
            }
        }
    }
    else if (geometry_type == 2) {
        printf("      geometry 2\n");  // a H-structure
        ipoint_t p;
        for (ip = 0; ip < F_ASIZE; ++ip) {
            allowed[ip] = 0;
        }
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
        // Elektroden
        pol1.k[0] = 0.0; pol1.k[1] = 1.0; pol1.k[2] = -3.0;
        pol2.k[0] = 0.0; pol2.k[1] =-1.0; pol2.k[2] = -3.0;
        emin = 0.5; 
        // Druck
        pressure_to_x_minpos = -dl.k[0]/2+1;           // accelerate particles between these two x-limits 
        pressure_to_x_maxpos = -dl.k[0]/2+1+0.1;       //
        pressure_diff_at_y = 0.0;                  // above this only pressure_ofs/2 applied (diff pressure for H-Structures)
    }
    else if (geometry_type == 3) {         // bell shape
        printf("      geometry 3\n");
        for (ip = 0; ip < F_ASIZE; ++ip) {
            allowed[ip] = 1;
        }
        save_frame(fi);
    }
    
    fclose(fi); 
    
    
    /* Generate random initial positions */
    printf("Init: random positions\n");
    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].r.k[kp]  = 0.0;
        }                                           
    }

    for (ip = 0; ip < num_particle; ++ip) {
    
        /* for the particles */
        do {
            if (geometry_type == 3) {
                liste[ip].r.k[0] = (r2s_() - 0.5) * dl.k[0]/3;
                liste[ip].r.k[1] = (r2s_() - 0.5) * dl.k[1]/3;
                liste[ip].r.k[2] = (r2s_()) * dl.k[2]/2/3;
            }
            else {
                if (liste[ip].left >= 0) {
                    float_t d_scale, r_dist;
                    fpoint_t dr_tmp;
                    int ip2 = liste[ip].left;
                    
                    r_dist = 0;
                    for (kp = 0; kp < n_dim; kp++) {
                        dr_tmp.k[kp] = r2s_() - 0.5;
                        r_dist += dr_tmp.k[kp]*dr_tmp.k[kp];
                    }

                    r_dist = sqrt(r_dist);
                    d_scale = 0.5 / r_dist;

                    for (kp = 0; kp < n_dim; kp++) {
                        liste[ip].r.k[kp]  = liste[ip2].r.k[kp] + dr_tmp.k[kp] * d_scale;
                    }
                }
                else {
                        for (kp = 0; kp < n_dim; kp++) {
                            liste[ip].r.k[kp] = (r2s_() - 0.5) * dl.k[kp];
                        }
                }
            }            
        } while (! ALLOWED(liste[ip].r) );
        if (ip % 1000 == 0) { printf("  i = %d\n", ip); }
    }

    /* Make sure that all particles are */
    /* inside the simulation box */
    foldb_();
    
    printf("Init: random velocities\n");

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
    
    
    /* Satisfy the condition that the */
    /* center-of-mass velocity vanishes */

    for (kp = 0; kp < n_dim; kp++) {
        vcm.k[kp] = vcm.k[kp] * dnum_particle_inv * pminv; 
    }

    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].v.k[kp] = liste[ip].v.k[kp] - vcm.k[kp];
        }
    }
    
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
    printf("Init: done!\n");
    return 0;
} /* genconf */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Shift (fold) particles back into primary simulation box. */
/*     NOTE: Calls of FOLDB are allowed only when the Verlet list */
/*     has to be reconstructed. */
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
/*     SUBROUTINE: */
/*     Integrate velocities using a time step of "hcd" for */
/*     conservative and dissipative forces, and a time step */
/*     of "hr" for random forces. The routine integrates */
/*     velocities over a time increment of dt/2, and therefore */
/*     integrates velocities in two parts. This approach is based */
/*     on the DPD-VV integration scheme (see the reference: */
/*     [G. Besold, I. Vattulainen, M. Karttunen, and J.M. Polson, */
/*     Phys. Rev. E Rapid Comm. vol. 62, R7611 (2000)]). */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int intv_(
    float_t pminv_, 
    float_t hcd_,     // dt/2
    float_t hr_)      // sqrt(dt)/2
{
    int ip, kp;

    for (ip = 0; ip < num_particle; ++ip) {
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
/*     SUBROUTINE: */
/*     Integrate the positions of the particles. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int intr_(
    float_t dt_)  // dt
{
    int ip,kp;
    const float_t t=0.2;
    
    static int first = 1;
    static float_t H, B, C;
    
    if (first) {   //   for the bell-shape
        first = 0;
        H = dl.k[2]/10*5;  // Z
        B = 0.01 * pow(M_PI,4.0/3) / pow(H,4.0/3);
        // B = 0.03*10/dl.k[1];  // Y X 
        C = (4*dl.k[1])*(4*dl.k[1]);
    }

    for (ip = 0; ip < num_particle; ++ip) {

        fpoint_t ri = liste[ip].r;
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].r.k[kp] = ri.k[kp] + liste[ip].v.k[kp] * dt_;
        }
        
        if (geometry_type == 3) {  //  for the bell-shape
            int ok = 0;
            if (liste[ip].r.k[2] >= 0.0) {
                float_t c2 = liste[ip].r.k[0]*liste[ip].r.k[0] + liste[ip].r.k[1]*liste[ip].r.k[1];
                if (c2 < C) {
                    float_t dr = H/(1+B*c2*c2);
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
            if (  ((liste[ip].r.k[0]>brick_xl)) 
                &&((liste[ip].r.k[0]<brick_xr))
                &&((liste[ip].r.k[1]<brick_yo))
                &&((liste[ip].r.k[1]>brick_yu)) ) {  // inside
                
                // and comes from
                if (ri.k[0] < brick_xl) {  // left
                    if (ri.k[1] > brick_yo) { // upper left
                        liste[ip].r.k[0] = brick_xl - (liste[ip].r.k[0] - brick_xl);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else if (ri.k[1] < brick_yu) { // lower left
                        liste[ip].r.k[0] = brick_xl - (liste[ip].r.k[0] - brick_xl);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else { // left
                        liste[ip].r.k[0] = brick_xl - (liste[ip].r.k[0] - brick_xl);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        if (wall_brakes) liste[ip].v.k[1] = 0;
                    }
                }
                else if (ri.k[0] > brick_xr) {  // right
                    if (ri.k[1] > brick_yo) { // upper right
                        liste[ip].r.k[0] = brick_xr + (brick_xr - liste[ip].r.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else if (ri.k[1] < brick_yu) { // lower right
                        liste[ip].r.k[0] = brick_xr + (brick_xr - liste[ip].r.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                    }
                    else  { // right
                        liste[ip].r.k[0] = brick_xr + (brick_xr - liste[ip].r.k[0]);
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                        if (wall_brakes) liste[ip].v.k[1] = 0;
                    }
                }
                else if (ri.k[1] > brick_yo) {  // top
                    if (ri.k[0] > brick_xr) { // top right
                        liste[ip].r.k[1] = brick_yo + (brick_yo - liste[ip].r.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }
                    else if (ri.k[0] < brick_xl) { // top left
                        liste[ip].r.k[1] = brick_yo + (brick_yo - liste[ip].r.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }                        
                    else { // top
                        liste[ip].r.k[1] = brick_yo + (brick_yo - liste[ip].r.k[1]);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        if (wall_brakes) liste[ip].v.k[0] = 0;
                    }
                }
                else if (ri.k[1] < brick_yu) {  // bottom
                    if (ri.k[0] > brick_xr) { // bottom right
                        liste[ip].r.k[1] = brick_yu - (liste[ip].r.k[1] - brick_yu);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }
                    else if (ri.k[0] < brick_xl) { // bottom left
                        liste[ip].r.k[1] = brick_yu - (liste[ip].r.k[1] - brick_yu);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        liste[ip].v.k[0] = - liste[ip].v.k[0];
                    }
                    else { // bottom left
                        liste[ip].r.k[1] = brick_yu - (liste[ip].r.k[1] - brick_yu);
                        liste[ip].v.k[1] = - liste[ip].v.k[1];
                        if (wall_brakes) liste[ip].v.k[0] = 0;
                    }
                }
            }
        }
        if (wall_at_y0) {
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
    }
    return 0;
} /* intr_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
// cell stack and chain handling
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void fill_cell_chains(int num_particle)
{
    int i,ci,ch, kp;

    for (i=0;i<F_ASIZE;i++) chainarray[i]=-1;
    
    for (i=0;i<num_particle;i++) {
        // calculate cell coords

        fpoint_t rr = liste[i].r;
        ipoint_t ii = {{0,0,0}};
        
        for (kp = 0; kp < n_dim; kp++) {
           rr.k[kp] -= Boundary(rr,kp);
           ii.k[kp] = (int)((rr.k[kp] + sys_size2+skin)/skin) % skin_size;  // +skin
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


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Calculate all forces: conservative, dissipative, and random */
/*     forces. */
/*     If needed (the int variable update = .true.), */
/*     then update the Verlet neighbor table. Otherwise use */
/*     the present table. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int fcdr_(
    int maxnab)
{
    /* Local variables */
    float_t  fcfac, fdfac, frfac;
    static int jpbeg;
    static int jpnab, jpend;
    int nlist, ip, jp, kp;
    float_t omega;
    float_t rijinv,rijsq, rrij;
    int nfi, nfj;
    fpoint_t ri, vi, fci, fdi, fri, eij, vij, rij, fcij, fdij, frij;
    int chaincell,nb,nbc;               // indices for cell chaining

    for (ip = 0; ip < num_particle; ++ip) {
        /* Initialize all forces */
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].fc.k[kp] = 0.;
            liste[ip].fd.k[kp] = 0.;
            liste[ip].fr.k[kp] = 0.;
        }
    }
    
    // printf(" %f %f %f %d\n", liste[0].r.k[0], liste[0].r.k[1],liste[0].r.k[2], liste[0].type);
    
    if (update) {
        /* Check if the Verlet list should */
        /* be updated */
        nvintv = 1;
        
        fill_cell_chains(num_particle);
        
        for (ip = 0; ip < num_particle; ++ip) {
            /* Ok, let us do it. */
            /* Save the particle positions. */
            liste[ip].r0 = liste[ip].r;
        }
        
        nlist = 0;
        for (ip = 0; ip < num_particle-1; ++ip) {
            /* Update the Verlet list and calculate */
            /* all forces */
            npoint[ip] = nlist + 1;
            
            ri = liste[ip].r;         /* Position of particle "ip"  */
            vi = liste[ip].v;         /* Velocity of particle "ip" */
            fci = liste[ip].fc;       /* Conservative forces acting on particle */
            fdi = liste[ip].fd;       /* Dissipative forces */
            fri = liste[ip].fr;       /* Random forces */
            nfi = liste[ip].type;     /* Type of particle "ip" */
        
            fpoint_t rr = liste[ip].r;
            ipoint_t ii = {{0,0,0}};

            for (kp = 0; kp < n_dim; kp++) {
               rr.k[kp] -= Boundary(rr,kp);
               ii.k[kp] = (int)(rr.k[kp] + sys_size2+skin)/skin;  // +skin
            } 

            for(nb=0; nb < (n_dim==2 ? 9:27); nb++) {
                int iix = (ii.k[0] + nb3Dx[nb]) % skin_size;
                int iiy = (ii.k[1] + nb3Dy[nb]) % skin_size;
                int iiz = (n_dim==2)? 0: (ii.k[2] + nb3Dz[nb]) % skin_size;

                nbc=(iiz*skin_size+iiy)*skin_size+iix;
                
                // error control
                if (nbc>=F_ASIZE)       {  
                     printf("Error GT in indexing  nb %d nbc %d  ip=%d\n",nb,nbc,ip);
                     break; 
                }
                else if (nbc<0) {  
                     printf("Error LT in indexing  nb %d nbc %d  ip=%d\n",nb,nbc,ip);
                     break;
                }
                chaincell=chainarray[nbc]; 
                //printf("%d:%4d(%d,%d,%d) ", nb, nbc,iix,iiy,iiz); 

                while (chaincell!=-1) {
                    jp=chaincell; 
                    chaincell=liste[chaincell].nextindex;
                                                                                                                                                                                                                             
                    if(jp>ip) {                                                                                                                                                                                                    
                    //for (jp = ip + 1; jp <= num_particle; ++jp) {
                        /* Go over all pairs within the Verlet list */

                        /* Vector from "jp" to "ip" */
                        rijsq = 0;                
                        for (kp = 0; kp < n_dim; kp++) {
                            rij.k[kp] = ri.k[kp] - liste[jp].r.k[kp];
                            rij.k[kp] -= Boundary(rij,kp);
                            rijsq += rij.k[kp]*rij.k[kp];
                        }
                        
                        
                        if (rijsq < skinsq) {
                            /* The particles are neigh */
                            ++nlist;
                            list_verl[nlist - 1] = jp;
                            if (nlist == maxnab) {
                                stopit("verlet-list too small");
                            }
                            if (rijsq < 1.) {
                                /* The two particles interaction */
                                
                                rrij = sqrt(rijsq);   /* Distance between the part */
                                rijinv = 1. / rrij;
                                omega = 1. - rrij;
                                
                                nfj = liste[jp].type;  /* Type of particle "jp"  */

                                fcfac = omega  - omegaofs;
                                fdfac = 0;
                                for (kp = 0; kp < n_dim; kp++) {
                                    eij.k[kp] = rij.k[kp] * rijinv;
                                
                                    /* Relative velocity */
                                    vij.k[kp] = vi.k[kp] - liste[jp].v.k[kp];
                                
                                    /* Calculate the pairwise conservative */
                                    fcij.k[kp] = fcfac * eij.k[kp] * ALPHA(nfi,nfj);
                                    
                                    fdfac += eij.k[kp] * vij.k[kp];
                                }
                                
                                fdfac = omega * omega * fdfac;
                                frfac = omega * rfac * (r2s_() * 2.0 - 1.0);
                                
                                for (kp = 0; kp < n_dim; kp++) {
                                    /* Calculate the pairwise dissipative force */
                                    fdij.k[kp] = fdfac * eij.k[kp];

                                    /* Calculate the pairwise random force */
                                    frij.k[kp] = frfac * eij.k[kp];
                    
                                    /* Update forces */
                                    fci.k[kp] += fcij.k[kp];
                                    fdi.k[kp] += fdij.k[kp];
                                    fri.k[kp] += frij.k[kp];

                                    liste[jp].fc.k[kp] = liste[jp].fc.k[kp] - fcij.k[kp];
                                    liste[jp].fd.k[kp] = liste[jp].fd.k[kp] - fdij.k[kp];
                                    liste[jp].fr.k[kp] = liste[jp].fr.k[kp] - frij.k[kp];
                                }
                            }
                        }
                    }
                }
            }
            //printf("\n");    
            for (kp = 0; kp < n_dim; kp++) {
                liste[ip].fc.k[kp] = fci.k[kp];
                liste[ip].fd.k[kp] = fdi.k[kp];
                liste[ip].fr.k[kp] = fri.k[kp];
            }
        }
        npoint[num_particle - 1] = nlist + 1;
        nvsize = nlist + 1;
    } 
    else {
        /* There is no need to update the Verlet */
        /* so let us use the current list and calc all forces */
        ++(nvintv);
        
        for (ip = 0; ip < num_particle-1; ++ip) {
            jpbeg = npoint[ip];
            jpend = npoint[ip+1] - 1;

            if (jpbeg <= jpend) {
                ri = liste[ip].r;
                vi = liste[ip].v;
                fci = liste[ip].fc;
                fdi = liste[ip].fd;
                fri = liste[ip].fr;
                nfi = liste[ip].type;
                
                for (jpnab = jpbeg; jpnab <= jpend; ++jpnab) {
                    jp = list_verl[jpnab - 1];
                    
                    /* Vector from "jp" to "ip" */
                    rijsq = 0;                
                    for (kp = 0; kp < n_dim; kp++) {
                        rij.k[kp] = ri.k[kp] - liste[jp].r.k[kp];
                        rij.k[kp] -= Boundary(rij,kp);
                        rijsq += rij.k[kp]*rij.k[kp];
                    }
                    
                    if (rijsq < 1.) {
                        rrij = sqrt(rijsq); /* Distance between two particles */
                        rijinv = 1. / rrij;
                        omega = 1. - rrij;

                        /* Type of particle "jp" ('1' or '2 */
                        nfj = liste[jp].type;

                        fcfac = omega - omegaofs;
                        fdfac = 0;
                        for (kp = 0; kp < n_dim; kp++) {
                            eij.k[kp] = rij.k[kp] * rijinv;
                        
                            /* Relative velocity */
                            vij.k[kp] = vi.k[kp] - liste[jp].v.k[kp];
                        
                            /* Calculate the pairwise conservative */
                            fcij.k[kp] = fcfac * eij.k[kp] * ALPHA(nfi,nfj);
                            
                            fdfac += eij.k[kp] * vij.k[kp];
                        }
                        
                        fdfac = omega * omega * fdfac;
                        frfac = omega * rfac * (r2s_() * 2.0 - 1.0);
                        
                        for (kp = 0; kp < n_dim; kp++) {
                            /* Calculate the pairwise dissipative for */
                            fdij.k[kp] = fdfac * eij.k[kp];

                            /* Calculate the pairwise random force */
                            frij.k[kp] = frfac * eij.k[kp];
            
                            /* Update forces */
                            fci.k[kp] += fcij.k[kp];
                            fdi.k[kp] += fdij.k[kp];
                            fri.k[kp] += frij.k[kp];

                            liste[jp].fc.k[kp] = liste[jp].fc.k[kp] - fcij.k[kp];
                            liste[jp].fd.k[kp] = liste[jp].fd.k[kp] - fdij.k[kp];
                            liste[jp].fr.k[kp] = liste[jp].fr.k[kp] - frij.k[kp];
                        }
                    }
                }
                for (kp = 0; kp < n_dim; kp++) {
                    liste[ip].fc.k[kp] = fci.k[kp];
                    liste[ip].fd.k[kp] = fdi.k[kp];
                    liste[ip].fr.k[kp] = fri.k[kp];
                }
            }
        }
    }
    for (ip = 0; ip < num_particle; ++ip) {
        /* Multiply remaining forces with prefactors */
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].fd.k[kp] = -(gamma_) * liste[ip].fd.k[kp];
            liste[ip].fr.k[kp] = sigma * liste[ip].fr.k[kp];
        }
    }

    /* Calculate conservative forces due to */
    /* harmonic springs between adjacent mono */
    
    for (jp = 0; jp < num_particle; ++jp) {
        if (liste[jp].right >= 0) {
            int jp1 = liste[jp].right;
            
            for (kp = 0; kp < n_dim; kp++) {
                rij.k[kp] = liste[jp].r.k[kp] - liste[jp1].r.k[kp];
                rij.k[kp] -= Boundary(rij,kp);

                liste[jp].fc.k[kp]   = liste[jp].fc.k[kp]  - d_bond__ * rij.k[kp];
                liste[jp1].fc.k[kp]  = liste[jp1].fc.k[kp] + d_bond__ * rij.k[kp];
            }
        }
    }

    if (efeld) {
        for (jp = 0; jp < num_particle; ++jp) {
            int typ = liste[jp].type;
            if (epsilon[typ] != 0.0) {
    
                rijsq = 0;            
                for (kp = 0; kp < n_dim; kp++) {
                    rij.k[kp] = liste[jp].r.k[kp] - pol1.k[kp];
                    rij.k[kp] -= Boundary(rij,kp);
                    rijsq += rij.k[kp]*rij.k[kp];
                }
                rrij = sqrt(rijsq);
                if (rrij > emin) {
                    rijinv = 1. / rrij;
                    omega = polarity * epsilonfac * epsilon[typ] * rijinv;

                    for (kp = 0; kp < n_dim; kp++) {
                        eij.k[kp] = rij.k[kp] * rijinv;
                        liste[jp].fc.k[kp]   = liste[jp].fc.k[kp] + eij.k[kp] * omega;
                    }
                }
                
                rijsq = 0;            
                for (kp = 0; kp < n_dim; kp++) {
                    rij.k[kp] = liste[jp].r.k[kp] - pol2.k[kp];
                    rij.k[kp] -= Boundary(rij,kp);
                    rijsq += rij.k[kp]*rij.k[kp];
                }
                rrij = sqrt(rijsq);
                if (rrij > emin) {
                    rijinv = 1. / rrij;
                    omega = -polarity * epsilonfac * epsilon[typ] * rijinv;

                    for (kp = 0; kp < n_dim; kp++) {
                        eij.k[kp] = rij.k[kp] * rijinv;
                        liste[jp].fc.k[kp]   = liste[jp].fc.k[kp] + eij.k[kp] * omega;
                    }
                }
            }
        }
    }
    return 0;
} /* fcdr_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Calculate dissipative forces only. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int fdr_()
{
    /* Local variables */
    fpoint_t  ri,rij,vi,fdi,fdij,eij,vij;
    float_t  fdfac;
    float_t omega;
    int jpbeg;
    int jpnab, jpend;
    int ip, jp, kp;
    float_t rijsq,rijinv, rrij;

    for (ip = 0; ip < num_particle; ++ip) {
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].fd.k[kp] = 0.;
        }
    }
    
    for (ip = 0; ip < num_particle-1; ++ip) {
        jpbeg = npoint[ip];
        jpend = npoint[ip+1] - 1;

        if (jpbeg <= jpend) {
            for (kp = 0; kp < n_dim; kp++) {
                ri.k[kp] = liste[ip].r.k[kp];
                vi.k[kp] = liste[ip].v.k[kp];
                fdi.k[kp] = liste[ip].fd.k[kp];
            }    
            for (jpnab = jpbeg; jpnab <= jpend; ++jpnab) {
                jp = list_verl[jpnab - 1];
                rijsq = 0;
                for (kp = 0; kp < n_dim; kp++) {
                    rij.k[kp] = ri.k[kp] - liste[jp].r.k[kp];
                    rij.k[kp] -= Boundary(rij,kp);
                    rijsq += rij.k[kp]*rij.k[kp]; 
                }
                if (rijsq < 1.) {
                    rrij = sqrt(rijsq);
                    rijinv = 1. / rrij;
                    omega = 1. - rrij - omegaofs;

                    fdfac = 0;
                    for (kp = 0; kp < n_dim; kp++) {
                        eij.k[kp] = rij.k[kp] * rijinv;
                        vij.k[kp] = vi.k[kp] - liste[jp].v.k[kp];
                        fdfac += eij.k[kp] * vij.k[kp];
                    }
                    fdfac = omega * omega * fdfac;

                    for (kp = 0; kp < n_dim; kp++) {
                        fdij.k[kp] = fdfac * eij.k[kp];
                        fdi.k[kp] += fdij.k[kp];
                        liste[jp].fd.k[kp] = liste[jp].fd.k[kp] - fdij.k[kp];
                    }
                }
            }
            for (kp = 0; kp < n_dim; kp++) {
                liste[ip].fd.k[kp] = fdi.k[kp];
            }
        }
    }
    
    for (ip = 0; ip < num_particle; ++ip) {
        /* Multiply forces with prefactors */
        for (kp = 0; kp < n_dim; kp++) {
            liste[ip].fd.k[kp] = -(gamma_) * liste[ip].fd.k[kp];
        }
    }
    return 0;
} /* fdr_ */



/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     SUBROUTINE: */
/*     Check if the Verlet neighbor table should be updated. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
int check_()
{
    int ip,kp;
    float_t dispmx = 0;

    for (ip = 0; ip < num_particle; ++ip) {
        /* Computing MAX */
        for (kp = 0; kp < n_dim; kp++) {
            dispmx = max(abs(liste[ip].r.k[kp] - liste[ip].r0.k[kp]),dispmx);
        }
    }
    dispmx = dispmx * sqrt(3.0) * 2.0;
    update = dispmx > skin - 1.;
    return 0;
} /* check_ */




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
    for (it = 1; it <= VACF_SIZE; ++it) {  /* the tracer diffusion coefficient */
        ivacf_loop__[it + VACF_SIZE] = it;  /* of the polymer motion. */
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
        if ((i>=0) && (i<VELO_ANZ)) {
            velo[i] += vi.k[0];
            veloanz[i] ++;
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
    float_t *pm, 
    float_t *pminv, 
    float_t *vcm, 
    float_t *temp)
{
    /* Local variables */
    static int ip,kp;
    fpoint_t vi, vcml = {{0,0,0}};
    float_t tmp;

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
    *vcm = dnum_particle_inv * *pminv * sqrt(tmp);
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
    float_t *dr_gyr__, 
    float_t *dr_ee__)
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
        for (kp=0; kp<n_dim; kp++) {
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
    float_t *dvacf)
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



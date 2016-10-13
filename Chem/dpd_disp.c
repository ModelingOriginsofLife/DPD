/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 20.07.2004  17:10h
  $Id: dpd_disp.c,v 1.22 2004/09/29 12:38:27 tmaeke Exp $
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

#include "SDL.h"
#include "SDL_etc.h"
#include "SDL_prim.h"

#include "dpd_vars.h"
#include "dpd_util.h"
#include "dpd_stat.h"
#include "dpd_file.h"
#include "dpd_disp.h"
#include "dpd_elec.h"
#include "dpd_geom.h"

/*Defines***********************************************************/
// graphics
//#define SHOW_ALL_BONDS	1	// controls wheher long bonds are drawn, looks crap if they are

#define MAX_X  256     // pixel a pixperdot size
#define MAX_Y  256     // MAX_X
#define PIXPERDOT 3
#define SDLRAND   0

#define MAXCOLOR 256

#define TYPOFS    10    // in color table
#define TYPSHADES 30    // number of shades per color per molecule=FClType
#define ZSHADES (TYPSHADES-16)
#define RGBSHADEFAC 6
#define TYP2COL(typ,shade) ((((typ)*TYPSHADES+(shade)) + TYPOFS)% MAXCOLOR)    

/*Types*************************************************************/

/*Vars**************************************************************/

    static bool_t xchg_yz = 0;                           // Display: exchange y with z 

    bool_t slice_mode = 0;                         // 0=all 1=slice from min_- to max_slice
    real_t min_slice = 0.0;
    real_t max_slice = 1.0;
    real_t delta_slice = 1.0;

    int snapsintv = 50;                     // snapshot interval
    char snapsfn[MAX_NAME] = "dpd_%06d.bmp"; // filename for snapshots
    int snaps_stp = 0;
    int snaps_start = 0;
    int snaps_end = 0;


    static Uint32 Color[MAXCOLOR];                       // color table
    static bool_t disp_pause = 0;
    static bool_t ena_disp[10] = {1,1,1,1,1,1,1,1,1,1};  // switch on/off displayed elements
   
    static int ofsx = 0;   
    static int ofsy = 0;
    static int ofsz = 0;
                
    static float sdl_xfac;
    static float sdl_yfac;
    static float sdl_zfac;
    static int   sdl_xofs;
    static int   sdl_yofs;
    static int   sdl_zofs;
    static int   sdl_xmod;
    static int   sdl_ymod;
    static int   sdl_zmod; 
    


/*Funcs*************************************************************/

/*******************************************************************/
#define toX(val) (int)(sdl_xfac *  (val) + sdl_xofs)
#define toY(val) (int)(sdl_yfac * -(val) + sdl_yofs)
#define toZ(val) (int)(sdl_zfac *  (val) + sdl_zofs)    // 0..100% 
    
#define otoX(val) ((toX(val)+ofsx)%sdl_xmod)
#define otoY(val) ((toY(val)+ofsy)%sdl_ymod)
#define otoZ(val) ((toZ(val)+ofsz)%sdl_zmod)


/*******************************************************************/
void snapshot()
{   
    if (snaps_stp>0  && (no_calc % snaps_stp == 0)) {
        if ((snaps_start != 0) && (no_calc<snaps_start)) return;
        if ((snaps_end != 0)   && (no_calc>snaps_end)) return;
        
/*        #ifdef __APPLE__
            SDL_SaveBMP(screen,snapsfn); 
            printf("Saved snapshot to %s\n", snapsfn);
        #else
				*/
            char fn[256], cmd[256];
            char tmpfile[] = "dpd.bmp";
            sprintf(fn, snapsfn, no_calc);
            sprintf(cmd, "convert %s %s", tmpfile, fn);
            SDL_SaveBMP(screen,tmpfile); 
            system (cmd);
            printf("Saved snapshot to %s and %s\n", tmpfile, fn);
            // convert -delay 20 dpd_*.png dpd.gif
       // #endif
        }
}

/*******************************************************************/
void drawline_ofs(
    real_t x0,
    real_t y0,
    real_t x1, 
    real_t y1,
    int colorno,
    Uint8 blend) 
{
    SDL_drawLine_Alpha(screen, otoX(x0), otoY(y0),  otoX(x1), otoY(y1), Color[colorno], blend); // Color[5], 255);
}    


/*******************************************************************/
void drawrect_ofs(
    real_t x0,
    real_t y0,
    real_t x1, 
    real_t y1,
    int colorno,
    Uint8 blend) 
{
    SDL_drawRect_Alpha(screen, otoX(x0), otoY(y0),  otoX(x1), otoY(y1), Color[colorno], blend);
}


/*******************************************************************/
void DrawIt()
{
    
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
        x0 = 300;
        y0 = 300;
    }
    
    //no_calc++;
    //if (++no_calc % display_intv == 0) { 

        ClearSurfaceRGB(screen,0,0,0);
        SDL_LockSurface(screen);

        // draw square of size 1.0x1.0
        SDL_drawRect_Alpha(screen, x0, y0, x0-dx, y0+dy, Color[5], 228);

        // cross in the center
        SDL_drawLine_Alpha(screen, otoX(0), otoY(-dl.k[1]/4), otoX(0), otoY(dl.k[1]/4), Color[5], 228);
        SDL_drawLine_Alpha(screen, otoX(-dl.k[0]/4),otoY(0),  otoX(dl.k[0]/4), otoY(0), Color[5], 228);
        
        
        if (with_brick) {
            SDL_drawRect_Alpha(screen, otoX(brick_ul.k[0]), otoY(brick_ul.k[1]), 
                                       otoX(brick_or.k[0]), otoY(brick_or.k[1]), Color[6], 150);
        }

        if (ena_disp[9]) { draw_polys(); }
                
        if (efeld) {  draw_elecs(); }
                
        int n;
        for (n=0; n<num_particle; n++) {
            if ((liste[n].inside)  && (liste[n].drawn == drawn)) {
                int col, SX,SY, SZ;

                if (ena_disp[liste[n].type]) {
                    if (slice_mode) {
                        if ((liste[n].r.k[2] >= min_slice) && (liste[n].r.k[2] <= max_slice)) {
                            SZ = otoZ(liste[n].r.k[2]);
                            col = Color[TYP2COL(liste[n].type, SZ)];
                            SX = otoX(liste[n].r.k[0]);
                            SY = otoY(liste[n].r.k[1]);

                            if (liste[n].left >= 0) {
								int bond_col = (liste[liste[n].left].type < liste[n].type) ? Color[TYP2COL(liste[liste[n].left].type, otoZ(liste[liste[n].left].r.k[2]))] : col;
                                int lastxx = otoX(liste[liste[n].left].r.k[0]);
                                int lastyy = otoY(liste[liste[n].left].r.k[1]);
// remember that bonds get drawn in three different pieces of this code depending on how the cookie crumbles
#ifdef SHOW_ALL_BONDS
								if ((abs(SX-lastxx)<MAX_X) && (abs(SY-lastyy)<MAX_Y))
#endif
								{
                                    SDL_drawLine_Alpha(screen,SX,SY,lastxx,lastyy, bond_col, 125);
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
                    }
                    else {
                        if (xchg_yz) {
                            SZ = otoZ(liste[n].r.k[1]);
                            col = Color[TYP2COL(liste[n].type, SZ)];
                            SX = otoX(liste[n].r.k[0]);
                            SY = otoY(liste[n].r.k[2]);

                            if (liste[n].left >= 0) {
								int bond_col = (liste[liste[n].left].type < liste[n].type) ? Color[TYP2COL(liste[liste[n].left].type, otoZ(liste[liste[n].left].r.k[2]))] : col;
                                int lastxx = otoX(liste[liste[n].left].r.k[0]);
                                int lastyy = otoY(liste[liste[n].left].r.k[2]);
#ifndef SHOW_ALL_BONDS
								if ((abs(SX-lastxx)<MAX_X) && (abs(SY-lastyy)<MAX_Y))
#endif
								{
                                    SDL_drawLine_Alpha(screen,SX,SY,lastxx,lastyy, bond_col,125);
                                }
                            }
                        }
                        else {                
                            SZ = otoZ(liste[n].r.k[2]);
                            col = Color[TYP2COL(liste[n].type, SZ)];
                            SX = otoX(liste[n].r.k[0]);
                            SY = otoY(liste[n].r.k[1]);

                            if (liste[n].left >= 0) {
								int bond_col = (liste[liste[n].left].type < liste[n].type) ? Color[TYP2COL(liste[liste[n].left].type, otoZ(liste[liste[n].left].r.k[2]))] : col;
                                int lastxx = otoX(liste[liste[n].left].r.k[0]);
                                int lastyy = otoY(liste[liste[n].left].r.k[1]);
#ifndef SHOW_ALL_BONDS
								if ((abs(SX-lastxx)<MAX_X) && (abs(SY-lastyy)<MAX_Y))
#endif
								{
                                    SDL_drawLine_Alpha(screen,SX,SY,lastxx,lastyy, bond_col,125);
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
        if (slice_mode) {
            sprintf(Str,"DPD %s %6d %6d  Z=%.3f-%.3f", ctrl_file_name, no_calc, outsiders, min_slice,max_slice);
        }
        else {
            sprintf(Str,"DPD %s %6d %6d", ctrl_file_name, no_calc, outsiders);
        }
        SDL_WM_SetCaption(Str,Str);
    //}
}


/*******************************************************************/
void slice_info()
{
    printf("\n   min=%5.3f  max=%5.3f   delta=%5.3f\n",
          min_slice, max_slice, delta_slice );
}
    
    


/*******************************************************************/
void helpinteract()
{
    printf("\n"
        "Program call:\n"
        "   dpd [option]\n"
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
        "   I:  show profile for X-velocity\n"
        "   J:  display state\n"
        "   P:  Create pressure to X at X=0  (vx += const)\n"
        "   E:  E-field\n"
        "   N:    Change polarity\n"
        "   Q:    Change field strength\n"
        "   K:  Change bonding force\n"
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
        "   z:  Exchange y-z for display (%d)\n"
        "   -:  toggle slicemode for display (%d)\n"
        "       slicemode: (arrow-keys)\n"
        "          up/down  move Z-coordinate\n"
        "          left/right  decr/increase delta-Z\n"
        "\n"
        "   S:   Save-Data\n"
        "   L:   Load-Data\n"
        "   O:   Load-new-design data\n"
        "   A:   Show actual parameter\n"
        "   V:   Save actual parameter in templatefile 'template.dat'\n"
        "   U:   Take a screenshot into 'dpd.bmp'\n"
        "\n"
        "   RET: Pause\n"
        "   H:   this help\n"
        "   ESC: End\n"
        "   G: End and design entry tool\n"
        "\n",  pressure_ofs, pressure_ofs, rfac, dt, xchg_yz, slice_mode);
}


/*******************************************************************/
int process_keysym(
    SDLKey sym,
    int is)
{            
    switch (sym) {
        case SDLK_MINUS   : slice_mode = 1-slice_mode; DrawIt(); slice_info(); break;  // !
        case SDLK_DOWN    : min_slice -= delta_slice/2; max_slice -= delta_slice/2; DrawIt(); slice_info(); break;
        case SDLK_UP      : min_slice += delta_slice/2; max_slice += delta_slice/2; DrawIt(); slice_info(); break;
        case SDLK_LEFT    : delta_slice /= 2.0; max_slice=min_slice+delta_slice; DrawIt(); slice_info(); break;
        case SDLK_RIGHT   : delta_slice *= 2.0; max_slice=min_slice+delta_slice; DrawIt(); slice_info(); break;
        
        case SDLK_a : template_info(stdout); break;
        case SDLK_b : wall_brakes = is; break;
        case SDLK_c : do_cons = is; break;
        case SDLK_d : do_diss = is; break;
        case SDLK_e : efeld = is; change_polar(is); break;
        case SDLK_f : change_rfac(is);  print_info(); break;
		case SDLK_g : break; // open
        case SDLK_h : helpinteract(); break;
        case SDLK_i : drawinfo = is; break;
        case SDLK_j : info_polys();  print_info(); break;
        case SDLK_k : change_bond(is);  print_info(); break;
        case SDLK_l : loaddata(1+is); DrawIt(); break;
        case SDLK_m : with_brick = is; break;
        case SDLK_n : change_polar(is); break;
        case SDLK_o : break; // open
        case SDLK_p : pressure_to_x = is; break;
        case SDLK_q : change_efield(is); break;
        case SDLK_r : do_rand = is; break;
        case SDLK_s : savedata(1+is); break;
        case SDLK_t : change_dt(is); print_info(); break;
        case SDLK_u : snapshot(); SDL_SaveBMP(screen,"dpd.bmp"); printf("Saved dpd.bmp\n"); break;
        case SDLK_v : save_template("template.dat"); break;
        case SDLK_w : wall_at_y0 = is; break;
        case SDLK_x : add_vx = 1; break;
        case SDLK_y : add_vy = 1; break;
        case SDLK_z : xchg_yz = 1-xchg_yz; DrawIt(); break;
        case SDLK_0 : ena_disp[0] = 1-ena_disp[0]; DrawIt(); break;
        case SDLK_1 : ena_disp[1] = 1-ena_disp[1]; DrawIt(); break;
        case SDLK_2 : ena_disp[2] = 1-ena_disp[2]; DrawIt(); break;
        case SDLK_3 : ena_disp[3] = 1-ena_disp[3]; DrawIt(); break;
        case SDLK_4 : ena_disp[4] = 1-ena_disp[4]; DrawIt(); break;
        case SDLK_5 : ena_disp[5] = 1-ena_disp[5]; DrawIt(); break;
        case SDLK_6 : ena_disp[6] = 1-ena_disp[6]; DrawIt(); break;
        case SDLK_7 : ena_disp[7] = 1-ena_disp[7]; DrawIt(); break;
        case SDLK_8 : ena_disp[8] = 1-ena_disp[8]; DrawIt(); break;
        case SDLK_9 : ena_disp[9] = 1-ena_disp[9]; DrawIt(); break;
        case SDLK_RETURN: {
                disp_pause = 1-disp_pause; 
                char Str[255];
                sprintf(Str,"DPD %d  %s", no_calc, disp_pause?"Pause":"");
                SDL_WM_SetCaption(Str,Str);
                break;
        }
        default: return -1; break;      
    }
    return 0;
}



/*******************************************************************/
int handleevent() 
{        
    int ret = 0;
    
    SDL_PollEvent(&event);
    if (Quit) ret = 1;  // exit(1);

    if (KeyDown) {
        int is = 1;
        //printf("Mod=%04X  Sym=%04X\n", MOD, SYM);
        if (MOD & KMOD_SHIFT) { 
           is = 0;
        }
        process_keysym(SYM,is);
        
        KeyDown = 0;
        //event.key.keysym.sym = 0;
        print_info2();
        
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
			 DrawIt();
    }
    return ret;
}    

/*******************************************************************/
int  wait_sdl() 
{   
    while ((disp_pause == 1) && !Quit) { 
        handleevent();
        // do not let the event cycle consume the whole processor
        RestCPU();  
        // SDL_Delay(200);
        // usleep(200000);
        // sleep(1);
    }
    return (Quit);
}    


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void display()
{
    DrawIt();
    handleevent();
    wait_sdl();
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
             int r = maxrgb(particle[typ].color[0]+ shade*RGBSHADEFAC);
             int g = maxrgb(particle[typ].color[1]+ shade*RGBSHADEFAC);
             int b = maxrgb(particle[typ].color[2]+ shade*RGBSHADEFAC);
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
    Color[9] = SDL_MapRGB(screen->format,0x60,0xe0,0xe0);  // 
    Color[10]= SDL_MapRGB(screen->format,0xe0,0x60,0x60);  // 
}


/*******************************************************************/
void sdl_set_size()
{
    sdl_xfac = PIXPERDOT*(MAX_X/(dl.k[0]+SDLRAND));
    sdl_yfac = PIXPERDOT*(MAX_Y/(dl.k[1]+SDLRAND));
    sdl_zfac = ZSHADES/dl.k[2];
    sdl_xofs = PIXPERDOT*MAX_X/2;
    sdl_yofs = PIXPERDOT*MAX_Y/2;
    sdl_zofs = TYPSHADES/2;
    sdl_xmod = MAX_X*PIXPERDOT;
    sdl_ymod = MAX_Y*PIXPERDOT;
    sdl_zmod = 1; 
}



/*******************************************************************/
void init_sdl(Uint32 flags)
{    
    FULLSCREEN=0;
    WIDTH=MAX_X*PIXPERDOT;
    HEIGHT=MAX_Y*PIXPERDOT;
    BPP=16;
    SOUND=0;

    SDL_INIT_ALL(flags);
    // set window caption
    SDL_WM_SetCaption("dpd","dpd");
    // clean screen
    ClearSurfaceRGB(screen,0,0,0);
    SDL_UpdateRect(screen,0,0,0,0);
    
    createcolors();
    sdl_set_size();
}

/*******************************************************************/
void end_sdl() 
{    
//    SDL_SaveBMP(screen,"dpd.bmp");   
    SDL_Quit();
}    

/*******************************************************************/

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 20.07.2004  17:10h
  $Id: dpd_disp.h,v 1.13 2004/09/29 12:38:27 tmaeke Exp $
********************************************************************/

#ifndef _dpd_disp_h_
#define _dpd_disp_h_

/*Include***********************************************************/
#include "SDL.h"

/*Defines***********************************************************/


/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

void snapshot();
void helpinteract();
int process_keysym(SDLKey sym, int is);
int handleevent();
void DrawIt();
int maxrgb(int);
void createcolors();
void init_sdl(Uint32);
int wait_sdl();
void sdl_set_size();
void end_sdl(); 
void display();

void drawline_ofs(real_t x0, real_t y0, real_t x1, real_t y1, int colorno, Uint8 blend);
void drawrect_ofs(real_t x0, real_t y0, real_t x1, real_t y1, int colorno, Uint8 blend);

#ifdef __cplusplus
}
#endif

/*Types*************************************************************/

/*Vars**************************************************************/

    extern char snapsfn[MAX_NAME];
    extern int snaps_stp;
    extern int snaps_start;
    extern int snaps_end;
    extern int snapsintv;             // snapshot interval

    extern bool_t slice_mode;     // 0=all 1=slice from min_- to max_slice
    extern real_t min_slice;
    extern real_t max_slice;
    extern real_t delta_slice;

/*******************************************************************/
#endif  /* _dpd_disp_h_ */

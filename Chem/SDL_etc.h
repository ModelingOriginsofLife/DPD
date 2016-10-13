/*
  Some stuff around SDL
*/

#ifndef _SDL_etc_h
#define _SDL_etc_h

#include "SDL.h"

#define GAME_TICK 40


/* Set up for C function definitions, even when using C++ */
#ifdef __cplusplus
extern "C" {
#endif

extern SDL_Surface *screen;
extern SDL_Event event;

extern Uint8	BPP;
extern char	TURBO;
extern char	FULLSCREEN;
extern char	RESIZABLE;
extern char	IconF[256];

extern SDL_Surface *IconS;

struct MouseStruct{
	int X,Y;
	char B;
};

extern struct MouseStruct Mouse;

extern char Quit,Left,Right,Up,Down,Fire,Alt,Ctrl,Space,Setup,KeyDown;
extern Uint16 UC;
extern SDLKey SYM;
extern SDLMod MOD;

extern unsigned int  WIDTH,HEIGHT;

extern int MVOLUME;
extern int SVOLUME;
extern char SOUND;

void        VIDEO_INIT(void);
SDL_Surface *NewSurface(int W,int H,Uint32 Fl);
void        ClearSurfaceRGB(SDL_Surface *Co, Uint8 R, Uint8 G, Uint8 B);
void        ClearSurface(SDL_Surface *Co);
void        DrawSpr(SDL_Surface *Co,SDL_Surface *Kam,int X,int Y);
void        R_Blit(SDL_Surface *Co,SDL_Surface *Kam,int SX, int SY, int SW, int SH, int X, int Y);
Uint32	    getpixel(SDL_Surface *surface, int x, int y);
void	    putpixel(SDL_Surface *surface, int x, int y, Uint32 pixel);
int         E_FILTER(const SDL_Event *event);
void        SDL_INIT_ALL(Uint32);
SDL_Cursor  *CreateCursor(const char *image[]);
void        RestCPU();


/* Ends C function definitions when using C++ */
#ifdef __cplusplus
};
#endif

#endif


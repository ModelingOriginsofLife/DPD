#include <stdlib.h>
#include <string.h>
#include "SDL_etc.h"

//#define DEBUGOUT(str)  fprintf(stderr,"SDL " str "\n") 
#define DEBUGOUT(str) 

Uint8	BPP=0;
char	TURBO=0;
char	FULLSCREEN=0;
char	RESIZABLE=0;
char	IconF[256]="";

SDL_Surface *screen;
SDL_Surface *IconS=NULL;

SDL_Event event;
char Quit=0,Left=0,Right=0,Up=0,Down=0,Fire=0,Alt=0,Ctrl=0,Space=0,Setup=0,KeyDown=0;
Uint16 UC=0;
SDLKey SYM=0;
SDLMod MOD=0;
struct MouseStruct Mouse;

unsigned int  WIDTH=320,HEIGHT=200;

Uint32 LastTicks=0;

int MVOLUME=64;
int SVOLUME=128;
char SOUND=1;
//Mix_Music *Music;

// --- Music ---
/*
void MusicLoopback(void){
//	printf("Module finished...\n");
	Mix_PlayMusic(Music,0);
	Mix_VolumeMusic(MVOLUME);
}

void PlaySample(Mix_Chunk *Snd){
	int Ch;
	if (SOUND)
		if (Mix_PlayChannel(-1,Snd,0)<0){
			Ch=Mix_GroupOldest(-1);
			Mix_HaltChannel(Ch);
			Mix_PlayChannel(Ch,Snd,0);
		}
}

void MusicStart(void){
	if (SOUND){
		Mix_HookMusicFinished(&MusicLoopback);
		Mix_VolumeMusic(MVOLUME);
		Mix_FadeInMusic(Music,0,3000);
	}
}

void MusicStop(void){
	if (SOUND){
		Mix_FadeOutMusic(1000);
		SDL_Delay(1000);
	}
}
*/

// --- Video ---

void VIDEO_INIT(void){
	Uint32 videoflags=0;

	if (BPP==0)
		videoflags = SDL_ANYFORMAT;
	if (FULLSCREEN) videoflags=videoflags | SDL_FULLSCREEN;
	if (RESIZABLE) videoflags=videoflags | SDL_RESIZABLE;
	if (TURBO) videoflags=videoflags | SDL_HWSURFACE;
	else videoflags=videoflags | SDL_SWSURFACE;

	if ( (screen=SDL_SetVideoMode(WIDTH,HEIGHT,BPP,videoflags)) == NULL ) {
		fprintf(stderr, "Couldn't set %ix%i video mode: %s\n",WIDTH,HEIGHT,SDL_GetError());
		exit(2);
	}
}


SDL_Surface *NewSurface(int W,int H,Uint32 Fl){
	return SDL_CreateRGBSurface(Fl,W,H,BPP,screen->format->Rmask,screen->format->Gmask,screen->format->Bmask,screen->format->Amask);
}


void ClearSurfaceRGB(SDL_Surface *Co, Uint8 R, Uint8 G, Uint8 B){
Uint32 color;
SDL_Rect clip;

	clip.x = 0;
	clip.y = 0;
	clip.w = Co->w;
	clip.h = Co->h;

	color=SDL_MapRGBA(Co->format,R,G,B,0);

	SDL_FillRect (Co, &clip, color);
}


void ClearSurface(SDL_Surface *Co){
	ClearSurfaceRGB(Co,255,0,255);
}


void DrawSpr(SDL_Surface *Co,SDL_Surface *Kam,int X,int Y){
SDL_Rect Dest;
	if (Co==NULL) return;
	Dest.x=X;
	Dest.y=Y;
	Dest.w=Co->w;
	Dest.h=Co->h;
	SDL_BlitSurface(Co,NULL,Kam,&Dest);
}


void R_Blit(SDL_Surface *Co,SDL_Surface *Kam,int SX, int SY, int SW, int SH, int X, int Y){
SDL_Rect Src,Dest;
	if (Co==NULL) return;
	Src.x=SX;
	Src.y=SY;
	Src.w=SW;
	Src.h=SH;
	Dest.x=X;
	Dest.y=Y;
	Dest.w=SW;
	Dest.h=SH;
	SDL_BlitSurface(Co,&Src,Kam,&Dest);
}


Uint32 getpixel(SDL_Surface *surface, int x, int y)
{
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to retrieve */
    Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;

    switch(bpp) {
    case 1:
        return *p;

    case 2:
        return *(Uint16 *)p;

    case 3:
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN)
            return p[0] << 16 | p[1] << 8 | p[2];
        else
            return p[0] | p[1] << 8 | p[2] << 16;

    case 4:
        return *(Uint32 *)p;

    default:
        return 0;       /* shouldn't happen, but avoids warnings */
    }
}


void putpixel(SDL_Surface *surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to set */
    Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;

//	if ((x<0)||(y<0)||(x>=surface->w)||(y>=surface->h))
//		return;

    switch(bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16 *)p = pixel;
        break;

    case 3:
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        } else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32 *)p = pixel;
        break;
    }
}

/*
char LoadSprite(SDL_Surface **Co, char *Soubor, Uint8 Trans){
	SDL_Surface *P;
	
	if ((P=IMG_Load(Soubor))==NULL){
		fprintf(stderr,"LoadSprite: Error loading '%s'!\n",Soubor);
		return 0;
	}else{
		fprintf(stderr,"LoadSprite: Loaded '%s'\n",Soubor);
		SDL_SetColorKey(P,SDL_SRCCOLORKEY|SDL_RLEACCEL,SDL_MapRGB(P->format,255,0,255));
		if (Trans)
			SDL_SetAlpha(P,SDL_SRCALPHA,Trans);
		*Co=SDL_DisplayFormat(P);
		SDL_FreeSurface(P);
	}
	return 1;
}
*/


// --- Event filter ---

int E_FILTER(const SDL_Event *event) {

	
	if ((RESIZABLE)&&(event->type == SDL_VIDEORESIZE)) {
		WIDTH=event->resize.w;
		HEIGHT=event->resize.h;
		VIDEO_INIT();
		return 0;
	}

	if (event->type==SDL_KEYDOWN){
		if (event->key.keysym.sym==SDLK_LEFT) Left=1;
		if (event->key.keysym.sym==SDLK_RIGHT) Right=1;
		if (event->key.keysym.sym==SDLK_UP) Up=1;
		if (event->key.keysym.sym==SDLK_DOWN) Down=1;
		if (event->key.keysym.sym==SDLK_RCTRL) Fire=1;
		if (event->key.keysym.sym==SDLK_LCTRL) Fire=1;
		if (event->key.keysym.sym==SDLK_SPACE) Space=1;
		if (event->key.keysym.sym==SDLK_ESCAPE) Quit=1;
		if (event->key.keysym.sym==SDLK_LALT) Alt=1;
		if (event->key.keysym.sym==SDLK_LCTRL) Ctrl=1;
		if (event->key.keysym.sym==SDLK_x && Alt) Quit=2;
		if ((event->key.keysym.sym==SDLK_s) && Alt) Setup=1;
		if ((event->key.keysym.sym==SDLK_q) && Alt) Quit=2;
		if ((event->key.keysym.sym==SDLK_F4) && Alt) Quit=2;
		UC=event->key.keysym.unicode;
		SYM=event->key.keysym.sym;
		MOD=event->key.keysym.mod;
		KeyDown=1;
		return(0);
	}

	if (event->type==SDL_KEYUP){
		if (event->key.keysym.sym==SDLK_LEFT) Left=0;
		if (event->key.keysym.sym==SDLK_RIGHT) Right=0;
		if (event->key.keysym.sym==SDLK_UP) Up=0;
		if (event->key.keysym.sym==SDLK_DOWN) Down=0;
		if (event->key.keysym.sym==SDLK_RCTRL) Fire=0;
		if (event->key.keysym.sym==SDLK_LCTRL) Fire=0;
		if (event->key.keysym.sym==SDLK_SPACE) Space=0;
		if (event->key.keysym.sym==SDLK_ESCAPE) Quit=0;
		if (event->key.keysym.sym==SDLK_LALT) Alt=0;
		if (event->key.keysym.sym==SDLK_LCTRL) Ctrl=0;
		KeyDown=0;
		return 0;
	}

	if ( event->type == SDL_MOUSEBUTTONDOWN ) {
		if (event->button.button==SDL_BUTTON_LEFT)
			Mouse.B=1;
		if (event->button.button==SDL_BUTTON_RIGHT)
			Mouse.B=2;
		if (event->button.button==SDL_BUTTON_MIDDLE)
			Mouse.B=3;
		return 0;
	}

	if ( event->type == SDL_MOUSEBUTTONUP ) {
		Mouse.B=0;
		return 0;
	}

	if ( event->type == SDL_MOUSEMOTION ) {
        Mouse.X=event->motion.x;
		Mouse.Y=event->motion.y;
		return(0);
    }

	if (event->type==SDL_QUIT) Quit=2;
	
	return(1);
}


void SDL_INIT_ALL(void){
	Uint32 InitFlags;

	InitFlags=SDL_INIT_VIDEO;
	if (SOUND)
	 InitFlags|=SDL_INIT_AUDIO;

	DEBUGOUT("1a");
	if ( SDL_Init(InitFlags) < 0 ) {
		fprintf(stderr, "Couldn't initialize SDL: %s\n",SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);

	DEBUGOUT("1c");

	if (SOUND)
/*		if ( Mix_OpenAudio(22050, MIX_DEFAULT_FORMAT, 2, 1024) < 0 ){// 22050, 8192
			fprintf(stderr,"Warning: Couldn't set 22050 Hz 16-bit audio\n- Reason: %s\n",SDL_GetError());
			SOUND=0;
		}
*/
//	if ((IconF!=NULL)&&(IconF[0]!='\0'))
//		SDL_WM_SetIcon(IMG_Load(IconF), NULL);
//	else
	if (IconS!=NULL)
		SDL_WM_SetIcon(IconS, NULL);

	DEBUGOUT("1b");

	SDL_ShowCursor(SDL_ENABLE);
	SDL_SetEventFilter(E_FILTER);
	SDL_EnableKeyRepeat(0,1);
	VIDEO_INIT();
	DEBUGOUT("1c");
//	if (SOUND)
//		Mix_Volume(-1,SVOLUME);
}



SDL_Cursor *CreateCursor(const char *image[])
{
  int i, row, col;
  Uint8 data[4*32];
  Uint8 mask[4*32];
  int hot_x, hot_y;

  i = -1;
  for ( row=0; row<32; ++row ) {
    for ( col=0; col<32; ++col ) {
      if ( col % 8 ) {
        data[i] <<= 1;
        mask[i] <<= 1;
      } else {
        ++i;
        data[i] = mask[i] = 0;
      }
      switch (image[4+row][col]) {
        case 'X':
          data[i] |= 0x01;
          mask[i] |= 0x01;
          break;
        case '.':
          mask[i] |= 0x01;
          break;
        case ' ':
          break;
      }
    }
  }
  sscanf(image[4+row], "%d,%d", &hot_x, &hot_y);
  return SDL_CreateCursor(data, mask, 32, 32, hot_x, hot_y);
}


void RestCPU(){
	Uint32 ActTicks=SDL_GetTicks();
	if ((ActTicks-LastTicks)<GAME_TICK)
		SDL_Delay(GAME_TICK-(ActTicks-LastTicks));
	LastTicks=ActTicks;
}

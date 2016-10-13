/*
 *  rgtest.h
 *  evoselfc
 *
 *  Created by John McCaskill on Mon Mar 15 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

struct Lat
{
	int type;
	float x;
	float y;
	float vx;
	float vy;
	float fx;
	float fy;
	int updated;
	struct Lat * next;
	struct Lat * prev;
};

typedef struct Lat Lat;
extern Lat **Lpl;

extern void latinit();
extern void latkine(int pos);
extern void latmove(int pos);
extern void lbndy(Lat **Plane);

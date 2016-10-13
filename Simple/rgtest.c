/*
 *  rgtest.c
 *  evoselfc
 *
 *  Created by John McCaskill on Mon Mar 15 2004.
 *  Edited into final dpd form by July 19th 2004
 *
 */

#include "evoself.h"

extern const float RAND_MAX_INV;
#define RANDF (rand()*RAND_MAX_INV - 0.5)

const float dt=0.01;


									// basic dpd momentum update of all particles at lattice site
void latkine(int pos)
{
	int i;
	float x1,y1,dx,dy,lx,ly;
	float r,r2,s,rinv,invcut;
	float fx,fy,fxcons,fycons,fxdiss,fydiss,fxrand,fyrand;
	Lat *tmp1,*tmp2;
	
	const float cutoff=1.0;
	const float cutoffsq=cutoff*cutoff;
	const float k1=-3.0;
	const float k2=-12.0;
	const float sqrtdt=sqrt(dt);
	const float g=0.1;
	const float kTdpd=1.0;

	invcut=1./cutoff;
	tmp1=Lpl[pos];
	while(tmp1!=NULL) {
		x1=tmp1->x;y1=tmp1->y;
		tmp1->updated=2;
		for(i=0;i<9;i++) {
			lx=i%3-1;
			ly=(i/3)%3-1;
			
			tmp2=Lpl[nbmoorepb(pos,i)];

			while(tmp2!=NULL) {
				if(tmp2!=tmp1) {
				 if(tmp2->updated!=2) {
				  dx=x1-lx-tmp2->x;
				  dy=y1-ly-tmp2->y;
				  if((r2=dx*dx+dy*dy)<cutoffsq) {
								// calculate force contributions to update momenta of particles
					r=sqrt(r2);
					s=g*r*invcut;
					rinv=1./r;
										// 1. conservative force
					if(tmp1->type==tmp2->type) {
						fxcons = -k1*dx;
						fycons = -k1*dy;
					}
					else {
						fxcons = -k2*dx;
						fycons = -k2*dy;					
					}
										// 2. dissipative force
					fxdiss = s*s*(tmp2->vx-tmp1->vx);
					fydiss = s*s*(tmp2->vy-tmp1->vy);
										// 3. random force 					
					fxrand = kTdpd*RANDF*s;
					fyrand = kTdpd*RANDF*s;
					fx=(fxcons+fxdiss)*dt+fxrand*sqrtdt;
					fy=(fycons+fydiss)*dt+fyrand*sqrtdt;
					tmp2->fx-=fx;
					tmp1->fx+=fx;
					tmp2->fy-=fy;
					tmp1->fy+=fy;
				   }
				  }
				}
				tmp2=tmp2->next;
			}
		}								

								// change in velocities
		tmp1->vx+=tmp1->fx;
		tmp1->vy+=tmp1->fy;
		tmp1=tmp1->next;
	}
			
}
									// basic dpd move of all particles at lattice site
void latmove(int pos)
{
	int ix,iy,newpos;
	Lat *tmp1,*tmp2;


	tmp1=Lpl[pos];
	while(tmp1!=NULL) {
	  tmp2=tmp1->next;
	  if(tmp1->updated!=1) {
		tmp1->x+=tmp1->vx*dt;
		tmp1->y+=tmp1->vy*dt;
		tmp1->fx=0.;
		tmp1->fy=0.;
		tmp1->updated=1;
							// calculate destination lattice offsets for moved particle
		if      (tmp1->x >= 0.5) {
			ix=1;tmp1->x-=1.0;
		}
		else if (tmp1->x < -0.5) {
			ix=-1;tmp1->x+=1.0;
		}
		else ix=0.;
		if      (tmp1->y >= 0.5) {
			iy=1;tmp1->y-=1.0;
		}
		else if (tmp1->y < -0.5) {
			iy=-1;tmp1->y+=1.0;
		}
		else iy=0.;
							// move particle to new neighboring lattice site
		if((ix!=0) || (iy!=0)) {
			newpos=nbmoorepb(pos,(iy+1)*3+ix+1);	// index of neighboring destination
													// cut out particle from doubly linked list
													// particle is top of list
			if(tmp1->prev==NULL) {
				if(tmp1->next!=NULL) tmp1->next->prev=NULL;
				Lpl[pos]=tmp1->next;
			}
			else {									// particle lower in list
				if(tmp1->next!=NULL) tmp1->next->prev=tmp1->prev;
				tmp1->prev->next=tmp1->next;
				tmp1->prev=NULL;
			}
													// insert particle at top of neighbor list
			if(Lpl[newpos]!=NULL) Lpl[newpos]->prev=tmp1;
			tmp1->next=Lpl[newpos];
			Lpl[newpos]=tmp1;
		}
	  }
	  tmp1=tmp2;
	}
			
}

void latinit()
{
	Lat * tmp, *tmp1;
	int i,j,xx,yy;
	
	mprintf1("nparticles=%d",nparticles);
	
	for(j=0;j<nparticles;j++) {
		xx=modrand(rowl-2)+1;
		yy=modrand(nrows-2)+1;
		i=yy*rowl+xx;
						// insert particle at top of list
		tmp1=Lpl[i];
		tmp = Lpl[i] = (Lat *) calloc(1,sizeof(Lat));
		tmp->next = tmp1;
		tmp->prev = NULL;
		if(tmp1!=NULL) tmp1->prev = tmp;

						// random coordinate offsets
		tmp->x=RANDF;
		tmp->y=RANDF;
		tmp->vx=RANDF*4.0;
		tmp->vy=RANDF*4.0;
		tmp->fx=0.;
		tmp->fy=0.;
		tmp->type=modrand(2);
		tmp->updated=0;
	}
	lbndy(Lpl);
}

//* DPD -VV 
        /* 1. Calculate all forces (conservative, dissipative, and random forces) */
            
        /* 2. Integrate velocities over dt/2 */
            
        /* 3. Calculate dissipative forces */
		
		/* 4. Integrate velocities over dt/2 */
		
		/* 5. Integrate particle positions over dt */


void lbndy(Lat **Plane)	// like cbndy but for Lat pointer structure
{
     register int i1,i2,itop,ibot;
     register int i,pk;
     
	 for(pk=0;pk<nsites;pk+=planel) {	  // depth loop        
							// x direction (along row)
		itop = rowl;
		ibot = (nrows-2)*rowl;	 
		i1=0;i2=ibot+rowl;
	 
		for(i=0;i<rowl;i++){
			Plane[i1++]=Plane[ibot++];
			Plane[i2++]=Plane[itop++];
		}
							// y direction
     
		i1=0;i2=rowl-1;
	 
		for(i=0;i<nrows;i++)
		{
			Plane[i1]=Plane[i2-1];
			Plane[i2]=Plane[i1+1];
			i1+=rowl;i2+=rowl;
		}
	} // pk for depth
}
	
 



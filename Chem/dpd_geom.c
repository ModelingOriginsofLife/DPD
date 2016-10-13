/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                                             Thomas Maeke  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 28.07.2004  11:50h
  $Id: dpd_geom.c,v 1.9 2004/09/28 16:25:43 tmaeke Exp $
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

/* #include <signal.h> */
/* #include <time.h>   */
/* #include <limits.h> */
/* #include <errno.h>  */

#include "dpd_vars.h"
#include "dpd_geom.h"
#include "dpd_disp.h"
#include "dpd_init.h"

/*Defines***********************************************************/

/*Types*************************************************************/

/*Vars**************************************************************/

const int magic_elem = 0x233445AB;
const int magic_poly = 0x34455667;
const int magic_elec = 0x89712ABC;

int magics[max_elem_type] = {
    0x12233445,
    0x23344556,
    0x34455667,
    0x45566778,
    0x56677889,
    0x67788990,
    0x78899001,
    0x89900112,
    0x9AABBCCD,
    0x10DEF012,
    0xA1122345,
    0xbb122345,
    0x12ADBC17
};

/* ****************************************************** */
/* ****************************************************** */
//elem_p coord_base = NULL;    
elem_p geom_channels = NULL;
elem_p geom_electrodes = NULL;
elem_p geom_sensors = NULL;

int draw_poly = 0;




/*Funcs*************************************************************/


/* ****************************************************** */
void indent_print(  
    int i,
    FILE* fi) 
{
    for (; i>0; i--) {
        fprintf(fi," ");
    }
}    
        

/* ****************************************************** */
/* Routines working on string_p                           */
/* ****************************************************** */

/* ****************************************************** */
string_p string_delete(
    string_p str)
{
    /*
        delete only the contents of s
    */
    if (str) {
        if ((str->len > 0) && (str->txt != NULL)  ) {
            free(str->txt);
        }
        str->txt = NULL;
        str->len = 0;
    }
    return str;
}    


/* ****************************************************** */
string_p string_repl(
    string_p str,
    char *src)
{    
    if (str==NULL) {
        str = string_malloc(0);
    }
    if (str && src) {
        str = string_delete(str);
        str->len = strlen(src);
        if (str->len > 0) {
            str->txt = (char*)malloc(str->len+1);
            if (str->txt && src) {
                strcpy(str->txt, src);
            }
        }
    }
    return str;
}


/* ****************************************************** */
string_p string_repl_s(
    string_p str,
    string_p src)
{    
    str = string_delete(str);
    if (src && src->txt) {
        str = string_repl(str, src->txt);
    }
    return str;
}


/* ****************************************************** */
string_p string_malloc(
    char * src)
{
    /*
       char * xx="hallo";
       string_p str;
       str = string_malloc(xx); 
    */
    
    string_p str = (string_p) malloc(sizeof(string_t));
    str->txt = NULL;
    str->len = 0;
    if (str) {
        str = string_repl(str, src);
    }
    return str;
}


/* ****************************************************** */
string_p string_free(
    string_p str)
{    
    /* 
        deletes the hole structure
        s is no longer valid
    */
    if (str) {
       string_delete(str);
       free (str); 
    }
    return NULL;
}


/* ****************************************************** */
void string_print(
    string_p str,
    int indent,
    FILE* fi)
{    
    /* 
        deletes the hole structure
        s is no longer valid
    */
    indent_print(indent,fi);
    if (str) {
        if (str->txt) {
           fprintf(fi, "%s\n", str->txt);
        }
        else fprintf(fi, "empty strg\n");
    }
    else fprintf(fi, "null strg_p\n");
}



/* ****************************************************** */
/* Routines working on poly_p                             */
/* ****************************************************** */

/* ****************************************************** */
poly_p polygon_malloc(
    void)
{
    poly_p poly = (poly_p)malloc(sizeof(poly_t));
    if (poly) {
        poly->magic = magic_poly;
        poly->len = 0;
        poly->colno = 1;
        poly->dat = (fpoint_t*) malloc(poly->len * sizeof(fpoint_t));
    }
    return poly;
}


/* ****************************************************** */
void * polygon_free (
    poly_p poly)
{
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (free)\n");
            return NULL;
        }
        if (poly->len > 0) {
            free (poly->dat);
            poly->magic = 0x555aaa5a;
            poly->len = 0;
            poly->colno = 0;
        }
        free (poly);
    }
    return NULL;
}

/* ****************************************************** */
poly_p polygon_add_point(
    poly_p poly,
    fpoint_t data)
{
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (add)\n");
            return NULL;
        }
        poly->len++;
        poly->dat = (fpoint_t*) realloc(poly->dat, poly->len*sizeof(fpoint_t));
        if (poly->dat) {
            poly->dat[poly->len-1] = data;
        }
    }
    return poly;
}    
    
/* ****************************************************** */
poly_p polygon_set_color(
    poly_p poly,
    int colno)
{
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (add)\n");
            return NULL;
        }
        poly->colno = colno;
    }
    return poly;
}    
    
/* ****************************************************** */
void polygon_print(
    poly_p thisp,
    int indent,
    FILE* fi)
{
    if (thisp && thisp->len>0) {
        if (thisp->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (print)\n");
            return;
        }
        int i;
        indent_print(indent,fi);
        fprintf(fi, "Polygon: col:%d\n", thisp->colno);
        for (i=0; i<thisp->len; i++) {
            if (i%2==0) indent_print(indent,fi);
            fprintf(fi, "  (%7.3f,%7.3f,%7.3f)", thisp->dat[i].k[0], thisp->dat[i].k[1],thisp->dat[i].k[2]);
            if (i%2==1 || i==thisp->len-1) fprintf(fi, "\n");
        }
    }
}        

/* ****************************************************** */
void polygon_draw(
    poly_p poly)
{
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (draw) %X\n", poly->magic);
            return;
        }
        int i;
        for (i=0; i<poly->len; i++) {
            drawline_ofs(poly->dat[i].k[0], 
                         poly->dat[i].k[1], 
                         poly->dat[(i+1)%poly->len].k[0], 
                         poly->dat[(i+1)%poly->len].k[1],
                         poly->colno,
                         255);
        }
    }
}


/* ****************************************************** */
#if 0
inline int isLeft( float_t x0, float_t y0,
                   float_t x1, float_t y1,
                   float_t x2, float_t y2)
{
    return ( (x1 - x0) * (y2 - y0)
           - (x2 - x0) * (y1 - y0) );
}

bool PointInPolygon(float x, float y)
{
    int    wn = 0;    // the winding number counter
                    
    // loop through all edges of the polygon
    int n,i;
    for (i=0; i<n-1; i++)  // edge from V[i] to V[i+1]
    {
        if  (poly[i]->y <= y) {            // start y <= pt->y
            if (poly[i+1]->y > y)          // an upward crossing
                if (isLeft( poly[i]->x,poly[i]->y,
                            poly[i+1]->x,poly[i+1]->y,
                            x,y) > 0)      // P left of edge
                    ++wn;                  // have a valid up intersect
        }
        else {                             // start y > P.y (no test needed)
            if (poly[i+1]->y <= y)         // a downward crossing
                if (isLeft( poly[i]->x,poly[i]->y,
                            poly[i+1]->x,poly[i+1]->y,
                            x,y) < 0)      // P right of edge
                    --wn;                  // have a valid down intersect
        }
     }
     if (wn==0)  return false;
     return true;
}
#endif
    




/* ****************************************************** */
bool_t polygon_inside_xy(
    poly_p poly,
    fpoint_t data)
{
    bool_t side = 0;
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (inside) %X\n", poly->magic);
            return -1;
        }
        int i,j,n;
        float_t x = data.k[0];
        float_t y = data.k[1];
        n = poly->len;
        
        /*-------------------------------------------------
        
             (Y[i]<=y and y<Y[j])  or  (Y[j]<=y and y<Y[i])
        
          and 
                  X[j]-X[i]
             x <  --------- * (y-Y[i]) + X[i]
                  Y[j]-Y[i]
                  
        --------------------------------------------------*/
        
        for (i=0, j=n-1; i<n; j=i++) {
            if (
                (
                 ( (poly->dat[i].k[1] <= y) && (y < poly->dat[j].k[1]) ) ||
                 ( (poly->dat[j].k[1] <= y) && (y < poly->dat[i].k[1]) )
                )
                &&
                ( x < 
                  (poly->dat[j].k[0] - poly->dat[i].k[0]) *   //    X[j]-X[i]
                  (y - poly->dat[i].k[1]) /                   // *  y-Y[i]
                  (poly->dat[j].k[1] - poly->dat[i].k[1])     // /  Y[j]-Y[i]
                  + poly->dat[i].k[0]                         // +  X[i]
                )
               ) 
            {
                side = ! side;
            }
        }
    }
    return side ? 1:0;
}


/* ****************************************************** */
void polygon_min_max(
    poly_p poly,
    fpoint_t * min,
    fpoint_t * max)
{
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (minmax) %X\n", poly->magic);
            return;
        }
        int i,k;
        for (i=0; i<poly->len; i++) {
            for (k=0; k<3; k++) {
                if (poly->dat[i].k[k] < min->k[k]) min->k[k]=poly->dat[i].k[k];
                if (poly->dat[i].k[k] > max->k[k]) max->k[k]=poly->dat[i].k[k];
            }
        }
    }
}

/* ****************************************************** */
void polygon_sub_mul(  // poly = (poly-sub)*mul
    poly_p poly,
    fpoint_t * sub,
    fpoint_t * mul)
{
    if (poly) {
        if (poly->magic != magic_poly) {
            fprintf(stderr,"Error magic poly (submul) %X\n", poly->magic);
            return;
        }
        int i,k;
        for (i=0; i<poly->len; i++) {
            for (k=0; k<3; k++) {
                poly->dat[i].k[k] = (poly->dat[i].k[k]-sub->k[k])*mul->k[k];
            }
        }
    }
}


/* ****************************************************** */
/* Routines working on electrode_p                        */
/* ****************************************************** */

/* ****************************************************** */
electrode_p electrode_malloc(
    void)
{
    electrode_p elec = (electrode_p)malloc(sizeof(electrode_t));
    if (elec) {
        elec->magic = magic_elec;
        elec->active = 0;     // 
        elec->position = (fpoint_t) {{0,0,0}};   // center of electrode
        elec->size     = (fpoint_t) {{1,1,1}};   // position +/- size/2
        elec->polarity = 1.0;   // -1 .. 0 .. 1
        elec->dutycycle = 100;  // 0 .. 100%
        elec->seq_ofs = 0;      // ofs in seq.list
        elec->poly = polygon_malloc();
        elec->ownpoly = 1;  // we have allocates this poly
    }
    return elec;
}


/* ****************************************************** */
void * electrode_free (
    electrode_p elec)
{
    if (elec) {
        if (elec->magic != magic_elec) {
            fprintf(stderr,"Error magic elec (free) %X\n",elec->magic);
            return NULL;
        }
        if (elec->ownpoly) polygon_free(elec->poly);
        elec->magic = 0x555aaa5a;
        free (elec);
    }
    return NULL;
}


/* ****************************************************** */
void electrode_print(
    electrode_p thisp,
    int indent,
    FILE* fi)
{
    if (thisp) {
        if (thisp->magic != magic_elec) {
            fprintf(stderr,"Error magic elec (print) %X\n",thisp->magic);
//            return;
        }
        indent_print(indent,fi);
        fprintf(fi,"Electrode: act:%d pol:%3.1f duty:%-3d seq:%d own:%d\n",
                    thisp->active, thisp->polarity,
                    thisp->dutycycle, thisp->seq_ofs, thisp->ownpoly);
        polygon_print(thisp->poly, indent+4,fi);
    }
}        


/* ****************************************************** */
void electrode_draw(
    electrode_p thisp)
{
    if (thisp) {
        if (thisp->magic != magic_elec) {
            fprintf(stderr,"Error magic elec (draw) %X\n",thisp->magic);
            return;
        }
        if (thisp->active) polygon_draw(thisp->poly);
    }
}        


/* ****************************************************** */
electrode_p electrode_set_color(
    electrode_p elec,
    int colno)
{
    if (elec) {
        if (elec->magic != magic_elec) {
            fprintf(stderr,"Error magic elec (setcol) %X\n",elec->magic);
            return NULL;
        }
        polygon_set_color(elec->poly, colno);
    }
    return elec;
}    
    
/* ****************************************************** */
electrode_p electrode_add_point(
    electrode_p elec,
    fpoint_t data)
{
    if (elec) {
        if (elec->magic != magic_elec) {
            fprintf(stderr,"Error magic elec (addpoint) %X\n",elec->magic);
            return NULL;
        }
        polygon_add_point(elec->poly, data);
    }
    return elec;
}    
    
/* ****************************************************** */
electrode_p electrode_add_polygon(
    electrode_p elec,
    poly_p poly)
{
    if (elec) {
        if (elec->magic != magic_elec) {
            fprintf(stderr,"Error magic elec (addpoly) %X\n",elec->magic);
            return NULL;
        }
        if (elec->poly) {
            elec->poly = (poly_p) polygon_free(elec->poly);
        }
        elec->poly = poly;
        elec->ownpoly = 0; // poly is allocated elsewhere
        
        fpoint_t min = {{ 1E10, 1E10, 1E10}};
        fpoint_t max = {{-1E10,-1E10,-1E10}};
        polygon_min_max(poly, &min, &max);
        int kp;
        for (kp=0; kp<3; kp++) {
            elec->position.k[kp] = (max.k[kp]+min.k[kp])/2;
            elec->size.k[kp] = (max.k[kp]-min.k[kp])/2;
        }
    }
    return elec;
}    
    

/* ****************************************************** */
/* Rooutines working on elem_p                            */
/* ****************************************************** */

/* ****************************************************** */
void elem_print(
    elem_p thisp,
    int indent,
    FILE* fi)
{
    indent_print(indent,fi);
    fprintf(fi, "elem: ");
    if (thisp) {
        if (thisp->magic != magic_elem) {
            fprintf(fi,"Error magic elem (print) ");
        }
        fprintf(fi, "%2d %08X %p %p %p  ",thisp->typ, thisp->magic, thisp, thisp->next, thisp->prev);
        string_print(thisp->name, 0, fi);
        if (thisp->load) {
            switch (thisp->typ) {
                case polygon: polygon_print((poly_p)thisp->load, indent+4, fi);  break;
                case electrode: electrode_print((electrode_p)thisp->load, indent+4, fi);  break;
                default: break;
            }
        }
    }
    else fprintf(fi, "null elem_p\n");
}    


/* ****************************************************** */
void elem_draw(
    elem_p thisp)
{
    if (thisp) {
        if (thisp->magic != magic_elem) {
            fprintf(stderr,"Error magic elem (draw) %X",thisp->magic);
        }
        if (thisp->load) {
            switch (thisp->typ) {
                case polygon: polygon_draw((poly_p)thisp->load);  break;
                case electrode: electrode_draw((electrode_p)thisp->load);  break;
                default: break;
            }
        }
    }
    else fprintf(stderr, "null elem_p\n");
}    


/* ****************************************************** */
elem_p elem_malloc(
    char * name,
    geom_type_t typ,
    elem_p liste)
{
    elem_p ele = (elem_p) malloc(sizeof(elem_t));
    if (ele) {
        ele->name = string_malloc(name);
        ele->prev = ele;
        ele->next = ele;
        ele->typ = typ;
        ele->load = NULL;
        if (liste) {
            //printf("insert\n");
            /*
              Liste -n-> <-p-ele-n-> <-p-eleN-n-> <-p-eleM
               \                                        /
                +---p-> <-n-eleX-p-> <-n-eleY-p-> <-n--+
            */
            elem_p tmpn = liste->next;
            elem_p tmpp = liste->next->prev;
            liste->next->prev = ele;
            liste->next = ele;
            ele->next = tmpn;
            ele->prev = tmpp;
        }
        ele->magic = magic_elem;
    }
    return ele;
}



/* ****************************************************** */
elem_p elem_free(
    elem_p thisp)
{
    if (thisp) {
        if (thisp->magic != magic_elem) {
            fprintf(stderr,"Error magic elem: %2d %08X %p\n", thisp->typ, thisp->magic, thisp);
            return thisp;
        }
        fprintf(stdout,"free "); elem_print(thisp,0,stdout);
        if (thisp->load) {
            switch (thisp->typ) {
                case polygon: thisp->load = polygon_free((poly_p)thisp->load); break;
                case electrode: thisp->load = electrode_free((electrode_p)thisp->load); break;
                default: break;
            }
        }
        if (thisp->load != NULL) {
            fprintf(stderr,"Error payload: %p\n", thisp->load);
            return thisp;
        }
        elem_p n = thisp->next;
        elem_p p = thisp->prev;
        if (thisp->prev) {
            thisp->prev->next = n;
        }
        if (thisp->next) {
            thisp->next->prev = p;
        }
        thisp->magic = 0x5555aaaa;
        thisp->typ = invalid;
        thisp->next = NULL;
        thisp->prev = NULL;
        string_free(thisp->name);
        free(thisp);

        if (n == p) { return NULL; }
        else { return n; }
    }
    return NULL;
}    


/* ****************************************************** */
elem_p elem_add(
    char * name,
    geom_type_t typ,
    elem_p liste,
    void * payload)
{
    elem_p ele = elem_malloc(name,typ,liste);
    if (ele) {
        ele->load = payload;
    }
    return ele;
}


/* ****************************************************** */
/* Routines: lists of elem_p                              */
/* ****************************************************** */

/* ****************************************************** */
void elem_list_print(
    elem_p thisp,
    int indent,
    FILE* fi)
{
    elem_p p = thisp;
    indent_print(indent,fi);
    fprintf(fi, "////////\n");
    while (p) {
        elem_print(p,indent+4, fi);
        p = p->next;
        if (p == thisp) break;
    }
    indent_print(indent,fi);
    fprintf(fi, "\\\\\\\\\\\\\\\\\n");
}       


/* ****************************************************** */
elem_p elem_list_free(
    elem_p thisp)
{
    elem_p p = thisp;
    while (p) {
        p = elem_free(p);
        //elem_list_print(p,8,stdout);
    }
    return NULL;
}       


/* ****************************************************** */
elem_p elem_list_findname(
    elem_p thisp,
    char * name)
{
    elem_p p = thisp;
    while (p) {
        if (!strcmp(p->name->txt, name)) {
            return p;
        }
        p = p->next;
        if (p == thisp) break;
    }
    return NULL;
}       


/* ****************************************************** */
elem_p elem_list_findname_qual(
    elem_p thisp,
    char * name,
    geom_type_t qual)
{
    elem_p p = thisp;
    while (p) {
        if ((p->typ == qual) && !strcmp(p->name->txt, name)) {
            return p;
        }
        p = p->next;
        if (p == thisp) break;
    }
    return NULL;
}       


/* ****************************************************** */
void elem_list_draw( 
    elem_p base,
    geom_type_t qual)
{
    elem_p p = base;
    //fprintf(stderr,"draw %s\n", p->name->txt);
    while (p) {
        if (p->typ==qual && p->load) {
            elem_draw(p);
            //fprintf(stderr,"   draw %s\n", p->name->txt);
        }
        p = p->next;
        if (p == base) break;
    }
}




/* ****************************************************** */
/* Routines working on an elem_p with poly_p load         */
/* ****************************************************** */

/* ****************************************************** */
elem_p poly_malloc(
    char * name,
    elem_p liste)
{
    poly_p poly = polygon_malloc();
    return elem_add (name,polygon,liste, (void*)poly);
}


/* ****************************************************** */
elem_p poly_setcolor(
    elem_p thisp, 
    int colno)
{
    if (thisp && thisp->typ==polygon && thisp->load) {
        polygon_set_color((poly_p)thisp->load, colno);
    }
    return thisp;
}    
    
    
/* ****************************************************** */
elem_p poly_add_point(
    elem_p thisp,
    fpoint_t data)
{
    if (thisp && thisp->typ==polygon && thisp->load) {
        polygon_add_point((poly_p)thisp->load, data);
    }
    return thisp;
}    
    
    

/* ****************************************************** */
void poly_min_max(
    elem_p elem,
    fpoint_t * min,
    fpoint_t * max)
{
    if (elem && elem->typ==polygon && elem->load) {
        polygon_min_max((poly_p)elem->load, min,max);
    }
}

/* ****************************************************** */
void poly_sub_mul(
    elem_p elem,
    fpoint_t * sub,
    fpoint_t * mul)
{
    if (elem && elem->typ==polygon && elem->load) {
        polygon_sub_mul((poly_p)elem->load, sub,mul);
    }
}


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void poly_allow_cube(
    elem_p polyparent)
{
    ipoint_t p;
    poly_p poly = (poly_p)polyparent->load;
    
    fprintf(stderr,"  %s\n", polyparent->name->txt);
    for (p.k[2] = -dli.k[2]/2-1; p.k[2] <= dli.k[2]/2+1; p.k[2]++) {
    
        for (p.k[1] = -dli.k[1]/2-1; p.k[1] <= dli.k[1]/2+1; p.k[1]++) { 
        
            for (p.k[0] = -dli.k[0]/2-1; p.k[0] <= dli.k[0]/2+1; p.k[0]++) {
            
                fpoint_t f = {{p.k[0], p.k[1], p.k[2]}};
                bool_t erg = polygon_inside_xy(poly, f);
                if (erg) {
                    ALLOWED(p) = 1;
                    fprintf(stderr,"%d ",erg);
                }
            }
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
    }
}    





/* ****************************************************** */
void poly_translate(
    real_t size)
{
    //int mode=1;
    /*
        translate coordinates from all polygons
        to fit into simulation area
    */
    fpoint_t min = {{ 1E10, 1E10, 1E10}};
    fpoint_t max = {{-1E10,-1E10,-1E10}};
    //fpoint_t one = {{    1,    1,    1}};
    
    elem_p p = geom_channels;
    while (p) {
        if (p->typ==polygon) poly_min_max(p, &min, &max);
        p = p->next;
        if (p == geom_channels) break;
    }
    
    fprintf(stdout, "min: %7.4f,%7.4f,%7.4f\n", min.k[0], min.k[1],min.k[2]);
    fprintf(stdout, "max: %7.4f,%7.4f,%7.4f\n", max.k[0], max.k[1],max.k[2]);

    /* rescale hole system
       to fit into +/- sys_size/2
    */
    int k,n=-1;
    fpoint_t delta;
    real_t maxk = -1E10;
    for (k=0; k<3; k++) {
        if (max.k[k] - min.k[k] > maxk) {
           n = k;
           maxk = max.k[k] - min.k[k];
        }
        min.k[k] = min.k[k] + (max.k[k]-min.k[k])/2;
    }    
    real_t fac = size / (maxk*1.1);
    delta = (fpoint_t){{ fac, fac, fac }};
    
    /* rescale coordinates */
    p = geom_channels;
    while (p) {
        if (p->typ==polygon) poly_sub_mul(p, &min, &delta);
        p = p->next;
        if (p == geom_channels) break;
    }
    p = geom_electrodes;
    while (p) {
        if (p->typ==polygon /*electrode*/) poly_sub_mul(p, &min, &delta);
        p = p->next;
        if (p == geom_electrodes) break;
    }
    p = geom_sensors;
    while (p) {
        if (p->typ==polygon) poly_sub_mul(p, &min, &delta);
        p = p->next;
        if (p == geom_sensors) break;
    }

    /* reduce simulation area */    
    fprintf(stderr," ------------ reduce simulation area -----------\n");
    min = (fpoint_t){{ 1E10, 1E10, 1E10}};
    max = (fpoint_t){{-1E10,-1E10,-1E10}};
    p = geom_channels;
    while (p) {
        if (p->typ==polygon) poly_min_max(p, &min, &max);
        p = p->next;
        if (p == geom_channels) break;
    }
    
    #if 0
    if (mode == 0) { // std
        geom_deny_all();  
        p = geom_channels;
        while (p) {
            if ((p->typ==polygon) && (p->name->txt[1] == 'C')) {   // channel 
                fpoint_t lmin = {{ 1E10, 1E10, 1E10}};
                fpoint_t lmax = {{-1E10,-1E10,-1E10}};
                poly_min_max(p, &lmin, &lmax);
                fprintf(stdout, "  deltas: %s  %7.4f,%7.4f,%7.4f\n", 
                                 p->name->txt, lmax.k[0]-lmin.k[0],
                                 lmax.k[1]-lmin.k[1], lmax.k[2]-lmin.k[2] );
                geom_allow_cube(&lmin, &lmax);
            }
            p = p->next;
            if (p == geom_channels) break;
        }
    }
    else {
    #endif
        geom_deny_all();  
        p = geom_channels;
        while (p) {
            if (p->typ==polygon) poly_allow_cube(p); // (poly_p)p->load);
            p = p->next;
            if (p == geom_channels) break;
        }
    //}
        
    fprintf(stdout, "after rescaling:\n");
    fprintf(stdout, "min: %7.4f,%7.4f,%7.4f\n", min.k[0], min.k[1],min.k[2]);
    fprintf(stdout, "max: %7.4f,%7.4f,%7.4f\n", max.k[0], max.k[1],max.k[2]);
    genconf();
    //set_positions_velocity();
    update = 1;
    draw_poly = 1;
}



/* ****************************************************** */
/* Routines working on an elem_p with electrode_p load    */
/* ****************************************************** */

/* ****************************************************** */
elem_p elec_malloc(
    char * name,
    elem_p liste)
{
    electrode_p elec = electrode_malloc();
    return elem_add (name,electrode,liste, (void*)elec);
}


/* ****************************************************** */
elem_p elec_new(
    char * name,
    elem_p liste,
    int seq,
    real_t pol,
    int duty,
    bool_t act)
{
    elem_p elem = elec_malloc(name, liste);
    electrode_p elec = (electrode_p) elem->load;
    if (elec) {
        elec->active = act;
        elec->polarity = pol;
        elec->dutycycle = duty;
        elec->seq_ofs = seq;
    }
    return elem;
}


/* ****************************************************** */
elem_p elec_set_color(
    elem_p thisp, 
    int colno)
{
    if (thisp && thisp->typ==electrode && thisp->load) {
        electrode_set_color((electrode_p)thisp->load, colno);
    }
    return thisp;
}    
    
    
/* ****************************************************** */
elem_p elec_add_point(
    elem_p thisp,
    fpoint_t data)
{
    if (thisp && thisp->typ==electrode && thisp->load) {
        electrode_add_point((electrode_p)thisp->load, data);
    }
    return thisp;
}    
    
/* ****************************************************** */
elem_p elec_add_polygon(
    elem_p thisp,
    poly_p poly)
{
    if (thisp && thisp->typ==electrode && thisp->load) {
        electrode_add_polygon((electrode_p)thisp->load, poly);
    }
    return thisp;
}    
    
    



/* ****************************************************** */
/* Routines working on all lists                          */
/* ****************************************************** */


/* ****************************************************** */
void info_polys()
{
    elem_list_print(geom_sensors,1,stdout);
    elem_list_print(geom_electrodes,1,stdout);
    elem_list_print(geom_channels,1,stdout);
}


/* ****************************************************** */
void draw_polys()
{
   static int first=1;
   if (first) {
       first = 0;
       elem_list_print(geom_electrodes,1,stdout);
       elem_list_print(geom_channels,1,stdout);
   }
   if (draw_poly) {
       elem_list_draw(geom_channels, polygon);
       elem_list_draw(geom_sensors, polygon);
       //elem_list_draw(geom_electrodes, polygon);
       //elem_list_draw(geom_electrodes, electrode);
   }
}


/* ****************************************************** */
void init_coords() 
{
//    coord_base = elem_malloc("base",nothing,0);
    geom_channels = elem_malloc("base_channels",nothing,0);
    geom_electrodes = elem_malloc("base_electrodes",nothing,0);;
    geom_sensors = elem_malloc("base_sensors",nothing,0);;
}

void done_coords()
{
    elem_list_free(geom_sensors);
    elem_list_free(geom_electrodes);
    elem_list_free(geom_channels);
}

    
/* ****************************************************** */
/* gcc dpd_geom.c && a.out                                */
/* ****************************************************** */
#if 0
int main() { 
 
  string_p s_p = string_malloc("hallo");
  string_print(s_p, 2,stdout);
  s_p = string_repl(s_p, "holla hier");
  string_print(s_p, 2,stdout);
  s_p = string_delete(s_p);  
  string_print(s_p, 4,stdout);
  s_p = string_free(s_p);
  string_print(s_p, 4,stdout);
  
  string_p s2 = string_malloc(NULL);
  string_print(s2, 2,stdout);
  s2 = string_repl_s(s_p, s2);
  string_print(s2, 2,stdout);
  
  string_p s3 = string_repl(0,"haha");
  string_print(s3, 3,stdout);
  s2 = string_repl_s(s2, s3);
  string_print(s2, 3,stdout);
  
  s3 = string_free(s3);
  string_print(s3, 4,stdout);
  
  elem_p liste = elem_malloc("liste",tag,0);
    elem_print(liste, 2, stdout);

  elem_p e1 =elem_add("haha", tag, liste, 0);
  elem_p e2 =elem_add("huha", tag, liste, 0);
  elem_p e3 =elem_add("hahu", electrode, liste, 0);
  elem_p e4 =elem_add("huhu", sensor, liste, 0);
  elem_p e5 =elem_add("holla", tag, liste, 0);
  
  elem_p fi = elem_list_findname(liste,"hahu");
  elem_print(fi,2,stdout);
      
  e1 = poly_malloc("poly1",liste);
  e2 = poly_malloc("poly2",liste);
  
  poly_add_point(e1, (fpoint_t){{1,2,9}} );
  poly_add_point(e2, (fpoint_t){{2,5,3}} );
  poly_add_point(e1, (fpoint_t){{1,6,7}} );
  poly_add_point(e2, (fpoint_t){{2,4,3}} );
      
    
  elem_list_print(liste,4, stdout);
  elem_free(e3);
  elem_list_print(liste,4, stdout);

  poly_min_max(liste);

#if 0    
    elem_print(liste, 2,stdout);
    elem_print(e1, 2,stdout);
    elem_print(e2, 2,stdout);
    elem_print(e3, 2,stdout);
    elem_print(e4, 2,stdout);
    elem_print(e5, 2,stdout);


  elem_list_free(liste);
  
    elem_print(liste, 2,stdout);
    elem_print(e1, 2,stdout);
    elem_print(e2, 2,stdout);
    elem_print(e3, 2,stdout);
    elem_print(e4, 2,stdout);
    elem_print(e5, 2,stdout);
#endif    
  liste = elem_free(liste);
    elem_print(liste, 2,stdout);
  liste = elem_list_free(liste);
    elem_print(liste, 2,stdout);
  
  return 0; 
}
#endif

/*******************************************************************/

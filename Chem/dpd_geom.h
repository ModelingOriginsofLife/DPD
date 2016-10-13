/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                                             Thomas Maeke  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 28.07.2004  11:50h
  $Id: dpd_geom.h,v 1.6 2004/09/10 07:50:10 tmaeke Exp $
********************************************************************/

#ifndef _dpd_geom_h_
#define _dpd_geom_h_

/*Include***********************************************************/

/*Defines***********************************************************/


/*Types*************************************************************/

typedef enum {
    invalid = -1,
    nothing = 0,
    area, body, channel, electrode, finput, 
    transaction, foutput, polygon, 
    sensor, tag, surface, vertex, 
    max_elem_type
} geom_type_t;

/* ****************************************************** */
typedef struct {
    char *  txt;
    int     len;
} string_t;
typedef string_t *string_p; 

/* ****************************************************** */
typedef struct elem_s elem_t;
struct elem_s {
    int            magic;
    string_p       name;
    geom_type_t    typ;
    elem_t*        prev;
    elem_t*        next;
    void*          load; // payload
};
typedef elem_t *elem_p;
    

/* ****************************************************** */
typedef struct poly_s poly_t;
struct poly_s {
    int       magic;
    int       len;
    int       colno;
    fpoint_t* dat;
};
typedef poly_t *poly_p;


/* ****************************************************** */
typedef struct electrode_s electrode_t;
struct electrode_s {
    int      magic;
    bool_t   active;     // 
    fpoint_t position;   // center of electrode
    fpoint_t size;       // position +/- size/2
    real_t   polarity;   // -1 .. 0 .. 1
    int      dutycycle;  // 0 .. 100%
    int      seq_ofs;    // ofs in seq.list
    poly_p   poly;
    bool_t   ownpoly;    // electrodes ownes poly
};
typedef electrode_t * electrode_p;


/* ****************************************************** */

/*Funcs*************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

// string_p
string_p string_delete(string_p str);
string_p string_repl(string_p str, char *src);
string_p string_repl_s(string_p str, string_p src);
string_p string_malloc(char * src);
string_p string_free(string_p str);

// poly_p
poly_p polygon_malloc(void);
void * polygon_free (poly_p poly);
poly_p polygon_add_point(poly_p poly, fpoint_t data);
poly_p polygon_set_color(poly_p poly, int colno);
void polygon_print(poly_p thisp, int indent, FILE* fi);
void polygon_draw(poly_p poly);
void polygon_min_max(poly_p poly, fpoint_t * min, fpoint_t * max);
void polygon_sub_mul(poly_p poly, fpoint_t * sub, fpoint_t * mul);
bool_t polygon_inside_xy(poly_p poly, fpoint_t data);

// electrode_p
electrode_p electrode_malloc(void);
void * electrode_free (electrode_p elec);
void electrode_print(electrode_p thisp, int indent, FILE* fi);
electrode_p electrode_set_color(electrode_p elec, int colno);
electrode_p electrode_add_point(electrode_p elec, fpoint_t data);
electrode_p electrode_add_polygon(electrode_p elec, poly_p poly);
    


// elem_p
void   elem_print(elem_p thisp, int indent, FILE* fi);
void elem_draw(elem_p thisp);
elem_p elem_malloc(char * name, geom_type_t typ, elem_p liste);
elem_p elem_free(elem_p thisp);
elem_p elem_add(char * name, geom_type_t typ, elem_p liste, void * payload);


// elem_p list
void   elem_list_print(elem_p thisp, int indent, FILE* fi);
elem_p elem_list_free(elem_p thisp);
elem_p elem_list_findname(elem_p thisp, char * name);
elem_p elem_list_findname_qual(elem_p thisp, char * name, geom_type_t qual);
void elem_list_draw(elem_p base, geom_type_t qual);


// elem_p + poly_p
elem_p poly_malloc(char * name, elem_p liste);
elem_p poly_free(elem_p thisp);
elem_p poly_add_point(elem_p thisp, fpoint_t data);
elem_p poly_setcolor(elem_p thisp, int colno);
void poly_min_max(elem_p elem, fpoint_t * min, fpoint_t * max);
void poly_sub_mul(elem_p elem, fpoint_t * min, fpoint_t * max);
void poly_translate(real_t size);
 

// elem_p + electorde_p
elem_p elec_malloc(char * name, elem_p liste);
elem_p elec_new(char * name, elem_p liste, int seq, real_t pol, int duty, bool_t act);
elem_p elec_setcolor(elem_p thisp, int colno);
elem_p elec_add_point(elem_p thisp, fpoint_t data);
elem_p elec_add_polygon(elem_p thisp, poly_p poly);


// all lists
void init_coords();
void done_coords();
void draw_polys();
void info_polys();


#ifdef __cplusplus
}
#endif

/*Vars**************************************************************/

//extern elem_p coord_base;
extern elem_p geom_channels;
extern elem_p geom_electrodes;
extern elem_p geom_sensors;

/*******************************************************************/
#endif  /* _dpd_geom_h_ */

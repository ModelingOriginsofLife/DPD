/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      Protolife                      */
/*                                                                     */
/*  (c) 2004                             Andrew Buchanan               */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

#ifndef _dpd_dynamic_h_
#define _dpd_dynamic_h_

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
	void do_census();
	void reset_census();
	void reset_reaction_list();
	void print_census();
	void print_reactions();

	bool_t particles_should_dynamically_bond(int particle_one, int particle_two);
	int dynamically_bond_particles(int particle_one, int particle_two);
	int dynamically_unbond_particles(int particle_one, int particle_two);
	int get_polymer(int p);
	int reverse_int(int target);
	int traverse_left(int start);
	bool_t joined(int ap, int bp); // return true if the particles are part of the same polymer
	
//	void set_poly_pop(int polymer, int count);
	// privatish functions
	void decrease_poly_pop(int polymer);
	void increase_poly_pop(int polymer);
	void increase_reax_pop(char reaction[]);
//	void decrease_reax_pop(char reaction[]);
	
#ifdef __cplusplus
}
#endif

/*******************************************************************/
#endif  /* _dpd_dynamic_h_ */

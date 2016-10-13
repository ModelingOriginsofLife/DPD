/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      Protolife                      */
/*                                                                     */
/*  (c) 2004                             Andrew Buchanan               */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "dpd_vars.h"
#include "dpd_dynamic.h"
#include "dpd_util.h"

void do_census() {
	// new and improved deadly accurate
	// now count again and this time store it in poly_pop_rec
	reset_census();
	int i;
	for (i=0; i<num_particle; i++) {
		if (liste[i].left < 0 && liste[i].right < 0) {
			// particle is free
			//printf("Increasing pop of %d\n", liste[i].type);
			increase_poly_pop(liste[i].type);
		} else if (liste[i].left >= 0 && liste[i].right < 0) {
			// particle is the end of a polymer
			//printf("Increasing pop of %d\n", get_polymer(i));
			increase_poly_pop(get_polymer(i));
		}
	}	
}

void reset_census() {
	poly_pop_rec *current = poly_pop_list_root;
	poly_pop_rec *last;
	while (current != NULL) {
		last = current;
		current = current->next;
		free(last);
	}
	poly_pop_list_root = NULL;
}

void reset_reaction_list() {
	reaction_record *current = reaction_list_root;
	reaction_record *last;
	while (current != NULL) {
		last=current;
		current = current->next;
		free(last);
	}
	reaction_list_root = NULL;
}

int dynamically_bond_particles(int p_one, int p_two) {
	// first figure out what polymers we're dealing with
	int poly_one = get_polymer(p_one);
	int poly_two = get_polymer(p_two);
	
	//create bond making sure to bond on correct side
	if ((liste[p_one].right == -1) && (liste[p_two].left == -1)) {
		liste[p_one].right = p_two;
		liste[p_two].left = p_one;
	}
	else if ((liste[p_one].left == -1) && (liste[p_two].right == -1)) {
		liste[p_one].left = p_two;
		liste[p_two].right = p_one;
	}
	else {
		return 1; // one of the particles already fully bonded
	}
	
	int final_poly = get_polymer(p_one);   // both p_one and p_two should yield the same now

	char reaction[MAX_REACTION]="";
	int greater = poly_one > poly_two ? poly_one : poly_two;
	int lesser = poly_one > poly_two ? poly_two : poly_one;
	sprintf(reaction, "%d + %d -> %d", greater, lesser, final_poly);
	// printf("just had this reaction: %s\n",reaction);
	increase_reax_pop(reaction);
	
	//update polymer count
	// take away one from each of current polymers
	//printf("decreasing %d and %d\n", poly_one, poly_two);
	decrease_poly_pop(poly_one);
	decrease_poly_pop(poly_two);
	// add one to newly formed polymer
	increase_poly_pop(final_poly);
	
	return 0;
}

void increase_reax_pop(char reaction[]) {
	reaction_record *current = reaction_list_root;
	reaction_record *last=current;
	while (current != NULL && strcmp(current->reaction,reaction)) {
		last=current;
		current = current->next;
	}
	if (current == NULL) {
		// this is the first time we're seeing this reaction
		if ((current = (reaction_record *) calloc(1, sizeof(reaction_record))) == NULL) {
			stopit("Error allocating memory");
		}
		strncpy(current->reaction, reaction,MAX_REACTION);
		current->count = 0;
		current->next = NULL;
		if (last == NULL) { 
			// this must be the first reaction
			reaction_list_root = current;
		} else {
			last->next = current; 
		}
	}
	current->count++;
}

bool_t joined(int ap, int bp) {
	// BUGBUG not loop safe
	int check = ap;
	// check left
	// BUGBUG maybe faster to check both ways a bit since we know they're close
	// enough to each other to be reacting?
	while (liste[check].left >= 0) {
		if (liste[check].left == bp) return TRUE;
		check = liste[check].left;
	}
	check = ap;
	while (liste[check].right >= 0) {
		if (liste[check].right == bp) return TRUE;
		check = liste[check].right;
	}
	
	return FALSE;	
}

int dynamically_unbond_particles(int p_one, int p_two) {
	int poly = get_polymer(p_one); // or p_two
	
	if (liste[p_one].right == p_two) {
		liste[p_one].right = -1;
		liste[p_two].left = -1;
	} else if (liste[p_one].left == p_two) {
		liste[p_one].left = -1;
		liste[p_two].right = -1;
	} else {
		return 1; // particles aren't bonded
	}

	int poly_one = get_polymer(p_one);
	int poly_two = get_polymer(p_two);

	char reaction[MAX_REACTION]="";
	int greater = poly_one > poly_two ? poly_one : poly_two;
	int lesser = poly_one > poly_two ? poly_two : poly_one;
	sprintf(reaction, "%d -> %d + %d", poly, greater, lesser);
	// printf("just had this reaction: %s\n",reaction);
	increase_reax_pop(reaction);
	
	//printf("decreasing %d\n", poly);	
	decrease_poly_pop(poly);
	increase_poly_pop(poly_one);
	increase_poly_pop(poly_two);
	
	return 0;
}

int traverse_left(int start) {
	int leftmost = start;
	while (liste[leftmost].left >= 0 && liste[leftmost].left != start) {
		leftmost = liste[leftmost].left;
	}
	return leftmost;
}

// given a particle, return the polymer of which it's a part of
int get_polymer(int p) {
	int index = traverse_left(p);
	int poly=0;
	// printf("getting polymer for %d, leftmost is %d\n", p, index);
	while (index >= 0) { //BUGBUG not looping polymer safe
		poly = 10 * poly + liste[index].type;
		// printf("found type %d, poly now %d\n",liste[index].type, poly);
		index = liste[index].right;
	}
	int poly_rev = reverse_int(poly); // check the reverse of poly and take whichever's greater
	if (poly < poly_rev) { poly = poly_rev; }
	return poly;
}

int reverse_int(int target) {
	int rev=0;
	int tmp=target;
	while (tmp > 0) {
		rev = 10 * rev + (tmp % 10);
		tmp = tmp / 10;
	}
	return rev;
}

void increase_poly_pop(int polymer) {
	poly_pop_rec *current=poly_pop_list_root;
	poly_pop_rec *last=current;
	while (current != NULL && current->poly != polymer) {
		last=current;
		current = current->next;
	}
	if (current == NULL) {
		if ((current = (poly_pop_rec *) calloc(1, sizeof(poly_pop_rec))) == NULL) {
			stopit("Error allocating memory");
		}
		current->poly = polymer;
		current->count = 0;
		current->next = NULL;
		if (last != NULL) {
			last->next = current;
		}
		else {
			// this must be the first time
			poly_pop_list_root = current;
			//printf("You should only see this message once\n");
			//stopit("Error increasing poly_pop (last == NULL)");
		}
	}
	current->count++;
}

void decrease_poly_pop(int polymer) {
	poly_pop_rec *current=poly_pop_list_root;
	while (current != NULL && current->poly != polymer) {
		current = current->next;
	}
	if (current == NULL || current->count == 0) {
		printf("Current polymer is %d\n", polymer);
		stopit("Error: decrease_poly_pop led to NULL or 0");
	}
	current->count--;
}

bool_t particles_should_dynamically_bond(int p_one, int p_two) {
	/* check that they have complementary open spots (i.e., left-right and not left-left) */
	if ((liste[p_one].right >= 0) && (liste[p_two].right >= 0)) { return FALSE; }
	if ((liste[p_one].left >= 0) && (liste[p_two].left >= 0)) { return FALSE; }
	
	/* check that neither particle is already fully bonded */
	if (((liste[p_one].right >= 0) && (liste[p_one].left >= 0)) || ((liste[p_two].right >= 0) && (liste[p_two].left >= 0))) { return FALSE; }
	
	/* check that they're not already bonded to each other.  No double bonds */
	if ((liste[p_one].right == p_two) || (liste[p_one].left == p_two) || (liste[p_two].right == p_one) || (liste[p_two].left == p_one)) { return FALSE; }
		
	/* check to see how long they are altogether */
	int count=0;
	int old_index_one=-1, old_index_two=-1;  
	int index = p_one;
	while ((index >= 0) && (count <= max_dyn_polylen)) {
		count++;
		old_index_one = index;
		index = liste[p_one].right >= 0 ? liste[index].right : liste[index].left;
		// traverse either left or right depending on which side p_one was already bonded on
	}    
	index = p_two;
	while ((index >= 0) && (count <= max_dyn_polylen)) {
		count++;
		old_index_two = index;
		index = liste[p_two].right >= 0 ? liste[index].right : liste[index].left;
	}
	if (count > max_dyn_polylen) {
		return FALSE;
	}
	
	/* prevent looping if DYN_LOOPS is set to FALSE, otherwise don't worry about it */
	if (dyn_loops == FALSE) {
		if (old_index_two == p_one || old_index_one == p_two) {
			if (!(old_index_two == p_one && old_index_one == p_two)) {
				printf("Error: strange bond found when examining particles %d and %d\n", p_one, p_two);
			}
			return FALSE;
		}
	}
	return TRUE;
}

void print_reactions() {
	reaction_record *current = reaction_list_root;
	while (current != NULL) {
		printf("%s\t%d\n", current->reaction, current->count);
		current = current->next;
	}
}

void print_census() {
	poly_pop_rec *current=poly_pop_list_root;
	while (current != NULL) {
		if (current->count > 0) {
			printf("%d = %d\n", current->poly, current->count);
		}
		current = current->next;
	}
}


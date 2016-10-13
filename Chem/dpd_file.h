/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*                                                                     */
/*  dpd-simulation                                                     */
/*                                      BioMIP Ruhr-Univ. Bochum       */
/*                                                                     */
/*  (c) 2004                             Thomas Maeke, John McCaskill  */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*******************************************************************
  tmm 20.07.2004  08:10h
  $Id: dpd_file.h,v 1.5 2004/08/02 15:11:07 tmaeke Exp $
********************************************************************/

#ifndef _dpd_file_h_
#define _dpd_file_h_

/*Include***********************************************************/

/*Defines***********************************************************/
#define MAX_FILE_LN 10000

/*Funcs*************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

int interpret_line(char * line);
void read_ctrl_file(char * fname);

void template_info(FILE * fi); 
void save_template(char * fname);

void helpinvoke();
void process_cmdline(int argc, char *argv[]);

void kill_bookkeeping_files();
void output_census();
void output_reactions();

#ifdef __cplusplus
}
#endif

            

/*Types*************************************************************/

/*Vars**************************************************************/

extern char * errors[];

/*******************************************************************/
#endif  /* _dpd_file_h_ */

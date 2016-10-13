/*******************************************************************
  tmm 03.08.2004  10:13h
  $Id: dpd_flow.h,v 1.1 2004/08/03 12:15:58 tmaeke Exp ${fn},v 1.10 2003/08/11 09:56:01 tmaeke Exp tmaeke $
********************************************************************/

#ifndef _dpd_flow_h_
#define _dpd_flow_h_

/*Include***********************************************************/

/*Defines***********************************************************/

#define MAX_IO_FLOW 3      // maxinum no of in-/output-flows

/*Types*************************************************************/
    
    typedef struct {
        fpoint_t  mi, ma;    // rectangular area 
        fpoint_t  mid;       // middle of mi,ma
        int       rate;      // inflow-rate, particles/dt
    } ioflow_t;
    

/*Vars**************************************************************/
    
    extern bool_t   inout_flow;
    extern int      inflows;
    extern int      outflows;
    extern ioflow_t inflow_list[MAX_IO_FLOW]; 
    extern ioflow_t outflow_list[MAX_IO_FLOW];
    extern index_t  outside_chain;
    

/*Funcs*************************************************************/

    void add_inflow(char * tag, int i, real_t rate);
    void add_outflow(char * tag, int i, real_t rate);

    bool_t test_outflow(fpoint_t * where);
    void flow_out(index_t ip);
    int flow_in(int inflow);
    bool_t flow_in_all();
    
    


/*******************************************************************/
#endif  /* _dpd_flow_h_ */

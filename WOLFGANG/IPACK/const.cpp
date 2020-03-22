#include "io_incl.h"
#include <assert.h>
#include <xtime.h>
#include "const.h"

void init_timer()
{
  timer -> set_name(    t_main,     "main");
  timer -> set_name(     t_vrr,      "vrr");
  timer -> set_name(    t_hrr1,     "hrr1");
  timer -> set_name(    t_hrr2,     "hrr2");
  timer -> set_name(t_contract, "contract");
  timer -> set_name(   t_addto,   "add_to");
  timer -> set_name( t_doublet,  "doublet");
  timer -> set_name(      t_sp,"single-pt");
  timer -> set_name(       t_f,"   f(T,m)");
  timer -> set_name(t_2p_total, "2P-total");
  timer -> set_name(t_2p_output, "2P-out");
  timer -> set_name(t_2p_disth, "2P-disth");
  timer -> set_name(t_2p_distv, "2P-distv");


};

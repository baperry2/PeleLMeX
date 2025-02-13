#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("T_0", prob_parm->T_0);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("T_in", prob_parm->T_in);
  pp.query("u_in", prob_parm->u_in);
  pp.query("u_0", prob_parm->u_0);
}

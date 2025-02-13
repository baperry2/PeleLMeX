#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  // Shared params
  amrex::Real T_0 = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T_in = 400.0_rt;
  amrex::Real u_in = 10.0;
  amrex::Real u_0 = 0.0;
};
#endif

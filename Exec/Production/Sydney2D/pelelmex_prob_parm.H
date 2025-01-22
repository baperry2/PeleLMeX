#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <PelePhysics.H>
//#include "mechanism.H"

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean    = 101325.0_rt;
  amrex::Real v_jet     = 57.0_rt;
  amrex::Real v_pilot   = 15.0_rt;
  amrex::Real v_coflow  = 10.0_rt;
  amrex::Real Do_jet    = 0.0075_rt;
  amrex::Real Di_pilot  = 0.0080_rt;
  amrex::Real Do_pilot  = 0.0100_rt;
  amrex::Real Di_coflow = 0.0110_rt;
  amrex::Real Y_jet[NUM_SPECIES] = {0.0};
  amrex::Real Y_pilot[NUM_SPECIES] = {0.0};
  amrex::Real Y_coflow[NUM_SPECIES] = {0.0};
  amrex::Real center_xy[2] = {0.0};
  int zstar_pos = 1;

  int turb_inflow_type = 0;  
  int nr;
  int nt;
  int nvar;
  int inflowNtime;
  amrex::Real inflowFreq;
  amrex::Real timeFromInflow;
  amrex::Real timeInflowMax;
  amrex::Real thetaMax;
  amrex::Real* d_timeInput;
  amrex::Real* d_rM;
  amrex::Real* d_thetaM;
  amrex::Real* d_Ut;
  amrex::Real* d_Ur;
  amrex::Real* d_Uz;
  amrex::Real* d_Zmix;
  pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type> const*
    eosparm;
};
#endif

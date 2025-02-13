#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  // Shared params
  amrex::Real T_mean = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real Y_fuel = 0.0445_rt;
  amrex::Real Y_o2 = 0.223_rt;
  amrex::Real T_hot = 400.0_rt;
  amrex::Real Twall = 300.0_rt;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real EBinflow_vel = 20.0;
  amrex::Real EBinflow_Yfuel = 0.5;
  amrex::Real EBinflow_T = 320.0;

  amrex::Real d_cyl = 7.542;
  amrex::Real l_nozzle = 4.789;
  amrex::Real d_nozzle = 1.780;
  amrex::Real d_hole = 0.5;
  amrex::GpuArray<amrex::Real, 3> l_holes = {
    5.184, 6.5 + 5.184, 4.846 + 5.184 + 6.5}; // 5.184, 11.684, 16.53
  amrex::GpuArray<int, 3> n_holes = {6, 12, 12};
  amrex::GpuArray<amrex::Real, 3> offset_angle_holes = {
    0.0, 0.0, 10 * M_PI / 180};
  amrex::Real u_hole = 0.0;
  amrex::Real T_hole = 300.0;
  amrex::Real l_combustor = 23.711;
  int n_jets = 6;
  amrex::Real d_jet = 0.35;
  amrex::Real umag_jet = 0.0;
  amrex::Real angle_jet = 0.0; // in degrees
  amrex::Real T_jet = 300.0;
  amrex::Real t_jet_scale = 1e-3;
  amrex::Real jet_ch4 = 1.0;
  amrex::Real d_cone = 2.2;
  amrex::Real l_cone = 1.118;
  amrex::Real l_exit = 1.0;
  amrex::Real d_swirler = 2.8;
  amrex::Real u_swirl_ax = 0.0;
  amrex::Real u_swirl_az = 0.0;
  amrex::Real T_swirl = 300.0;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real l_annulus = 0.0;
  amrex::Real r_annulus_inner = 0.0;
  amrex::Real r_annulus_outer = 0.0;
  int outlet_type = 1;
  amrex::Real p_outlet = p0;
  bool do_inlet_cone = false;
  bool reinit_dilution_holes = false;
  bool reinit_mixture_fractions = false;
};
#endif

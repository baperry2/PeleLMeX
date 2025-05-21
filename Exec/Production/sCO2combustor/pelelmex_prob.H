#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <PMFData.H>
#include <PelePhysics.H>
#include <PeleLMeX_ProblemSpecificFunctions.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#define NROWS 3

struct GeomParms
{
  amrex::Real d_cyl = 0.07542;
  amrex::Real l_nozzle = 0.04789;
  amrex::Real d_nozzle = 0.025;
  amrex::Real d_hole = 0.005;
  amrex::GpuArray<amrex::Real, NROWS> l_holes = {
    0.05184, 0.065 + 0.05184,
    0.04846 + 0.05184 + 0.065}; // 5.184, 11.684, 16.53
  amrex::GpuArray<int, NROWS> n_holes = {6, 12, 12};
  amrex::GpuArray<amrex::Real, NROWS> offset_angle_holes = {
    0.0, 0.0, 10 * M_PI / 180};
  amrex::Real l_combustor = 0.23711;
  amrex::Real l_exit = 0.010;
  int do_dilution_holes = 4;
};

struct MyProbParm : public ProbParmDefault
{
  // Shared params
  amrex::Real T_0 = 298.0;
  amrex::Real P_mean = 101325.0;
  amrex::Real T_in = 400.0;
  amrex::Real u_in = 10.0;
  amrex::Real u_0 = 0.0;
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata(
    int i,
    int j,
    int k,
    int is_incompressible,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();

    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    massfrac[O2_ID] = 0.23;
    massfrac[N2_ID] = 1.0 - massfrac[O2_ID];

    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.u_0;
                 , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);

    if (!(is_incompressible == 0)) {
      return;
    }

    state(i, j, k, TEMP) = prob_parm.T_0;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
    state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
    }
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    if (sgn == 1) {
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      massfrac[O2_ID] = 0.23;
      massfrac[N2_ID] = 1.0 - massfrac[O2_ID];

      auto eos = pele::physics::PhysicsType::eos();

      AMREX_D_TERM(s_ext[VELX] = prob_parm.u_in;, s_ext[VELY] = 0.0;
                   , s_ext[VELZ] = 0.0);

      s_ext[TEMP] = prob_parm.T_in;

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] =
        RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
      }
    }
  }

  static void parse_params(GeomParms* prob_parm_device)
  {
    // Parse params
    amrex::ParmParse pp("prob");
    pp.query("d_cyl", prob_parm_device->d_cyl);
    pp.query("l_nozzle", prob_parm_device->l_nozzle);
    pp.query("d_nozzle", prob_parm_device->d_nozzle);
    amrex::Vector<amrex::Real> l_holes(NROWS, 0.0);
    for (int i = 0; i < NROWS; i++) {
      l_holes[i] = prob_parm_device->l_holes[i];
    }
    pp.queryarr("l_holes", l_holes, 0, NROWS);
    for (int i = 0; i < NROWS; i++) {
      prob_parm_device->l_holes[i] = l_holes[i];
    }
    pp.query("d_hole", prob_parm_device->d_hole);
    pp.query("l_combustor", prob_parm_device->l_combustor);
    pp.query("l_exit", prob_parm_device->l_exit);
    pp.query("do_dilution_holes", prob_parm_device->do_dilution_holes);
  }

  static void EBUserDefined(
    const amrex::Geometry& geom,
    const int required_coarsening_level,
    const int max_coarsening_level)
  {
    // ParmParse your geometry parameters

    // Build geometry pieces using EB2::* methods

    // Build your geometry shop using EB2::makeShop

    // Build geom using EB2::Build

    // We shoulnd't be here, copy this file in your run folder
    // and implement your geometry

    GeomParms gparms;
    parse_params(&gparms);
    GeomParms const* pp = &gparms;

    amrex::EB2::CylinderIF main(
      0.5 * pp->d_cyl, 2.0 * pp->l_combustor, 0,
      {AMREX_D_DECL(0.5 * pp->l_combustor, 0, 0)}, true);

    amrex::RealArray center_hole = {AMREX_D_DECL(0.0, 0.5 * pp->d_cyl, 0)};
    amrex::EB2::CylinderIF hole(
      0.5 * pp->d_hole, 2.0 * pp->d_hole, 1, center_hole, true);

    amrex::Vector<amrex::Real> angles_holes[NROWS];
    for (int i = 0; i < NROWS; i++) {
      angles_holes[i].resize(pp->n_holes[i]);
      for (int n = 0; n < pp->n_holes[i]; n++) {
        angles_holes[i][n] =
          n * 360.0 / pp->n_holes[i] * M_PI / 180.0 + pp->offset_angle_holes[i];
      }
    }

    amrex::RealArray offset_holes[NROWS] = {
      {AMREX_D_DECL(pp->l_holes[0], 0, 0)},
      {AMREX_D_DECL(pp->l_holes[1], 0, 0)},
      {AMREX_D_DECL(pp->l_holes[2], 0, 0)}};

    auto first_hole = amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][0], 0), offset_holes[0]);
    auto first_row = amrex::EB2::makeIntersection(
      first_hole,
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[0][1], 0), offset_holes[0]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[0][2], 0), offset_holes[0]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[0][3], 0), offset_holes[0]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[0][4], 0), offset_holes[0]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[0][5], 0), offset_holes[0]));
    auto all_holes = amrex::EB2::makeIntersection(
      first_row,
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][0], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][1], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][2], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][3], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][4], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][5], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][6], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][7], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][8], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][9], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][10], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[1][11], 0), offset_holes[1]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][0], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][1], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][2], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][3], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][4], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][5], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][6], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][7], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][8], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][9], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][10], 0), offset_holes[2]),
      amrex::EB2::translate(
        amrex::EB2::rotate(hole, angles_holes[2][11], 0), offset_holes[2]));

    amrex::Real slope_nozzle =
      (0.5 * pp->d_cyl - 0.5 * pp->d_nozzle) / pp->l_nozzle;
    amrex::Real norm = -1.0 / slope_nozzle;
    amrex::Real nmag = std::sqrt(1 + 1 / (norm * norm));
#if AMREX_SPACEDIM == 2
    amrex::EB2::PlaneIF nozzle_plane(
      {AMREX_D_DECL(0, 0, 0)},
      {AMREX_D_DECL(slope_nozzle / nmag, 1 / nmag, 0.0)}, true);
    auto nozzle = amrex::EB2::translate(
      nozzle_plane,
      {AMREX_D_DECL(
        pp->l_combustor + 0.5 * pp->d_nozzle / slope_nozzle, 0, 0)});
#elif AMREX_SPACEDIM == 3
    amrex::EB2::PlaneIF nozzle_plane(
      {AMREX_D_DECL(0, 0, 0)},
      {AMREX_D_DECL(1 / nmag, slope_nozzle / nmag, 0.0)}, true);
    auto nozzle = amrex::EB2::translate(
      amrex::EB2::rotate(amrex::EB2::lathe(nozzle_plane), 90 * M_PI / 180, 1),
      {AMREX_D_DECL(
        pp->l_combustor + 0.5 * pp->d_nozzle / slope_nozzle, 0, 0)});
#endif

    amrex::EB2::CylinderIF exit(
      0.5 * pp->d_nozzle, 4 * pp->l_exit, 0,
      {AMREX_D_DECL(pp->l_combustor, 0, 0)}, true);
    auto nozzle_exit = amrex::EB2::makeIntersection(nozzle, exit);

    if (pp->do_dilution_holes == 0) {
      auto polys = amrex::EB2::makeUnion(main, nozzle_exit);
      auto gshop = amrex::EB2::makeShop(polys);
      amrex::EB2::Build(
        gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
    } else if (pp->do_dilution_holes == 1) {
      auto comb = amrex::EB2::makeIntersection(main, first_hole);
      auto polys = amrex::EB2::makeUnion(comb, nozzle_exit);
      auto gshop = amrex::EB2::makeShop(polys);
      amrex::EB2::Build(
        gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
    } else if (pp->do_dilution_holes == 2) {
      auto comb = amrex::EB2::makeIntersection(main, first_row);
      auto polys = amrex::EB2::makeUnion(comb, nozzle_exit);
      auto gshop = amrex::EB2::makeShop(polys);
      amrex::EB2::Build(
        gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
    } else if (pp->do_dilution_holes == 3) {
      auto comb = amrex::EB2::makeIntersection(main, all_holes);
      auto polys = amrex::EB2::makeUnion(comb, nozzle_exit);
      auto gshop = amrex::EB2::makeShop(polys);
      amrex::EB2::Build(
        gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
    } else {
      amrex::Abort("Invalid input option");
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;

#endif

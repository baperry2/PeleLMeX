#ifndef EBUSERDEFINED_H
#define EBUSERDEFINED_H

#define NROWS 3

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(
  const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
  amrex::Real& EBflagType,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  EBflagType = 0.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(
  const amrex::Real xEBface[AMREX_SPACEDIM],
  amrex::Real s_ext[NVAR],
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  amrex::Abort("Should not be here");
}

struct GeomParms
{
  amrex::Real d_cyl = 7.542;
  amrex::Real l_nozzle = 4.789;
  amrex::Real d_nozzle = 1.780;
  amrex::Real d_hole = 0.5;
  amrex::GpuArray<amrex::Real, NROWS> l_holes = {
    5.184, 6.5 + 5.184, 4.846 + 5.184 + 6.5}; // 5.184, 11.684, 16.53
  amrex::GpuArray<int, NROWS> n_holes = {6, 12, 12};
  amrex::GpuArray<amrex::Real, NROWS> offset_angle_holes = {
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

  int turb_inflow_type = 0;
  int inflowNtime = 0;
  int nr = 0;
  int nt = 0;
  amrex::Real inflowFreq = 0.0;
  amrex::Real timeFromInflow = 0.0;
  amrex::Real timeInflowMax = 0.0;
  amrex::Real rMax = 0.0;
  amrex::Real thetaMax = 0.0;
  amrex::Real vInflowBulk = 30.0; // m/s
  amrex::Real* d_timeInput = nullptr;
  amrex::Real* d_rM = nullptr;
  amrex::Real* d_thetaM = nullptr;
  amrex::Real* d_Uz = nullptr;
  amrex::Real* d_Ur = nullptr;
  amrex::Real* d_Ut = nullptr;
  amrex::Real* d_jet_time_offsets = nullptr;
};

void
parse_params(GeomParms* prob_parm_device)
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
  pp.query("u_hole", prob_parm_device->u_hole);
  pp.query("T_hole", prob_parm_device->T_hole);
  pp.query("l_combustor", prob_parm_device->l_combustor);
  pp.query("n_jets", prob_parm_device->n_jets);
  pp.query("d_jet", prob_parm_device->d_jet);
  pp.query("umag_jet", prob_parm_device->umag_jet);
  pp.query("angle_jet", prob_parm_device->angle_jet);
  pp.query("T_jet", prob_parm_device->T_jet);
  pp.query("t_jet_scale", prob_parm_device->t_jet_scale);
  pp.query("jet_ch4", prob_parm_device->jet_ch4);
  pp.query("d_cone", prob_parm_device->d_cone);
  pp.query("l_cone", prob_parm_device->l_cone);
  pp.query("l_exit", prob_parm_device->l_exit);
  pp.query("d_swirler", prob_parm_device->d_swirler);
  pp.query("u_swirl_ax", prob_parm_device->u_swirl_ax);
  pp.query("u_swirl_az", prob_parm_device->u_swirl_az);
  pp.query("T_swirl", prob_parm_device->T_swirl);
  pp.query("p0", prob_parm_device->p0);
  pp.query("T0", prob_parm_device->T0);
  pp.query("do_inlet_cone", prob_parm_device->do_inlet_cone);
  pp.query("reinit_dilution_holes", prob_parm_device->reinit_dilution_holes);
  pp.query(
    "reinit_mixture_fractions", prob_parm_device->reinit_mixture_fractions);
  pp.query("outlet_type", prob_parm_device->outlet_type);
  pp.query("p_outlet", prob_parm_device->p_outlet);

  // Info for annulus cutting into cone
  const amrex::Real slope_cone =
    0.5 * prob_parm_device->d_cone / prob_parm_device->l_cone;
  const amrex::Real p = 0.5 * prob_parm_device->d_cone;
  const amrex::Real area_jets = 6.0 * M_PI * (0.5 * prob_parm_device->d_jet) *
                                (0.5 * prob_parm_device->d_jet);
  // prob_parm_device->l_annulus = 0.5 * prob_parm_device->d_jet;
  prob_parm_device->l_annulus =
    area_jets /
    (0.5 * M_PI * slope_cone * p -
     0.25 * M_PI * prob_parm_device->l_cone * slope_cone * slope_cone);
  prob_parm_device->r_annulus_inner =
    -slope_cone *
      (0.5 * prob_parm_device->l_cone + 0.5 * prob_parm_device->l_annulus) +
    p;
  prob_parm_device->r_annulus_outer =
    -slope_cone *
      (0.5 * prob_parm_device->l_cone - 0.5 * prob_parm_device->l_annulus) +
    p;
}

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(
  const Geometry& geom,
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

  auto comb = amrex::EB2::makeIntersection(
    main,
    amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][0], 0), offset_holes[0]),
    amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][1], 0), offset_holes[0]),
    amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][2], 0), offset_holes[0]),
    amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][3], 0), offset_holes[0]),
    amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][4], 0), offset_holes[0]),
    amrex::EB2::translate(
      amrex::EB2::rotate(hole, angles_holes[0][5], 0), offset_holes[0]),
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

  const amrex::Real slope_cone = 0.5 * pp->d_cone / pp->l_cone;
  amrex::Real norm = -1.0 / slope_cone;
  amrex::Real nmag = std::sqrt(1 + 1 / (norm * norm));
  amrex::EB2::PlaneIF cone_plane(
    {AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(1 / nmag, slope_cone / nmag, 0.0)},
    false);
  auto cone = amrex::EB2::translate(
    amrex::EB2::rotate(amrex::EB2::lathe(cone_plane), 90 * M_PI / 180, 1),
    {AMREX_D_DECL(pp->l_cone, 0, 0)});

  // Jet at angle to cone, fails due to multi-cut cells
  // amrex::EB2::CylinderIF jet(0.5 * pp->d_jet, 2.0 * pp->d_jet, 1, {0.0, 0,
  // 0}, true); const amrex::Real angle_jet = -std::atan2(0.5 * pp->d_cone,
  // pp->l_cone); auto polys = amrex::EB2::makeIntersection(
  //   cone, amrex::EB2::translate(amrex::EB2::rotate(jet, angle_jet, 2), {0,
  //   0.25*pp->d_cone,0}));

  // Annulus cutting into cone
  amrex::EB2::CylinderIF inner_annulus(
    pp->r_annulus_inner, pp->l_annulus, 0,
    {AMREX_D_DECL(0.5 * pp->l_cone, 0.0, 0.0)}, false);
  amrex::EB2::CylinderIF outer_annulus(
    pp->r_annulus_outer, pp->l_annulus, 0,
    {AMREX_D_DECL(0.5 * pp->l_cone, 0.0, 0.0)}, true);
  auto annulus = amrex::EB2::makeUnion(inner_annulus, outer_annulus);

  // Swirler before cone
  amrex::EB2::CylinderIF swirler(
    0.5 * pp->d_cone, 10.0, 0, {AMREX_D_DECL(0, 0, 0)}, false);

  // Wall above swirler section
  const amrex::Real* problo = geom.ProbLo();
  amrex::EB2::CylinderIF inner(
    0.5 * pp->d_swirler, 2.0 * (0.0 - problo[0]), 0,
    {AMREX_D_DECL(problo[0], 0, 0)}, true);
  amrex::EB2::CylinderIF outer(
    0.5 * pp->d_cyl, 2.0 * (0.0 - problo[0]), 0,
    {AMREX_D_DECL(problo[0], 0, 0)}, false);
  auto left_wall = amrex::EB2::makeIntersection(inner, outer);

  // Full inlet
  auto inlet = amrex::EB2::makeUnion(
    amrex::EB2::makeIntersection(cone, annulus, swirler), left_wall);

  amrex::Real slope_nozzle =
    (0.5 * pp->d_cyl - 0.5 * pp->d_nozzle) / pp->l_nozzle;
  norm = -1.0 / slope_nozzle;
  nmag = std::sqrt(1 + 1 / (norm * norm));
#if AMREX_SPACEDIM == 2
  amrex::EB2::PlaneIF nozzle_plane(
    {AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(slope_nozzle / nmag, 1 / nmag, 0.0)},
    true);
  auto nozzle = amrex::EB2::translate(
    nozzle_plane,
    {AMREX_D_DECL(pp->l_combustor + 0.5 * pp->d_nozzle / slope_nozzle, 0, 0)});
#elif AMREX_SPACEDIM == 3
  amrex::EB2::PlaneIF nozzle_plane(
    {AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(1 / nmag, slope_nozzle / nmag, 0.0)},
    true);
  auto nozzle = amrex::EB2::translate(
    amrex::EB2::rotate(amrex::EB2::lathe(nozzle_plane), 90 * M_PI / 180, 1),
    {AMREX_D_DECL(pp->l_combustor + 0.5 * pp->d_nozzle / slope_nozzle, 0, 0)});
#endif

  amrex::EB2::CylinderIF exit(
    0.5 * pp->d_nozzle, 4 * pp->l_exit, 0,
    {AMREX_D_DECL(pp->l_combustor, 0, 0)}, true);
  auto nozzle_exit = amrex::EB2::makeIntersection(nozzle, exit);

  if (pp->do_inlet_cone) {
    auto polys = amrex::EB2::makeUnion(comb, inlet, nozzle_exit);
    auto gshop = amrex::EB2::makeShop(polys);
    amrex::EB2::Build(
      gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
  } else {
    auto polys = amrex::EB2::makeUnion(comb, nozzle_exit);
    auto gshop = amrex::EB2::makeShop(polys);
    amrex::EB2::Build(
      gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
  }
}
#endif
#endif

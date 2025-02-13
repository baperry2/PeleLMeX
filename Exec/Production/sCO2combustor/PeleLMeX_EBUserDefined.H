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
  amrex::Real l_combustor = 23.711;
  amrex::Real l_exit = 1.0;
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
  pp.query("l_combustor", prob_parm_device->l_combustor);
  pp.query("l_exit", prob_parm_device->l_exit);
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

  amrex::Real slope_nozzle =
    (0.5 * pp->d_cyl - 0.5 * pp->d_nozzle) / pp->l_nozzle;
  amrex::Real norm = -1.0 / slope_nozzle;
  amrex::Real nmag = std::sqrt(1 + 1 / (norm * norm));
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

  auto polys = amrex::EB2::makeUnion(comb, nozzle_exit);
  auto gshop = amrex::EB2::makeShop(polys);
  amrex::EB2::Build(
    gshop, geom, required_coarsening_level, max_coarsening_level, 4, false);
}
#endif
#endif

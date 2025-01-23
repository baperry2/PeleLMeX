#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

// ----------------------------------------------------------- 
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

AMREX_FORCE_INLINE
void
read_binary_int(std::ifstream& infile, int* store_int, int n_int=1 )
{
  infile.read(reinterpret_cast<char*>(store_int), n_int*sizeof(int));
}

AMREX_FORCE_INLINE
void
read_binary_double(std::ifstream& infile, double* store_double, int n_double=1)
{
  infile.read(reinterpret_cast<char*>(store_double), n_double*sizeof(double));
}

AMREX_FORCE_INLINE
void
read_binary_strshort(std::ifstream& infile)
{
  char tmp[8];
  infile.read(reinterpret_cast<char*>(&tmp), 8*sizeof(char));
  amrex::Print() << tmp << std::endl;
  return;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
InterpolateVelInflow(
  const amrex::Real radius,
  amrex::Real& theta,
  amrex::Real& u,
  amrex::Real& v,
  amrex::Real& w,
  amrex::Real& Zmix,
  const amrex::Real timeInp,
  const ProbParm& prob_parm)
{
  // Assume radius passed is in m
  amrex::Real radSI = radius;
  // Get radius location for interpolation
  int idR = 0;
  locate(prob_parm.d_rM, prob_parm.nr, radSI, idR);
  int idRp1 = idR + 1;

  theta = std::fmod(theta, prob_parm.thetaMax);
  int idT = 0;
  locate(prob_parm.d_thetaM, prob_parm.nt, theta, idT);
  int idTp1 = idT + 1;

  amrex::Real timeInflow = std::fmod(timeInp, prob_parm.timeInflowMax);
  int indxTime = 0;
  locate(prob_parm.d_timeInput, prob_parm.inflowNtime, timeInflow, indxTime);
  const int indxTimeP1 = indxTime + 1;

  // Interpolate in space
  amrex::Real fR1 = 0.0, fR0 = 0.0;
  if (radSI <= prob_parm.d_rM[0]) {
    fR1 = 0.0;
    fR0 = 1.0 - fR1;
  } else if (radSI >= prob_parm.d_rM[prob_parm.nr - 1]) {
    fR1 = 1.0;
    fR0 = 0.0;
    idRp1 = prob_parm.nr - 1;
    idR = prob_parm.nr - 2;
  } else {
    fR1 = (radSI - prob_parm.d_rM[idR]) /
          (prob_parm.d_rM[idRp1] - prob_parm.d_rM[idR]);
    fR0 = 1.0 - fR1;
  }

  amrex::Real fT1 = 0.0, fT0 = 0.0;
  // theta grid is periodic with period thetaMax, should already be constrained to that period
  AMREX_ASSERT(theta >= 0.0);
  AMREX_ASSERT(theta <= prob_parm.thetaMax);
  if (theta <= prob_parm.d_thetaM[0] or theta >= prob_parm.d_thetaM[prob_parm.nt -1]) {
    idT = 0;
    idTp1 = prob_parm.nt - 1;
    amrex::Real dtheta = prob_parm.d_thetaM[0] - prob_parm.d_thetaM[prob_parm.nt -1] + prob_parm.thetaMax;
    if (theta <= prob_parm.d_thetaM[0]) {
      fT1 = (prob_parm.d_thetaM[0] - theta)/dtheta;
      fT0 = 1.0 - fT1;
    } else {
      fT0 = (theta - prob_parm.d_thetaM[prob_parm.nt -1])/dtheta;
      fT1 = 1.0 - fT0;
    }
  } else {
    fT1 = (theta - prob_parm.d_thetaM[idT]) /
      (prob_parm.d_thetaM[idTp1] - prob_parm.d_thetaM[idT]);
    fT0 = 1.0 - fT1;
  }

  // Interpolate in time
  amrex::Real fS1 = 0.0, fS0 = 0.0;
  if (timeInflow <= prob_parm.d_timeInput[0]) {
    fS1 = 0.0;
    fS0 = 1.0 - fS1;
  } else {
    fS1 = (timeInflow - prob_parm.d_timeInput[indxTime]) /
          (prob_parm.d_timeInput[indxTimeP1] - prob_parm.d_timeInput[indxTime]);
    fS0 = 1.0 - fS1;
  }

  // Now Interpolate
  const amrex::Real ur_dns =
    fR0 *
      (fT0 *
         prob_parm.d_Ur[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 * prob_parm
               .d_Ur[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Ur[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 * prob_parm
               .d_Ur[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);
  const amrex::Real ut_dns =
    fR0 *
      (fT0 *
         prob_parm.d_Ut[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 * prob_parm
               .d_Ut[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Ut[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 * prob_parm
               .d_Ut[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);
  const amrex::Real uz_dns =
    fR0 *
      (fT0 *
         prob_parm.d_Uz[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 * prob_parm
               .d_Uz[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Uz[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 * prob_parm
               .d_Uz[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);
  const amrex::Real  zmix_dns =
    fR0 *
      (fT0 *
       prob_parm.d_Zmix[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 * prob_parm
               .d_Zmix[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Zmix[(indxTime * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 * prob_parm
               .d_Zmix[(indxTime * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);

  const amrex::Real ur_dns_1 =
    fR0 *
      (fT0 * prob_parm
               .d_Ur[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 *
         prob_parm
           .d_Ur[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Ur[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 *
         prob_parm
           .d_Zmix[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);
  const amrex::Real ut_dns_1 =
    fR0 *
      (fT0 * prob_parm
               .d_Ut[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 *
         prob_parm
           .d_Ut[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Ut[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 *
         prob_parm
           .d_Ut[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);
  const amrex::Real uz_dns_1 =
    fR0 *
      (fT0 * prob_parm
               .d_Uz[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 *
         prob_parm
           .d_Uz[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Uz[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 *
         prob_parm
           .d_Uz[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);
  const amrex::Real zmix_dns_1 =
    fR0 *
      (fT0 * prob_parm
               .d_Zmix[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idR] +
       fT1 *
         prob_parm
           .d_Zmix[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idR]) +
    fR1 *
      (fT0 * prob_parm
               .d_Zmix[(indxTimeP1 * prob_parm.nt + idT) * prob_parm.nr + idRp1] +
       fT1 *
         prob_parm
           .d_Zmix[(indxTimeP1 * prob_parm.nt + idTp1) * prob_parm.nr + idRp1]);

  const amrex::Real uzInterp = fS0 * uz_dns + fS1 * uz_dns_1;
  const amrex::Real urInterp = fS0 * ur_dns + fS1 * ur_dns_1;
  const amrex::Real utInterp = fS0 * ut_dns + fS1 * ut_dns_1;
  const amrex::Real zmixInterp = fS0 * zmix_dns + fS1 * zmix_dns_1;

  // Transform to cartesian coordinates
  u = (urInterp * cos(theta) - radSI * utInterp * sin(theta));
  v = (urInterp * sin(theta) + radSI * utInterp * cos(theta));
  w = uzInterp;
  Zmix = zmixInterp;
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata (int i, int j, int k,
                 int is_incompressible,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& aux,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    auto eos = pele::physics::PhysicsType::eos(prob_parm.eosparm);

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    for (int n = 0; n < NUM_SPECIES; n++)
    {
        massfrac[n] = prob_parm.Y_coflow[n];
    }

    state(i,j,k,VELX) = 0.0;
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = prob_parm.v_coflow;

    amrex::Real rho_cgs;
    amrex::Real dummy = 0.0;

    eos.PYT2R(dummy, massfrac, dummy, rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.HY2T(dummy, massfrac, state(i,j,k,TEMP));
    state(i,j,k,RHOH) = 0.0;   // No H for manifold

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real *x, /*x[AMREX_SPACEDIM]*/
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  if (sgn == 1 && idir == 2) {
    const amrex::Real dummy = 0.0;
    amrex::Real rho_cgs;
    auto eos = pele::physics::PhysicsType::eos(prob_parm.eosparm);

    // Jet function
    amrex::Real radius = std::sqrt(std::pow(x[0]-prob_parm.center_xy[0],2) + std::pow(x[1]-prob_parm.center_xy[1],2));
    amrex::Real angle = 
          3.1415926535897932384626433 - atan2(x[1] - prob_parm.center_xy[1], -(x[0] - prob_parm.center_xy[0]));
    if ( radius < prob_parm.Do_jet*0.5 ) {
      // Jet
      if (prob_parm.turb_inflow_type == 1) {
        amrex::Real Y_jet[NUM_SPECIES];
        InterpolateVelInflow(radius, angle, s_ext[VELX], s_ext[VELY], s_ext[VELZ], Y_jet[prob_parm.zstar_pos], time, prob_parm);
	// Ignore Z from inflow data
        for (int n = 0; n < NUM_SPECIES; n++) {
          Y_jet[n] = prob_parm.Y_jet[n];
        }
        eos.PYT2R(dummy, Y_jet, dummy, rho_cgs);
        s_ext[DENSITY] = rho_cgs * 1.0e3;
        eos.HY2T(dummy, Y_jet, s_ext[TEMP]);
        s_ext[RHOH] = 0.0;   // No H for manifold
        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[FIRSTSPEC+n] = Y_jet[n] * s_ext[DENSITY];
        }
      } else {
        s_ext[VELX] = 0.0;
        s_ext[VELY] = 0.0;
        s_ext[VELZ] = prob_parm.v_jet;
        eos.PYT2R(dummy, prob_parm.Y_jet, dummy, rho_cgs);
        s_ext[DENSITY] = rho_cgs * 1.0e3;
        eos.HY2T(dummy, prob_parm.Y_jet, s_ext[TEMP]);
        s_ext[RHOH] = 0.0;   // No H for manifold
        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[FIRSTSPEC+n] = prob_parm.Y_jet[n] * s_ext[DENSITY];
        }
      }
    } else if (radius >= 0.5*prob_parm.Di_pilot && radius < 0.5*prob_parm.Do_pilot) {
      // Pilot
      s_ext[VELX] = 0.0;
      s_ext[VELY] = 0.0;
      s_ext[VELZ] = prob_parm.v_pilot;
      eos.PYT2R(dummy, prob_parm.Y_pilot, dummy, rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;
      eos.HY2T(dummy, prob_parm.Y_pilot, s_ext[TEMP]);
      s_ext[RHOH] = 0.0;   // No H for manifold
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC+n] = prob_parm.Y_pilot[n] * s_ext[DENSITY];
      }
    } else if (radius >= 0.5*prob_parm.Di_coflow) {
      // Coflow
      s_ext[VELX] = 0.0;
      s_ext[VELY] = 0.0;
      s_ext[VELZ] = prob_parm.v_coflow;
      eos.PYT2R(dummy, prob_parm.Y_coflow, dummy, rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;
      eos.HY2T(dummy, prob_parm.Y_coflow, s_ext[TEMP]);
      s_ext[RHOH] = 0.0;   // No H for manifold
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC+n] = prob_parm.Y_coflow[n] * s_ext[DENSITY];
      }
    } else {
      // Wall
      s_ext[VELX] = 0.0;
      s_ext[VELY] = 0.0;
      s_ext[VELZ] = 0.0;
      eos.PYT2R(dummy, prob_parm.Y_coflow, dummy, rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;
      eos.HY2T(dummy, prob_parm.Y_coflow, s_ext[TEMP]);
      s_ext[RHOH] = 0.0;   // No H for manifold
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC+n] = prob_parm.Y_coflow[n] * s_ext[DENSITY];
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)

     // TODO: zero species and thermal diffusivities at wall
}
#endif

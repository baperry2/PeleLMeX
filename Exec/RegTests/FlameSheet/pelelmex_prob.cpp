#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  pp.query("pertmag", PeleLM::prob_parm->pertmag);
  pp.query("pertlength", PeleLM::prob_parm->pertlength);
  pp.query("ref_frame", PeleLM::prob_parm->ref_frame);
  if (PeleLM::prob_parm->ref_frame == PmfReferenceFrame::fixed_velocity) {
    pp.get("moving_frame_velocity", PeleLM::prob_parm->moving_frame_velocity);
  }

  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();
  PeleLM::pmf_data.initialize();

  PeleLM::prob_parm->tabfunc_parms = new pele::physics::TabulatedFunctionParams;
  std::string bc_data_prefix;
  pp.get("bc_data_prefix", bc_data_prefix);
  PeleLM::prob_parm->tabfunc_parms->host_only_parm().parm_parse_prefix = bc_data_prefix;
  PeleLM::prob_parm->tabfunc_parms->initialize();
  PeleLM::prob_parm->tabfunc_data_d =  static_cast<const pele::physics::TabulatedFunctionData*>(PeleLM::prob_parm->tabfunc_parms->device_parm());
  const bool required = true;
  const auto host_parm = &PeleLM::prob_parm->tabfunc_parms->host_parm();
  PeleLM::prob_parm->idx_out_zmix = get_var_index("mean-f", host_parm, required);
  PeleLM::prob_parm->idx_out_xvel = get_var_index("mean-x-velocity", host_parm, required);
  PeleLM::prob_parm->idx_out_yvel = get_var_index("mean-y-velocity", host_parm, required);
  PeleLM::prob_parm->idx_out_zvel = get_var_index("mean-z-velocity", host_parm, required); 
}

void
PeleLM::freeProbParm()
{
  PeleLM::pmf_data.deallocate();
  PeleLM::prob_parm->tabfunc_parms->deallocate();
  delete PeleLM::prob_parm->tabfunc_parms;
}

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

  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();
  PeleLM::pmf_data.initialize();
  PeleLM::prob_parm->tabfunc_par = new pele::physics::TabulatedFunctionParams();
  PeleLM::prob_parm->tabfunc_par->host_only_parm().parm_parse_prefix = "test";
  PeleLM::prob_parm->tabfunc_par->initialize();
  PeleLM::prob_parm->tabfunc_dat = static_cast<const pele::physics::TabulatedFunctionData*>(PeleLM::prob_parm->tabfunc_par->device_parm());
  const std::string varname = "T";
  PeleLM::prob_parm->idxT = pele::physics::get_var_index(varname.c_str(), &PeleLM::prob_parm->tabfunc_par->host_parm());
}

void
PeleLM::freeProbParm()
{
  PeleLM::pmf_data.deallocate();
  PeleLM::prob_parm->tabfunc_par->deallocate();
  delete PeleLM::prob_parm->tabfunc_par;
}

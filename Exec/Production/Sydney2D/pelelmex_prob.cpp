#include <PeleLMeX.H>
#include <pelelmex_prob.H>

void PeleLM::readProbParm()
{
  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();

  amrex::ParmParse pp("prob");

   pp.query("v_jet",  prob_parm->v_jet);
   pp.query("v_pilot",  prob_parm->v_pilot);
   pp.query("v_coflow",  prob_parm->v_coflow);
   pp.query("Do_jet",   prob_parm->Do_jet);
   pp.query("Di_pilot",   prob_parm->Di_pilot);
   pp.query("Do_pilot",   prob_parm->Do_pilot);
   pp.query("Di_coflow",   prob_parm->Di_coflow);
   pp.query("P_mean",   prob_parm->P_mean);
   pp.query("zstar_pos", prob_parm->zstar_pos);

   std::vector<amrex::Real> Y_jet_temp(NUM_SPECIES), Y_pilot_temp(NUM_SPECIES), Y_coflow_temp(NUM_SPECIES);
   pp.getarr("Y_jet",  Y_jet_temp);
   pp.getarr("Y_pilot",  Y_pilot_temp);
   pp.getarr("Y_coflow",  Y_coflow_temp); 
   for (int n = 0; n < NUM_SPECIES; n++) {
     prob_parm->Y_jet[n] = Y_jet_temp[n];
     prob_parm->Y_pilot[n] = Y_pilot_temp[n];
     prob_parm->Y_coflow[n] = Y_coflow_temp[n];
   }

   AMREX_ALWAYS_ASSERT( prob_parm->Di_pilot >=  prob_parm->Do_jet);
   AMREX_ALWAYS_ASSERT( prob_parm->Do_pilot >=  prob_parm->Di_pilot);
   AMREX_ALWAYS_ASSERT( prob_parm->Di_coflow >=  prob_parm->Do_pilot);

   auto problo = geom[0].ProbLo();
   auto probhi = geom[0].ProbHi();
   prob_parm->center_xy[0] = 0.5 * (probhi[0] + problo[0]);
   prob_parm->center_xy[1] = 0.5 * (probhi[1] + problo[1]);

   // stuff for turbulent inflow for jet
   pp.query("turb_inflow_type",prob_parm->turb_inflow_type);
   if (prob_parm->turb_inflow_type == 1) {

     // Open file
     std::string iname;
     pp.get("turb_inflow_file", iname);
     std::ifstream infile(iname, std::ios::in | std::ios::binary);
     if (not infile.is_open()) {
       amrex::Abort("Unable to open input file " + iname);
     }

     // Read
     read_binary_int(infile, &prob_parm->inflowNtime);
     prob_parm->inflowNtime = 2000;
     read_binary_int(infile, &prob_parm->nr);
     read_binary_int(infile, &prob_parm->nt);
     read_binary_int(infile, &prob_parm->nvar);
     amrex::Print() << "inflow file: nr=" << prob_parm->nr
                    << " ntheta=" <<  prob_parm->nt
                    << " nvar=" << prob_parm->nvar << std::endl;
     read_binary_double(infile, &prob_parm->inflowFreq);
     read_binary_double(infile, &prob_parm->timeFromInflow);
     const int ny = prob_parm->nr;
     const int nz = prob_parm->nt;
     const int ncell = prob_parm->nr*prob_parm->nt;
     const int ntot = prob_parm->inflowNtime*ncell;

     // Read and discard variable names
     amrex::Print() << "Variables: " << std::endl;
     for (int i = 0; i < prob_parm->nvar; ++i) {
       read_binary_strshort(infile);
     }
     int dummy; read_binary_int(infile, &dummy);

     // Create temporary CPU vectors to store data
     std::vector<amrex::Real> timeInput (prob_parm->inflowNtime ,0.0);
     std::vector<amrex::Real> rM (prob_parm->nr, 0.0);
     std::vector<amrex::Real> thetaM (prob_parm->nt, 0.0);
     std::vector<amrex::Real> Uz (ntot, 0.0);
     std::vector<amrex::Real> Ur (ntot, 0.0);
     std::vector<amrex::Real> Ut (ntot, 0.0);
     std::vector<amrex::Real> Ur_temp (ncell, 0.0);
     std::vector<amrex::Real> Ut_temp (ncell, 0.0);
     std::vector<amrex::Real> Zmix (ntot, 0.0);
     std::vector<amrex::Real> Zmix2 (ntot, 0.0);
     
     // Read data into arrays

     // time 
     for (int i = 0; i < prob_parm->inflowNtime ; ++i) {
       timeInput[i] = i*prob_parm->inflowFreq;
     }
     prob_parm->timeInflowMax = timeInput[prob_parm->inflowNtime-1];

     // radial coordinate: note values correspond to midpoints between file coordinates
     amrex::Real prev_val;
     read_binary_double(infile, &prev_val);
     for (int i = 0; i < prob_parm->nr ; ++i) {
       rM[i] = prev_val;
       read_binary_double(infile, &prev_val);
       rM[i] += prev_val;
       rM[i] *= 0.5;
     }

     // theta coordinate: note values correspond to midpoints between file coordinates
     read_binary_double(infile, &prev_val);
     for (int i = 0; i < prob_parm->nt ; ++i) {
       thetaM[i] = prev_val;
       read_binary_double(infile, &prev_val);
       thetaM[i] += prev_val;
       thetaM[i] *= 0.5;
     }
     prob_parm->thetaMax = prev_val;

     // Read the data
     for (int i = 0; i < prob_parm->inflowNtime; i++) {
       read_binary_double(infile, &Uz[i*ncell], ncell);
       read_binary_double(infile, Ur_temp.data(), ncell);
       read_binary_double(infile, Ut_temp.data(), ncell);
       read_binary_double(infile, &Zmix[i*ncell], ncell);
       read_binary_double(infile, &Zmix2[i*ncell], ncell);
       
       // Interpolate Ur and Utheta, which are stored at faces rather than centers
       for (int k = 0; k < nz; k++) {
         for (int j = 0; j < ny-1; j++) {
           Ur[i*ncell+j+k*ny] = 0.5*(Ur_temp[j+k*ny]+Ur_temp[j+k*ny+1]);
         }
         Ur[i*ncell+ny-1+k*ny] = 0.5*Ur_temp[ny-1+k*ny];
       }
       for (int j = 0; j < ny; j++) {
         for (int k = 0; k < nz-1; k++) {
           Ut[i*ncell+j+k*ny] = 0.5*(Ut_temp[j+k*ny]+Ut_temp[j+(k+1)*ny]);
         }
         Ut[i*ncell+j+(nz-1)*ny] = 0.5*(Ut_temp[j+(nz-1)*ny]+Ut_temp[j]);
       }
     }

     // Allocate device data storage and move to device
     prob_parm->d_timeInput = (amrex::Real*) amrex::The_Arena()->alloc(prob_parm->inflowNtime*sizeof(amrex::Real));
     prob_parm->d_rM = (amrex::Real*) amrex::The_Arena()->alloc(prob_parm->nr*sizeof(amrex::Real));
     prob_parm->d_thetaM = (amrex::Real*) amrex::The_Arena()->alloc(prob_parm->nt*sizeof(amrex::Real));
     prob_parm->d_Uz = (amrex::Real*) amrex::The_Arena()->alloc(ntot*sizeof(amrex::Real));
     prob_parm->d_Ur = (amrex::Real*) amrex::The_Arena()->alloc(ntot*sizeof(amrex::Real));
     prob_parm->d_Ut = (amrex::Real*) amrex::The_Arena()->alloc(ntot*sizeof(amrex::Real));
     prob_parm->d_Zmix = (amrex::Real*) amrex::The_Arena()->alloc(ntot*sizeof(amrex::Real));
     // Zmix2 not saved on device because we don't use it (yet)

     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      timeInput.begin(),
                      timeInput.end(),
                      prob_parm->d_timeInput);
     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      rM.begin(),
                      rM.end(),
                      prob_parm->d_rM);
     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      thetaM.begin(),
                      thetaM.end(),
                      prob_parm->d_thetaM);
     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      Uz.begin(),
                      Uz.end(),
                      prob_parm->d_Uz);
     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      Ur.begin(),
                      Ur.end(),
                      prob_parm->d_Ur);
     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      Ut.begin(),
                      Ut.end(),
                      prob_parm->d_Ut);
     amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                      Zmix.begin(),
                      Zmix.end(),
                      prob_parm->d_Zmix);

   } // end of turbinflow stuff

   
}

## sCO2 Combustor

Nonreacting version of sCO2 combustor case. Very stripped down to give linear
solver issues without other complicating factors. Basically just flow through
a tube with a converging nozzle at the end. The edges of the tube have small
bumps, which is what causes linear solver issues for the MAC projection.
The number of the bumps is controlled with `prob.do_dilution_holes`:

   prob.do_dilution_holes = 0 -> no bumps, converges easily
   prob.do_dilution_holes = 1 -> one bump, still converges easily
   prob.do_dilution_holes = 2 -> one row of 6 bumps, converges but slowly
   prob.do_dilution_holes = 3 -> three rows of bumps (30 total), bad convergence


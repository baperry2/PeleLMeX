## FlameSheet
A 2D (or 3D) harmonically perturbed flame sheet, initial solution from a Cantera simulation provided for 3 mechanisms
(drm19, dodecane\_lu and dodecane\_lu\_qss). This is the basis for [weak scaling studies](https://amrex-combustion.github.io/PeleLMeX/manual/html/Performances.html) in PeleLMeX and tests all the
reactive pieces of the algorithm as well as transport options (Unity Lewis number, Soret effect, ...).
More details on the case setup and step-by-step instructions can be found in this [tutorial](https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials_FlameSheet.html).

Comparisons of PeleLMeX results for a methane/air flame against Cantera at several resolutions are reported [here](https://amrex-combustion.github.io/PeleLMeX/manual/html/Validation.html#laminar-premixed-flame).

ANISOTROPIC GRIDS

Compare with the following inputs:

  flamesheet-drm19-2d.inp (isotropic) - works well
  flamesheet-drm19-2d-2x.inp (dx = 2*dy) - MAC projection converges, but requires more iterations
  flamesheet-drm19-2d-3x.inp (dx = 3*dt) - MAC projection fails

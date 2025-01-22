# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = TRUE
USE_CUDA        = FALSE
USE_HIP         = FALSE
USE_SYCL        = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PelePhysics
Eos_Model       = Manifold
Chemistry_Model = Null
Transport_Model = Manifold
Manifold_Dim = 2
Manifold_Type = Table

# PeleLMeX
USE_EFIELD      = FALSE

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX

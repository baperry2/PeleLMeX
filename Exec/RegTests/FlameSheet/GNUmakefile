TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro


# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE
FUEGO_GAS       = TRUE
USE_LIBTORCH    = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
ifeq ($(FUEGO_GAS), TRUE)
	Chemistry_Model = drm19
	Eos_Model		= Fuego
	Transport_Model = Simple
else
	Eos_Model = Manifold
	Chemistry_Model = Manifold
	Transport_Model = Manifold
endif

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX

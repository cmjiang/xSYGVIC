#
#  Compiler and flags
#

# ifort compiler
# --------------
# FC = ifort

# gfortran compiler
# -----------------
FC = gfortran

# flags
# -----
FC_FLAGS = -O3

# Linking to MKL
# --------------
# LAPACK_DIR = /opt/intel/mkl/lib/intel64
# LAPACK_LIB = -L$(LAPACK_DIR) -lmkl_intel_lp64 -lmkl_intel_thread \
#              -lmkl_core -liomp5 -lpthread

# Linking to BLAS and LAPACK shipped with this package (not optimized)
# --------------------------------------------------------------------
LAPACK_LIB = $(MAIN_DIR)/liblapack.a

# GSEP library
# ------------
GSEP_LIB = libgsep.a

# Archiver and flags
# ------------------
ARCH       = ar
ARFLAG     = cr
RANLIB     = ranlib

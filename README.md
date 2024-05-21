# Proto_Bisicles
* AMR Ice sheet code ported to the device using the Proto infrastructure
* BISICLES is an adaptive mesh ice sheet model.
* This particular variant is being used to develop a proto-based ice sheet solver
and, as such, will be a frozen version of the ice sheet solver except for a few,
well-identified locations.

# Changes from the real BISICLES:
* BISICLES has a very nice system for version control that uses svn and python that  has been taken out of this repository since I did not want to reprogram it for git.
* The fortran-based FAS solver in BISICLES works with AMRFAS which is an early implementation
  of AMRFASMultiGrid (which lives in Chombo). FASIceSolver is therefore being rewritten to use the more recent interface.

# Changed source files:
* AmrIce.{H,cpp}

# New Source files:
* Proto_FAS_IceSolver.H
* Proto_FAS_IceSolver.cpp

# New Directories:
* code/exec2D/fas_fortran
* code/exec2D/_fas_proto
* code/test_harness


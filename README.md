# Proto_Bisicles
* AMR Ice sheet code ported to the device using the Proto infrastructure
* BISICLES is an adaptive mesh ice sheet model.
* This particular variant is being used to develop a proto-based ice sheet solver
and, as such, will be a frozen version of the ice sheet solver except for a few,
well-identified locations.
* I gutted out of the make system its svn and python dependencies because they were getting in
the way.    All changes to the make system should not be propogated as the real BISICLES has
a system that holds together for its users. 

# Changed source files:
* AmrIce.cpp  


# New Source files:
* Proto_FAS_IceSolver.H
* Proto_FAS_IceSolver.cpp

# New Directories:
* exec_fortran
* exec_proto
* src_proto
* test_harness


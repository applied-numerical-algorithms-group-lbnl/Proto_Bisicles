# Proto_Bisicles
* AMR Ice sheet code ported to the device using the Proto infrastructure
* BISICLES is an adaptive mesh ice sheet model.
* This particular variant is being used to develop a proto-based ice sheet solver
and, as such, will be a frozen version of the ice sheet solver except for a few,
well-identified locations.


# Branches
* master holds standard BISICLES
* dev     is the devleopment branch. I will try to leave dev a in stable state.
* dtg_dev is for my compulsive code saves.

# Notes:
* BISICLES has a very nice system for version control that uses svn and python.
* This system has been taken out of this repository since I did not want to reprogram it for git.
* Also, the BISICLES make system is a bit odd.   I just use the Chombo one so there  is no mk directory here.
* The old fortran-based FAS solver in BISICLES works with AMRFAS, which lives in a separate svn repository.
* It did not converge for me so...
* I wrote ChF_FAS_Ice_Solver based on AMRFASMultiGrid (which lives in Chombo).
* The next task is to write Proto_FAS_Ice_Solver (also based on AMRFASMultiGrid) but using  Proto.
* The Proto infrastructure is for for performance portability.   It will allow the velocity solve to run on the device.


# Changed source files:
* AmrIce.{H,cpp}

# New Source files:
*   CHF_FAS_Ice_Solver.H
* Proto_FAS_Ice_Solver.H (to come)

# New Directories:
* code/exec2D/_fas_fortran holds fortran-based FAS solver input templates
* code/exec2D/_fas_proto holds proto-based FAS solver input templates
* code/test_harness will be a Chombo-style test harness for simulation campaigns.


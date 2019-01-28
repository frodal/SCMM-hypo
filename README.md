# SIMLab Crystal Mechanics model - Hypoelastic formulation (SCMM-Hypo)
User defined material model for ABAQUS Standard/Explicit

## Prerequisites
Before compiling the user material subroutine for the FEM code, install and check:

- A Fortran compiler compatible with the FE code (e.g. Intel Visual Fortran)
- Check the license of the FE solver, and the version of the FE solver

## Compilation

1. First case, run a simulation with a local compiled library (the compilation will be done each time a simulation is run)
  - Clone the repository: `SCMM-hypo` from the SIMLab project on www.code.sintef.no
  - Change the current directory to the `SCMM-hypo` folder.
  - Copy the simulation input (e.g. `MySim.inp`) to the current directory.
  - Run the simulation using the command: `abaqus double job=MySim user=HypoImp int`
  - Use the option `user=HypoImp` for the implicit solver (ABAQUS/Standard).
  - Use the option `user=HypoExp` for the explicit solver (ABAQUS/Explicit).

2. Second case, compile a library and share it:
  - If needed, make a directory where to place the ABAQUS library (e.g. library-dir = `/home/username/bin/abaqus` for Linux or `C:\Users\username\abaqus` for Windows)
  - Clone the repository: `SCMM-hypo` from the SIMLab project on www.code.sintef.no
  - Change the current directory to the `SCMM-hypo` folder.
  - Run the command: `abaqus make library=HypoImp directory=library-dir`
  - Use the option `library=HypoImp` for the implicit solver (ABAQUS/Standard).
  - Use the option `library=HypoExp` for the explicit solver (ABAQUS/Explicit).
  - Change the environment file to specify the path of the new library:
      - If needed, create an environment file, `abaqus_v6.env`, in your home directory and/or the current directory. Settings in the home directory file will be applied to all jobs that you run. Settings in the current directory file will be applied only to jobs run from the current directory.
	  - Add the entry: `usub_lib_dir='library-dir'`
	  
	  
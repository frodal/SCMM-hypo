# SIMLab Crystal Mechanics model - Hypoelastic formulation (SCMM-Hypo)

User defined material model for Abaqus/Standard and Abaqus/Explicit

## Prerequisites

Before compiling the user material subroutine for the FEM code, install and check:

- A Fortran compiler compatible with the FE code (e.g. Intel Visual Fortran)
- Check the license of the FE solver, and the version of the FE solver

## Compilation

Follow point 1 or 2 below to compile

1. First case, run a simulation with a local compiled library (the compilation will be done each time a simulation is run)
    - Clone the repository: `SCMM-hypo` from the SIMLab project on www.code.sintef.no or gitte.kt.ntnu.no
    - Change the current directory to the `SCMM-hypo` folder
    - Copy the simulation input (e.g. `MySim.inp`) to the current directory
    - Run the simulation using the command: `abaqus double job=MySim user=HypoImp int`
    - Use the option `user=HypoImp` for the implicit solver (Abaqus/Standard)
    - Use the option `user=HypoExp` for the explicit solver (Abaqus/Explicit)
2. Second case, compile a library:
    - If needed, make a directory where to place the Abaqus library (e.g. library-dir = `/home/username/bin/abaqus` for Linux or `C:\Users\username\abaqus` for Windows)
    - Clone the repository: `SCMM-hypo` from the SIMLab project on www.code.sintef.no or gitte.kt.ntnu.no
    - Change the current directory to the `SCMM-hypo` folder
    - Run the command: `abaqus make library=HypoImp directory=library-dir`
    - Use the option `library=HypoImp` for the implicit solver (Abaqus/Standard)
    - Use the option `library=HypoExp` for the explicit solver (Abaqus/Explicit)
    - Change the environment file to specify the path of the new library:
      - If needed, create an environment file, `abaqus_v6.env`, in your home directory and/or the current directory. Settings in the home directory file will be applied to all jobs that you run. Settings in the current directory file will be applied only to jobs run from the current directory
      - Add the entry: `usub_lib_dir='library-dir'`

## Tests

To run the tests:

1. Compile the subroutine using method nr. 2 above, both for Abaqus/Explicit and Abaqus/Standard
2. Edit the `abaqus_v6.env` file in the `./Tests/Abaqus/` folder so that it points to the folder where the compiled libraries are located, i.e., edit the entry `usub_lib_dir='library-dir'`
3. Make sure that you have Abaqus and Python 3 installed with the necessary Python libraries (See `./Tests/Test.py`)
4. Change the current directory to the `SCMM-hypo` folder
5. Run the command: `python3 ./Tests/Test.py run --location=0` to run the tests (To run the tests on Snurre use `--location=1` instead, note that the Python script must be run from Snurre in this case)
6. When the Abaqus jobs have finished, run the command: `python3 ./Tests/Test.py post` to post-process the tests
7. To clean the test directory, run the command: `python3 ./Tests/Test.py clean`

## Contributing

To contribute:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Contact

Bjørn Håkon Frodal - [@frodal](https://github.com/frodal) - bjorn.h.frodal@ntnu.no

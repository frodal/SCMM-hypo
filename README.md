# SIMLab Crystal Mechanics model - Hypoelastic formulation (SCMM-Hypo)

User defined material model for Abaqus/Standard and Abaqus/Explicit

## Cite

Please cite the articles listed below if you use this model

```bibtex
@article{Frodal.et.al.2019,
title = {Modelling and simulation of ductile failure in textured aluminium alloys subjected to compression-tension loading},
journal = {International Journal of Plasticity},
volume = {118},
pages = {36-69},
year = {2019},
issn = {0749-6419},
doi = {https://doi.org/10.1016/j.ijplas.2019.01.008},
url = {https://www.sciencedirect.com/science/article/pii/S0749641918305904},
author = {Bjørn Håkon Frodal and Lars Edvard Blystad Dæhli and Tore Børvik and Odd Sture Hopperstad}
}

@article{Frodal.et.al.2021,
title = {On the coupling of damage and single crystal plasticity for ductile polycrystalline materials},
journal = {International Journal of Plasticity},
pages = {102996},
year = {2021},
issn = {0749-6419},
doi = {https://doi.org/10.1016/j.ijplas.2021.102996},
url = {https://www.sciencedirect.com/science/article/pii/S0749641921000711},
author = {Bjørn Håkon Frodal and Susanne Thomesen and Tore Børvik and Odd Sture Hopperstad}
}
```

## Subroutine input

| Property number | Material input                                                                     |
|:---------------:|:---------------------------------------------------------------------------------- |
|               1 | Elastic constant *c*<sub>11</sub>                                                  |
|               2 | Elastic constant *c*<sub>12</sub>                                                  |
|               3 | Elastic constant *c*<sub>44</sub>                                                  |
|               4 | Reference slip rate, *γ*<sub>0</sub>                                               |
|               5 | Instantaneous strain rate sensitivity, *m*                                         |
|               6 | Initial critical resolved shear stress, *τ*<sub>*c*0</sub>                         |
|               7 | Latent hardening coefficient, *q*                                                  |
|               8 | Texture flag (1=Euler angles from material card, 2=Euler angles from history card) |
|               9 | Initial Euler angle, *ϕ*<sub>1</sub> in degree                                     |
|              10 | Initial Euler angle, Φ in degree                                                   |
|              11 | Initial Euler angle, *ϕ*<sub>2</sub> in degree                                     |
|              12 | Hardening flag (1 for Voce or 2 for Kalidindi et al. (1992))                       |
|              13 | Hardening parameter, *h*<sub>0</sub> or *θ*<sub>1</sub> depending on hflag         |
|              14 | Hardening parameter, *τ*<sub>s</sub> or *τ*<sub>1</sub> depending on hflag         |
|              15 | Hardening parameter, *a* or *θ*<sub>2</sub> depending on hflag                     |
|              16 | Hardening parameter, 0.0 or *τ*<sub>2</sub>                                        |
|              17 | Tangent operator flag (1=Elastic tangent operator, 2=Consistent tangent operator)  |
|              18 | Initial damage, *f*<sub>0</sub>                                                    |
|              19 | Critical damage, *f*<sub>*c*</sub>                                                 |
|              20 | Damage evolution parameter, *q*<sub>1</sub>                                        |
|              21 | Damage evolution parameter, *q*<sub>2</sub>                                        |

***Warning: Do not use a local coordinate system (CSYS) or a material orientation with this subroutine in Abaqus Standard. This will break the co-rotational formulation.***

Note that the Tangent operator flag, property number 17, is only used for Abaqus Standard,
and the Quasi-Newton solution technique should be used when the Elastic tangent operator is
selected. The Quasi-Newton solution technique is located under the step settings in the Abaqus
CAE.

Note also that the last Hardening parameter, property number 16, should be put to 0.0 when
the Kalidindi et al. (1992) hardening model is used with Abaqus Standard.

## Subroutine output

| Variable number | Solution dependent variable                                         |
|:---------------:|:------------------------------------------------------------------- |
|             1-3 | Euler angles, *ϕ*<sub>1</sub>, Φ, *ϕ*<sub>2</sub>                   |
|            4-12 | Components of the rotation tensor ***R***                           |
|           13-24 | Critical resolved shear stresses *τ*<sub>*c*</sub><sup>(α)</sup>    |
|              25 | Accumulated plastic shear strain *Γ*                                |
|              26 | Equivalent von Mises stress *σ*<sub>eq</sub>                        |
|              27 | Equivalent von Mises plastic strain *ε*<sup>p</sup><sub>eq</sub>    |
|              28 | Number of sub-steps in the current time step n<sub>sub</sub>        |
|              29 | Damage variable, *f*                                                |
|              30 | Status variable used for element deletion in Abaqus/Explicit        |

## Prerequisites

Before compiling the user material subroutine for the FEM code, install and check:

- A Fortran compiler compatible with the FE code (e.g. Intel Visual Fortran)
- Check the license of the FE solver, and the version of the FE solver

## Compilation

Follow point 1 or 2 below to compile

1. First case, run a simulation with a local compiled library (the compilation will be done each time a simulation is run)
    - Clone the repository: `SCMM-hypo` from the SIMLab project on www.code.sintef.no or from github.com
    - Change the current directory to the `SCMM-hypo` folder
    - Copy the simulation input (e.g. `MySim.inp`) to the current directory
    - Run the simulation using the command: `abaqus double job=MySim user=HypoImp int`
    - Use the option `user=HypoImp` for the implicit solver (Abaqus/Standard)
    - Use the option `user=HypoExp` for the explicit solver (Abaqus/Explicit)
2. Second case, compile a library:
    - If needed, make a directory where to place the Abaqus library (e.g. library-dir = `/home/username/bin/abaqus` for Linux or `C:\\Users\\username\\abaqus` for Windows)
    - Clone the repository: `SCMM-hypo` from the SIMLab project on www.code.sintef.no or from github.com
    - Change the current directory to the `SCMM-hypo` folder
    - Run the command: `abaqus make library=HypoImp directory=library-dir`
    - Use the option `library=HypoImp` for the implicit solver (Abaqus/Standard)
    - Use the option `library=HypoExp` for the explicit solver (Abaqus/Explicit)
    - Change the environment file to specify the path of the new library:
      - If needed, create an environment file, `abaqus_v6.env`, in your home directory and/or the current directory. Settings in the home directory file will be applied to all jobs that you run. Settings in the current directory file will be applied only to jobs run from the current directory
      - Add the entry: `usub_lib_dir='library-dir'`

### Compiler directives

Note that the subroutines use preprocessor directives to determine if it is compiled for Abaqus/Standard, Abaqus/Explicit or neither. Specific choices can also be made at compile time instead of at runtime, e.g., the choice of which hardening model to use. For further information on the availible preprocessor definitions see the `Definitions.f` file. These preprocessor directives can be used to minimize the computational time by including only necessary model features, and removing unecessary conditional statements.

If the subroutines are to be included in another file, use the preprocessor include directive (e.g., `#include 'HypoImp.f'`) instead of the Fortran include statement (e.g., `include 'HypoImp.f'`).

## Tests

To run the tests:

1. Compile the subroutine using method nr. 2 above, both for Abaqus/Explicit and Abaqus/Standard
2. Edit the `abaqus_v6.env` file in the `./Tests/Abaqus/` folder so that it points to the folder where the compiled libraries are located, i.e., edit the entry `usub_lib_dir='library-dir'`
3. Make sure that you have Abaqus and Python 3 installed with the necessary Python libraries (See `./Tests/Test.py`)
4. Change the current directory to the `SCMM-hypo` folder
5. Run the command: `python ./Tests/Test.py run` to run the tests (The flag `--interactive_off` can be appended this command to run all the Abaqus tests simultaneously)
6. When the Abaqus jobs have finished, run the command: `python ./Tests/Test.py post` to post-process the tests (The flag `--plot` can be appended this command to plot the results)
7. To clean the test directory, run the command: `python ./Tests/Test.py clean`

## Contributing

To contribute:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Contact

Bjørn Håkon Frodal - [@frodal](https://github.com/frodal) - bjorn.h.frodal@ntnu.no

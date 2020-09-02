<p align="center">
<img width="250" height="154" src="Common/doc/logoSU2small.png">
</p>

# feature_DDRANS branch of SU2
This feature branch includes capability for Data-Driven RANS simulations, whereby the turbulent Reynolds stress field is fixed to a given field. The given field would typically arise from data-driven (supervised machine learning) predictions e.g. see [here](https://github.com/ascillitoe/mondrian_turbulence). The field could also be a Reynolds stress field arising directly from an LES or DNS simulation. 

### How to use
1. Run a standard RANS simulation until convergence. 
2. Restart from this baseline RANS solution i.e. `RESTART_SOL=YES` and `READ_BINARY_RESTART=NO`, with `USING_DD=YES`. 
3. The solver will run with the Reynolds stress field set perturbed from the original RANS field by the amount set in the .csv restart file (see notes below). 

### Notes
- Currently only works with the incompressible solver (`SOLVER=INC_RANS`) and SST model (`KIND_TURB_MODEL= SST`).
- The Reynolds stress field is set by perturbing the original RANS stresses by the perturbations defined in the restart file. Six additional fields must be added to the .csv restart file after the k and omega columns (binary not currently implemented):
   1. delta_x: x perturbation to eigenvalue (in Cartesian coords).
   2. delta_y: y perturbation to eigenvalue (in Cartesian coords).
   3. h2: 2nd element of unit quarternion defining eigenvector rotation.
   4. h3: 3rd element of unit quarternion defining eigenvector rotation.
   5. h4: 4th element of unit quarternion defining eigenvector rotation.
   6. delta_logk: Log10 of perturbation to turbulent kinetic energy.
   
See the file `convert.py` for an example helper script to read the perturbations from a .vtk file and write to an existing .csv restart file. h1 is calculated within the solver from h2,h3 and h4. Setting h2=h3=h4=0 leads to no eigenvector perturbations. (see [here](https://arxiv.org/abs/1709.05683) for more info on this perturbation strategy).

### Settings
`USING_DD=YES/NO`: Turn data-driven RANS on/off. 

`WRITE_DD=YES/NO`: Write Aij_ML (the anisotropy tensor field after applying perturbations) to solution files for debugging.

`DD_GAMMA_MAX`: Maximum blending factor to blend original RANS and new Aij fields together i.e. `Aij_new = (1-gamma_max)*Aij_orig + gamma_max*Aij_DD`

`DD_N_MAX`: The blending factors gamma and alpha are ramped up until iteration `n=n_max`, where `gamma=gamma_max` and `alpha=alpha_max`.

`DD_ALPHA_MAX`: Blend the turb. kinetic energy (TKE) field between the exact one i.e. with the k production set as the exact production term corresponding to the perturbed Reynolds stress field i.e. Pk=tau:grad(U), and the original RANS TKE field perturbed by delta_logk i.e. `TKE=(1-alpha)*TKE_exact + alpha*TKE_perturbed`. `DD_ALPHA_MAX=0` gives a strategy like the "continuation" solver described by [Kaandorp and Dwight](https://arxiv.org/abs/1810.08794), `DD_ALPHA_MAX=1` gives a strategy similar to [Wang et al.](https://arxiv.org/abs/1606.07987). Ideally, if the DNS anisotropy tensor is fed in (by perturbing the eigen-values and vectors), using the exact TKE production (`DD_ALPHA_MAX=0`) should result in convergence to something near to the DNS flow field (assuming no significant non-local turbulence effects). However, in practice, this appears to be sensitive to the initial condition and mesh quality (perhaps because both of these can lead to under prediction of grad(U)). Therefore, blending in some of the perturbed TKE field (at least at first) can be useful. 

`WRITE_DIST=YES/NO`: Write the wall distance field to the solution file (can be useful as an input feature for machine learning approaches).


# SU2 (ver. 7.0.6 "Blackbird"): The Open-Source CFD Code

Computational analysis tools have revolutionized the way we design engineering systems, but most established codes are proprietary, unavailable, or prohibitively expensive for many users. The SU2 team is changing this, making multiphysics analysis and design optimization freely available as open-source software and involving everyone in its creation and development. 

For an overview of the technical details in SU2, please see the following AIAA Journal article:

"SU2: An open-source suite for multiphysics simulation and design," AIAA Journal, 54(3):828-846, 2016. http://arc.aiaa.org/doi/10.2514/1.J053813

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

Continuous Integration:<br/>
[![Regression Testing](https://github.com/su2code/SU2/workflows/Regression%20Testing/badge.svg?branch=develop)](https://github.com/su2code/SU2/actions)
[![Release](https://github.com/su2code/SU2/workflows/Release%20Management/badge.svg?branch=develop)](https://github.com/su2code/SU2/actions)

Code Quality:<br/>
[![CodeFactor](https://www.codefactor.io/repository/github/su2code/su2/badge)](https://www.codefactor.io/repository/github/su2code/su2)

# SU2 Introduction

SU2 is a suite of open-source software tools written in C++ for the numerical solution of partial differential equations (PDE) and performing PDE constrained optimization. 

The primary applications are computational fluid dynamics and aerodynamic shape optimization, but has been extended to treat more general equations such as electrodynamics and chemically reacting flows. 

You will find more information and the latest news in:
   - SU2 Home Page: https://su2code.github.io
   - GitHub repository: https://github.com/su2code
   - CFD Online: http://www.cfd-online.com/Forums/su2/
   - Twitter: https://twitter.com/su2code
   - Facebook: https://www.facebook.com/su2code


# SU2 Installation

## Precompiled binaries for Linux, MacOS, Windows

You can find precompiled binaries of the latest version on our [download page](https://su2code.github.io/download/) or under [releases](https://github.com/su2code/SU2/releases).

## Build SU2
The build system of SU2 is based on a combination of [meson](http://mesonbuild.com/) (as the front-end) and [ninja](https://ninja-build.org/) (as the back-end). Meson is an open source build system meant to be both extremely fast, and, even more importantly, as user friendly as possible. Ninja is a small low-level build system with a focus on speed. 

Short summary of the minimal requirements:

- C/C++ compiler
- Python 3

**Note:** all other necessary build tools and dependencies are shipped with the source code or are downloaded automatically.

If you have these tools installed, you can create a configuration using the `meson.py` found in the root source code folder:
```
./meson.py build
```
Use `ninja` to compile and install the code

```
./ninja -C build install
```

For more information on how to install and build SU2 on Linux, MacOS or Windows, have a look at the [documentation](https://su2code.github.io/docs_v7/).

##  SU2 Path setup

When installation is complete, please be sure to add the `$SU2_HOME` and `$SU2_RUN` environment variables, and update your `$PATH` with `$SU2_RUN`. 

For example, add these lines to your `.bashrc` file:
```
export SU2_RUN="your_prefix/bin"
export SU2_HOME="/path/to/SU2vX.X.X/"
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$SU2_RUN:$PYTHONPATH
```

`$SU2_RUN` should point to the folder where all binaries and python scripts were installed. This is the prefix you set with the --prefix option to meson. Note that the bin/ directory is automatically added to your prefix path.

`$SU2_HOME` should point to the root directory of the source code distribution, i.e., `/path/to/SU2vX.X.X/`.

Thanks for building, and happy optimizing!

- The SU2 Development Team


# SU2 Developers


We follow the popular "GitFlow" branching model for scalable development. In the SU2 repository, the master branch represents the latest stable major or minor release (7.0, 6.2.0, etc.), it should only be modified during version releases. Work that is staged for release is put into the develop branch via Pull Requests on GitHub from various "feature" branches where folks do their day-to-day work on the code. At release time, the work that has been merged into the develop branch is pushed to the master branch and tagged as a release.

SU2 is being developed by individuals and organized teams all around the world. 

A list of current contributors can be found in the AUTHORS.md file.

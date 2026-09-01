# Introduction

This repository supplements the paper *Safe-by-Design: Approximate Nonlinear Model Predictive Control with Real-time Feasibility*. A [preprint](https://arxiv.org/abs/2509.22422) is available
It provides more details on the implementation and a short user guide to either reproduce
the figures of the paper or re-run the MATLAB scripts.

# Overview

The numerical result section is split into two parts. The first part is
to compared the proposed approach with other state-of-the-art nonlinear
control methods. The second part demonstrates only the proposed approach
in a more challenging scenario. After providing the basic *installation*
instructions and what is needed, a short user guide for the two
aforementioned parts is provided.

## Needed Software and Setup

- MATLAB 2023b or newer is required (newer versions should also work,
  but are not tested!)

- For the pre-computation step we make use of
  CaΣoS [1] v1.0.0-rc. The used version can be found
  in the repository. After cloning this repository, get the submodule
  CaΣoS with
  ```
  git clone https://github.com/iFR-ACSO/TAC25-Inf-MPC.git TAC-Inf-MPC
  cd TAC-Inf-MPC
  git submodule update --init --recursive
  ```
  or in one go
  ```
  git clone --recurse-submodules https://github.com/iFR-ACSO/TAC25-Inf-MPC.git TAC-Inf-MPC

  ```


- An alternative is to use the provided `.zip` in the [DARUS repository](https://doi.org/10.18419/DARUS-5297). Download the data set and unzip it. A copy of CaΣoS submodule is provided.
  Follow the instructions below. 

- MOSEK [2] v11.0.4 is used as the underlying SDP
  solver for the pre-computation step. Academic licenses can be obtained
  from <https://www.mosek.com/products/academic-licenses/>. Follow the
  installation instructions from MOSEK.

- For online optimizations, we make use of
  CasADi [3] v3.6.7 to setup the problems and can
  be obtained from <https://web.casadi.org/get/>. We used
  qrqp [4] for the (S)QP and
  IPOPT [5] to solve the discrete-time OCP.
  The solvers are included in CasADi so no additional installation
  required.

- For reproduction of the plots in the [papert](https://arxiv.org/abs/2509.22422), large `.mat` files (>50 MB) are needed. They are stored in the [DARUS repository](https://doi.org/10.18419/DARUS-5297).

After installing all the above needed software:

1. Open Matlab and navigate to the main folder. You should see three folder, one with the initialized submodule (or Copy)
2. Add the CaΣoS main folder to the MATLAB path.
3. You can either re-run everything (see below) or if you just want to reproduce the data (i.e., plots and tables), you might need to copy paste some large `.mat` files from the [DARUS repository](https://doi.org/10.18419/DARUS-5297).


## Running Examples and Reproducing Results

In the following, it is assumed everything is installed and setup as
explained before.

## Part I: Comparative Study

### Overview Files

For the CBF-CLF-QP and infinetesimal MPC scheme the user finds a synthesis script
to compute the terminal conditions. These use an initial guess, which is
loaded from a `.mat` file. Once the synthesis is done, the terminal
conditions are also stored in a `.mat` file. This `.mat` files are the
loaded into the workspace in the simulation scripts.

`.mex` functions for the simulation are generated to improve simulation
time. To compile the `.mex` functions a C/C++ compiler is required.

For the full-horizon NMPC formulations (Ipopt, RTI), the terminal set is
also pre-computed, i.e., the maximum stable level set, which can be
found in `full_MPC_IPOPT/termIngredient_full.m`. The level set is manually set in the simulation script. After the simulations
run, plots from all runs are generated and `.mat ` file of the whole workspace is stored. Due to the larger size, the `.mat` files are post-processed. You can find the full workspace `.mat` in the [DARUS repository](https://doi.org/10.18419/DARUS-5297). The post-processed `.mat`
results are stored in the main folder.

### Folder Structure

```text
Comparison_singleAxis/
├── CBF_CLF_QP         # CBF-CLF synthesis and simulation
├── full_MPC_alpaqa    # Full-horizon NMPC simulation using alpaqa solver (or fatrop)
├── full_MPC_IPOPT     # Full-horizon NMPC simulation using IPOPT solver
├── helperFunc         # Folder that contains helper functions (e.g. MRP→Euler)
├── horizon-one        # Simulation of the horizon-one MPC
├── inf_MPC            # Synthesis and simulation of inf.MPC
├── poly_controlLaw    # Simulation of the poly. control law from proposed approach
└── RTI_full_MPC       # Full-horizon NMPC simulation using custom RTI solver
```


### Reproduction

We provide `.mat` files for each individual approach in the main folder, of the
actual simulation results. Run `comparison_MultipleRuns.m`, to reproduce
the table from the paper and to get the plot for the single axis
rotation. Due to the large amount of data, the full-horizon NMPC
formulations (RTI and IPOPT) have post-processing scripts in their
folders (see above). Once the actual simulation ran, the post-processing script
reduces the data to the comparison minimum. The data of the complete
workspace for the full-horizon formulations is provided in the [DARUS repository](https://doi.org/10.18419/DARUS-5297). This pre-computation results can be found in the corresponding folder of each method.

### Re-running 

The user is welcome to run the scripts and functions to reproduce the results. It should be noted
that, for example, the full-horizon NMPC formulation might take a
significant amount of time for the simulation in MATLAB.

## Part II: Performance Test

### Overview Files

The second part only considers the infinetesimal MPC scheme. The user
finds a synthesis script to compute the terminal conditions. These use an
initial guess, which is loaded from a `.mat` file. Once the synthesis is
done, the terminal conditions are also stored in a `.mat` file. 
This`.mat` files are the loaded into the workspace in the simulation
scripts.

`.mex` functions for the simulation are generated to improve simulation
time. To compile the `.mex` functions a C/C++ compiler is required.

### Folder Structure

Not all MATLAB scripts and `.mat` files are listed. Only the most important
functions.
```text
Three-Axis_constrained_quartic_compFun/
├── helperFunction/               # Folder that contains helper functions (e.g. MRP→Euler)
├── innerApprox_constraintSet.m   # Compute inner-approximation of constraint set
├── synthesis_CBF_CLF.m           # Script for the synthesis of the terminal conditions
└── inf_MPC_simulation.m          # Script to run the Monte-Carlo simulations
```


### Reproduction

Run `evaluation.m` to get all plots from the paper and additional once. Due to the large size of the `.mat` file, you have to copy and paste it from the [DARUS repository](https://doi.org/10.18419/DARUS-5297) in the corresponding folder. Or, you re-run the simulation. However, you might receive different results because a new uniform distribution is calculated for the initial conditions.

### Re-running

The inner-approximation of the constraint set (i.e., the constraint set)
for this scenario can be pre-computed using `innerApprox_constraintSet.m`.
The approximated constraint set is manually inserted into the synthesis
script`synthesis_CBF_CLF.m`. Important to know is that the simulations
script `inf_MPC_simulation.m` computes a new uniform distribution if
re-run. Thus, different results to the paper are expected!


### Final remarks
In case of problems, questions or remarks, please contact the corresponding authors (see below). 
Jan Olucak: jan.olucak@ifr.uni-stuttgart.de
Torbjørn Cunis: torbjoern.cunis@ifr.uni-stuttgart.de
Arthur Castello B. de Oliveira: castello.a@northeastern.edu


### Citation
Please cite the paper as 
```
@misc{olucak2025safebydesignapproximatenonlinearmodel,
      title={Safe-by-Design: Approximate Nonlinear Model Predictive Control with Real Time Feasibility}, 
      author={Jan Olucak and Arthur Castello B. de Oliveira and Torbjørn Cunis},
      year={2025},
      eprint={2509.22422},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2509.22422}, 
}
```

The supplemantary material can be cited with
```
@data{DARUS-5297_2025,
author = {Olucak, Jan and Castello Branco de Oliveira, Arthur and Cunis, Torbjørn},
publisher = {DaRUS},
title = {{Supplementary Material for: Safe-by-Design Approximate Nonlinear Model Predictive Control with Realtime Feasibility}},
year = {2025},
version = {DRAFT VERSION},
doi = {10.18419/DARUS-5297},
url = {https://doi.org/10.18419/DARUS-5297}
}

```

### References


[1] T. Cunis and J. Olucak, “CaΣoS: A nonlinear sum-of-squares optimization suite,” in 2025 Amer-
ican Control Conference, (Boulder, CO), pp. 1659–1666, 2025.

[2] E. D. Andersen and K. D. Andersen, “The Mosek Interior Point Optimizer for Linear Program-
ming: An Implementation of the Homogeneous Algorithm,” in High Performance Optimization
(H. Frenk, K. Roos, T. Terlaky, and S. Zhang, eds.), 2000.

[3] J. A. E. Andersson, J. Gillis, G. Horn, J. B. Rawlings, and M. Diehl, “CasADi: a software frame-
work for nonlinear optimization and optimal control,” Mathematical Programming Computation,
vol. 11, no. 1, pp. 1–36, 2019.

[4] J. A. Andersson and J. B. Rawlings, “Sensitivity analysis for nonlinear programming in casadi,”
IFAC-PapersOnLine, vol. 51, no. 20, pp. 331–336, 2018. 6th IFAC Conference on Nonlinear
Model Predictive Control NMPC 2018.

[5] A. Wächter and L. T. Biegler, “On the implementation of an interior-point filter line-search
algorithm for large-scale nonlinear programming,” Mathematical Programming, vol. 106, pp. 25–
57, Mar. 2006.


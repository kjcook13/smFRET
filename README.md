# R Code for Simulating smFRET 
This repository contains scripts to simulate smFRET model in the paper entitled "Time-Heterogeneity of the Förster Radius from Dipole Orientational Dynamics Impacts Single-Molecule FRET Experiments".
The preprint can be found here: https://arxiv.org/abs/2404.09883 

## Repository Contents
README.md: Documentation for understanding and using simulation codes.

## General Code Information
Each code is well commented with definitions of each parameter. Each code outputs a .rds file containing 4 columns of information extracted from the simulations.
That is Bin_FRET = FRET, Bin_life_D = Lifetime, Var_life = Variance of Lifetime, ds = Dynamic Shift, and df = Donor fluorescence rate.
The user can decide which values they are interested in plotting.
To plot the simulations, 

## Isotropic_Final.R
Script to simulate the anisotropic spring model in Section IIA and Figure 4A.

## Anisotropic_Final.R
Script to simulate the anisotropic spring model in Section IIA and Figure 4C.

## Elastic_Pendulum_Final.R
Script to simulate the elastic pendulum model in Section IIB and Figure 4E.

## Elastic_Pendulum_WithKappa_Final.R
Script to simulate the elastic pendulum model with kappa, orientational dynamics, in Section IIC and Figure 4G. 





# WSIDM
 A boltzmann code in C to compute the cosmic evolution of interacting dark matter

 Daniel Egana-Ugrinovic, 2021.

maincode.c is the core of the code.
Physics constants are defined in background.h
Derivatives.h defines ODE derivatives as macros
integrals.h are functions for Euler integration
RKF.h contains the ODE for the cosmic perturbations
structures.h defines the perturbations to be solved in structures
thermodynamics.h defines basic thermodynamics functions, such as sound speed

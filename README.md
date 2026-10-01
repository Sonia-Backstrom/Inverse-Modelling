# Inverse Modelling of an Ice Sheet

## Overview

This project investigates an inverse modelling problem for an ice sheet. 
The aim is to estimate the underlying bedrock from observed surface 
elevation and horizontal surface velocity.

The bedrock is represented by a prescribed two-bump model. The shapes and 
positions of the two bumps are known, while their amplitudes are unknown.

The project consists of four main parts:

1. Calculate ice thickness from the observed surface elevation and velocity.
2. Reconstruct the bedrock.
3. Estimate the two bedrock bump amplitudes using a least-squares fit.
4. Investigate how random measurement noise affects the estimated amplitudes.

A 200-year forward simulation is also performed using the inferred bedrock.





## Repository Structure

```text
.
├── README.md
├── icemodel.py
├── RunIceSheetSimulation.py
├── inverse_model.ipynb
├── coordinates.npy
├── surface.npy
├── ux.npy
└── uz.npy

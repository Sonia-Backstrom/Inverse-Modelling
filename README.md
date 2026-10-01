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

---

## Method

The supplied ice-sheet model uses the shallow ice approximation (SIA).

The observed data consist of:

- Surface elevation, \(h(x)\)
- Horizontal surface velocity, \(u_x(x)\)

### 1. Surface gradient

The surface gradient is calculated from the observed surface elevation 
using finite differences. For interior grid points, the central difference is

\[
\frac{dh}{dx}
\approx
\frac{h_{j+1}-h_{j-1}}{2\Delta x}.
\]

### 2. Ice thickness

The horizontal velocity relation used in the inverse model is

\[
u_x =
-\frac{2A(\rho g)^n}{n+1}
\left(\frac{dh}{dx}\right)^n
H^{n+1}.
\]

Rearranging gives

\[
H =
\left[
\frac{u_x(n+1)}
{-2A(\rho g)^n
\left(\frac{dh}{dx}\right)^n}
\right]^{1/(n+1)}.
\]

This allows the ice thickness \(H\) to be estimated from the observed 
surface elevation and velocity.

### 3. Bedrock reconstruction

Ice thickness is related to surface elevation and bedrock elevation by

\[
H=h-b.
\]

Therefore,

\[
\boxed{b=h-H}.
\]

The reconstructed bedrock is then fitted using the prescribed two-bump model.

### 4. Two-bump bedrock model

The bedrock is represented as

\[
b(x)
=
b_{\mathrm{base}}(x)
+
a_1B_1(x)
+
a_2B_2(x),
\]

where:

- \(b_{\mathrm{base}}(x)\) is the base bedrock slope,
- \(B_1(x)\) is the known first bump shape,
- \(B_2(x)\) is the known second bump shape,
- \(a_1\) and \(a_2\) are the unknown amplitudes.

The amplitudes are estimated using a least-squares fit:

\[
\min_{a_1,a_2}
\sum_j
\left[
b_{\mathrm{rec}}(x_j)
-
\left(
b_{\mathrm{base}}(x_j)
+a_1B_1(x_j)
+a_2B_2(x_j)
\right)
\right]^2.
\]

---

## Results

The inverse calculation gives the following bedrock amplitudes:

\[
\boxed{a_1 \approx 0.8000}
\]

\[
\boxed{a_2 \approx -0.4000}
\]

Therefore, the fitted bedrock is

\[
b_{\mathrm{fit}}(x)
=
b_{\mathrm{base}}(x)
+
0.8B_1(x)
-
0.4B_2(x).
\]

The fitted two-bump profile closely follows the reconstructed bedrock.

The calculated ice thickness ranges approximately from

\[
473.3\ \text{m}
\]

to

\[
954.5\ \text{m}.
\]

---

## 200-Year Forward Simulation

The inferred bedrock amplitudes

\[
a_1=0.8,\qquad a_2=-0.4
\]

are used as input to the supplied forward ice-sheet model.

The model is then simulated for 200 years. The resulting surface elevation 
and bedrock are plotted to demonstrate the ice-sheet evolution produced by 
the inferred bedrock.

The forward simulation is separate from the inverse calculation because 
the inverse problem does not use the unknown observation time.

---

## Noise Experiment

To investigate the effect of measurement uncertainty, random noise is 
added to both:

- surface elevation \(h\)
- horizontal velocity \(u_x\)

The **same direct inverse method** is then applied at each noise level.

The following noise levels are tested:

\[
0\%,\quad 1\%,\quad 2\%,\quad 5\%,\quad 10\%.
\]

The same underlying random noise pattern is used and its magnitude is 
scaled for each noise level.

### Noise results

| Noise level | \(a_1\) | \(a_2\) |
|-------------|--------:|--------:|
| 0%          | 0.8000  | -0.4000 |
| 1%          | 0.1953  | -0.4932 |
| 2%          | 0.6118  | -0.5201 |
| 5%          | 1.5751  | 0.8293 |
| 10%         | 1.6844  | 1.3866 |

The results show that measurement noise can substantially change the 
estimated bedrock amplitudes. The changes are not monotonic with noise 
level.

At higher noise levels, the estimated amplitudes can differ considerably 
from the noise-free values.

---

## Model Parameters

The main model parameters are:

| Parameter | Value |
|-----------|-------|
| Domain length \(L\) | \(10^5\) |
| \(A\) | 300 |
| \(\rho\) | \(9.138\times10^{-19}\) |
| \(g\) | \(9.7692\times10^{15}\) |
| Flow-law exponent \(n\) | 3 |
| Grid points | 51 |
| Initial ice thickness | 800 m |
| Base slope angle | \(0.2^\circ\) |

---

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

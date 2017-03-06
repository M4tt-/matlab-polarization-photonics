# matlab-polarization-photonics
A small collection of .m's that I've used in photonics experiments.

Each file is independent of one another. Ensure your MATLAB working directory can see each of these files if you seek to use them.

PlotPS.m - Plots a single Stokes vector on the Poincare Sphere.
StokesPolarimetry.m - Perform Stokes polarimetry to characterize a general polarization state. 
calcFidelity.m - Calculates the quantum fidelity between two density matrices.
calcPurity.m - Calculates the quantum purity of a density matrix.
calcSNR_CCD.m - Calculates the signal-to-noise ratio from intensity data. Ideal for CCD cameras.
gaussFit.m - Fits gaussian parameters mu (mean) and sigma (std) to a normal-like distibution.
jVec2sVec.m - Computes a Stokes vector from a Jones vector.
sVec2ell.m - Compute elliptical parameters (orientation, eccentricity, semi-major-minor axes lengths) from Stokes vector.

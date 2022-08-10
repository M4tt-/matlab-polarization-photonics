# matlab-polarization-photonics

A small collection of .m's that I've used in photonics experiments.

## Getting Started

Pull the master branch of this repository and ensure your MATLAB working directory can see each of these files if you seek to use them.

## Generally Useful Content

These files could be generally useful for engineers working with polarization optics.

- `gaussFit.m` - Fits gaussian parameters mu (mean) and sigma (std) to a normal-like distibution.
- `calcFidelity.m` - Calculates the quantum fidelity between two density matrices.
- `calcPurity.m` - Calculates the quantum purity of a density matrix.
- `PlotPS.m` - Plots a single Stokes vector on the Poincare Sphere.
- `jVec2sVec.m` - Computes a Stokes vector from a Jones vector.
- `sVec2ell.m` - Compute elliptical parameters (orientation, eccentricity, semi-major-minor axes lengths) from Stokes vector.
- `stokes2DensityMat.m` - Covert a Stokes vector to a density matrix

## Full Content

- `PlotPS.m` - Plots a single Stokes vector on the Poincare Sphere.
- `StokesPolarimetry.m` - Perform Stokes polarimetry to characterize a general polarization state. 
- `calcFidelity.m` - Calculates the quantum fidelity between two density matrices.
- `calcPurity.m` - Calculates the quantum purity of a density matrix.
- `calcSNR_CCD.m` - Calculates the signal-to-noise ratio from intensity data. Ideal for CCD cameras.
- `gaussFit.m` - Fits gaussian parameters mu (mean) and sigma (std) to a normal-like distibution.
- `jVec2sVec.m` - Computes a Stokes vector from a Jones vector.
- `sVec2ell.m` - Compute elliptical parameters (orientation, eccentricity, semi-major-minor axes lengths) from Stokes vector.
- `stokes2DensityMat.m` - Covert a Stokes vector to a density matrix
- `calibratePhase_CCD.m` - A calibration script for the phase-grayscale relationship of a spatial light modulator. Uses a CCD camera.

## Author

Matt Runyon

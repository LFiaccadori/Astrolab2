Overview
This repository is organized around two main analytical tasks:

Bias Analysis

Understand and process bias frames to evaluate the readout noise and electronic offset in CCD detectors.
Compute the median bias frame for use in calibration.
Analyze the statistical properties of the bias and validate assumptions about its uniformity across the frame.
Flat Field Analysis

Process flat frames to correct for non-uniform pixel sensitivity, dust, and optical effects like vignetting.
Normalize flat frames and compute the median flat frame for calibration.
Study variations in flat illumination and their sources (e.g., lamp fluctuations).
Key Features
Data Preprocessing:
Includes tools to read and organize astronomical FITS datasets using Python libraries like astropy, numpy, and matplotlib.
Statistical Analysis:
Analyze readout noise, pixel response variations, and errors associated with both bias and flat frames.
Error Propagation:
Evaluate the total error in calibration frames, including contributions from readout noise, photon noise, and temporal variations in flat illumination.
Visualization:
Generate diagnostic plots to validate assumptions and visualize calibration frames.
Dependencies
To run the code in this repository, make sure you have the following Python libraries installed:

numpy
matplotlib
astropy
pickle
How to Use
Download the Dataset

Obtain the FITS files for your group and move them to the appropriate directory (e.g., TASTE_analysis folder).
Bias Analysis

Extract information from bias FITS files.
Compute the median bias frame.
Perform statistical analysis on the bias data to evaluate the readout noise and uniformity.
Flat Field Analysis

Normalize flat frames using the computed bias.
Analyze and correct for non-uniform illumination across the detector.
Compute the median flat field for final calibration.
Run Jupyter Notebooks

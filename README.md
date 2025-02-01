# TASTE Analysis Pipeline

## Overview
This repository contains a complete pipeline for processing astronomical data from the TASTE project. The analysis follows a structured sequence, from bias and flat corrections to aperture and differential photometry. The project is implemented in Python using Jupyter notebooks and relies on the `astropy`, `numpy`, `matplotlib`, and `pickle` libraries for data handling and visualization.

## Project Structure
The analysis is divided into the following steps:

1. **Bias Analysis**  
   - Loads and organizes bias frames.
   - Extracts and analyzes FITS metadata.
   - Computes the median bias frame and its statistical properties.
   - Saves processed data for further use.

2. **Flat Analysis**  
   - Reads and processes flat-field frames.
   - Identifies and removes overscan regions.
   - Normalizes the flat frames for uniform response calibration.
   - Saves the median flat frame for later corrections.

3. **Science Frame Analysis**  
   - Applies gain correction, bias subtraction, and flat-field correction.
   - Estimates error propagation through the corrections.
   - Saves the processed science images.

4. **Centroid Measurement**  
   - Determines the centroid of the target star in each frame.
   - Implements iterative refinement for precise centroid positioning.

5. **Aperture Photometry**  
   - Implements sky background subtraction.
   - Computes aperture photometry for target and reference stars.
   - Tests multiple aperture sizes to optimize signal-to-noise ratio.

6. **Aperture Photometry Class**  
   - Defines a Python class for handling aperture photometry in a structured manner.
   - Enables easy reusability and integration across multiple datasets.

7. **Differential Photometry**  
   - Normalizes and calibrates the photometry using reference stars.
   - Tracks environmental and instrumental variations.
   - Selects the best combination of reference stars to minimize noise.

## Installation
Ensure you have Python 3.x installed along with the required dependencies:

```sh
pip install numpy matplotlib astropy pickle
```

Alternatively, you can create a virtual environment:

```sh
python -m venv taste_env
source taste_env/bin/activate  # On Windows, use taste_env\Scripts\activate
pip install -r requirements.txt
```

## Usage
The workflow is designed to be executed in Jupyter notebooks:

1. Download and organize the dataset in the `TASTE_analysis` directory.
2. Run the bias, flat, and science correction notebooks sequentially.
3. Process centroid measurement and perform aperture photometry.
4. Use the `AperturePhotometry` class to extract photometric data.
5. Perform differential photometry and analyze the final light curves.

## Output
Processed data is saved in pickle files for easy retrieval. The final results include:
- Median bias and flat frames.
- Corrected science frames.
- Aperture photometry measurements.
- Differential photometry light curves.

## Data Attribution
The observational data used in this project is provided by the **University of Padova**. The data is not included in this repository and is subject to the terms and conditions set by the University.

## License
The code and analysis in this repository are developed by the project owner and released under the **MIT License**. See `LICENSE` for details. Note that the license applies only to the code and methodology, not to the observational data, which remains the property of the **University of Padova**.

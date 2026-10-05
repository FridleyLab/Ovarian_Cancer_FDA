
# Ovarian_Cancer_FDA

Code associated with the manuscript:

"Functional Data Analysis of Spatial Clustering Identifies Prognostic T Cell Patterns in Ovarian Cancer."

## Files

### FPCA_Simulation_Study.R

Functions and analysis code used to:

- Generate simulated mulitplex immunofluorescence (mIF) ROIs
- Estimate spatial clustering trajectories
- Perform functional principal component analysis (FPCA)
- Fit survival models for both fixed radii and FPCA methods

### FPCA_Survival_Analysis.R

Functions and analysis code used to:

- Generate spatial clustering trajectories
- Perform FPCA
- Fit survival models
- Conduct study-specific analyses

### FPCA_vs_Fixed_Radius_Analysis.R

Functions and analysis code used to:

- Generate spatial clustering trajectories
- Perform FPCA
- Fit survival models for both fixed radii and FPCA methods

### Meta_Analysis_Plots.Rmd

R Markdown workflow used to:

- Compile study-specific results
- Perform random-effects meta-analysis

### Major R Packages Utilized

***`scSpatialSim`***, ***`spatialTIME`***, ***`mxfda`***, ***`survival`***, ***`meta`***


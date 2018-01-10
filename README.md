# pygs
Library dedicated to Geospatial computations. Individual functions can be understood using their docstrings  
Modules include:  
gslinalg.py - Computations related to regression and correlation for 3D datasets (Spatio temporal datasets)  
gsncl.py - Porting NCL like functions to Python like area averaging, area sum and volume computation  
gsplots.py - Simplifying geospatial plotting and hovmoller diagrams. Other functions will be added as and when required  
statmisc.py - Small modele with statistical functions related to significance level calculations, time series, and so on... 

### Changelog v0.12
1. Function to compute spatial autocorrelation
2. Separate function padding functions for nd arrays
3. Extend options for contour map and cyclic points
4. MIT Licensing

### Changelog v0.11
1. netCDF4 variables can now be coordinate subscripted
2. Masking in wgt_area_average fixed
3. Added thresholding to Contour maps (bad idea)
4. Added distutils packaging with setup.py. Package can now be installed with setup.py or pip

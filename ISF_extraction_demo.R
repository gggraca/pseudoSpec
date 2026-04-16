#### Script to obtain a in-source spectra from an LC-MS chromatogram #####
#### Goncalo Graca, Imperial College London, 2026 ####

# This script requires two R packages: 'MSnbase', 'ggplot2' and 'gridExtra'
# The following section checks if the packages are installed and loads them
# otherwise, these will be installed and loaded
if(!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if(!require("MSnbase", quietly = TRUE))
    BiocManager::install("MSnbase")
if(!require("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if(!require("gridExtra", quietly = TRUE))
    install.packages("gridExtra")

# The R function 'extractISF.R' is needed and should be stored in the working directory
# it then needs to be loaded:
source("extractISF.R")

# In this example we will extract an in-source fragment spectrum from a
# lipidomics human serum LC-MS chromatogram that can be downloaded from:
# https://zenodo.org/records/17408169/files/Lipid_Positive_QC.mzML
# and stored into the working directory
# We'll extract the in-source fragment spectrum for the feature 585.2691 m/z and 72.8 s
# Run the 'extractISF' function with the default options: 
# RT tolerance of 10s, mz tolerance of 0.01 and correlation threshold of 0.8
# By default an image file containing the correlated ions EICs and ISF pseudo-spectrum
# will be saved in the working directory
# Optionally, the feature EIC and apex spectrum can also be saved into the working directory as .pdf
# as well as the The resulting ISF spectrum (m/z and intensity) and correlation coefficient values as .csv
# if plotResults and  saveResults are set to TRUE, respectively:
extractISF(dataPath="Lipid_Positive_QC.mzML", featureMz=585.2691, 
            featureRT=72.8, corThresh=0.8, RTtol=10, mztol=0.01, 
            plotResults=TRUE, saveResults=TRUE)

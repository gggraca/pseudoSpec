# Extract in-source spectra from an LC-MS chromatograms 

Code to obtain in-source fragmention pseudo-MS mass spectra.

### Demo script to obtain a in-source spectra from an LC-MS chromatogram 

The demo example presented below makes use of the function 'extractISF.R', and also depends on the R packages 'MSnbase', 'ggplot2' and 'gridExtra'.

Firstly, check that the packages are installed:
```
if(!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if(!require("MSnbase", quietly = TRUE))
    BiocManager::install("MSnbase")
if(!require("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if(!require("gridExtra", quietly = TRUE))
    install.packages("gridExtra")
```

Also, the 'extractISF.R' function must be downloaded into the working directory and loaded:
```
source('extractISF.R')
```

This example shows how to extract an in-source fragment spectrum from a lipidomics human serum LC-MS chromatogram that can be downloaded from:
[https://zenodo.org/records/17408169/files/Lipid_Positive_QC.mzML]

The file should be stored into the working directory.

To extract the in-source fragment spectrum for the feature 585.2691 m/z and 72.8 s, the following command should be run:
```
extractISF(dataPath="Lipid_Positive_QC.mzML", featureMz=585.2691, 
            featureRT=72.8, corThresh=0.8, RTtol=10, mztol=0.01, 
            plotResults=TRUE, saveResults=TRUE)
```
By default an image file containing the correlated ions EICs and ISF pseudo-spectrum will be saved in the working directory:
[pseudoMS.png]
Optionally, the feature EIC and apex spectrum can also be saved into the working directory as .pdf, as well as the The resulting ISF 
spectrum (m/z and intensity) and correlation coefficient values as .csv, if plotResults and  saveResults are set to TRUE, respectively:
[feature_selection.png]


# Extract in-source fragments from MS1 spectra given a feature
# Using PCA on EICs to select EICs with score close to the feature of interest
# Might not work exactly as well as the example below if there is significant co-elution

# Read data and packages
library(MSnbase)
library(ggplot2)
library(pcaMethods)
rawData <- readMSData('Lipid_Positive_QC.mzML', mode='onDisk', centroided. = TRUE)

# Filter MS1 spectra. This is important particularly when the data has more than one type of scans, which is the case.
msFilteredData <- filterMsLevel(rawData, msLevel = 1)

# Read the retention times of the filtered MS1 scans.
RTs <- rtime(msFilteredData)

# Extract the EIC for the feature of interest (585.2691 m/z and 72.8 s) and locate its apex.
# The retention time window here is set to approximately +/- 20 s, but more accurate values could be obtained from the peakPantheR ROIs.
# Same applies to the mz values, which here an interval of +/- 0.01 was considered.

scans <- which(RTs >= 60 & RTs <= 90) # scans corresponding to the feature of interest

rtFilteredData <- msFilteredData[scans] # only scans at the RTs of interest

spectraData <- spectra(rtFilteredData) # spectra for corresponding to the scans of interest 

feic <- filterMz(rtFilteredData, mz = c(585.2591, 585.2791)) # EIC object

fint <- as.numeric(intensity(feic)) # intensity for the EIC

frt <- rtime(feic)

apex <- which.max(fint) # apex for the EIC

# Gather the m/z and intensity values for the spectrum corresponding to the apex of the EIC of the feature of interest. 
# The most intense features can be used, for instance above 10000 of intensity.
apexSpec <- data.frame(mz = mz(spectraData[[apex]]), intensity = intensity(spectraData[[apex]])) # spectrum to get ISFs from
idx <- which(apexSpec$intensity > 10000)
mzSelected <- apexSpec$mz[idx]

# Generate the EICs for the selected peaks. 
# function to get all EICs for the selected peaks from the apex spectrum
# mz tolerance set to 0.01
getEIC <- function(mz, rtFilteredData, tol = 0.01){
  eic <- filterMz(rtFilteredData, mz = c(mz-tol, mz+tol))
  rt <- as.numeric(rtime(eic))
  int <- intensity(eic)
  tmp <- lapply(int, function(x) which(length(x) > 1 | length(x) == 0)) # some features contain 2 intensity values?
  idx <- which(tmp == 1) # map features with two intensity values
  idx <- as.numeric(idx)
  int[idx] <- NA # replace features with two intensity values by NA
  int <- unlist(int)
  df <- data.frame(rt =rt, intensity = int)
  return(df)
}

r <- lapply(mzSelected, function(x) getEIC(x, rtFilteredData, tol = 0.01))

# create a matrix for the EICs including the feature of interest
m <- data.frame(intensity = fint)
colnames(m) <- as.character(585.2691)
for(i in 1:length(mzSelected)){
  tmp <- round(mzSelected[i], 4)
  m <-  cbind(m, r[[i]][,2])
  colnames(m)[i+1] <- tmp
}

# make one EIC per row 
m <- t(m)
# PCA with NIPALS to account for NAs
p <- pca(m, method="nipals", nPcs = 5, scale = "uv", center = TRUE)

# features more similar to the feature of interest should be closer to it in PC1
# calculate distance to determine the other points
d <- abs(p@scores[1,1] - p@scores[,1])
# select the closest < 1
idx <- which(d < 1 & d > 0)
idx <- idx-1 # -1 one to account for the feature of interest removal 
# select features with the highest correlation coefficient
isf <- mzSelected[idx]

# obtain the pseudo-in source fragment (ISF) spectrum
idx <- sapply(isf, function(x) grep(x, apexSpec$mz))
isfSpec <- apexSpec[idx,]



#' Function to obtain pseudo-MS in-source fragment spectrum 
#' for a feature of interest of an LC-MS experiment
#'
#' @author  Goncalo Graca (Imperial College London)
#' 11 June 2025
#'
#' @param rawData OnDiskMSnExp object from MSnbase
#' @param rtMin minimum retention time of the feature of interest
#' @param rtMax maximum retention time of the feature of interest
#' @param mzMin minimum m/z value of the feature of interest
#' @param mzMax maximum m/z value of the feature of interest
#' @param mztol m/z tolerance for fragment search
#' @param itol intensity tolerance for feature mass spectrum peak-picking 
#' @param cthr Pearson correlation coefficient threshold for EIC profile correlation
#' @return the in-source fragment pseudo-MS spectrum as a 'Spectra1' class objects from the MSnbase package
#' @export

getISF <- function(rawData,
                   rtMin,	rtMax,	
                   mzMin,	mzMax, 
                   mztol = 0.01, 
                   itol = 10000, 
                   cthr = 0.95){
  
  # Filter MS1 spectra. This is important particularly when the data has more than one type of scans, which is the case.
  msFilteredData <- filterMsLevel(rawData, msLevel = 1)
  
  # Read the retention times of the filtered MS1 scans.
  RTs <- rtime(msFilteredData)
  
  # get feature EIC and apex
  scans <- which(RTs >= rtMin & RTs <= rtMax) # scans corresponding to the feature of interest
  rtFilteredData <- msFilteredData[scans] # only scans at the RTs of interest
  spectraData <-spectra(rtFilteredData) # spectra for corresponding to the scans of interest 
  feic <- filterMz(rtFilteredData, mz = c(mzMin, mzMax)) # EIC object
  fint <- as.numeric(intensity(feic)) # intensity for the EIC
  frt <-rtime(feic)
  apex <- which.max(fint) # apex for the EIC
  
  # Gather the m/z and intensity values for the spectrum corresponding to the apex of the EIC of the feature of interest
  # spectrum to get ISFs from:
  apexSpec <- data.frame(mz = mz(spectraData[[apex]]), intensity = intensity(spectraData[[apex]])) 
  idx <- which(apexSpec$intensity > itol)
  mzSelected <- apexSpec$mz[idx]
  
  # Generate the EICs for the selected peaks from the apex spectrum
  # function to get all EICs for the selected peaks from the apex spectrum
  # mz tolerance set to 0.01
  getEIC <- function(mz, rtFilteredData, mztol){
    eic <- filterMz(rtFilteredData, mz = c(mz-mztol, mz+mztol))
    rt <- as.numeric(rtime(eic))
    int <- intensity(eic)
    tmp <- lapply(int, function(x) which(length(x) > 1 | length(x) == 0))
    idx <- which(tmp == 1)
    idx <- as.numeric(idx)
    int[idx] <- NA
    int <- unlist(int)
    df <- data.frame(rt =rt, intensity = int)
    return(df)
  }
  
  r <- lapply(mzSelected, function(x) getEIC(x, rtFilteredData, mztol))
  
  # Obtain the in-source fragment pseudo-spectrum  from the EICs of the features
  # that best correlate with the EIC of the feature of interest
  
  # correlate the feature EIC with the other features EICs
  c <- lapply(r, function(x) cor(fint, x[,2], use = "pairwise.complete.obs"))
  c <- unlist(c)
  
  # select features with the highest correlation coefficient
  isf <- mzSelected[which(c > cthr)]
  
  # obtain the pseudo-in source fragment (ISF) spectrum
  idx <- sapply(isf, function(x) grep(x, apexSpec$mz))
  isfSpec <- apexSpec[idx,]
  
  # convert the feature and the ISF spectra into "Spectrum1" class objects
  # to use spectra comparison functionality from the MSnbase package
  isfSpec <- new("Spectrum1", mz = isfSpec$mz, intensity = isfSpec$intensity, 
                 centroided = TRUE,
                 polarity = unique(polarity(rawData)),
                 rt = as.numeric(frt[apex]))
  
  return(isfSpec)
}

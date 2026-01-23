library(MSnbase)
library(ggplot2)

# set working directory
setwd("/Users/GG/Documents/")

# read data
rawData <- readMSData("Lipid_Positive_QC.mzML", mode="onDisk", centroided.=TRUE)

# filter MS1 scans
msFilteredData <- filterMsLevel(rawData, msLevel=1)

# read RTs of all the scans
RTs <- rtime(msFilteredData)

# Extract the EIC for the feature of interest (585.2691 m/z and 72.8 s) 
# and locate its apex.
# The retention time window here is set to approximately +/- 20 s, 
# but more accurate values could be obtained from the peakPantheR ROIs. 
# Same applies to the mz values, which here an interval of +/- 0.01 was considered.

scans <- which(RTs >= 60 & RTs <= 90) # scans corresponding to the feature of interest
rtFilteredData <- msFilteredData[scans] # only scans at the RTs of interest
spectraData <- spectra(rtFilteredData) # spectra for corresponding to the scans of interest

# generate EIC object for the feature of interest
feic <- filterMz(rtFilteredData, mz = c(585.2591, 585.2791)) 

# separate intensity and RT of the EIC
fint <- as.numeric(intensity(feic)) 
frt <- rtime(feic)
# apex for the EIC
apex <- which.max(fint) 

plot(frt, fint, type="b", pch=19, xlab="Retention time (s)", 
     ylab="Intensity (a.u.)", main="585.2691 m/z, 72.8s")
abline(v = frt[apex], col = "red", lty = 3)

# Gather the m/z and intensity values for the spectrum corresponding to the apex 
# of the EIC of the feature of interest.
# spectrum to get ISFs from:
apexSpec <- data.frame(mz=mz(spectraData[[apex]]), 
                       intensity=intensity(spectraData[[apex]])) 

# estimate noise intensity as 10x the upper hinge of the intensity boxplot
# this is empirically defined...there might be a better way of doing this:
bplt <- boxplot(apexSpec$intensity, plot=FALSE) 
thres <- bplt$stats[5]*10
idx <- which(apexSpec$intensity > thres)
mzSelected <- apexSpec$mz[idx] # selection

# plot spectrum with RT threshold
plot(apexSpec[,1],apexSpec[,2], type="h", xlab="m/z", ylab="Intensity", 
     main=paste("Retention time", frt[apex], "s"))
abline(h = thres, col = "red", lty = 3)

# Generate EICs
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

# correlate the feature EIC with the other features EICs
c <- lapply(r, function(x) cor(fint, x[,2], use = "pairwise.complete.obs"))
c <- unlist(c)

# select features with the highest correlation coefficient
mzISF <- mzSelected[which(c > 0.95)]

# obtain the pseudo-in source fragment (ISF) spectrum
idx <- sapply(mzISF, function(x) grep(x, apexSpec$mz))
isfSpec <- apexSpec[idx,]

# Plot the overlap between EICs
selEICs <- r[which(c > 0.95)]

for(i in seq_along(selEICs)){
    selEICs[[i]][,"mz"] <- round(isfSpec[i], 4)
}

df <- do.call(rbind,selEICs)
p1 <- ggplot2::ggplot(df[!is.na(df$intensity),],
                      ggplot2::aes(x=rt, y=intensity, colour=factor(mz))) +
    ggplot2::geom_line() +
    ggplot2::labs(x="RT (s)", y = "Intensity (a.u.)", 
                  colour="fragments") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))

p1

# function to plot MS pseudo-spectrum:
plotPseudoMS <- function(spectrum){
    df <- as.data.frame(spectrum)
    p <- ggplot2::ggplot(df,
                          ggplot2::aes(x=mz, y=intensity, label=round(mz, 3))) +
        ggplot2::geom_segment(ggplot2::aes(xend=mz, yend=0),
                              color="red", lwd=0.5) +
        ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ylim(0, max(df[,2]) + 0.1*max(df[,2])) +
        ggplot2::xlim(min(df[,1])-50, max(df[,1])+50) +
        ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")
    p
}

plotPseudoMS(isfSpec)
# put the correlation coefficient in the table

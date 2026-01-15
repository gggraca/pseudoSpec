#' Extracts the in-source fragment (ISF) spectrum of a given feature from LC-MS data
#'  
#' @description
#' Function to extract the in-source fragment spectrum for one features of interest
#' Save the result to workspace and plots and table to file (optionally). Requires
#' functionality from "MSnbase", "ggplot2" and "gridExta" packages. 
#' 
#' @author Goncalo Graca
#' 
#' @param dataPath Either the file name (.mzML) if present in the working directory,
#' or the full path to its location.
#' @param featureMz The feature of interest m/z.
#' @param featureRT The feature of interest RT in seconds.
#' @param corThresh Correlation coefficient threshold (default 0.8).
#' @param RTtol RT tolerance in seconds, to define the RT window where to search the 
#' feature of interest and related ions.
#' @param mztol m/z tolerance for searching the feature of interest.
#' @param plotResults To specify if feature EIC and apex spectrum, the correlated 
#' EICs and ISF plots should be saved to file (pdf). The default value is TRUE. 
#' @param saveResults The resulting ISF spectrum (m/z and intensity) 
#' and correlation coefficient values are saved to .csv. The default is TRUE.
#' 
#' @importFrom MSnbase readMSdata
#' @importFrom gridExtra grid.arrange
#' 
#' @return Stores the extracted ISF spectrum to the workspace environment,
#' and saves feature and ISF plots and table to disk if "plotResults" and "saveResults"
#' are set to "TRUE".
#'  
extractISF <- function(dataPath, featureMz, featureRT, 
						corThresh=0.8, RTtol=10, mztol=0.01, 
						plotResults=TRUE, saveResults=TRUE){
	# read data from the datapath
	rawData <- readMSData("Lipid_Positive_QC.mzML", mode="onDisk", centroided.=TRUE)
	msFilteredData <- filterMsLevel(rawData, msLevel=1)
	
	# read RTs of all the scans
	RTs <- rtime(msFilteredData)
	scans <- which(RTs >= featureRT-RTtol & RTs <= featureRT+RTtol) # scans corresponding to the feature of interest
	rtFilteredData <- msFilteredData[scans] # only scans at the RTs of interest
	spectraData <- spectra(rtFilteredData) # spectra for corresponding to the scans of interest
	
	# generate EIC object for the feature of interest
	featureEIC <- filterMz(rtFilteredData, mz=c(featureMz-mztol, featureMz+mztol))
	# separate intensity and RT of the EIC
	fint <- as.numeric(intensity(featureEIC)) 
	frt <- rtime(featureEIC)
	# apex (index) of the EIC
	apex <- which.max(fint)
	# apex spectrum to get ISFs from:
	apexSpec <- data.frame(mz=mz(spectraData[[apex]]), 
						   intensity=intensity(spectraData[[apex]]))
	# estimate noise intensity as 10x the upper hinge of the intensity boxplot
	# this is empirically defined...there might be a better way of doing this:
	bplt <- boxplot(apexSpec$intensity, plot=FALSE) 
	thres <- bplt$stats[5]*10
	idx <- which(apexSpec$intensity > thres)
	mzSelected <- apexSpec$mz[idx] # selection
	
	# Generate EICs
	r <- lapply(mzSelected, function(x) getEIC(x, rtFilteredData, tol=mztol))
	
	# Correlate the feature EIC with the other features EICs
	c <- lapply(r, function(x) cor(fint, x[,2], use = "pairwise.complete.obs"))
	c <- unlist(c)
	
	# select features with the highest correlation coefficient
	mzISF <- mzSelected[which(c > corThresh)]
	
	# obtain the pseudo-in source fragment (ISF) spectrum
	idx <- sapply(mzISF, function(x) grep(x, apexSpec$mz))
	isfSpec <- apexSpec[idx,]
	isfSpec$correlation <- c[which(c > corThresh)] # add correlation coefficient
	
	# Plot the overlap between EICs and ISF
	if(plotResults){
		# plot feature EIC and spectrum
		pdf(file=paste("MS1_feature_",featureMz, "mz_", featureRT, "s", ".pdf", sep=""),
			width=10, height=10)
		par(mfrow=c(2,1))
		plot(frt, fint, type="b", pch=19, xlab="Retention time (s)", 
			 ylab="Intensity (a.u.)", 
			 main=paste("EIC of feature", featureMz, "m/z,", featureRT, "s"))
		abline(v = frt[apex], col = "red", lty = 3)
		# plot spectrum with RT threshold
		plot(apexSpec[,1],apexSpec[,2], type="h", xlab="m/z", ylab="Intensity", 
			 main=paste("MS1 Spectrum at", frt[apex], "s"))
		abline(h = thres, col = "red", lty = 3)
		dev.off()
		# plot correlated EICs
		selEICs <- r[which(c > corThresh)]
		for(i in seq_along(selEICs)){
			selEICs[[i]][,"mz"] <- round(isfSpec[i,1], 4)
		}
		df1 <- do.call(rbind,selEICs)
		p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
							  ggplot2::aes(x=rt, y=intensity, colour=factor(mz))) +
			ggplot2::geom_line() +
			ggplot2::labs(x="RT (s)", y = "Intensity (a.u.)", 
						  colour="fragments") +
			ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
			ggtitle(paste("EICs correlated at", "C >", corThresh, "with feature", 
						  featureMz, "m/z,", featureRT, "s"))
		# plot ISF
		df2 <- as.data.frame(isfSpec)
		p2 <- ggplot2::ggplot(df2,
							  ggplot2::aes(x=mz, y=intensity, label=round(mz, 3))) +
			ggplot2::geom_segment(ggplot2::aes(xend=mz, yend=0), lwd=0.5) +
			ggplot2::geom_text(ggplot2::aes(color=correlation),
							   size=3, angle=45, hjust=0, vjust=0) +
			ggplot2::theme_minimal() +
			ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
			ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
			ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
			ggplot2::labs(x = "m/z", y = "Intensity (a.u.)") +
			ggplot2::scale_color_gradientn(colors=rainbow(nrow(isfSpec))) +
			ggtitle(paste("Pseudo-MS/MS spectrum of feature", featureMz, 
						  "m/z,", featureRT, "s"))
		pdf(file=paste("LCMS_in-source_",featureMz, "mz_", featureRT, "s", ".pdf", sep=""),
			width=10, height=10)
		gridExtra::grid.arrange(p1, p2, nrow=2)
		dev.off()
	}
	
	if(saveResults){
		write.csv(isfSpec, 
				  paste("LCMS_in-source_spectrum_",featureMz, "mz_", featureRT, "s", ".csv", sep=""),
				  row.names=FALSE)
	}
	
	assign("isfSpec", isfSpec, envir = .GlobalEnv)
}	
	
	
## helper functions-------------------------------------------------------

# function to get all EICs for the selected peaks from the apex spectrum
# mz tolerance set to 0.01
getEIC <- function(mz, rtFilteredData, tol=0.01){
	eic <- filterMz(rtFilteredData, mz=c(mz-tol, mz+tol))
	rt <- as.numeric(rtime(eic))
	int <- intensity(eic)
	tmp <- lapply(int, function(x) which(length(x) > 1 | length(x) == 0))
	# some features contain 2 intensity values?
	idx <- which(tmp == 1) # map features with two intensity values
	idx <- as.numeric(idx)
	int[idx] <- NA # replace features with two intensity values by NA
	int <- unlist(int)
	df <- data.frame(rt =rt, intensity=int)
	return(df)
}

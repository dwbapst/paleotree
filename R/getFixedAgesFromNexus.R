


# @param originalNexusFile Filename (and possibly path too) to the original NEXUS file for this analysis.



# this function will look for, scan, and parse an
	# associated NEXUS file. 
	
	
# Ignoring any commented lines (ie. anything between "[   ]" ), commands
# for fixing taxa will be identified, parsed and returned to the user, either as a message
# pinted to the R console if output is read to a file, or as a attribute named 'fixed ages'
# if output as an R object (formatted as a two-column table of OTU names and their respective fixed ages).
# The search for the NEXUS file is controlled with argument \code{originalNexusFile}

# Please note: this has a while() loop in it for removing nested series of
# square brackets (i.e. treated as comments in NEXUS files) then files with
# ridicuously nested series of brackets may cause this code to take a while
# to complete, or may even cause it to hang.


## Example for testing
# origNexusFile <- "D:\\dave\\research\\0 devonian terebrat tip dating\\terebratDev_FAD-LAD_05-08-17.nex"
# getMrBFixedAgesFromNexus(origNexusFile)



getMrBFixedAgesFromNexus <- function(origNexusFile, shouldEndWithNex = TRUE){
	#		
	#
	##############################
	if(shouldEndWithNex){
		# does it end with .nex already?
		endNex <- grepl("\\.nex$", 
			origNexusFile , ignore.case=TRUE)
		# if the nexus file does not end with .nex, add it
		if(!endNex){
			origNexusFile <- paste0(origNexusFile,".nex")
			}
		}
	#############################
	# get the nexus file
	origNexus <- scan(file = origNexusFile,
		what = character(), #sep = "\n",
		quiet = TRUE, 
		comment.char = "[", 
		strip.white = TRUE)
	#
	origNexus <- readChar(origNexusFile, 
		file.info(origNexusFile)$size)
	# remove white space
	origNexus <- gsub(origNexus,pattern = " ",replacement = "")
	origNexus <- gsub(origNexus,pattern = "\t",replacement = "")
	# remove brackets and everything between brackets	
	while(grepl(origNexus,pattern = "\\[|\\]")){
		origNexus <- gsub(origNexus,pattern = "\\[[^\\[\\]]*\\]",replacement = "",
			perl = TRUE)
		}
	# split at new lines
	origNexus <- unlist(strsplit(x = origNexus,split = "\r",fixed = TRUE))
	origNexus <- unlist(strsplit(x = origNexus,split = "\n",fixed = TRUE))
	# remove ""
	origNexus <- origNexus[origNexus != ""]
	# find lines with calibrate
	hasCalibrate <- grepl(pattern = "calibrate",x = origNexus,ignore.case = TRUE)
	hasCalibrate <- origNexus[hasCalibrate]
	# remove calibrate
	hasCalibrate <- gsub(x = hasCalibrate,pattern = "calibrate",replacement = "")
	# find calibrate lines with fixed
	hasFixed <- grepl(pattern = "=fixed\\(",x = hasCalibrate,ignore.case = TRUE)
	hasFixed <- hasCalibrate[hasFixed]
	# remove fixed
	hasFixed <- gsub(x = hasFixed,pattern = "fixed\\(",replacement = "")
	# remove ;]
	hasFixed <- gsub(x = hasFixed,pattern = "\\);",replacement = "")
	fixedList <- strsplit(hasFixed,split = "=")
	fixedMatrix <- matrix(sapply(fixedList,function(x) x),,2,byrow = TRUE)
	fixedTable <- as.data.frame(fixedMatrix)
	fixedTable[,2] <- as.numeric(fixedMatrix[,2])
	if(nrow(fixedTable)==0){
		warning("No OTUs with fixed ages found in original NEXUS file")
		}
	colnames(fixedTable) <- c("OTUname","fixedAge")
	return(fixedTable)
	}


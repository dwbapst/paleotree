
	# @param getFixedAgesNexus If \code{TRUE}, this function will also look for, scan, and parse an
	# associated NEXUS file. Ignoring any commented lines (ie. anything between "[   ]" ), commands
	# for fixing taxa will be identifiedd, parsed and returned to the user, either as a message
	# pinted to the R console if output is read to a file, or as a attribute named 'fixed ages'
	# if output as an R object (formatted as a two-column table of OTU names and their respective fixed ages).
	# The search for the NEXUS file is controlled with argument \code{originalNexusFile}



getMrBFixedAgesFromNexus<-function(origNexusFile){

			
		origNexusFile<-"D:\\dave\\research\\0 devonian terebrat tip dating\\terebratDev_FAD-LAD_05-08-17.nex"


	# get the nexus file
	origNexus<-scan(file=origNexusFile,what = character(), #sep = "\n",
		quiet = TRUE, comment.char = "[", strip.white = TRUE)


	origNexus <- readChar(origNexusFile, file.info(origNexusFile)$size)
	# remove white space
	origNexus<-gsub(origNexus,pattern=" ",replacement="")
	origNexus<-gsub(origNexus,pattern="\t",replacement="")
	# remove brackets and everything between brackets	
	while(grepl(origNexus,pattern="\\[|\\]")){
		origNexus<-gsub(origNexus,pattern="\\[[^\\[\\]]*\\]",replacement="",
			perl=TRUE)
		}

	# split at new lines
	origNexus<-unlist(strsplit(x=origNexus,split="\r",fixed=TRUE))
	origNexus<-unlist(strsplit(x=origNexus,split="\n",fixed=TRUE))
	# remove ""
	origNexus<-origNexus[origNexus!=""]
	# find lines with calibrate
	hasCalibrate<-grepl(pattern="calibrate",x=origNexus,ignore.case=TRUE)
	hasCalibrate<-origNexus[hasCalibrate]
	# remove calibrate
	hasCalibrate<-gsub(x=hasCalibrate,pattern="calibrate",replacement="")
	# find calibrate lines with fixed
	hasFixed<-grepl(pattern="=fixed\\(",x=hasCalibrate,ignore.case=TRUE)
	hasFixed<-hasCalibrate[hasFixed]
	# remove fixed
	hasFixed<-gsub(x=hasFixed,pattern="fixed\\(",replacement="")
	# remove ;]
	hasFixed<-gsub(x=hasFixed,pattern="\\);",replacement="")
	fixedTable<-strsplit(hasFixed,split="=")
	sapply(fixedTable,function(x) x)
	
	
	fixedTable<-
		
	return(fixedTable)
	}

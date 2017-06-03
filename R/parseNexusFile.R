
parseNexusFile<-function(origNexusFile=origNexusFile,asIs=TRUE){
	# if asIs, then the morph Nexus never gets broken down or parsed beyond being scanned
	morphNexusAsIs<-readLines(con=origNexusFile,warn=FALSE)
	if(!asIs){
		# cleaner form based on ape's read.nexus.data
		morphNexus<-scan(file=origNexusFile,what = character(), sep = "\n",
			quiet = TRUE, comment.char = "[", strip.white = TRUE)
		#
		# find the NTAX line
		ntaxLine<-grepl(morphNexus,pattern="\\bNTAX",ignore.case=TRUE)
		if(sum(ntaxLine)>1){
			stop("More than one line containing 'NTAX' found in the provided NEXUS file")}
		if(sum(ntaxLine)<1){
			stop("No line containing 'NTAX' found in the provided NEXUS file")}
		# get number of taxa (more regexp borrowed from read.nexus.data)
		oldNtax<-as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
			"\\3", morphNexus[ntaxLine], perl = TRUE, ignore.case = TRUE))
		# get other pieces of ntax line
		ntaxLineFirst<-sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
			"\\1\\2", morphNexus[ntaxLine], perl = TRUE, ignore.case = TRUE)
		ntaxLineLast<-sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
			"\\4", morphNexus[ntaxLine], perl = TRUE, ignore.case = TRUE)
		#
		# find the matrix line
		matrixLine<-grepl(morphNexus,pattern="\\bMATRIX",ignore.case=TRUE)
		if(sum(matrixLine)>1){
			stop("More than one line containing 'MATRIX' found in the provided NEXUS file")
			}
		if(sum(matrixLine)<1){
			stop("No line containing 'MATRIX' found in the provided NEXUS file")
			}
		#
		# find the next semi-colon line
		matrixEnd<-(which(matrixLine):length(morphNexus))[
			grepl(morphNexus[which(matrixLine):length(morphNexus)],
				pattern=";")][1]
		if(length(matrixEnd)!=1){
			stop("Cannot find semicolon on line following 'MATRIX' line in NEXUS file")
			}
		############################################
		# get Header 1 - start to line right before ntax
		headerOne<-morphNexus[1:(which(ntaxLine)-1)]
		# get header 2 - line right after ntax to MATRIX
		headerTwo<-morphNexus[(which(ntaxLine)+1):which(matrixLine)]
		# get footer - new semicolon line, + semicolon+1 to end
		footer<-c(";",morphNexus[(matrixEnd+1):length(morphNexus)])
		#####################################################################
		# get taxon data
		#
		# isolate the lines after matrix, to semicolon
		taxonLines<-morphNexus[(which(matrixLine)+1):matrixEnd]
		# remove semicolon only lines
		taxonLines<-gsub(";","",taxonLines)
		taxonLines<-taxonLines[taxonLines!=""]
		taxonLines<-gsub("\t"," ",taxonLines)
		# NOW it should be equal to number of taxa 
			# (note if character across multiple lines, would be problematic...)
		if(length(taxonLines)!=oldNtax){
			stop("The number of apparent lines in the matrix of the input NEXUS file doesn't match the given NTAX in the header of that file")
			}
		# 
		# now to get the taxon names and character codings
		taxonData<-regmatches(taxonLines, regexpr(" ", taxonLines), invert = TRUE)
		# clean of whitespace - actually this really isn't necessary
		#taxonData<-lapply(taxonData,sapply,function(z) gsub(pattern=" ",replacement="",z))
		#taxonData<-lapply(taxonData,sapply,function(z) gsub(pattern="\\t",replacement="",z))
		# get names
		taxonNames<-sapply(taxonData,function(z) z[[1]])
		charData<-sapply(taxonData,function(z) z[[2]])
		#
		#########################################
		# make a function
		#
		remakeDataBlockFun<-function(newTaxaTable,taxonNames=taxonNames,charData=charData,
					ntaxLineFirst=ntaxLineFirst,ntaxLineLast=ntaxLineLast,
				headerOne=headerOne,headerTwo=headerTwo,footer=footer){
			# given data on new taxa (with old taxa), rebuild NEXUS block
			# input: a matrix with column 1 = new taxon names
				# column 2 = old taxon names
			# identify all old taxa
			origMatch<-sapply(newTaxaTable[,2],function(x) which(x==taxonNames))
			if(is.list(origMatch)){
				stop("more than one match of taxon names when remaking data NEXUS block")
				}
			# test that all old taxa are in taxonNames
			if(any(is.na(origMatch))){
				stop("No match for original taxon name when remaking data NEXUS block")
				}
			###############################
			# make the new character matrix lines
			newData<-cbind(newTaxaTable[,1],charData[origMatch])
			newTaxonBlock<-apply(newData,1,function(x) paste0(x[1],"  ",x[2]))
			##################
			# make the new block
			#
			# get new ntax line with new ntax
			newNtaxValue<-nrow(newTaxaTable)
			newNtaxLine<-paste0(ntaxLineFirst," ",newNtaxValue," ",ntaxLineLast)
			# make new block
			newBlock<-c(headerOne,newNtaxLine,headerTwo,
				"",newTaxonBlock,"",footer)
			return(newBlock)
			}
		#
		#if NOT asIs, need to output names, and a function for rebuilding NEXUS file
		res<-list(taxonNames=taxonNames,remakeDataBlockFun=remakeDataBlockFun,morphNexusAsIs=morphNexusAsIs)
	}else{
		#if just asIs
		res<-morphNexusAsIs
		}
	return(res)
	}


# origNexusFile<-"D:\\dave\\workspace\\mrbayes\\mat.nex"
# parseNexusFile(origNexusFile=origNexusFile,asIs=TRUE)
# parseNexusFile(origNexusFile,asIs=FALSE)

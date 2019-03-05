# hidden utility functions for getting phylopics


getPhyloPicUIDdataFromPBDB <- function(picIDs){
	imgNumAPIurl <- "https://paleobiodb.org/data1.2/taxa/thumb.txt?id="
	URLwithNums <- paste0(imgNumAPIurl,picIDs)	
	names(URLwithNums) <- names(picIDs)
	res <- do.call(rbind,
          lapply(URLwithNums, read.csv, stringsAsFactors=FALSE))
	return(res)
	}
	
getPhyloPicFromPhyloPic <- function(picUID){
	# get image info
	picInfoURL <- paste0("http://phylopic.org/api/a/image/",
		picUID,"?options=credit+licenseURL+pngFiles")
	if(RCurl::url.exists(picInfoURL)){
		picInfo <- jsonlite::fromJSON(picInfoURLs[i],  
			simplifyVector = FALSE)
		#################
		#
		if(length(picInfo$result$pngFiles)>0){
			#need to check if there is a URL at all
			picPNGurl <- picInfo$result$pngFiles[[
				length(picInfo$result$pngFiles)
				]]$url
			picPNGurl <- paste0("http://phylopic.org",
				picPNGurl)
			##############
			if(RCurl::url.exists(picPNGurl)){
				#need to check if there is a PNG at that URL
				pngContent <- RCurl::getURLContent(picPNGurl)
				pngContent  <-  png::readPNG(pngContent)
			}else{
				pngContent <- NULL
				}
		}else{
			pngContent <- NULL
			}
	}else{
		pngContent <- NULL
		}
	return(pngContent)
	}

	
getPhyloPicPNG_PBDB<-function(
		picID_PBDB, 
		noiseThreshold = 0.1,
		rescalePNG = TRUE, 
		trimPNG = TRUE,
		makeMonochrome = FALSE,
		plotComparison = FALSE){
	############################################
	#	
	# require(png);require(RCurl)
	# png::readPNG RCurl::getURLContent
	#
	# pause 1 second
	Sys.sleep(runif(1,1,2))
	#
	# GET IMAGE
	# get the URL address for the pic via API
	apiPicURL <- paste0(
		"http://paleobiodb.org/data1.2/taxa/thumb.png?id=",
		picID_PBDB)
	#
	# get picPNG		
	picPNG <-  png::readPNG(RCurl::getURLContent(apiPicURL))
	############
	# phylopic PNG is on the fourth slice
	#  image(picPNG [,,4])
	########################
	return(picPNG)
	}
	
getPhyloPicIDNumFromPBDB <- function(taxaData, tree){
	# check or obtain the PBDB phylopic IDs from PBDB
		# get the phylopic-specific IDs  as well
	#
	if(is.null(taxaData)){
		# get image ID numbers using PBDB API calls for each
			# tip taxon in the tree using the tip labels
		tiptaxa <- paste0(tree$tip.label, 
			collapse = ",")
		apiAddressTaxa <- paste0(
			"http://paleobiodb.org/data1.2/taxa/list.txt?name=",
			tiptaxa, "&rel=exact&show=img"
			)	
		# call PBDB API
		tiptaxaData <- read.csv(apiAddressTaxa,
			stringsAsFactors = FALSE)
		# get the image IDs
		phylopicIDsPBDB<- tiptaxaData $image_no[
			match(tree$tip.label, tiptaxaData $taxon_name)]
		names(phylopicIDsPBDB) <- tree$tip.label
	}else{
		#
		# get phylo pic IDs and label with tip labels
		phylopicIDsPBDB <- as.character(taxaData$image_no[
			match(tree$tip.label, taxaData$taxon_name)])
		names(phylopicIDsPBDB) <- tree$tip.label
		#
		# CHECKS
		# does the provided number of IDs equal the number
			# of tips?
		if(length(phylopicIDsPBDB) != Ntip(tree)){
			#print(tree)
			#print(phylopicIDsPBDB)
			stop(paste0(
				"If provided, phylopicIDsPBDB must be have the same",
				" length as the number of tips on the input tree"))
			}
		# make sure its a character vector
		if(!is.character(phylopicIDsPBDB)){
			phylopicIDsPBDB <- as.character(phylopicIDsPBDB)
			if(!is.character(phylopicIDsPBDB)){		
				stop("Cannot coerce phylopicIDsPBDB to be type character")
				}
			}
		# Are there NAs or blanks?
		if(any(is.na(phylopicIDsPBDB))){
			stop("If provided, phylopicIDsPBDB cannot contain any NA values")
			}
		if(any(phylopicIDsPBDB == "")){
			stop("If provided, phylopicIDsPBDB cannot contain any blank values")
			}
		# make sure user provided IDs are labeled
		if(length(names(phylopicIDsPBDB)) != Ntip(tree)){
			print(length(names(phylopicIDsPBDB)))
			print(phylopicIDsPBDB)
			stop(paste0("If provided, phylopicIDsPBDB must have",
				" labels matching the tree's tip labels"))	
			}
		}
	# make sure phylopicIDsPBDB is sorted relative to tip order
	phylopicIDsPBDB <- phylopicIDsPBDB[
		match(tree$tip.label, names(phylopicIDsPBDB))
		]
	return(phylopicIDsPBDB)
	}
	
prepPhyloPic<-function(
		picPNG, 
		noiseThreshold = 0.1,
		rescalePNG = TRUE, 
		trimPNG = TRUE,
		makeMonochrome = FALSE,
		plotComparison = FALSE){
	############################################
	# RESCALE PALETTE
	sliceOriginal <- picPNG[,,4]
	if(rescalePNG){
		#rescale pic so that min is 0 and max is 1
		picPNG[,,4] <- picPNG[,,4]-abs(min(picPNG[,,4]))
		picPNG <- picPNG/max(picPNG)
		#picPNG <- picPNG^0.75
		}
	#
	###################
	# TRIM THE PHYLOPIC
		# lots of phylopics have contiguous whitespace at the top/bottom
	sliceContrasted <- picPNG[,,4]
	if(trimPNG){
		sliceContrasted[sliceContrasted  < noiseThreshold] <- 0
		# find all rows of the PNG from the top AND bottom
			# that are non-contiguous whitespace
			# these are to be SAVED
			# use 'contrasted' version
		saveRows <- ((cumsum(apply(sliceContrasted , 1, sum)) > 0)
			 & rev(cumsum(rev(apply(sliceContrasted , 1, sum))) > 0))
		#
		# turns out lots phylopics also have whitespace on their left/right
			# need to trim this too
		# find all COLUMNS of the PNG from the RIGHT AND LEFT
			# that are non-contiguous whitespace
			# these are to be SAVED
		saveCols <- ((cumsum(apply(sliceContrasted , 2, sum)) > 0)
			 & rev(cumsum(rev(apply(sliceContrasted , 2, sum))) > 0))
		#
		# remove whitespace from all slices of the array
		picPNG <- picPNG[saveRows,saveCols,]
		}
	##############
	if(makeMonochrome){
		picPNG <- picPNG^0.001
		}
	##############
	#
	if(plotComparison){
		# plots a comparison of three images
			# as a diagnostic mode
		layout(1:3)
		par(mar=c(0,0,0,0))
		graphics::image(sliceOriginal)
		graphics::image(sliceContrasted)
		graphics::image(picPNG [,,4])
		}
	return(picPNG)
	}
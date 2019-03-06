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
		picID_PBDB){
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
	
matchTaxaColor <- function(
			taxaColorOld, 
			taxaNames,
			transparency = 1
			){
		###########################################
		taxaColorNew <- rep(NULL, length(taxaNames))
		if(!is.null(taxaColorOld)){
			# if taxaColorOld is NULL, all are 'black'
				# just make it black
				# make all taxonColor NULL
			# don't have to anything because its already null!
			##
			##
			if(!is.character(taxaColorOld)){
				# else, taxaColorOld must be type character
					# if not, FAIL
				stop("taxaColor must be either NULL or type character")
				}
			if(length(taxaColorOld)==1){
				# if its length 1
					# if its a value that matches a tip label, color that taxon "red"
				if(any(taxaColorOld == taxaNames)){
					#this code will make a single 'focal' taxon a bright red	
						# red is "#FF0000FF"
					taxaColorNew[taxaNames == taxaColorOld] <- "red"
				}else{
					# if its a value that does not match a tip label, coerce to a color
						# if not colors, FAIL (will check later when converting to hex values)
						# if colors, all taxa will be in that color
					taxaColorNew <- rep(taxaColorOld ,length(taxaNames))
					}
			}else{
				if(length(taxaColorOld) != length(taxaNames)){
					# if its not length 1, it must be same length as number of tips
						# if not same length as number of tips, FAIL
					stop("If taxaColor is not null or length 1, it must be the same length as number of tips")
					}
				# if colors, each taxa will be in that color, in same order as tip.labels
					# if right length, value is expected to be a color
				taxaColorNew <- taxaColorOld
				# if not colors, FAIL (will check later when converting to hex values)
				}	
		}
	if(any(is.null(taxaColorNew))){
		#
		if(length(transparency) == 1){
			transparency <- rep(transparency, length(taxaNames))
		}else{
			if(length(transparency) == length(taxaNames)){
				stop("If transparency is not length 1, it must be the same length as number of tips")
				}
			}
		#
		for(i in which(!is.null(taxaColorNew))){
			# coerce all non null values to a hex value
				# red should be "#FF0000FF"
			# if not a color, FAIL
			taxaColorNew[i] <- convertColor2Hex(
				colorName = taxaColorNew[i]
				)
			}
		}
	return(taxaColorNew)
	}

convertColor2Hex <- function(
		colorName,
		transparency = 1
		){
	######################
	# check that transparency is between 0 and 1
	if(transparency<0 | transparency>1){
		stop("transparency values must be between 0 and 1")
		}
	alpha <- transparency * 255
	############################
	# coerce all non null values to a hex value
	newColor <- t(col2rgb(colorName))
	#
	newColor <- rgb(
		newColor,
		alpha = alpha,
		maxColorValue = 255
		)
	# red should be "#FF0000FF"
			# yep!
	# if not a color, FAIL	
	return(newColor)
	}

## Cache phylopic images for a tree
## 
## This function takes a tree and saves phylopic thumbnails locally,
## in order to dramatically speed up later rendering.
## 

## @param tree A phylogeny of class \code{phylo} which will be
## plotted, with the terminal tip taxa replaced by silhouettes.
## The tree will be plotted with edge lengths.

## @param taxaDataPBDB  A \code{data.frame} of taxonomic data from
## the Paleobiology Database containing an \code{$image_no} variable,
## as returned when \code{show = "img"} is used. See \emph{Details}.

## @param cacheDir Where to save the output PNGs


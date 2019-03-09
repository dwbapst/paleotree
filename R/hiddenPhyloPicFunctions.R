# hidden utility functions for getting phylopics


getPhyloPicUIDsTableFromPBDB <- function(picIDs){
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
		picInfo <- jsonlite::fromJSON(picInfoURL,  
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
		taxaData <- read.csv(apiAddressTaxa,
			stringsAsFactors = FALSE)
		# get the PBDB image IDs and label with tip labels
		phylopicIDsPBDB<- taxaData$image_no[
			match(tree$tip.label, taxaData$taxon_name)
			]
		names(phylopicIDsPBDB) <- tree$tip.label
	}else{
		#
		# get the PBDB image IDs and label with tip labels
		phylopicIDsPBDB <- taxaData$image_no[
			match(tree$tip.label, taxaData$taxon_name)
			]
		phylopicIDsPBDB <- as.character(phylopicIDsPBDB)
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
	
matchTaxaColor <- function(
		taxaColorOld, 
		taxaNames,
		transparency = 1
		){
	###########################################
	taxaColorNew <- rep(NA, length(taxaNames))
	#print(taxaNames)
	#
	if(!is.null(taxaColorOld)){
		# if taxaColorOld is NULL, all are 'black'
			# just make it black
			# make all taxonColor NA
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
				stop("If taxaColor is not NULL or length 1, it must be the same length as number of tips")
				}
			# if colors, each taxa will be in that color, in same order as tip.labels
				# if right length, value is expected to be a color
			taxaColorNew <- taxaColorOld
			# if not colors, FAIL (will check later when converting to hex values)
			}	
		}
	if(any(!is.na(taxaColorNew))){
		#
		if(length(transparency) == 1){
			transparency <- rep(transparency, length(taxaNames))
		}else{
			if(length(transparency) == length(taxaNames)){
				stop("If transparency is not length 1, it must be the same length as number of tips")
				}
			}
		#
		#print(taxaColorNew)
		#print(which(!is.na(taxaColorNew)))
		#
		for(i in which(!is.na(taxaColorNew))){
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
	newColor <- t(grDevices::col2rgb(colorName))
	#
	newColor <- grDevices::rgb(
		newColor,
		alpha = alpha,
		maxColorValue = 255
		)
	#print(newColor)
	#
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


getCoordsPhyloPic <- function(
		xx,
		yy,
		sizeScale,
		plotAspRatio,
		picAspRatio,
		orientation
		){
	############################################
	# modify offset and size adjustment based on orientation
	###
	if(orientation == "rightwards"){
		# adjustment of sizeScale
			# need to modify sizeScale relative to aspect ratio
		if(picAspRatio < 1){
			# its flatish so correct it by aspect ratio
			picSize <- sizeScale * (picAspRatio^0.7)
		}else{
			# then its skinny, not flat
			# don't do anything
			picSize <- sizeScale 
			}
		###########################################
		# GET THE COORDINATES
		#
		# offset is sizeScale/2 by default
		offset <- sizeScale * 0.9 * plotAspRatio
		#
		x<-xx+offset
		y<-yy
		#adjust the position of the sides for the image
		#
		xAdj <- (picSize/2) * (plotAspRatio/picAspRatio) 
		yAdj <- picSize/2
		}
	###############
	###
	if(orientation == "upwards"){
		# adjustment of sizeScale
			# need to modify sizeScale relative to aspect ratio
		if(picAspRatio > 100){
			# its skinny so correct it by aspect ratio
			picSize <- sizeScale * (picAspRatio^1)
		}else{
			# then its flattish, don't do anything
			picSize <- sizeScale 
			}
		###########################################
		# GET THE COORDINATES
		#
		# offset is sizeScale/2 by default
		offset <- sizeScale * 0.5 * 1/plotAspRatio
		#
		x<-xx
		y<-yy+offset
		#adjust the position of the sides for the image
		#
		xAdj <- picSize/2 * (plotAspRatio/picAspRatio) 
		yAdj <- (picSize/2) 
		}
	######
	#		
	finalCoords <- list(
		xleft = x - xAdj ,
		ybottom = y - yAdj ,
		xright = x + xAdj ,
		ytop = y + yAdj
		)
	return(finalCoords)
	}	

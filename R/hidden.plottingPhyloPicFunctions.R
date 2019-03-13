# hidden functions related to plotting individual phylopics
	
plotSinglePhyloPic <- function(
		picPNG,
		whichTip,
		lastPP,
		sizeScale,
		offsetPic,
		orientation,
		plotAspRatio,
		taxonColor 
		){
		#
	#########################################
	#get aspect ratio
		# ratio of # of pixel dimensions
	picAspRatio <- dim(picPNG)[1]/dim(picPNG)[2]
	#
	#points(lastPP$xx,lastPP$yy)	
	#
	xyCoords <- getCoordsPhyloPic(
		xx = lastPP$xx[whichTip],
		yy = lastPP$yy[whichTip],
		sizeScale = sizeScale,
		offsetPic = offsetPic,
		plotAspRatio = plotAspRatio,
		picAspRatio = picAspRatio,
		orientation = orientation			
		)
	#
	##################################################
	# plot the picPNG using graphics::rasterImage
	picPNG_raster <- grDevices::as.raster(picPNG)
	# want color?
		# replaces all black tiles with the color in taxonColor
	if(!is.na(taxonColor)){
		#print(taxonColor)
		picPNG_raster[which(picPNG_raster=="#000000FF")] <- taxonColor
		#graphics::rasterImage(picPNG_raster [,,4])
		#stop()
		}
	# now plot the phylopic
	graphics::rasterImage(
		image = picPNG_raster,	
		xleft = xyCoords$xleft,
		ybottom = xyCoords$ybottom,
		xright = xyCoords$xright,
		ytop = xyCoords$ytop,
		interpolate = TRUE
		)
	# cool	
	}		

getCoordsPhyloPic <- function(
		xx,
		yy,
		sizeScale,
		offsetPic,
		plotAspRatio,
		picAspRatio,
		orientation
		){
	#############
	# modify offsetPic and size adjustment based on orientation
	#
	if(orientation == "rightwards"){
	############################################
	# adjustment of sizeScale
		# need to modify sizeScale relative to aspect ratio
	if(picAspRatio < 1){
		# its flatish so correct it by aspect ratio
		picSize <- sizeScale * (picAspRatio)
	}else{
		# then its skinny, not flat
		# don't do anything
		picSize <- sizeScale 
		}
		#
		###########################################
		# GET THE COORDINATES
		# add offsetPic calculated based on the plot limits
		x <- xx + offsetPic
		y <- yy
		#
		#adjust the position of the sides for the image
		xAdj <- (picSize/2) * (plotAspRatio/picAspRatio) 
		yAdj <- (picSize/2)
		}
	###############
	###
	if(orientation == "upwards"){
		############################################
	# adjustment of sizeScale
		# need to modify sizeScale relative to aspect ratio
	if(picAspRatio < 1){
		# its flatish so correct it by aspect ratio
		picSize <- sizeScale 
	}else{
		# then its skinny, not flat
		# don't do anything
		picSize <- sizeScale / (picAspRatio)
		}
		#
		# GET THE COORDINATES
		# add offsetPic calculated based on the plot limits
		x <- xx
		y <- yy + offsetPic
		#
		#adjust the position of the sides for the image
		xAdj <- (picSize/2) #* (plotAspRatio/picAspRatio) 
		yAdj <- (picSize/2) / (plotAspRatio/picAspRatio) 
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

	
prepPhyloPic<-function(
		picPNG, 
		noiseThreshold = 0.1,
		rescalePNG = TRUE, 
		trimPNG = TRUE,
		colorGradient = NULL,
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
	if(colorGradient == "trueMonochrome"){
		picPNG <- picPNG^0.001
		picPNG[picPNG < 0.2] <- 0
		picPNG[picPNG >= 0.2] <- 1
		}
	if(colorGradient == "increaseDisparity"){
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




getSortedMinMaxStratRanges <- function(timeTree,rangesFourDate ){
	#rangesFourDate <- timeTree$tipTaxonFourDateRanges	
	#print(rangesFourDate)
	#########################
	# get the ranges
	rangesMinMax <- rangesFourDate[
		,c("firstapp_max_ma","lastapp_min_ma")
		]
	rangesMinMax <- as.matrix(rangesMinMax)	
	rownames(rangesMinMax) <- rownames(rangesFourDate)
	# transform rangesMinMax for
	   # same backwards timescale as tree
	#rangesMinMax <- timeTree$root.time - rangesMinMax
	rangesMinMaxSorted <- matrix(,Ntip(timeTree),2)
	for(i in 1:Ntip(timeTree)){
		whichRanges <- rownames(rangesMinMax) == (timeTree$tip.label[i])
		#print(rangesMinMaxSorted);print(rangesMinMax[whichRanges,])
		rangesMinMaxSorted[i,] <- rangesMinMax[whichRanges,]
		}
	colnames(rangesMinMaxSorted) <- c("firstMax","lastMin")
	rownames(rangesMinMaxSorted) <- timeTree$tip.label
	return(rangesMinMaxSorted)
	}
		
		
plotTaxonStratDurations <- function(
		rangesMinMax,
		orientation,
		XX, YY,
		boxWidth = 0.7,
		boxCol 
		){
	########################
	#
	boxCol[is.na(boxCol)] <- "#000000FF"
	#
	if(orientation == "rightwards"){	
		for(i in 1:nrow(rangesMinMax)){
			yCent <- YY[i]
			ageMax <- rangesMinMax[i,"firstMax"]
			ageMin <- rangesMinMax[i,"lastMin"]
			#print(c(ageMax,ageMin,yCent))
			#points(ageMax,yCent,cex=2)
			#
			graphics::rect(	
				xleft = ageMax, 
				xright = ageMin,
				ytop = yCent + (boxWidth/2),
				ybottom = yCent - (boxWidth/2),
				col = boxCol[i],
				border = boxCol[i]
				)
			# get new xx and yy for end of ranges
			XX[i] <- ageMin
			YY[i] <- yCent 
			}	
		}
	if(orientation == "upwards"){
		for(i in 1:nrow(rangesMinMax)){
			xCent <- XX[i]
			ageMax <- rangesMinMax[i,"firstMax"]
			ageMin <- rangesMinMax[i,"lastMin"]
			graphics::rect(	
				ytop = ageMin, 
				ybottom = ageMax,
				xleft = xCent - (boxWidth/2),
				xright = xCent + (boxWidth/2),
				col = boxCol[i],
				border = boxCol[i]
				)
			# get new xx and yy for end of ranges
			XX[i] <- xCent
			YY[i] <- ageMin
			}		
		}		
	# return newXX,YY
	newXY <- list(XX = XX, YY = YY)
	return(newXY)
	}

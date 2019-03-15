
getSortedMinMaxStratRanges <- function(timeTree,
		rangesFourDate = timeTree$tipTaxonFourDateRanges
		){
	#########################
	rangesFourDate <- timeTree$tipTaxonFourDateRanges	
	# get the ranges
	rangesMinMax <- rangesFourDate[
		,c("firstapp_max_ma","lastapp_min_ma")
		]
	colnames(rangesMinMax) <- c("firstMax","lastMin")
	# transform rangesMinMax for
	   # same backwards timescale as tree
	rangesMinMax <- timeTree$root.time - rangesMinMax
	rangesMinMaxSorted <- 
	for(i in 1:Ntip(timeTree)){
		whichRanges <- rownames(rangesMinMax) == timeTree$tip.label[i]
		rangesMinMaxSorted[i,] <- rangesMinMax[whichRanges,]
		}
	return(rangesMinMaxSorted)
	}
		
		
plotTaxonStratDurations <- function(
		rangesMinMax,
		orientation,
		XX <- lastPP$xx,
		YY <- lastPP$yy,
		boxWidth = 0.7,
		boxCol = "black"
		){
	########################
	#
	if(orientation == "rightwards"){
		for(i in 1:nrow(rangesMinMax)){
			yCent <- YY[i]
			ageMax <- rangesMinMax[i,"firstMax"]
			ageMin <- rangesMinMax[i,"lastMin"]
			rect(	
				xleft = ageMax, 
				xright = ageMin,
				ytop = yCent + (boxWidth/2),
				ybottom = yCent - (boxWidth/2),
				col = boxCol
				)
			# get new xx and yy for end of ranges
			XX[i] <- ageMin
			YY[i] <- yCent 
			}	
		}
	if(orientation == "upwards"){
		for(i in 1:nrow(rangesMinMax)){
			XCent <- XX[i]
			ageMax <- rangesMinMax[i,"firstMax"]
			ageMin <- rangesMinMax[i,"lastMin"]
			rect(	
				ytop = ageMin, 
				ybottom = ageMax,
				xleft = xCent - (boxWidth/2),
				xright = xCent + (boxWidth/2),
				col = boxCol
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

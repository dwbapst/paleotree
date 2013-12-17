plotMultiDiv<-function(results,plotLogRich=FALSE,timelims=NULL,plotMultCurves=FALSE,
		multRainbow=TRUE,divPalette=NULL){
	#plots the median diversity curve for a multiDiv() result
	int.start<-results[[1]][,1]
	int.end<-results[[1]][,2]
	times1<-c(int.start,(int.end+((int.start-int.end)/100)))
	if(plotMultCurves){
		divs<-results[[2]]	#here's my div information
		divs1<-rbind(divs,divs)[order(times1),]
		times1<-sort(times1)
		#set up the general plotting window
		if(plotLogRich){
			y_lim<-c(min(divs1[divs1>=1]),max(divs1[divs1>=1]))
			plot(times1[divs1[,1]>0],divs1[divs1[,1]>0,1],type="n",ylim=y_lim,log="y",
					xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Log Lineage/Taxic Richness",
				main=paste("Multiple Diversity Curves"))
		}else{
			y_lim<-c(min(divs1),max(divs1))
			plot(times1,divs1[,1],type="n",ylim=y_lim,
					xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Lineage/Taxic Richness",
				main=paste("Multiple Diversity Curves"))
			}
		if(is.null(divPalette)){
			if(multRainbow){divPalette<-sample(rainbow(ncol(divs1)))
				}else{divPalette<-rep(1,ncol(divs1))}
			}
		for(i in 1:ncol(divs1)){	#plot each line
			lines(times1,divs1[,i],lwd=3,col=divPalette[i])
			}
	}else{
		mdiv<-results[[3]]
		mdiv1<-rbind(mdiv,mdiv)[order(times1),]
		times1<-sort(times1)
		if(plotLogRich){
			mdiv1[mdiv1[,2]<1,2]<-1;mdiv1[mdiv1[,3]<1,3]<-1
			y_lim<-c(min(mdiv1[mdiv1>=1]),max(mdiv1[mdiv1>=1]))
			plot(times1[mdiv1[,3]>0],mdiv1[mdiv1[,3]>0,3],type="n",ylim=y_lim,log="y",
					xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Log Lineage/Taxic Richness",
				main=paste("Median Diversity Curve"))
		}else{
			y_lim<-c(min(mdiv1),max(mdiv1))
			plot(times1,mdiv1[,3],type="n",ylim=y_lim,
					xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Lineage/Taxic Richness",
				main=paste("Median Diversity Curve"))
			}
		polygon(c(times1,rev(times1)),c(mdiv1[,3],rev(mdiv1[,2])),col="gray",border=NA)
		lines(times1,mdiv1[,1],lwd=3)
		}
	}
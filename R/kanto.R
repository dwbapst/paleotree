#' Example Species Abundances Tables
#'
#' A totally fictional example of species abundance data, for testing functions
#' that require a site-by-taxon table of community ecology data.

#' @name kanto
#' @rdname kanto
#' @aliases kanto

#' @details
#' A classic dataset of ecological data collected by Satoshi and Okido, consisting of
#' individual counts for 54 terrestrial faunal and floral species,
#' fron 23 sites across the mainland Kanto region.
#'
#' Different ontogenetic stages were compounded and recorded by the common name for the
#' first ontogenetic stage, with some inconsistency for species whose earliest stage have
#' only been recently recognized. When separate names are commonly applied to sexual
#' dimorphic forms, these were also combined and a single common name was used.
#'
#' \emph{Note: This data is totally made-up, and a satirical homage to
#' a well-known video game series, and thus should constitute fair-use.}

#' @format 
#' A table of type integer, representing terrestrial fauna and flora abundance counts.

#' @source 
#' Pokemon And All Respective Names are Trademark and Copyright of Nintendo 1996-2015.

#' @keywords datasets

#' @docType data

#' @examples
#'
#' data(kanto)
#'
#' #visualize site abundances as barplots
#' barplotAbund<-function(x){
#' 	x<-x[,colSums(x)>0]
#' 	layout(1:(nrow(x)+1))
#' 	xpar<-par(mar=c(0,7,2,0))
#' 	for(i in 1:(nrow(x)-1)){
#' 		barplot(x[i,],ylab=rownames(x)[i],
#' 			names.arg="")
#' 		}
#' 	barplot(x[nrow(x),],
#' 		ylab=rownames(x)[nrow(x)],las=3)
#' 	par(xpar)
#' 	layout(1)
#' 	mtext("Abundances",side=2,line=3,adj=0.8)
#' 	}
#' 
#' #first five sites
#' kanto5<-kanto[1:5,]
#' barplotAbund(kanto5)
#'
#' donttest{
#'
#' require(vegan)
#' bcDist<-vegdist(abundances,method="bray")
#'
#'
#' }



#'
NULL
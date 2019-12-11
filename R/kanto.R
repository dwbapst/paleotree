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
#' from 23 sites across the mainland Kanto region.
#' 
#' Different ontogenetic stages were compounded and recorded by the common name for the
#' first ontogenetic stage, with some inconsistency for species whose earliest stage have
#' only been recently recognized. When separate names are commonly applied to sexual
#' dimorphic forms, these were also combined and a single common name was used.
#' 
#' \emph{Note: This data is a totally made-up, satirical homage to
#' a well-known video game series (thus constituting fair-use).}

#' @format 
#' A table of type integer, representing terrestrial fauna and flora abundance counts.

#' @source 
#' Pokemon And All Respective Names are Trademark and Copyright of Nintendo 1996-2015.

#' @seealso
#' \code{\link{twoWayEcologyCluster}}, \code{\link{communityEcology}}

#' @keywords datasets

#' @docType data

#' @examples
#' 
#' data(kanto)
#' 
#' #visualize site abundances as barplots
#' barplotAbund <- function(x){
#' 	x <- x[,colSums(x)>0]
#' 	layout(1:(nrow(x)+1))
#' 	xpar <- par(mar = c(0,7,2,0))
#' 	for(i in 1:(nrow(x)-1)){
#' 		barplot(x[i,],ylab = rownames(x)[i],
#' 			names.arg = "")
#' 		}
#' 	barplot(x[nrow(x),],
#' 		ylab = rownames(x)[nrow(x)],las = 3)
#' 	par(xpar)
#' 	layout(1)
#' 	mtext("Abundances",side = 2,line = 3,adj = 0.8)
#' 	}
#' 
#' #first five sites
#' kanto5 <- kanto[1:5,]
#' barplotAbund(kanto5)
#' 
#' #get pairwise Spearman rho coefficients
#' rhoCoeff <- pairwiseSpearmanRho(kanto,dropAbsent = "bothAbsent")
#' 
#' #what are the nearest-neighbor rhos (largest rho correlations)?
#' diag(rhoCoeff) <- NA
#' rhoNearest <- apply(rhoCoeff,1,max,na.rm = TRUE)
#' rhoNearest
#' 
#' # We can see the power plant sample is extremely different from the rest
#' 
#' # measure evenness: Hurlbert's PIE
#' 
#' kantoPIE <- HurlbertPIE(kanto)
#' 
#' # compare to dominance (relative abundance of most abundant taxon)
#' 
#' dominance <- apply(kanto,1,function(x) max(x)/sum(x) )
#' 
#' plot(kantoPIE,dominance)
#' 
#' # relatively strong relationship!
#' 
#' 
#' \dontrun{
#' #########################################
#' #################################################
#' #########################################################
#' # Some Cool Ecology Stuff With Other Packages
#' 
#' # basically all the analyses & visualizations
#' 		#for ecology in R that I think are awesome
#' 
#' 
#' ##########################################
#' ###########################
#' #Ordination (PCO, DCA)
#' 
#' #get bray-curtis distances
#' library(vegan)
#' bcDist <- vegdist(kanto,method = "bray")
#' 
#' # do a PCO on the bray-curtis distances
#' pcoRes <- pcoa(bcDist,correction = "lingoes")
#' scores <- pcoRes$vectors
#' # plot the PCO
#' plot(scores,type = "n")
#' text(labels = rownames(kanto),scores[,1],scores[,2],cex = 0.5)
#' 
#' # the way the power plant and the pokemon tower converge
#' 	# is very suspicious: may be distortion due to a long gradient
#' 
#' # do a DCA instead with vegan's decorana
#' dcaRes <- decorana(kanto)
#' # plot using native vegan functions
#' 	   #will show species scores in red
#' plot(dcaRes,cex = 0.5)
#' #kind of messy
#' 
#' #show just the sites scores
#' plot(dcaRes,cex = 0.5,display = "sites")
#' 
#' #show just the species scores
#' plot(dcaRes,cex = 0.5,display = "species")
#' 
#' #well, that's pretty cool
#' 
#' #######################
#' #get the nearest neighbor for each site
#'     # based on pair-wise rho coefficients
#' rhoNeighbor <- apply(rhoCoeff,1,function(x)
#' 	    rownames(kanto)[tail(order(x,na.last = NA),1)])
#' 
#' #let's plot the nearest neighbor connections with igraph
#' NNtable <- cbind(rownames(kanto),rhoNeighbor)
#' 
#' # now plot with igraph
#' library(igraph)
#' NNlist <- graph.data.frame(NNtable)
#' plot(NNlist)
#' 
#' #arrows point at the nearest neighbor of each sample
#' 	    # based on maximum Spearman rho correlation
#' 
#' #########################################
#' #######################################################
#' # Two Way Cluster With Heatmap
#' 
#' # This example based on code provided by Max Christie
#' 
#' # load pheatmap library for this example
#' library(pheatmap) 
#' 
#' # get distance matrices for sites and taxa
#' 	# based on bray-curtis dist
#' 	# standardized to total abundance
#' 
#' # standardize site matrix to relative abundance
#' siteStand <- decostand(kanto, method = "total")
#' # site distance matrix (Bray-Curtis)
#' siteDist <- vegdist(siteStand, "bray", diag = TRUE)
#' 
#' # standardize taxa matrix to relative abundance
#' taxaStand <- decostand(t(kanto), method = "total")
#' # taxa distance matrix (Bray-Curtis)
#' taxaDist <- vegdist(taxaStand, "bray", diag = TRUE)
#' 
#' ### Need to set graphic parameters for table
#' 
#' # Check out range of values for relative abundance
#' # hist(myStand) # none get very high...
#' 
#' # number of breaks: number of colors for heatmap
#' nBreaks <- 15
#' 
#' # set underValue
#' 	# anything below this counts as not appearing
#' 	# at that site for visualization purposes
#' underValue <-  min(siteStand[siteStand>0])-min(siteStand[siteStand>0])/10
#' # set overValue (max relative abundance)
#' overValue <- max(siteStand)
#' # you can set your breaks to any sequence you want
#' 	# and they don't have to be the same length.  
#' 	# You can do this manually too.
#' # here we added a 0 to 'underValue' bin to 
#' 	# the heatmap, making this bin essentially 0.
#' colorBreaks <- c(0,seq(underValue,max(siteStand), 
#' 	by = max(siteStand)/(nBreaks-1)))
#' # here we used the function rainbow to create a vector of colors.  
#' 	# You can set these colors yourself too.  
#' # It is important that this vector is one element 
#' 	# less than the myBreaks vector
#' rainColors <- rainbow(nBreaks) 
#' # now we can add "white" onto the vector, 
#' 	# this will be the first color bin, 
#' 	# which we're going to set to be (essentially) 0.  
#' rainColors <- c("white", rainColors) 
#' # If you don't add white, taxa at 0 abundance get colored in
#' 
#' ### Plot the 2-Way Cluster
#' 
#' # heatmap, with user-set colors
#' # feed the function a distance matrix we wanted to use.  
#' 	#siteDist and taxaDist made above by vegdist (bray-curtis)
#' # scale is the relative abundance, let's label it as such
#' 
#' dev.new(width = 10)
#'
#' #for some reason, mtext() doesn't recognize pheatmap as plot.new
#' plot.new(width = 7) 
#'
#' pheatmap(
#'    siteStand, 
#' 	  clustering_method = "ward.D", 
#' 	  clustering_distance_rows = siteDist, 
#' 	  clustering_distance_cols = taxaDist,
#'    color = rainColors, 
#'    breaks = colorBreaks
#'    )
#' mtext("Relative Abundance",
#'    side = 4, line = -1.4, adj = 0.95)
#' 
#' # pretty cool looking!
#' 
#' ########################
#' # even better: 
#'    # twoWayEcologyCluster in paleotree
#' 
#' dev.new(width=10)	
#' 
#' twoWayEcologyCluster(
#'     xDist = siteDist,
#'     yDist = taxaDist,
#'     propAbund = siteStandKanto,
#'     cex.axisLabels = 0.8
#'     )
#' 
#' #########################################
#' #########################################################
#' ## Testing for differences between groups of sites
#' 
#' #is there a difference between routes and non-routes
#' groups <- rep(0, nrow(kanto))
#' groups[grep(rownames(kanto), pattern = "Route")] <- 1
#' 
#' #anosim (in vegan)
#' 	#are distances within groups smaller than distances between?
#' library(vegan)
#' anosim(dat = kanto, grouping = groups)
#' 
#' # we could also use PERMANOVA instead
#'     # this is generally considered more robust than ANOSIM
#  # use function 'adonis' from vegan
#'       # note that group needs to be factor for PERMANOVA
#' groupsAsFactor <- factor(groups)
#' adonis(kanto ~ groupsAsFactor)
#' 
#' # both analyses are very significant
#' 
#' ####################################################################
#' # SIMPER analysis (SIMalarity PERcentages) in Vegan
#' # which taxa contribute most to the difference between groups?
#'      # this might be 'index' taxa for different communities
#' # beware: it might also be the taxa that vary most within groups
#' 
#' simperResult <- simper(comm = kanto, group = groupsAsFactor)
#' simperResult
#' 
#' # these are the species that account for at least 70% of
#' # differences between groups, based on Bray-Curtis distances
#' 
#' # can see % contribtion for all species with summary()
#'     # as well as more detail in general...
#' summary(simperResult)
#' 
#' # other analyses to look into:
#'    # SimProf to test clusters from a cluster analysis...
#' 
#' #########################################################
#' # alternative for differentiating groups:
#'    # using multivariate GLMs in mvabund
#' 
#' library(mvabund)
#' 
#' ft <- manyglm(formula = kanto ~ groupsAsFactor)
#' anova(ft)
#' 
#' # also highly significant!
#' # note that this method though uses absolute abundances
#' # it will not accepted
#'     # which are usually impossible to get 
#' 
#' }
#' 




#' 
NULL
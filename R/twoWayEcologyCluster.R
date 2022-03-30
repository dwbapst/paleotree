#' R-Mode vs Q-Mode Two-Way Cluster Analyses and Abundance Plot for Community Ecology Data
#' 
#' This mode plots both R-mode (across sites) and Q-mode (across taxa) dendrograms
#' for a community ecology dataset, with branches aligned with a grid of dots
#' representing the relative abundance of taxa at each site in the dataset.

#' @details
#' You might be able to apply this to datasets that aren't community ecology datasets
#' of proportional abundance, but I can't guarantee or even predict what will happen.
#' 

#' @param xDist The pair-wise distance matrix for the cluster diagram drawn along the
#' horizontal axis of the graphic. Should be a distance matrix, or a matrix that can
#' be coerced to a distance matrix, for the same number of units as rows in \code{propAbund}.

#' @param yDist The pair-wise distance matrix for the cluster diagram drawn along the
#' vertical axis of the graphic. Should be a distance matrix, or a matrix that can
#' be coerced to a distance matrix, for the same number of units as columns in \code{propAbund}.

#' @param propAbund A matrix of abundance data, preferably relative
#' abundance scaled as proportions of the total number of individuals
#' at each site. This data determines the size scale of the taxon/site dots.

#' @param clustMethod The agglomerative clustering method used, as with
#' argument \code{method} with function \code{hclust}. 
#' \code{clustMethod} must be one of \code{"average"} 
#' (the default method for this function, 
#' also known as average-linkage or as UPGMA), 
#' \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"},
#' \code{"mcquitty"} (also known as WPGMA), \code{"median"}
#' (also known as WPGMC) or \code{"centroid"} (also known as UPGMC).

#' @param marginBetween Argument controlling space placed between
#' the cluster diagrams and the abundance plot. Default is 0.1.
 
#' @param abundExpansion An argument that is a multiplier controlling the size
#' of dots plotted for reflecting relative abundance.

#' @param cex.axisLabels Character expansion parameter for controlling the plotting
#' of axis labels on the abundance dot-grid only.

#' @param trimChar How many characters should the axis labels be trimmed to?
#' Default is 5, which means only the first five letters of each taxon/site
#' label will be shown on the dot-abundance plot.

#' @param xAxisLabel The label placed on the horizontal axis of the plot.

#' @param yAxisLabel The label placed on the vertical axis of the plot.

#' @return
#' This function creates a plot, and returns nothing, not even invisible output.

#' @seealso
#' Several other functions for community ecology data in paleotree are described
#' at the \code{\link{communityEcology}} help file. 
#' Also see the example dataset, \code{\link{kanto}}.

#' @author David W. Bapst

#' @references
#' The function here was designed to emulate previous published 'two-way'
#' cluster diagrams, particularly the one in Miller, 1988:
#' 
#' Miller, A. I. 1988. Spatial Resolution in Subfossil Molluscan
#' Remains: Implications for Paleobiological Analyses. \emph{Paleobiology} 14(1):91-103.
#' 

#' @examples
#' set.seed(1)
#' 
#' # generate random community ecology data
#'     # using a Poisson distribution
#' data<-matrix(rpois(5*7,1),5,7)
#' 
#' # get relative abundance, distance matrices
#' propAbundMat<-t(apply(data,1,function(x) x/sum(x)))
#' rownames(propAbundMat)<-paste0("site ", 1:nrow(propAbundMat))
#' colnames(propAbundMat)<-paste0("taxon ", 1:ncol(propAbundMat))
#' 
#' # for simplicity, let's calculate
#'     # the pairwise square chord distance
#'     # between sites and taxa
#' 
#' squareChordDist<-function(mat){
#'     res<-apply(mat,1,function(x)
#'         apply(mat,1,function(y)
#'             sum((sqrt(x)-sqrt(y))^2)
#'             )
#'         )
#'     #
#'     res<-as.dist(res)
#'     return(res)
#'     }
#' 
#' # its not a very popular distance metric
#'     # but it will do
#'     # quite popular in palynology
#' 
#' siteDist<-squareChordDist(propAbundMat)
#' taxaDist<-squareChordDist(t(propAbundMat))
#' 
#' dev.new(width=10)    
#' 
#' twoWayEcologyCluster(
#'     xDist = siteDist, 
#'     yDist = taxaDist,
#'     propAbund = propAbundMat
#'     )
#' 
#' \dontrun{
#' 
#' # now let's try an example with the example kanto dataset
#' # and use bray-curtis distance from vegan
#' 
#' library(vegan)
#' 
#' data(kanto)
#' 
#' # get distance matrices for sites and taxa
#'     # based on bray-curtis dist
#'     # standardized to total abundance
#' 
#' # standardize site matrix to relative abundance
#' siteStandKanto <- decostand(kanto, method = "total")
#' 
#' # calculate site distance matrix (Bray-Curtis)
#' siteDistKanto <- vegdist(siteStandKanto, "bray")
#' 
#' # calculate taxa distance matrix (Bray-Curtis)
#'     # from transposed standardized site matrix 
#' taxaDistKanto <- vegdist(t(siteStandKanto), "bray")
#' 
#' dev.new(width=10)    
#' 
#' twoWayEcologyCluster(
#'     xDist = siteDistKanto,
#'     yDist = taxaDistKanto,
#'     propAbund = siteStandKanto,
#'     cex.axisLabels = 0.8
#'     )
#' 
#' }



#' @name twoWayEcologyCluster
#' @rdname twoWayEcologyCluster
#' @export
twoWayEcologyCluster<-function(
        xDist, yDist, propAbund,
        # method
        # the agglomeration method to be used. This should be (an unambiguous abbreviation
        # of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
        # "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
        clustMethod = "average",
        #
        # space to place between cluster diagrams and abundance plot
        marginBetween = 0.1,
        # multiplier applies for plotting dots reflecting abundance
        abundExpansion = 3, 
        # axis label cex
        cex.axisLabels = 1,
        # how many characters to trim the axis labels to
        trimChar = 5,
        xAxisLabel = "Across Sites",
        yAxisLabel = "Across Taxa"
        ){
    #
    #
    xClust<-hclust(xDist, method=clustMethod)
    yClust<-hclust(yDist, method=clustMethod)
    #
    # convert cluster diagrams to phylo format
        # and double edge lengths
    xPhylo <- as.phylo(xClust)
    xPhylo$edge.length <- xPhylo$edge.length*2
    yPhylo <- as.phylo(yClust)
    yPhylo$edge.length <- yPhylo$edge.length*2
    # 
    # reorder propAbund
    propAbund <- propAbund[xClust$order, yClust$order]
    #
    # rescale propAbund
    propScale <- (propAbund/max(propAbund))*abundExpansion
    # alter names of propAbund
    #rownames(propAbund)<-paste0(
    #            substr(rownames(propAbund),
    #                start=1,stop=5)
    #            ,"...")
    rownames(propAbund) <- strtrim(rownames(propAbund), width = trimChar)
    colnames(propAbund) <- strtrim(colnames(propAbund), width = trimChar)
    #
    ####################################
    ## PLOTTING IT
    ##
    # make panes
    oldPar<-par(no.readonly = TRUE)
    layout(mat=rbind(3:4,1:2),
        heights=c(1,5),
        widths=c(5,1)
        )
    #
    # (lower left)
    ## major abund dot plot
    #
    par(
        mar=c(5,5,
        marginBetween,
        marginBetween)
        )
    iSet<-1:nrow(propAbund)
    jSet<-1:ncol(propAbund)
    #
    plot( 
        x=c(1,nrow(propAbund)),
        y=c(1,ncol(propAbund)),
        main="",
        xlab=xAxisLabel,
        ylab=yAxisLabel,
        type="n",
        axes=FALSE)
    #
    for(i in iSet){
        for(j in jSet){
            points(
                i,j,
                pch=16,col="black",
                cex=propScale[i,j]
                )
            }
        }
    #
    axis(
        side=1,
        at=iSet,
        las=2,
        lwd.ticks=0,
        #padj=-1,
        mgp=c(3,0.1,0),
        labels = rownames(propAbund),
        cex.axis = cex.axisLabels
        )
    axis(
        side=2,
        at=jSet,
        las=2,
        lwd.ticks=0,
        mgp=c(3,0.1,0),
        labels = colnames(propAbund),
        cex.axis = cex.axisLabels
        )
    #
    #################
    # lower right
    par(mar=c(5,0,marginBetween,1))
    plot(
        yPhylo,
        direction = "leftwards",
        show.tip.label=FALSE
        )
    axis(
        side=1,
        at=axTicks(side=1)[-1]
        )
    mtext(
        "Height",
        cex=0.8,
        side=1,
        line=2.5
        )
    #
    #################
    
    # upper left
    par(
        mar=c(0,5,1,marginBetween)
        )
    plot(xPhylo,
        direction="downwards",    
        show.tip.label=FALSE
        )
    axis(side=2,
        at=axTicks(side=2)[-1])
    mtext("Height",cex=0.8,
        side=2,line=2.5)
    #
    ################################
    par(mar=c(0.5, 0.5, 0.5, 0.5))
    plot(1, 1, type="n", axes=FALSE)
    legendScale <- min(propAbund[propAbund!=0])
    legendScale <- c(legendScale,
                        ((max(propAbund)-legendScale)/2) + legendScale,
                        max(propAbund)
                        )
    legend(x="center",
        legend = round(legendScale,2),
        pt.cex = abundExpansion * (legendScale/max(legendScale)),
        pch = 16)
    #
    layout(1)
    suppressWarnings(par(oldPar))
    #
    }
    

#' Cladogram and Range Data for the Retiolitinae
#' 
#' The majority rule consensus cladogram for 22 genera from the Retiolitinae, a
#' clade of Silurian retiolitids, along with discrete time interval data, all
#' taken from the same publication (Bates et al., 2005).
#' 
#' @details Interval dates were taken from Sadler et al. (2009). These zones were not a
#' 1-1 match to those in Bates et al., so it took some merging and splitting by
#' the package author, so buyer beware.
#' 
#' @name retiolitinae

#' @rdname retiolitinae

#' @aliases retiolitinae retioRanges retioTree

#' @docType data

#' @format This dataset is composed of two objects, \code{retioTree} (an ape 'phylo'
#' object containing the consensus cladogram) and \code{retioRanges}, a list
#' containing two matrices. The first matrix describes the first and last
#' interval times for 20 Silurian graptolite zones and the second matrix
#' describes when the various genera on the cladogram first and last appear in
#' those graptolite zones. (In other words, \code{retioRanges} has the 'timeList'
#' format called by some paleotree functions).

#' @source 

#' Source for cladogram and zonal ranges for genera: 
#'
#' Bates, D. E. B., A. Kozlowska, and A. C. Lenz. 2005. Silurian retiolitid graptolites:
#' Morphology and evolution. \emph{Acta Palaeontologica Polonica} 50(4):705-720.
#' 
#' Source for interval dates for graptolite zones: 
#'
#' Sadler, P. M., R. A. Cooper, and M. Melchin. 2009. High-resolution, early Paleozoic (Ordovician-Silurian)
#' time scales. \emph{Geological Society of America Bulletin} 121(5-6):887-906.

#' @seealso For more example graptolite datasets, see \code{\link{graptDisparity}}

#' @keywords datasets

#' @examples
#' 
#' #load data
#' data(retiolitinae)
#' 
#' #Can plot discrete time interval diversity curve with retioRanges
#' taxicDivDisc(retioRanges)
#' 
#' #Can plot the unscaled cladogram
#' plot(retioTree)
#' 
#' #Use basic time-scaling (terminal branches only go to FADs)
#' ttree<-bin_timePaleoPhy(tree=retioTree,timeList=retioRanges,type="basic",
#' 	ntrees=1,plot=TRUE)
#' 
#' #Note that this function creates stochastic time-scaled trees...
#' 	#A sample of 1 is not representative!
#' 
#' #phylogenetic diversity curve
#' phyloDiv(ttree)
#' 
#' 
NULL
#' Scale Tree to Unit-Length
#' 
#' Rescales all edges of a phylogeny to be equal to a single unit (1, or "unit-length").
#' 
#' @details Probably not a good way to scale a tree for comparative studies.
#' What does it mean to scale every edge of the phylogeny to the same length?
#'
#' This is not a rhetorical question. First, consider that on a 'reconstructed' tree with only
#' extant taxa, it would mean assuming the time between births of new
#' lineages that survive to the modern is extremely constant over evolutionary
#' history (because the unit-length wouldn't change, unlike the birth-death model,
#' which assumes lineages that survive to the modern 
#' accumulate at an accelerating exponential rate, even with constant birth and death rates).
#' 
#' A paleontological tree (say, under the Fossilized Birth-Death Model) treated with 
#' this 'unit-length' approach would assuming constancy and rigid homogeneity of 
#' the timing between the birth (origination events) of new lineages that 
#' (a) survive to the modern day, or 
#' (b) are sampled at some future point in the fossil record. 
#' We should assume even with constant extinction and fossilization rates
#' that such lineages should occur more frequently as we approach the
#' present-day.
#' 
#' Note that in neither of those cases, the 'unit-length' branch-scaling approach
#' does not produce trees whose edge lengths somehow represent 
#' the 'speciational' model, where evolutionary change is 
#' entirely 'cladogenetic' (ala punctuated equilibrium) 
#' and associated only with branching events. This would only be true
#' on the true, perfectly sampled tree, which isn't what anyone has.
#' 
#' Thus, overall, the value of the 'unit-length' approach is rather questionable.

#' @param tree A phylogeny as an object of class \code{phylo}.
#' 
#' @return Returns the modified phylogeny as an object of class \code{phylo}. Any
#' \code{$root.time} element is removed.

#' @seealso As an alternative to using \code{unitLengthTree} in comparative studies,
#' see \code{\link{timePaleoPhy}}. Or nearly anything, really...
#' 
#' See also \code{speciationalTree} in the package geiger, which does
#' essentially the same thing as \code{unitLengthTree}.

#' @examples
#' 
#' set.seed(444)
#' tree <- rtree(10)
#' 
#' layout(1:2)
#' plot(tree)
#' plot(unitLengthTree(tree))
#' layout(1)
#' 
#' @export unitLengthTree
unitLengthTree <- function(tree){
	if(!inherits(tree, "phylo")){
		stop("tree is not of class phylo")
		}
	tree$edge.length <- rep(1,Nedge(tree))
	tree$root.time <- NULL
	return(tree)
	}

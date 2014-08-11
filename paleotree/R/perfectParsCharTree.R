#' Simulate a Set of Parsimony-Informative Characters for a Phylogeny
#'
#' Creates a simulated set of parsimony-informative characters for a given rooted phylogeny,
#' with characters shared out equally across nodes in the phylogeny, with any remaining characters
#' assigned randomly to nodes.

#' @details
#' This function takes some given number of characters code and places them 
#'
#' This function assumes, like almost every function in paleotree, that the tree given is rooted, even if the
#' most basal node is a polytomy.

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso 

#' @author David W. Bapst 

#' @references

#' @examples
#' 

#' @export perfectParsCharTree
perfectParsCharTree<-function(tree,nchar){
	#simulate a perfect character dataset (parsimony informative binary chars) for a given tree
	charMat<-matrix(0,Ntip(tree),nchar)
	rownames(charMat)<-tree$tip.label
	desc<-sapply(prop.part(tree),function(x) tree$tip.label[x])
	desc<-desc[sapply(desc,length)!=Ntip(tree)]	        #get rid of root node that contains all taxa (not pars informative!
	nnode<-length(desc)
	if(nchar>nnode){
		#repeat desc if nchar multiple of nnode
		if((nchar %/% nnode) >1){
			for(i in 1:((nchar %/% nnode)-1) ){
				desc[(length(desc)+1):(length(desc)+nnode)]<-desc[1:nnode]
				}
			}
		if((nchar%%length(desc)) != 0){
			nrand<-nchar-length(desc)
			desc[(length(desc)+1):nchar]<-desc[sample(1:nnode,nrand,replace=TRUE)]
			message(paste("Randomly sampling nodes for",nrand,"extra characters"))
			}
	}else{
		if(nnode>nchar){stop("nchar needs to be larger than the number of nodes")}
		}
	for(i in 1:nchar){
		charMat[desc[[i]],i]<-1
		}		
	return(charMat)
	}
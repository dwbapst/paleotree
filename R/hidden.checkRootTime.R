# not exported

checkRootTime <- function(tree, stopIfFail = FALSE){	
	# check that the tree and its root age makes sense
	if(exists(tree$root.time) & exists(tree$edge.length)){
		if(tree$root.time < max(node.depth.edgelength(tree))){
			if(stopIfFail){
				stop(paste0("Total tree depth (",
					max(node.depth.edgelength(tree)),
					") is greater than root age(",
					tree$root.time,
					"), so tips are in the future",
					"\nSomething has probably gone very wrong"))
				}else{
				warning(paste0("Total tree depth (",
					max(node.depth.edgelength(tree)),
					") is greater than root age(",
					tree$root.time,
					"), so tips are in the future",
					"\nSomething has probably gone very wrong"))
				}
			res <- FALSE
		}else{
			res <- TRUE
			}
		}
	}
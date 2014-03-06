perfectCharTree<-function(tree,nchar){
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
			desc[(length(desc)+1):nchar]<-desc[sample(1:nnode,nchar-length(desc),replace=TRUE)]
			message("Randomly sampling nodes for extra character changes")
			}
	}else{
		if(nnode>nchar){stop("nchar needs to be larger than the number of nodes")}
		}
	for(i in 1:nchar){
		charMat[desc[[i]],i]<-1
		}		
	return(charMat)
	}
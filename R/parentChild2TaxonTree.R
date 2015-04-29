#' Create a Taxonomy-Based Phylogeny ('Taxon Tree') from a Table of Parent-Child Taxon Relationships
#'
#' This function takes a two-column matrix of taxon names,
#' indicating a set of binary parent-taxon:child-taxon 
#' paired relationships with a common root, and returns
#' a phylogeny object of class 'phylo'.

#' @details
#' All taxa listed must be traceble via their parent-child relationships to a single,
#' common ancestor which will act as the root node for output phylogeny.

#' @param parentChild A two-column matrix of type \code{character} where
#' each element is a taxon name. Each row represents a parent-child relationship
#' with first the parent (column 1) taxon name and then the child (column 2).

#' @param tipSet This argument controls which taxa are selected as tip taxa for the
#' output tree. The default \code{tipSet="nonParents"} selects all child taxa which
#' are not listed as parents in \code{parentChild}. Alternatively, \code{tipSet="all"}
#' will add a tip to every internal node with the parent-taxon name encapsulated in
#' parentheses.

#' @inheritParams makePBDBtaxontree

#' @return
#' A phylogeny of class 'phylo', with tip taxa as controlled by argument \code{tipSet}.
#' The output tree is returned with no edge lengths.

#' @seealso \code{\link{makePBDBtaxontree}}, \code{\link{taxonTable2TaxonTree}}

#' @author David W. Bapst

#' @examples
#' 
#' #let's create a small, really cheesy example
#' pokexample<-rbind(cbind("Squirtadae",c("Squirtle","Blastoise","Wartortle")),
#' 	c("Shelloidea","Lapras"),c("Shelloidea","Squirtadae"),
#' 	c("Pokezooa","Shelloidea"),c("Pokezooa","Parasect"),
#' 	c("Rodentapokemorpha","Linoone"),c("Rodentapokemorpha","Sandshrew"),
#' 	c("Rodentapokemorpha","Pikachu"),c("Hirsutamona","Ursaring"),
#' 	c("Hirsutamona","Rodentapokemorpha"),c("Pokezooa","Hirsutamona"))
#' 
#' #Default: tipSet='nonParents'
#' pokeTree<-parentChild2TaxonTree(pokexample, tipSet="notParents")
#' plot(pokeTree);nodelabels(pokeTree$node.label)
#'
#' #Get ALL taxa as tips with tipSet='all'
#' pokeTree<-parentChild2TaxonTree(pokexample, tipSet="all")
#' plot(pokeTree);nodelabels(pokeTree$node.label)
#'
#' 
#' \dontrun{
#' 
#' # let's try a dataset where not all the taxon relationships lead to a common root
#' 
#' pokexample_bad<-rbind(cbind("Squirtadae",c("Squirtle","Blastoise","Wartortle")),
#' 	c("Shelloidea","Lapras"),c("Shelloidea","Squirtadae"),
#' 	c("Pokezooa","Shelloidea"),c("Pokezooa","Parasect"),
#' 	c("Rodentapokemorpha","Linoone"),c("Rodentapokemorpha","Sandshrew"),
#' 	c("Rodentapokemorpha","Pikachu"),c("Hirsutamona","Ursaring"),
#' 	c("Hirsutamona","Rodentapokemorpha"),c("Pokezooa","Hirsutamona"),
#' 	c("Umbrarcheota","Gengar"))
#' 
#' #this should return an error, as Gengar doesn't share common root
#' pokeTree<-parentChild2TaxonTree(pokexample_bad)
#' 
#' }
#' 

#' @name parentChild2TaxonTree
#' @rdname parentChild2TaxonTree
#' @export
parentChild2TaxonTree<-function(parentChild,tipSet="nonParents",cleanTree=TRUE){
	#takes a two column matrix of character class taxon names
		#each row is a relationship: parent, then child
	#CHECKS
	if(length(tipSet)!=1 | !is.character(tipSet)){stop("tipSet must be a single character element")}
	if(!is.character(parentChild)){
		message("parentChild isn't of class character, attempting to coerce")
		parentChild<-as.character(parentChild)}
	if(length(dim(parentChild))!=2){
		stop("parentChild must be a matrix of class character with two columns and multiple rows")
	}else{
		if(!(dim(parentChild)[2]==2 & dim(parentChild)[1]>1)){
			stop("parentChild must be a matrix of class character with two columns and multiple rows")
			}
		}
	#first, get nodeNames, with root name first
	nodeNames<-unique(parentChild[,1])
	whichRoot<-which(sapply(nodeNames,function(x) !any(x==parentChild[,2])))
	#check that there isn't more than one root
	if(length(whichRoot)>1){stop(paste("Not all taxable are traceable to a single common root \n",
		length(whichRoot),"possible roots found:",paste0(nodeNames[whichRoot],collapse=", ")))}
	#now resort nodeNames
	nodeNames<-c(nodeNames[whichRoot],nodeNames[-whichRoot])
	if(tipSet!="nonParents"){
		if(tipSet=="all"){
			parentChild<-rbind(parentChild,cbind(nodeNames,paste0("(",nodeNames,")")))
		}else{stop("tipSet must be one of either 'nonParents' or 'all'")}}
	#identify tip taxa, this will be all taxa who are not-parents
	notParents<-sapply(parentChild[,2],function(x) !any(x==parentChild[,1]))
	tipNames<-parentChild[notParents,2]
	#now convert parentChild matrix to edge matrix
	edgeMat<-matrix(,nrow(parentChild),ncol(parentChild))
	taxonNames<-c(tipNames,nodeNames)
	#convert internal nodes to Ntip+nodeNames ID
	edgeMat[,1]<-sapply(parentChild[,1],function(x) which(x==taxonNames))
	edgeMat[,2]<-sapply(parentChild[,2],function(x) which(x==taxonNames))
	#reorder edge
	edge<-edgeMat[order(edgeMat[,1],edgeMat[,2]),]
	#make the tree
	tree<-list(edge=edge,tip.label=tipNames,edge.length=NULL, #edge.length=rep(1,nrow(edge))
		Nnode=length(nodeNames),node.label=nodeNames)
	class(tree)<-"phylo"
	if(cleanTree){ #make it a good tree
		#collapse singles
		tree<-collapse.singles(tree)
		#check it
		if(!testEdgeMat(tree)){stop("Edge matrix has inconsistencies")}
		tree<-reorder(tree,"cladewise") 	#REORDER IT
		tree<-read.tree(text=write.tree(tree))
		if(Ntip(tree)!=length(tipNames)){stop("Taxa number changed while cleaning tree")}
		tree<-ladderize(tree)
		}
	#plot(tree);nodelabels(tree$node.label)
	return(tree)
	}
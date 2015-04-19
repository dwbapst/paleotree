#' Estimating the Minimum Number of Character Transitions Using Maximum Parsimony
#'
#' \code{minCharChange} is a function which takes a cladogram and a discrete trait and finds the
#' solutions of inferred character states for ancestral nodes that minimizes the number of
#' character state transitions (either gains or losses/reversals) for a given topology and a set of
#' discrete character data. \code{minCharChange} relies on \code{ancPropStateMat}, which is a wrapper
#' for \code{phangorn}'s function \code{ancestral.pars}.

#' @param trait A vector of trait values for a discrete character, preferably named with taxon names
#' identical to the tip labels on the input tree.

#' @param tree A cladogram of type 'phylo'. Any branch lengths are ignored.

#' @param randomMax The maximum number of cladograms examined when searching a large number of solutions
#' consistent with the reconstructed ancestral states from \code{ancestral.pars} with the minimum number
#' of character state transitions. If the number of potential solutions is less than randomMax, then
#' solutions are exhaustively searched.

#' @param maxParsimony If maxParsimony is TRUE (the default) then only solutions with the smallest
#' number of total transitions examined will be returned. Note that since solutions are stochastically
#' 'guessed' at, and the number of possible solutions may not be exhaustively searched, there may have
#' been solutions not examined with a lower number of transitions even if maxParsimony=TRUE. Regardless,
#' one may want to do maxParsimony=FALSE if one is interested in whether there are solutions with a
#' smaller number of gains or losses and thus wants to return all solutions.

#' @param printMinResult If TRUE (the default), a summary of the results is printed to the terminal. The
#' information in this summary may be more detailed if the results of the analysis are simpler (i.e. 
#' fewer unique solutions).

#' @param type The parsimony algorithm applied by \code{ancestral.pars}, which can apply one of two:
#' "MPR" (the default) is a relatively fast algorithm developed by Hamazawa et al. (1995) and Narushima
#' and Hanazawa (1997), which relies on reconstructing the states at each internal node by re-rooting at
#' that node.  "ACCTRAN", the 'accelerated transitions' algorithm (Swofford and Maddison, 1987), favors
#' character reversal over independent gains when there is ambiguity. The "ACCTRAN" option in
#' ancestral.pars avoids repeated rerooting of the tree to search for a smaller set of maximum-parsimony
#' solutions that satisfy the ACCTRAN algorithm, but does so by assigning edge weights.
#' As of phangorn v1.99-12, both of these algorithms apply
#' the Sankoff parsimony algorithm, which allows multifurcations (polytomies).
 
#' @param cost A matrix of the cost to change between states of the input character trait. If NULL (the
#' default), the character is assumed to be unordered with equal cost to change from any state to another.
#' Cost matrices only impact the "MPR" algorithm; if a cost matrix is given but 'type = "ACCTRAN"', an error
#' is issued.

#' @details
#' The wrapper function \code{ancPropStateMat} simply automates the application of functions
#' \code{ancestral.pars} and \code{phyDat} from \code{phangorn}, along with several additional checks
#' and code to present the result as a matrix, rather than a specialized list. 
#' 
#' Note that although the default \code{cost} argument assumes that multistate characters are unordered,
#' the results of character change will be reported as gains and losses relative to the numbering of the
#' states in the output \code{transitionSumChanges}, exactly as if they had been ordered. In the case
#' where the character is actually ordered, this may be
#' considered a conservative approach, as using a parsimony algorithm for unordered character states allows fewer
#' gains or losses to be counted on branches where multiple gains and losses are reported. If the character is
#' presumably unordered, however, then the gains and losses division is arbitrary nonsense and should be combined to
#' to obtain the total number of character changes.

#' @return
#' A list is invisibly returned containing the following elements:
#'
#' \describe{

#'  \item{\code{message}}{Describes the performance of \code{minCharChange} at searching for a minimum solution.}

#' \item{\code{sumTransitions}}{A vector recording the total number of necessary transitions (sum total of gains
#' and losses/reversal) for each solution; effectively the parsimony cost of each solution.}

#' \item{\code{minTransitions}}{A symmetrical matrix with number of rows and columns equal to the number of
#' character states, with values in each cell indicating the minimum number of transitions from one ancestral
#' state (i.e. the rows) to a descendant state (i.e. the columns), taken across the set of kept solutions
#' (dependent on which are kept as decided by argument \code{maxParsimony}).  Generally guaranteed not to
#' add up to the number of edges contained within the input tree, and thus may not represent any realistic
#' evolutionary scenario but does represent a conservative approach for asking 'what is the smallest possible
#' number of transitions from 0 to 1' or 'smallest possible number of transitions from 1 to 0', independently
#' of each other.}

#' \item{\code{solutionArray}}{A three-dimensional array, where for each solution, we have a matrix with edges
#' for rows and two columns indicating the ancestral and child nodes of that edge, with values indicating the
#' states inferred for those nodes in a particular solution.}

#' \item{\code{transitionArray}}{A labelled three-dimensional array where for each solution we have a symmetrical
#' matrix with number of rows and columns equal to the number of character states, with values in each cell
#' indicating the total number of transitions from one ancestral state (i.e. the rows) to a descendant state
#' (i.e. the columns).}

#' \item{\code{transitionSumChanges}}{Which is a three column matrix with a row for every solution, with the
#' values in the three columns measuring the number of edges (branches) inferred to respectively have gains,
#' no change or losses (i.e. reversals), as calculated relative to the order of character states.}

#' }  

#' @aliases ancPropStateMat

#' @seealso
#' The functions described here are effectively wrapers of \code{phangorn}'s function
#' \code{\link{ancestral.pars}}.

#' @author David W. Bapst

#' @references
#' Hanazawa, M., H. Narushima, and N. Minaka. 1995. Generating most parsimonious reconstructions on
#' a tree: A generalization of the Farris-Swofford-Maddison method. Discrete Applied Mathematics
#' 56(2-3):245-265.
#' 
#' Narushima, H., and M. Hanazawa. 1997. A more efficient algorithm for MPR problems in phylogeny.
#' Discrete Applied Mathematics 80(2-3):231-238.
#'
#' Schliep, K. P. 2011. phangorn: phylogenetic analysis in R. \emph{Bioinformatics} 27(4):592-593.
#'
#' Swofford, D. L., and W. P. Maddison. 1987. Reconstructing ancestral character states under
#' Wagner parsimony. Mathematical Biosciences 87(2):199-229.


#' @examples
#' # let's write a quick & dirty ancestral trait plotting function
#' 
#' quickAncPlot<-function(tree,ancData,cex=cex){
#' 	plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="upwards")
#' 	tiplabels(pch=16,col=(as.numeric(char[tree$tip.label])+1))
#' 	nodelabels(pie=ancData[-(1:Ntip(tree)),],cex=cex,piecol=2:5)	
#' 	}
#'
#' # example with retiolitid graptolite data
#' 
#' data(retiolitinae)
#' 
#' ancMPR<-ancPropStateMat(retioTree,trait=retioChar[,2],type="MPR")
#' ancACCTRAN<-ancPropStateMat(retioTree,trait=retioChar[,2],type="ACCTRAN")
#' 
#' #let's compare MPR versus ACCTRAN results
#' layout(1:2)
#' quickAncPlot(tree,ancMPR,cex=0.3)
#' text(x=8,y=15,"type='MPR'",cex=1.5)
#' quickAncPlot(tree,ancACCTRAN,cex=0.3)
#' text(x=9,y=15,"type='ACCTRAN'",cex=1.5)
#' 
#' minCharChange(retioTree,trait=retioChar[,2],type="MPR")
#' minCharChange(retioTree,trait=retioChar[,2],type="ACCTRAN")
#' 
#' # with simulated data
#' 
#' set.seed(444)
#' tree<-rtree(50)
#' #simulate under a likelihood model
#' char<-rTraitDisc(tree,k=3,rate=0.7)
#' tree$edge.length<-NULL
#' tree<-ladderize(tree)
#' 
#' ancMPR<-ancPropStateMat(tree,trait=char,type="MPR")
#' ancACCTRAN<-ancPropStateMat(tree,trait=char,type="ACCTRAN")
#' 
#' #let's compare MPR versus ACCTRAN results
#' layout(1:2)
#' quickAncPlot(tree,ancMPR,cex=0.3)
#' text(x=8,y=15,"type='MPR'",cex=1.5)
#' quickAncPlot(tree,ancACCTRAN,cex=0.3)
#' text(x=9,y=15,"type='ACCTRAN'",cex=1.5)
#' #MPR has much more uncertainty in node estimates
#' 	#but that doesn't mean ACCTRAN is preferable
#' 
#' \donttest{
#' # what ancPropStateMat automates (with lots of checks):
#' char1<-matrix(char,,1)
#' rownames(char1)<-names(char)
#' #translate into something for phangorn to read
#' char1<-phyDat(char1,type="USER",levels=sort(unique(char1)))
#' x<-ancestral.pars(tree,char1,type="MPR")
#' y<-ancestral.pars(tree,char1,type="ACCTRAN")
#' }
#' 
#' #estimating minimum number of transitions with MPR 
#' minCharChange(tree,trait=char,type="MPR")
#'
#' #and now with ACCTRAN
#' minCharChange(tree,trait=char,type="ACCTRAN")

#' @name minCharChange
#' @rdname minCharChange
#' @export
minCharChange<-function(trait, tree, randomMax=10000, maxParsimony=TRUE,
		 type="MPR", cost=NULL, printMinResult=TRUE){
	#randomMax=100;maxParsimony=TRUE;printMinResult=TRUE;type="MPR";cost=NULL
	#print result gives back a reasonable 
	ancMat<-ancPropStateMat(trait, tree, type=type, cost=cost)
	#num of potential solutions
	taxSol<-apply(ancMat,1,function(x) sum(x>0))	#taxSol = solution length of each taxon
	nSol<-prod(taxSol)
	#supposedly charN (my trait vector to be sampled) can be character, its fine
	charN<-colnames(ancMat)
	if(nSol>randomMax){
		solMat<-t(apply(ancMat,1,function(x) sample(charN[x>0],randomMax,replace=T)))
	}else{	
		#exhaustive search needed
		#first, build matrix of non-changing taxa
		noChange<-which	(taxSol==1)
		solMat<-matrix(sapply(noChange,function(x) charN[ancMat[x,]>0]),,1)
		rownames(solMat)<-noChange
		if(nSol>1){
			for(i in 2:max(taxSol)){
				changers<-which(taxSol==i)
				for(j in changers){
					solMat2<-lapply(charN[ancMat[j,]>0],function(x) rbind(solMat,x))
					solMat1<-solMat2[[1]]
					for(k in 2:length(solMat2)){solMat1<-cbind(solMat1,solMat2[[k]])}
					colnames(solMat1)<-NULL
					rownames(solMat1)<-c(rownames(solMat),j)
					solMat<-solMat1
					}
				}
			}
		solMat<-solMat[order(as.numeric(rownames(solMat))),,drop=FALSE]
		}
	#are all solMats unique? (yes, if TRUE)
	solUnq<-all(!sapply(1:ncol(solMat),function(x) 
		any(sapply((1:ncol(solMat))[-x],function(y) identical(solMat[,x],solMat[,y])))))
	#do I need to stop if not all solutions are unique???
	if(!solUnq){
		if(nSol>randomMax){
			#if random, then okay, I guess you might have non unique solutions
			solDup<-c(FALSE,sapply(2:ncol(solMat),function(x) 
				any(sapply((1:ncol(solMat))[1:(x-1)],function(y) identical(solMat[,x],solMat[,y])))))
			solMat<-solMat[,!solDup]
		}else{
			#if not random, then stop cause something is wrong!
			stop("Error: Not all solutions are unique, as calculated, despite random permutations not used. Please investigate or contact Dave Bapst.")
			}
		}
	#edgeSol is a 3D array, where for each solution, we have a matrix with edges for
		# rows and two columns indicating the ancestral node of that edge
		# and the child node of that edge, with values indicating the states
		# inferred for those nodes in a particular solution
	edgeSol<-array(,dim=c(Nedge(tree),2,ncol(solMat)))
	for(i in 1:ncol(solMat)){
		xSol<-solMat[,i]
		#rearrange as an edge matrix of transitions
		edgeSol[,,i]<-cbind(xSol[sapply(tree$edge[,1],function(x) which(x==names(xSol)))],
			xSol[sapply(tree$edge[,2],function(x) which(x==names(xSol)))])
		}
	#tranMat is a 3D array where for each solution we have a symmetrical matrix
		#equal to the number of character states, with values indicating the total
		#number of transitions from one ancestral state (given as the rows) to
		#a descendant state (given as columns)
	tranMat<-array(,dim=c(length(charN),length(charN),ncol(solMat)))
	rownames(tranMat)<-paste("anc.",colnames(ancMat),sep="")
	colnames(tranMat)<-paste("desc.",colnames(ancMat),sep="")
	#sumTran is the parsimony cost: number of gains+losses
	sumTran<-numeric()	
	for(i in 1:ncol(solMat)){
		edgeTran<-edgeSol[,,i]
		#turn into transition matrix
		tranMat1<-t(sapply(charN,function(x) sapply(charN,function(y) 
			sum(edgeTran[,1]==x & edgeTran[,2]==y))))
		tranMat[,,i]<-tranMat1
		diag(tranMat1)<-0
		sumTran[i]<-sum(tranMat1)
		#rows are the ancestor state, columns are the desc state
		}
	#are all tranMats unique? generally not
	#	unqTran<-sapply(1:length(tranMat),function(x) 
	#		any(sapply((1:length(tranMat))[-x],function(y) identical(tranMat[,,x],tranMat[,,y]))))	
	#hist(sumTran)
	if(nSol==1){
		maxPars<-1
		tranMat<-tranMat[,,1,drop=FALSE]
		edgeSol<-edgeSol[,,1,drop=FALSE]
	}else{
		maxPars<-which(sumTran==min(sumTran))
		if(maxParsimony){
			#select only most parsimonious solutions
			solMat<-solMat[,maxPars,drop=FALSE]
			tranMat<-tranMat[,,maxPars,drop=FALSE]
			sumTran<-sumTran[maxPars]
			}
		}
	#get # of gains and # of losses and # of no-change for each transition matrix
	tranSumChange<-t(sapply(lapply(1:dim(tranMat)[3],function(y) tranMat[,,y]),function(x) 
		c(sum(x[upper.tri(x)]),sum(diag(x)),sum(x[lower.tri(x)]))))
	colnames(tranSumChange)<-c("Gains","NoChanges","Losses")
	#get the minimum solution
	minTran<-apply(tranMat,c(1,2),min)
	#
	funcMess<-c(paste(nSol,"potential solutions under",type,",",length(maxPars),"most parsimonious solutions found"),
		ifelse(nSol>randomMax,"Solutions sampled stochastically","Solutions exhaustively checked"))
	if(printMinResult){
		if(length(maxPars)<6){
			print(list(message=funcMess,sumTransitions=sumTran,
				transitionArray=tranMat,minTransitions=minTran))
		}else{
			print(list(message=funcMess,sumTransitions=sumTran,minTransitions=minTran))
			}
		}
	return(invisible(list(message=funcMess,sumTransitions=sumTran,minTransitions=minTran,
		solutionArray=edgeSol,transitionArray=tranMat,transitionSumChanges=tranSumChange))) #
	}

#' @rdname minCharChange
#' @export
ancPropStateMat<-function(trait, tree, type="MPR", cost=NULL){
	#wrapper for phangorn's ancestral.pars that returns a fully labeled matrix indicating
		#the relative frequency of a node being reconstructed under a given state
	#require(phangorn)
		#return error if cost is not null and type=ACCTRAN
	if(type=="ACCTRAN" & !is.null(cost)){
		stop("cost matrix is inapplicable if ACCTRAN algorithm is used")}
	#check names
	if(is.null(names(trait))){
		if(Ntip(tree)!=length(trait)){
			stop("names(trait) missing and length(trait) isn't same as number of tips in tree!")
		}else{
			names(trait)<-tree$tip.label
			message("names(trait) missing \n","trait values will be assigned to taxa exactly as in tree$tip.label")
			}
		}
	char1<-matrix(trait,,1)
	rownames(char1)<-names(trait)
	#translate into something for phangorn to read
	char1<-phyDat(char1,type="USER",levels=sort(unique(char1)))
	#get anc states
	anc1<-ancestral.pars(tree,char1,type=type,cost=cost)
	#turn into a col-per-state matrix with each row a node or tip, numbered as in edge
	anc2<-matrix(unlist(anc1),,length(unique(char1)),byrow=T)
	#based on conversation with Klaus on 04-17-15
		#will treat output as if it was always ordered exactly as tips and nodes
		#are numbered in $edge; should be as basic as numbering 1:nrow
	rownames(anc2)<-1:nrow(anc2)
	#does that make sense for the tree
	if(nrow(anc2) != (Nnode(tree)+Ntip(tree))){
		stop("ancestral state matrix has wrong number of rows??")}
	#and now name the columns by the levels
	colnames(anc2)<-attributes(anc1)$levels
	return(anc2)
	}
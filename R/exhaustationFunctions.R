#' Analyses of the exhaustion of Character States Over Evolutionary History
#' 
#' The following functions are for measuring and fitting various
#' distributions to the gradual exhaustion of unexpressed
#' character states, as originally described by Wagner (2000,
#' Evolution). 
#' 
#' \code{accioExhaustionCurve} uses a Sankoff parsimony ancestral-reconstruction
#' algorithm (written by PJ Wagner, \emph{not} the one from \code{phangorn} used
#' elsewhere in paleotree) to calculate character changes across each branch
#' (internode edge) of a tree, and then reports the counts of character state
#" changes, new state changes, etc. 
#' 
#' \code{accioBestAcquisitionModel} takes output from \code{accioExhaustionCurve},
#' calculates one of two character change indices, and then fits a series of user-selected
#' models to the dataset, returning details pertaining to the best-fit model.
#' 
#' \code{charExhaustPlot} is a wrapper for \code{accioBestAcquisitionModel} that
#' produces a plot of the observed character change data against the
#' expectation under the best-fit model.
#'

#' @details
#' The functions \code{accioBestAcquisitionModel}  and \code{charExhaustPlot} offer
#' users two different options for examining character change: \code{totalAcc}
#' fits models to the total accumulated number of state changes over the phylogeny,
#' thus using exhaustion to explor the size and distribution of character space. The
#' other option \code{charAlt} fits models to the number of character that alter from
#' primitive to derived over phylogeny, thus reflecting the size and distribution of state space.
#' 
#' \code{accioExhaustionCurve} can order its reconstruction of change by stratigraphic order of first appearances. It is unclear what this means.

#' @param charData A \code{data.frame} of morphological character
#' codings (a morphological 'matrix'), with taxon units
#' as rows and characters as columns.

#' @param charTypes A vector of length equal to
#' the number of chaacters in \code{charData}, with elements indicating
#' whether the corresponding character in \code{charData}
#' is \code{"unordered"} or \code{"ordered"}. However,
#' if \code{length(charTypes) = 1}, then it is repeated for all taxa.
#' The default value for this argument is \code{"unordered"}.

#' @param outgroup A string matching to one
#' of the tip labels as given by \code{tip.label}, 

#' @param firstAppearances A vector, with length equal to the
#' same number of taxa (rows) as in \code{charData}, in the same corresponding order.

#' @param missingCharValue The string value indicating a missing
#' character coding value, by default \code{"?"}.

#' @param inapplicableValue The string value indicating an
#' inapplicable character coding value, by default \code{"-"}.

# @param changes A vector counting character change; the extraction of this is automated if you use .

#' @param phyloTree A phylogenetic tree of class \code{phylo} as used by package \code{ape}.

#' @param models A vector of type \code{character} naming models to be fit.
#' Default is \code{c("exponential","gamma","lognormal","zipf")}.

#' @param exhaustion_info The list of results output
#' from function \code{accioExhaustionCurve}.

#' @param changesType A single character value, indicating the character change data
#' to be assessed from the result of the character
#' exhaustion analysis, must be one of either 'totalAcc' (to the total number of
#' accumulated character changes, ideal for modeling the size and distribution of
#' \emph{state} space) or 'charAlt' (to plot the total number of character alterations,
#' ideal for modeling the size and distribution of \emph{character} space).

#' @param xlab Label for the X axis; \code{"Total Characters"} by default.

#' @param ylab Label for the Y axis. If not provided by the user,
#' a label based on the \code{changesType} argument will be used.

#' @param main Main title label for the plot. If not provided by
#' the user, a label based on the \code{changesType} argument will be used.

#' @param xsize Parameter controling size of the axes,
#' which are forced to be symmetric.

#' @return
#' \code{accioExhaustionCurve} outputs a list containing two objects: first,
#' a matrix named \code{exhaustion} consisting of three columns: \code{"Steps"},
#' \code{"Novel_States"}, and \code{"Novel_Characters"}, respectively giving
#' the counts of these respective values for each branch (internode edge).
#' The second element of this list is named \code{State_Derivations} and is
#' a count of how often each state, across all characters, was derived relative
#' to the primitive position along each internode edge.
#' 
#' The output of \code{accioBestAcquisitionModel} is a list object containing
#' information on the best-fit model, the parameters of that model, the calculated
#' probabilition distribution for that model at the same intervals, for use in quantile plots.	
#' 
#' \code{charExhaustPlot} produces a plot, and outputs no data.
	


#' @note
#' This family of functios presented here were originally written 
#' by Peter J. Wagner, and then modified and adapted by David W.
#' Bapst for wider release  in a CRAN-distributed
#' package: \code{paleotree}. This makes the code presented here
#' a very different beast than typical paleotree code, for
#' example, there are fewer tests of correct input type, length, etc.

#' @seealso
#' Also see \code{paleotree} functions \code{\link{minCharChange}} and
#' \code{\link{ancPropStateMat}}, the latter of which is a wrapper
#' for \code{phangorn}'s function \code{ancestral.pars}.

#' @author 
#' Initially written by Peter J. Wagner, with modification
#' and documentation by David W. Bapst.

#' @references
#' Wagner, P. J. 2000. Exhaustion of morphologic character
#' states among fossil taxa. \emph{Evolution} 54(2):365-386.

#' @examples
#' 
#' # get data
#' data(SongZhangDicrano)
#' 
#' dicranoTree<-cladogramDicranoX13
#' 
#' # modify char data
#' charMat<-data.matrix(charMatDicrano)
#' charMat[is.na(charMatDicrano)]<-0
#' charMat<-(charMat-1)
#' colnames(charMat)<-NULL
#' # replace missing values
#' charMat[is.na(charMatDicrano)]<-"?"
#' 
#' # the 'outgroup' is Exigraptus, first taxon listed in the matrix
#' exhaustionResults <- accioExhaustionCurve(phyloTree=dicranoTree,
#'    charData=charMat, charTypes="unordered",
#'    outgroup="Exigraptus_uniformis")
#'
#' # fits models to exhaustion for total accumulation
#' accioBestAcquisitionModel(exhaustion_info=exhaustionResults,
#'  changesType="totalAcc", 	
#'  models=c("exponential","gamma","lognormal","zipf")) 
#' 
#' # plot of exhausation of total accumulation of character states
#' charExhaustPlot(exhaustion_info=exhaustionResults,
#' 	   changesType="totalAcc")
#' 
#' # plot of exhausation of character alterations
#' charExhaustPlot(exhaustion_info=exhaustionResults,
#' 	   changesType="charAlt")
#' 	



# examples trash

#	# first get count state derivations to count total accumulation
# changesDicrano <- sort(rowSums(exhaustionResults$State_Derivations),decreasing=TRUE) 

# # stratigraphic age data
#   # Not clear where this taken from (collected by PJW)
# strat_data<-cbind(
#     c(0, 4, 7, 7, 6, 6, 6, 2, 2, 1, 1, 6, 7),
#     c(0, 5, 7, 7, 7, 7, 7, 2, 2, 1, 1, 7, 7))
# # get just first appearances from strat data
# FAs <- strat_data[,1]
#  
# 	  firstAppearances=FAs, 

# # get character types
# nchar <- dim(charMat)[2]
# types <- rep("unordered",nchar)	# set character types to unordered (1)



#' @name exhaustionFunctions
#' @rdname exhaustionFunctions
#' @export
accioExhaustionCurve <- function(phyloTree,charData,
	charTypes="unordered",outgroup=NULL,firstAppearances=NULL,
	missingCharValue="?",inapplicableValue="-")	{
	#
	# sorry pete, but paleotree is camelCase... ;)
	#
	# rewrite ape's phylo format tree
	### vtree is a vector tree, where each cell gives the ancestral HTU of the OTU or HTU
	###		HTU's are nodes, with the basal HTU = nTips+1 (i.e., 13 for 12 taxa)
	ntaxa<- dim(charData)[1]
	if(ntaxa!=Ntip(phyloTree)){
		stop("Number of tips on phylotree doesn't match number of taxa in charData")
		}
	mtree<-phyloTree
	mtree$tip.label<-match(phyloTree$tip.label,rownames(charData))
	mtree <- accio_matrix_tree_from_ape_tree(mtree)	
	#
	# based on routine used in Wagner (2000).  Uses simple parsimony optimization
	#	and (if available) stratigraphic data to order branches and show how
	#	many novel states appear with how many changes along each branch working
	#	up the tree.
	#	MODIFY to show number of characters triggered as well as # of states
	nchar <- dim(charData)[2]
	# test charTypes
	if(nchar != length(charTypes)){
		if(length(charTypes)==1){
			charTypes<-rep(charTypes,nchar)
		}else{
			stop("length of charTypes must be 1 or equal to number of characters")
			}
		}
	# convert charData row names to numbers
	rownames(charData)<-1:nrow(charData)
	#
	#get outgroupID
	if(is.null(outgroup)){
		outgroupID<-1
	}else{
		if(length(outgroup)!=1){
			stop("if given, outgroup must be of length 1")
			}
		outgroupID<-which(phyloTree$tip.label==outgroup)
		if(length(outgroupID) != 1){
			stop("There is no single match to the outgroup given among tip labels")
			}
		}
	#
	char_history <- accio_minimum_steps_history(mtree=mtree,charData=charData,charTypes=charTypes,
		missingCharValue=missingCharValue,inapplicableValue=inapplicableValue,outgroupID=outgroupID)
	#
	if (is.null(firstAppearances))	{
		# order branches from base of the tree
		branch_order <- accio_patristic_distance_from_base(mtree)	
		# for ties, put nodes before otus
		used_order <- branch_order + (1-(1:length(branch_order))/100)
		# list branches in order	
		exhaust_order <- order(used_order,decreasing=FALSE)		
	}else{
		if(ntaxa!=length(firstAppearances)){
			stop("length of firstAppearances doesn't match number of taxa in charData")
			}
		exhaust_order <- used_order <- accio_branch_order_with_stratigraphy(mtree,firstAppearances)
		### insert ordering with stratigraphy here.	
		}
	#
	branches <- max(mtree)
	branch_derivs <- vector(length=branches)
	for (b in 1:branches){
		branch_derivs[b] <- sum(char_history$Changes[,1]==b)
		}
	# get relevant branches (non-zero length) in order of appearance
	rel_branches <- exhaust_order[branch_derivs[exhaust_order]>0]
	rb <- length(rel_branches)
	#
	char_deriv <- vector(length=nchar)
	exhaust_matrix <- matrix(0,max(char_history$Changes[,2]),max(char_history$Changes[,3]))
	found <- total_steps <- novel_chars <- novel_states <- 0
	for (b in 1:rb)	{
		br <- rel_branches[b]
		#br_ch <- sum(char_history$Changes[,1]==br)
		#if (br_ch>0)	{
		local_change <- subset(char_history$Changes,char_history$Changes[,1]==br)
		char_deriv[local_change[,2]] <- char_deriv[local_change[,2]]+1
		novel_chars <- novel_chars + sum(char_deriv[local_change[,2]]==1)
		for (d in 1:branch_derivs[br])	{
			if (branch_derivs[br]>1){
				ch <- local_change[d,2]
				st <- local_change[d,3]
			}else{
				ch <- local_change[2]
				st <- local_change[3]
				}
			exhaust_matrix[ch,st] <- exhaust_matrix[ch,st]+1
			if (exhaust_matrix[ch,st]==1)	novel_states <- novel_states+1
			}
		total_steps <- total_steps+branch_derivs[br]
		if (exhaust_order[b]<10){
			if (max(mtree)<10){
				brnm <- paste("br_",exhaust_order[b],sep="")
			} else{
				if (max(mtree<100))	{
					brnm <- paste("br_0",exhaust_order[b],sep="")	
				}else{
					brnm <- paste("br_00",exhaust_order[b],sep="")	
					}
				}
		}else{
			if (exhaust_order[b]<100){
				if (max(mtree)<100){
					brnm <- paste("br_",exhaust_order[b],sep="")	
				}else{
					brnm <- paste("br_0",exhaust_order[b],sep="")	
					}				
				}
			}
		if (b==1){
			exhaustion <- c(total_steps,novel_states,novel_chars)
			branch_names <- brnm
		}else{
			exhaustion <- rbind(exhaustion,c(total_steps,novel_states,novel_chars))
			branch_names <- c(branch_names,brnm)
			}
		}
	#
	colnames(exhaustion) <- c("Steps","Novel_States","Novel_Characters")
	rownames(exhaustion) <- branch_names
	exhaust_output <- list(exhaustion,exhaust_matrix)
	names(exhaust_output) <- c("Exhaustion","State_Derivations")
	return(exhaust_output)
	}

	
	
#' @rdname exhaustionFunctions
#' @export
accioBestAcquisitionModel <- function(exhaustion_info,changesType,
		models=c("exponential","gamma","lognormal","zipf"))	{
	# check
	if(all(changesType!=c("totalAcc","charAlt"))){
		stop("changesType must be one of either 'totalAcc' or 'charAlt'")
		}
	# get changes
	if(changesType=="totalAcc"){
		# get the best model for the size and distribution of character space
		changes <- sort(rowSums(exhaustion_info$State_Derivations),decreasing=TRUE)	#char_changes
		}
	if(changesType=="charAlt"){
		# get the best model for the size and distribution of state space
		changes <- sort(array(exhaustion_info$State_Derivations)[array(exhaustion_info$State_Derivations)>0],decreasing=TRUE)	#state_derivs
		}	
	# get the best model based on change/derivation distributions
	best_H <- best_uniform <- optimize_uniform_abundance(counts=changes)
	best_model <- "uniform"
	best_dist <- rep(1/as.numeric(best_H[1]),as.numeric(best_H[1]))
	if (!is.na(match("exponential",models)))	{
		best_exponential <- optimize_geometric_abundance(counts=changes)
		if (best_exponential[length(best_exponential)] < best_H[length(best_H)])	{
			best_H <- best_exponential
			best_model <- "exponential"
			best_dist <- geometric_distribution(decay=best_exponential[1])
			}
		}
	if (!is.na(match("gamma",models)))	{
		best_gamma1 <- optimize_gamma_one_abundance(counts=changes)
		if (best_gamma1[length(best_gamma1)] < best_H[length(best_H)])	{
			best_H <- best_gamma1
			best_model <- "gamma"
			best_dist <- gamma_distribution(a=best_gamma1[1],b=best_gamma1[1],S=best_gamma1[2])
			}
		}
	if (!is.na(match("lognormal",models)))	{
		best_lognormal <- optimize_lognormal_abundance(counts=changes)
		if (best_lognormal[length(best_lognormal)] < best_H[length(best_H)])	{
			best_H <- best_lognormal
			best_model <- "lognormal"
			best_dist <- lognormal_distribution(mag=best_lognormal[1],S=best_lognormal[2])
			}
		}
	if (!is.na(match("zipf",models)))	{
		best_zipf <- optimize_zipf_abundance(counts=changes)
		if (best_zipf[length(best_zipf)] < best_H[length(best_H)])	{
			best_H <- best_zipf
			best_model <- "zipf"
			best_dist <- zipf_distribution(zm=best_zipf[1],S=best_zipf[2])
			}
		}
	output <- list(best_model,best_H,best_dist)
	names(output) <- c("Best_Model","Best_Model_Parameters","Best_Distribution")
	return(output)
	}

#' @rdname exhaustionFunctions
#' @export
charExhaustPlot<-function(exhaustion_info,changesType,xlab="Total Characters",ylab=NULL,main=NULL,xsize=3){
	if(all(changesType!=c("totalAcc","charAlt"))){
		stop("changesType must be one of either 'totalAcc' or 'charAlt'")
		}
	exhaustion <- exhaustion_info$Exhaustion
	nstep <- max(exhaustion_info$Exhaustion[,1])
	best_exhaustion <- accioBestAcquisitionModel(exhaustion_info=exhaustion_info,
			changesType=changesType, models=c("exponential","gamma","lognormal","zipf"))
	# get the best model for the size and distribution of character/state space
	ba <- best_exhaustion$Best_Distribution
	after3 <- exhaustion[2,1]
	hba <- expected_exhaustion_curve(rel_ab_dist=ba, nstep, after3, length(ba))
	#
	mxx <- 10*ceiling(max(exhaustion[,1])/10)
	if (mxx<100)	{
		mxx <- 5*ceiling(mxx/5)
	}else{
		mxx <- 10*ceiling(mxx/10)
		}
	#
	# plot activation of characters
	#xsize <- 3
	par(pin=c(xsize,xsize))
	#
	if(changesType=="totalAcc"){
		if(is.null(ylab)){
			ordinate <- "Novel States"
		}else{
			ordinate<-ylab
			}	
		if(is.null(main)){
			main_title <- "State Accumulation"	
		}else{
			main_title <- main
			}	
		plot(NA,type='n',axes=FALSE,main=main_title,
			xlab=abcissa,ylab=ordinate,
			xlim=c(0,mxx),ylim=c(0,mxx))
		draw_symmetric_axes(mxx,xsize)
		#
		lines((1:nstep),hba,lwd=2,lty=3)
		points(exhaustion[,1],exhaustion[,3],pch=21,bg="green",cex=1.25)
		}
	#
	if(changesType=="charAlt"){
		if(is.null(ylab)){
			ordinate <- "Altered Characters"
		}else{
			ordinate<-ylab
			}	
		if(is.null(main)){
			main_title <- "Character Alterations"
		}else{
			main_title <- main
			}	
		plot(NA,type='n',axes=FALSE,main=main_title,
			xlab=abcissa,ylab=ordinate,
			xlim=c(0,mxx),ylim=c(0,mxx))
		draw_symmetric_axes(mxx,xsize)
		#
		lines((1:nstep),hba,lwd=2,lty=2)
		points(exhaustion[,1],exhaustion[,2],pch=21,bg="skyblue",cex=1.25)
		}
	#
	#layout(1)
	}
	
	


	
##############################################################


### ROUTINES TO RECONSTRUCT CHARACTER CHANGE

# get history of character evolution implied by minimum steps
accio_minimum_steps_history <- function(mtree,charData,charTypes,missingCharValue="?",
		inapplicableValue="-",outgroupID=outgroupID)	{
	# mtree: matrix tree, where each row gives the branches stemming from a node
	# charData: character data
	# charTypes: 1 for unordered, 0 for ordered
	# missingCharValue: numeric representation for "?"
	# inapplicableValue: numeric representation for inapplicable or gap
	# outgroupID: the tip ID # of the outgroupID taxon
	otus <- dim(charData)[1]
	Nnodes <- dim(mtree)[1]
	nchar <- dim(charData)[2]
	steps <- vector(length=nchar)
	for (c in 1:nchar)	{
		cvector <- charData[,c]
		type <- charTypes[c]
		char_evolution <- Sankoff_character(mtree=mtree,cvector=cvector,
			type=type,missingCharValue=missingCharValue,
			inapplicableValue=inapplicableValue,outgroupID=outgroupID)
		steps[c] <- char_evolution$Steps
		if (c==1)	{
			full_matrix <- char_evolution$States
			changes_matrix <- char_evolution$Derivation
			}	else	{
			full_matrix <- cbind(full_matrix,char_evolution$States)
			changes_matrix <- cbind(changes_matrix,char_evolution$Derivation)
			}
		}
	char_labels <- vector(length=nchar)
	for (c in 1:nchar)	{
		if (c<10)	{
			if (nchar<10)	{
				char_labels[c] <- paste("ch_",c,sep="")
				}	else if (nchar>9 && nchar<100)	{
				char_labels[c] <- paste("ch_0",c,sep="")
				}	else	{
				char_labels[c] <- paste("ch_00",c,sep="")
				}
			}	else if (c<100)	{
				if (nchar<100)	{
				char_labels[c] <- paste("ch_",c,sep="")
				}	else	{
				char_labels[c] <- paste("ch_0",c,sep="")
				}
			}
		}

	ttl_br <- otus+Nnodes
	nonzero <- 1
	rbr <- 1
	for (br in 1:ttl_br)	{
		nonzero <- sum(changes_matrix[br,]>0)
		if (nonzero>0)	{
			dch <- (1:nchar)[changes_matrix[br,]>0]
			dst <- changes_matrix[br,changes_matrix[br,]>0]
			dbr <- rep(br,length(dch))
			if (rbr==1)	{
				branch_changes <- cbind(dbr,dch,dst)
				} else	{
				branch_changes <- rbind(branch_changes,cbind(dbr,dch,dst))
				}
			rbr <- rbr+1
			}
		}
	rownames(branch_changes) <- rep("",dim(branch_changes)[1])
	colnames(full_matrix) <- colnames(changes_matrix) <- char_labels
	output <- list(steps,full_matrix[((otus+1):(otus+Nnodes)),],branch_changes)
	names(output) <- c("Steps","Ancestral_Reconstructions","Changes")
	return(output)
	}

# get Sankoff matrix ancestral reconstructions
Sankoff_character <- function(mtree,cvector,type=1,missingCharValue="?",inapplicableValue="-",outgroupID=outgroupID)	{
	# method for reconstructing ancestral conditions.
	obs_states <- sort(unique(cvector[cvector>=0]))
	ttl_states <- length(obs_states)
	Nnodes <- dim(mtree)[1]
	otus <- length(cvector)
	scored <- (1:otus)
	sankoff_matrix <- matrix(1,otus+Nnodes,(length=ttl_states))
	for (s in 1:otus)	{
		if (cvector[s] >= 0)	{
			st <- match(cvector[s],obs_states)
			sankoff_matrix[s,st] <- 0
			}	else	sankoff_matrix[s,] <- 0
		}
	cvector <- c(cvector,rep(missingCharValue,Nnodes))
	node_rich <- vector(length=Nnodes)
	for (n in Nnodes:1)	{
		missing <- gap <- 0
		node_rich[n] <- f1 <- length(mtree[n,mtree[n,]>0])
		# if all taxa are scored
		# mtree[n,]
	#	if (sum(mtree[n,(1:f1)] %in% scored)==f1)	{
		ht <- n+otus
		sankoff_matrix[ht,] <- 0
		for (s in 1:f1)	{
			#### add something to deal with inapplicables here.
			# list states that demand more than minimal change
			sp <- mtree[n,s]
			if (cvector[sp]!=inapplicableValue && cvector[sp]!=missingCharValue)	{
				#	we will add a step to each of these because if sankoff is:
				#		0 1 1 for states 0, 1 & 2
				#	then we need 1 step from either 1 or 2
				n_s <- obs_states[sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
				a_s <- obs_states[!sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
				no_step <- match(n_s,obs_states)
				add_step <- match(a_s,obs_states)
				sankoff_matrix[ht,no_step] <- sankoff_matrix[ht,no_step]+min(sankoff_matrix[sp,])
				if (length(no_step) < length(obs_states))	{
					sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+min(sankoff_matrix[sp,])+1
		#			sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+sankoff_matrix[sp,add_step]
		#			sankoff_matrix[ht,]
					}
				}	else if (cvector[sp]==inapplicableValue)	{
				gap <- gap+1
				}	else if (cvector[sp]==missingCharValue)	{
				missing <- missing+1
				}
			}
		if (missing==f1)	{
			cvector[ht] <- missingCharValue
			scored <- c(scored,ht)
			}	else if (gap==f1)	{
			cvector[ht] <- inapplicableValue
			scored <- c(scored,ht)
			}	else if (sum(sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,]))==1)	{
			# see if we can fix a state
			ast <- match(min(sankoff_matrix[ht,]),sankoff_matrix[ht,])
			cvector[ht] <- obs_states[ast]
			scored <- c(scored,ht)
			}	else	{
			#print(obs_states[sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,])])
			cvector[ht] <- ravel_polymorph(obs_states[sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,])])
			}
		}

	base <- otus+1
	ch_steps <- min(sankoff_matrix[base,])

	if (cvector[base]<0 && (cvector[base]!=missingCharValue && cvector[base]!=inapplicableValue))	{
		#poss_starts <- (1:ttl_states)[sankoff_matrix[base,] %in% min(sankoff_matrix[base,])]
		#print(cvector[base])
		poss_starts <- unravel_polymorph(cvector[base])
		# IF the designated outgroupID is attached to the first node AND if it 
		#	has one of the possible nodal states, then assign that to basal node
		if (!is.na(match(outgroupID,mtree[1,])) && !is.na(match(cvector[outgroupID],poss_starts)))	{
			cvector[base] <- cvector[outgroupID]
			}	else	{
			# otherwise, just grab one of them at random: it really doesn't matter
			grab <- ceiling(runif(1)/(1/length(poss_starts)))
			cvector[base] <- poss_starts[grab]
			}
		}

	### start here: work up the tree, using cvector[htu] to set the state
	changes_above <- vector(length=Nnodes)
	for (n in 1:Nnodes)	{
		ht <- n+otus
		if (cvector[ht]!=inapplicableValue && cvector[ht]!=missingCharValue)	{
			f1 <- node_rich[n]
			# if all taxa are scored in mtree[n,]
			if (sum(mtree[n,(1:f1)] %in% scored)<f1)	{
				anc_st <- cvector[ht]
				unscored <- mtree[n,!mtree[n,(1:f1)] %in% scored]
				for (u in 1:length(unscored))	{
					if (length(unscored)>1)	{
						f2 <- unscored[u]
						}	else	{
						f2 <- unscored
						}
					# if parent value matches one of the possible values for the daughter node
					#	then go with that value
					if(anc_st>=0 && (!is.na(match(anc_st,obs_states[sankoff_matrix[f2,] %in% min(sankoff_matrix[f2,])]))))	{
						cvector[f2] <- anc_st
						}	else if (anc_st<0 && (anc_st!=missingCharValue && anc_st!=inapplicableValue))	{
						anc_poss <- unravel_polymorph(anc_st)
						f2_poss <- unravel_polymorph(cvector[f2])
						reduced_poss <- f2_poss[f2_poss %in% anc_poss]
						if (length(reduced_poss)==1)	{
							cvector[ht] <- cvector[f2] <- reduced_poss
							}	else	{
							#cvector[ht] <- cvector[f2] <- ravel_polymorph(reduced_poss)
							# if it is still up in the air at this point, then you might as
							#	well roll the dice!
							cvector[f2] <- reduced_poss[ceiling(runif(1)/(1/length(reduced_poss)))]
							}
						}	else	{
					# it really does not matter at this point what state we pick: it will all be the same.
						f2_poss <- unravel_polymorph(cvector[f2])
						cvector[f2] <- f2_poss[ceiling(runif(1)/(1/length(f2_poss)))]
						}
						# end case where we can assign ancestral condition to daughter node
					}	# end search of unscored nodes
				}	# I could add a routine here to compare sets of equally parsimonious states
			changes_above[n] <- min(sankoff_matrix[n+otus,])
			}		# if node has (0,1) and daughter node has (1,2), then go with 1
		}
	cchanges <- rep(0,length(cvector))
	for (n in 1:Nnodes)	{
		ht <- otus+n
		if (cvector[ht]!=missingCharValue && cvector[ht]!=inapplicableValue)	{
			f1 <- node_rich[n]
			for (f in 1:f1)	{
				f2 <- mtree[n,f]
				if (cvector[f2] != cvector[ht])
					if ((cvector[f2]!=missingCharValue 
					&& cvector[f2]!=inapplicableValue) 
					&& (cvector[ht]!=missingCharValue 
					&& cvector[ht]!=inapplicableValue))
						# if both taxa are resolved, then this is easy
						if (length(cvector[f2])==1 && length(cvector[ht])==1)
							cchanges[f2] <- match(cvector[f2],obs_states)
				}
			}
		}
	output <- list(ch_steps,cvector,cchanges)
	names(output) <- c("Steps","States","Derivation")
	return(output)
	}


#############################################################

#### ROUTINES TO READ DATA
# read character matrix

accio_data_from_nexus_file <- function(nexus_file_name, polymorphs=TRUE,
		 missingCharValue="?", inapplicableValue="-",label_char_cols=FALSE)	{
	# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
	# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
	# missingCharValue: value substituting for "?"
	# inapplicableValue: value substituting for gap ("-")
	nexus <- scan(file=nexus_file_name,what=character(),sep="\n")
	gp_scr <- "-"
	unkn_scr <- "?"
	if (!is.na(match("BEGIN DATA;",nexus)))	{
		# for old MacClade files
		ln <- match("BEGIN DATA;",nexus)
	#	test <- "	DIMENSIONS  NTAX=17 NCHAR=25;"
	#	txt <- strsplit(test,split="",fixed=TRUE)[[1]]
		txt <- strsplit(nexus[ln],split="",fixed=TRUE)[[1]]
		dd <- (1:length(txt))[txt %in% "="]
		notu <- 0
		i <- dd[1]+1
		while (txt[i]!=" ")	{
			notu <- (10*notu)+as.numeric(txt[i])
			i <- i+1
			}
		nchars <- 0
		i <- dd[2]+1
		while (txt[i]!=";")	{
			nchars <- (10*nchars)+as.numeric(txt[i])
			i <- i+1
			}
		ln <- ln+1
		txt <- strsplit(nexus[ln],split="",fixed=TRUE)[[1]]
		while (txt[2]!="F" && txt[3]!="O" && txt[9]!="S" && txt[10]!="Y")	{
			ln <- ln+1
			txt <- strsplit(nexus[ln],split="",fixed=TRUE)[[1]]
			}
		dd <- (1:length(txt))[txt %in% "="]
		gp_scr <- txt[dd[3]+1]
		unkn_scr <- txt[dd[2]+1]
		} else	{
		# routine for Mesquite files
		ln <- taxa_line <- match("BEGIN TAXA;",nexus)
		found <- FALSE
		while(!found)	{
			ln <- ln+1
			txt <- strsplit(nexus[ln],split="",fixed=TRUE)[[1]]
			found <- txt[2]=="D" && txt[3]=="I" && txt[4]=="M"
			}
		notu <- 0
		for (i in (match("=",txt)+1):(length(txt)-1))	notu <- (10*notu)+as.numeric(txt[i])

		ln <- taxa_line <- match("BEGIN CHARACTERS;",nexus)
		found <- FALSE
		while(!found)	{
			ln <- ln+1
			txt <- strsplit(nexus[ln],split="",fixed=TRUE)[[1]]
			found <- txt[2]=="D" && txt[3]=="I" && txt[4]=="M"
			}
		nchars <- 0
		for (i in (match("=",txt)+1):(length(txt)-1))	nchars <- (10*nchars)+as.numeric(txt[i])
		found <- FALSE
		while(!found)	{
			ln <- ln+1
			txt <- strsplit(nexus[ln],split="",fixed=TRUE)[[1]]
			found <- txt[2]=="F" && txt[3]=="O" && txt[4]=="R"
			}
		dd <- (1:length(txt))[txt %in% "="]
		gap_scr <- txt[dd[2]+2]
		unkn_scr <- txt[dd[3]+2]
		}

	# find where matrix starts
	nexus_alt <- gsub("\t","",nexus)
	matrix_starts <- match("MATRIX",nexus_alt)
	ms <- matrix_starts+1
	## MacClade files will have "[" above the actual data
	while (strsplit(nexus_alt[ms],split="",fixed=TRUE)[1]=="[")	ms <- ms+1
	taxa <- vector(length=notu)
	chmatrix <- matrix(0,notu,nchars)
	for (n in 1:notu)	{
		ms <- matrix_starts+n
		j <- strsplit(nexus[ms],split="",fixed=TRUE)[[1]]
		while (j[1]=="\t")	j<- j[2:length(j)]
		name <- ""
		if (j[1]=="'")	{
			i <- 2
			while (j[i]!="'")	{
				name <- paste(name,j[i],sep="")
				i <- i+1
				}
			}	else	{
			i <- 1
			while (j[i]!=" " && j[i]!="\t")	{
				name <- paste(name,j[i],sep="")
				i <- i+1
				}
			}
		taxa[n] <- name
		i <- i+1
		while (j[i]==" " || j[i]=="\t")	i <- i+1
		for (ch in 1:nchars)	{
			if (j[i]==unkn_scr)	{
				chmatrix[n,ch] <- missingCharValue
				}	else if (j[i]==gap_scr)	{
				chmatrix[n,ch] <- inapplicableValue
				}	else if (j[i]=="(" || j[i]=="{" || j[i]=="[")	{
				# polymorphic!
				i <- i+1
				k <- 0
	#			while (j[k]!=")" && j[k]!="}" && j[k]!="]")	k <- k+1
				while (j[i]!=")" && j[i]!="}" && j[i]!="]")	{
					while (j[i]==" " || j[i]==",")	i <- i+1
					if (j[i]>9)	{
						chmatrix[n,ch] <- chmatrix[n,ch]-((10^k)*convert_letter_state_to_numeric(j[i]))
						} else	{
						chmatrix[n,ch] <- chmatrix[n,ch]-((10^k)*as.numeric(j[i]))
						}
					k <- k+1
					i <- i+1
					}
				}	else if (j[i]>9)	{
				chmatrix[n,ch] <-  convert_letter_state_to_numeric(j[i])
				}	else	{
				chmatrix[n,ch] <- as.numeric(j[i])
				}
			i<-i+1
			}
		}
	## now get the number of states per character
	nstates <- vector(length=nchars)
	for (ch in 1:nchars)	{
		states <- sort(unique(chmatrix[,ch]))
		states <- states[states!=missingCharValue]
		states <- states[states!=inapplicableValue]
		if (min(states)<0)	{
			polys <- states[states<0]
			for (p in 1:length(polys))	{
				states <- c(states,unravel_polymorph(polys[p]))
				}
			states <- sort(unique(states[states>=0]))
			}
		nstates[ch] <- length(states)
		}

	# look for outgroup designation	
	outgroup <- -1
	if (!is.na(match("BEGIN SETS;",nexus)))	{
		tx_pt <- 1+match("BEGIN SETS;",nexus)	# look at taxon partitions
		splitNexusTxPt<-strsplit(nexus[tx_pt],"")[[1]]
		# below line refered 'nexusfile[tx_pt,]' changed to splitNexusTxPt
		while (outgroup==-1 && splitNexusTxPt[1]=="\t" 
			&& splitNexusTxPt[2]=="T" && splitNexusTxPt[3]=="A" 
			&& splitNexusTxPt[4]=="X" && splitNexusTxPt[5]=="P" 
			&& splitNexusTxPt[6]=="A")	{
				j <- strsplit(nexus[tx_pt],split=" ",fixed=TRUE)[[1]]
				if (!is.na(match("Outgroup",j)))	{
					otg <- match("Outgroup",j)+1
					while (j[otg]==":" || j[otg]==" " || j[otg]==""){
						otg <- otg+1
						}
					out <- strsplit(j[otg],split="",fixed=TRUE)[[1]]
					a <- 1
					outgroup <- 0
					while (out[a]!="," && out[a]!=";")	{
						outgroup <- (10*outgroup)+as.numeric(out[a])
						a <- a+1
						}
					}
				tx_pt <- tx_pt+1
				splitNexusTxPt<-strsplit(nexus[tx_pt],"")[[1]]
				}
		}

	if (!is.na(match("CHARSTATELABELS ",nexus_alt)))	{
		char_names <- vector(length=nchars)
		state_names <- matrix("",nchars,max(nstates)+1)
		ln <- match("CHARSTATELABELS ",nexus_alt)
		nexus_alt[ln+1] <- gsub(" ; ","",nexus_alt[ln+1])
		char_info <- strsplit(nexus_alt[ln+1],split=", ")[[1]]
		for (ch in 1:nchars)	{
			txt <- strsplit(char_info[ch],split="",fixed=TRUE)[[1]]
			i <- 1
			while (txt[i]!=" ")	i <- i+1
			i <- i+1
			if (txt[i]=="'")	{
				j <- i+1
				while (txt[j]!="'")	j <- j+1
				char_names[ch] <- ""
				for (k in i:j)
					char_names[ch] <- paste(char_names[ch],txt[k],sep="")
				}	else	{
				j <- i+1
				while (txt[j] != " " && j<length(txt))	j <- j+1
				char_names[ch] <- ""
				if (txt[j]==" ")	{
					h<-j-1
					} else	h <- j
				for (k in i:h)
					char_names[ch] <- paste(char_names[ch],txt[k],sep="")
				if (j<length(txt))	{
					while (txt[j]==" " || txt[j]=="/")	j <- j+1
					statenames <- ""
					for (k in j:length(txt))
						statenames <- paste(statenames,txt[k],sep="")
					sn <- strsplit(statenames,split=" ")[[1]]
					state_names[ch,1:length(sn)] <- sn
					}
				}
			}
		} else	{
		char_names <- state_names <- FALSE
		}

		#x <- list(taxa,chmatrix,strat_ranges,geography)
		#return (list(taxa,chmatrix,strat_ranges,geography))
	rownames(chmatrix) <- taxa
	if (length(char_names)>1 && label_char_cols)	{
		colnames(chmatrix) <- char_names
		}

	if (length(char_names)>1 && outgroup!=-1)	{
		output <- (list(taxa=taxa,chmatrix=chmatrix,nstates=nstates,char_names=char_names,state_names=state_names,outgroup=outgroup))
		names(output)<- c("OTUs","Matrix","States","Characters_Labels","States_Labels","Outgroup")
		} else if (length(char_names)>1){
		output <- (list(taxa=taxa,chmatrix=chmatrix,nstates=nstates,char_names=char_names,state_names=state_names))
		names(output)<- c("OTUs","Matrix","States","Characters_Labels","States_Labels")
		} else if (outgroup!=-1)	{
		output <- (list(taxa=taxa,chmatrix=chmatrix,nstates=nstates,outgroup=outgroup))
		names(output)<- c("OTUs","Matrix","States","Outgroup")
		}	else	{
		output <- (list(taxa=taxa,chmatrix=chmatrix,nstates=nstates))
		names(output)<- c("OTUs","Matrix","States")
		}
	return(output)
	}

# if A is used for state 10, then this will convert A -> 10
convert_letter_state_to_numeric <- function(letterstate)  {
	if (letterstate=="A") {	state <- 10
		} else if (letterstate=="B") {	state <- 11
		} else if (letterstate=="C") {	state <- 12
		} else if (letterstate=="D") {	state <- 13
		} else if (letterstate=="E") {	state <- 14
		} else if (letterstate=="F") {	state <- 15
		} else if (letterstate=="G") {	state <- 16
		} else if (letterstate=="H") {	state <- 17
		} else if (letterstate=="J") {	state <- 18
		} else if (letterstate=="K") {	state <- 19
		} else if (letterstate=="L") {	state <- 20
		} else if (letterstate=="M") {	state <- 21
		} else if (letterstate=="N") {	state <- 22
		} else if (letterstate=="P") {	state <- 23
		} else if (letterstate=="Q") {	state <- 24
		} else if (letterstate=="R") {	state <- 25
		} else if (letterstate=="S") {	state <- 26
		} else if (letterstate=="T") {	state <- 27
		} else if (letterstate=="U") {	state <- 28
		} else if (letterstate=="V") {	state <- 29
		} else if (letterstate=="W") {	state <- 30
		} else if (letterstate=="X") {	state <- 31
		} else if (letterstate=="Y") {	state <- 32
		} else if (letterstate=="Z") {	state <- 33
		}
	return(state)
	}

## why no I or O? (PJW)
	# I have no idea, Pete (DWB)




#poly <- -321
unravel_polymorph <- function(poly)	{
	poly<-as.numeric(poly)
	combo <- -1*poly
	sts <- 1+floor(log10(abs(combo)))
	polymorphics <- vector(length=sts)
	#
	base <- 10^(sts-1)
	for (s in 1:sts)	{
		polymorphics[s] <- floor(abs(combo)/base)
		combo <- combo%%base
		base <- base/10
		}
	return (polymorphics)
	}

# routine to read a text file with a nexus tree
read_nexus_tree <- function(nexustreefile) {
	# suppose tree in text file looks like: "(1,(2,((3,4),(5,(6,(7,(8,(9,10))))))))"
	# this returns:
	#	[1] 11 12 14 14 15 16 17 18 19 19 -1 11 12 13 13 15 16 17 18	
	nexustree <- scan(file=nexustreefile,what=character(),sep="\n")
	nexus_string <- strsplit(nexustree,split="",fixed=TRUE)[[1]]
	nodes <- 0
	for (i in 1:length(nexus_string))		if (nexus_string[i]=="(")	nodes <- nodes+1
	# get clades
	clades <- vector(length=nodes)
	for (c in 1:nodes)	clades[c] <- c
	# get taxa
	notu <- p <- 0
	for (i in 1:length(nexus_string))	{
		if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
			otu <- as.numeric(nexus_string[i])+(otu * (10^p))
			p <- p+1
			if (otu>notu)	notu <- otu
			} else {
			p <- otu <- 0
			}
		}
	vector_tree <- vector(length=notu+max(clades))
	for (c in 1:nodes)	clades[c] <- -1
	cl <- c <- 0
	i <- 1
	for (i in 1:length(nexus_string))	{
		if (nexus_string[i]=="(")	{
			sp <- p <- 0
			cl <- cl+1
			if (cl>1)	{
				vector_tree[notu+cl] <- clades[c]+notu
				} else vector_tree[notu+1] <- -1
			c <- c+1
			clades[c] <- cl
			} else if (nexus_string[i]==")")	{
			c <- c-1
			sp <- p <- 0
			} else if (nexus_string[i]==",")	{
			sp <- p <- 0
			} else if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
			sp <- as.numeric(nexus_string[i])+(sp*(10^p))
			p <- p+1
			if (nexus_string[i+1]<"0" || nexus_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
			}
		}

	return(vector_tree)
	}

# routine to convert vector tree to a matrix tree
# 	vector tree gives ancestral node of each otu & htu
# 	matrix tree give one row per node with each otu & htu in that node
#	NOTE: the basal node is htu = otu+1
accio_matrix_tree_from_vector_tree <- function(vector_tree)	{
	node_rosetta <- sort(unique(vector_tree[vector_tree>0]))
	Nnodes <- length(node_rosetta)
	maxtomy <- 2
	if (length(vector_tree)<((2*Nnodes)+1))	{
		for (n in 1:length(node_rosetta))	{
			if (sum(vector_tree %in% node_rosetta[n]) > maxtomy)	maxtomy <- sum(vector_tree %in% node_rosetta[n])
			}	
		}
	#order(vector_tree)[2:length(vector_tree)]
	node_rich <- vector(length=Nnodes)
	matrix_tree <- matrix(0,Nnodes,maxtomy)
	for (i in 1:length(vector_tree))	{
		if (vector_tree[i]>=node_rosetta[1])	{
			node <- match(vector_tree[i],node_rosetta)
			node_rich[node] <- node_rich[node]+1
			matrix_tree[node,node_rich[node]] <- i
			}
		}
	return(matrix_tree)
	}



# convert string of states into a polymorphic
ravel_polymorph <- function(polies)	{
	polies<-as.numeric(polies)
	poly <- 0
	for (i in 1:length(polies))	poly <- poly + (polies[i]*10^(i-1))
	return(-1*poly)
	}

### ROUTINES TO ORDER BRANCHES ON TREE
# get order of branches from base of tree up, accounting for stratigraphy
accio_branch_order_with_stratigraphy <- function(mtree,firstAppearances)	{
	branch_order <- accio_patristic_distance_from_base(mtree)	# order branches from base of the tree
	used_order <- branch_order + (1-(1:length(branch_order))/100)	# for ties, put nodes before otus
	strat_order <- date_taxa_on_tree_simple(mtree,firstAppearances)
	tosort <- cbind(strat_order,used_order)
	return(order(tosort[,1],tosort[,2],decreasing=FALSE))
	}

# get order of branches from base of tree up
accio_patristic_distance_from_base <- function(atree)	{
	ttus <- max(atree)
	pat_dist_from_base <- vector(length=ttus)
	if (length(dim(as.array(atree)))==1)	{
		base <- min(atree[atree>0])
		Nnodes <- 1+max(atree)-base
		pat_dist_from_base[base] <- 0
		for (n in (base+1):ttus)	pat_dist_from_base[n] <- 1+pat_dist_from_base[atree[n]]
		for (s in 1:(base-1))		pat_dist_from_base[s] <- 1+pat_dist_from_base[atree[s]]
		}	else 	{
		Nnodes <- dim(atree)[1]
		otus <- ttus-Nnodes
		base <- otus+1
		pat_dist_from_base[base] <- 0
		for (t in 1:Nnodes)	{
			f1 <- sum(atree[t,]>0)
			ht <- t+otus
			for (f in 1:f1)	{
				f2 <- atree[t,f]
				pat_dist_from_base[f2] <- pat_dist_from_base[ht]+1
				}
			}
		}
	return(pat_dist_from_base)
	}

# get simple dating of taxa on trees
date_taxa_on_tree_simple<-function(atree,firstAppearances)	{
	# warning - original code has mysterious variable 'vtree'
		# was this a typo of 'atree'?? replaced
	# vtree means 'vector tree', so I guess 'atree' should be a 'vtree' type object
	### vtree is a vector tree, where each cell gives the ancestral HTU of the OTU or HTU
	###		HTU's are nodes, with the basal HTU = notu+1 (i.e., 13 for 12 taxa)
	#
	ttus <- max(atree)	# total htus+otus
	dates <- vector(length=ttus)
	notu <- length(firstAppearances)	# number of otus
	if (length(dim(as.array(atree)))==1)	{
		# routine if tree given as a single vector with ancestral
		base <- notu+1
		end<-max(firstAppearances)
		for (nd in base:ttus)	dates[nd]<-end
		for (sp in 1:notu)	{
			if (atree[sp]>0)	{ 	#vtree?
				dates[sp]<-firstAppearances[sp]
				anc<-atree[sp] 	#vtree?
				if (dates[anc]>dates[sp])	dates[anc]<-dates[sp]
				}	# this is to make sure that species excluded from the tree do not mess up analysis
			}
		for (nd in ttus:(base+1))	{
			if (atree[nd]>0)	{ 	#vtree?
				anc<-atree[nd] 	#vtree?
				if (dates[anc]>dates[nd])	dates[anc]<-dates[nd]
				}	# this is to make sure that species excluded from the tree do not mess up analysis
			}
		} else	{
		Nnodes <- dim(atree)[1]
		dates[1:notu] <- firstAppearances
		for (n in Nnodes:1)	{
			f1 <- sum(atree[n,]>0)
			ht <- n+notu
			dates[ht] <- min(dates[atree[n,(1:f1)]])
			}
		}
	return(dates)
	}

### FUNCTIONS FOR EVALUATING DISTRIBUTIONS OF CHANGES
# second-order (modified) AIC
# second-order (modified) AIC
modified_AIC <- function(lnL,k,n)	{
	# L: log-likelihood; # k: parameters; # n: data points
	#if (n==k)	n <- n+2
	if (is.na(n / (n - k - 1)) || (n / (n - k - 1)<0) || (n / (n - k - 1))==Inf)	{
		aic_c <- AIC(lnL,k)
		}	else	aic_c <- (-2*lnL) + (2*k)*(n / (n - k - 1))
	return(aic_c)
	}

# convert abundances to "Fisher plot" giving # taxa with 1…N Finds
fisher_plot <- function(finds)	{
	return(hist(finds,breaks=(0:max(finds)),plot=FALSE)$counts)
	}

create_quantile_vector<-function(N)	{
	return(seq(1/(N+1),N/(N+1),by=1/(N+1)))
	}

# find expectations of this gamma at this sample size
# exp will give the expected number of taxa found 1…ncoll times
expected_abundances <- function(rel_ab_dist, nspec, S)	{
	expected <- vector(length=nspec)
	for (t in 1:S)	{
		for (i in 1:nspec)	{
			expected[i] <- expected[i]+dbinom(i,nspec,rel_ab_dist[t])
			}
		}
	return(expected)
	}

# expected sampling curve
expected_sampling_curve <- function(rel_ab_dist, nspec, S)	{
	exp_find <- c()
	for (i in 1:nspec)	{
		exp_find <- c(exp_find,expected_sampled_richness(rel_ab_dist, i, S))
		}
	return(exp_find)
	}

# expected exhaustion curve
expected_exhaustion_curve <- function(rel_ab_dist, nstep, after3, C)	{
	# like sampling curve, but first three branches are free!
	xx <- expected_sampling_curve(rel_ab_dist,nspec=nstep,S=C)
	d <- sum(xx<after3)
	xfinds <- xx[d]
	while (xfinds < after3)	{
		d <- d+0.01
		xfinds <- expected_sampled_richness(rel_ab_dist, d, S=C)
		}
	dd <- d-after3
	exp_states <- (1:after3)
	for (i in (after3+1):nstep)	{
		exp_states <- c(exp_states,expected_sampled_richness(rel_ab_dist, i+dd, S=C))
	#	exp_states <- c(exp_states,expected_sampled_richness(rel_ab_dist, i, S))
		}
	return(exp_states)
	}

# expected sampled
# get the likelihood of observed numbers of taxa with 1…N finds given expected numbers of taxa with 1…N finds
expected_sampled_richness <- function(rel_ab_dist, nspec, S, MINNO = 5e-324)	{
	exp_find <- 0
	for (t in 1:S)
	#	exp_find <- c(exp_find,(1-dbinom(0,nspec,rel_ab_dist[t])))
	#	exp_find <- exp_find+(1-dbinom(0,nspec,rel_ab_dist[t]))
		exp_find <- exp_find+(1-(1-rel_ab_dist[t])^nspec)
	return(exp_find)
	}


# get multinomial likelihood of distribution
distribution_loglikelihood_mul <- function(observed,expected,oS,hS,MINNO = 5e-324)	{
	#print(c(oS,hS))		# for debugging
	mxfind <- length(observed)							# maximum finds observed
	prop_expected <- expected[1:mxfind]/sum(expected)		# convert expected species to proportions
	prop_expected[prop_expected==0] <- MINNO
	lnlo <- observed*log(prop_expected)	# exact probability of observing # taxa with 1, 2, 3, etc. finds
	lnlo[is.na(lnlo)] <- 0
	#for (i in 1:mxfind)	if (is.na(lnlo[i]))	lnlo[i] <- 0
	sobs <- sum(lnlo)
	# log probability of observing X taxa from hypothesized hS taxa given expected # taxa with 1, 2, 3, etc. finds
	eS <- sum(expected)									# get expected sampled species
	while (eS==hS)	eS <- 0.99999*eS
		# this is a kluge to get around rounding error of eS->hS

	if (hS>=oS && round(eS,5)<round(hS,5))	{
		# get the probability of observing oS of hS species given that we expected to observe eS of hS species
		lnls <- lfactorial(hS)-(lfactorial(hS-oS)+lfactorial(oS))+(oS*log(eS/hS))+((hS-oS)*log((hS-eS)/hS))
		# this should be the norm: more species than observed hypothesized
		}	else	{
		lnls <- oS*log(MINNO)
		# if hypothesized true richness is less than observed, then this is impossible
		}
	return(sobs+lnls)
	}

# find best model with just one frequency
optimize_uniform_abundance <- function(counts, MAXNO = 1.797693e+308)	{
	# written 2017-01-28
	minS <- length(counts)
	nspec <- sum(counts)
	oS <- hS <- minS
	observed <- fisher_plot(finds=counts)
	mxlnl <- lnl <- -1*MAXNO
	while (lnl == mxlnl)	{
		rel_ab_dist <- rep(1/hS,hS)
		raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
		lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
		if (mxlnl < lnl)	mxlnl <- lnl
		hS <- hS+1
		}
	bH <- c(hS-1,round(mxlnl,3),round(modified_AIC(mxlnl,1,nspec),3))
	names(bH) <- c("Uniform_S","Uniform_log-likelihood","Uniform_AICc")
	return(bH)
	}

# some basic distributions
geometric_distribution <- function(decay, MINEXPN = 10^-10)	{
	#	rescale geometric to base rate r
	if (decay>1)	decay <- (1/decay)
	S <- round(1+((log(MINEXPN)-log(1-decay))/log(decay)))
	ranks <- (1:S)
	prop <- (1-decay)*(decay^(ranks-1))
	return(prop)
	}

# get log-likelihood of geometric given decay
loglikelihood_geometric_rad <- function(decay,nspec,oS,observed, MINNO = 5e-324)	{
	# p0[1]: r	# p0[2]: decay	# p0[3]: S
	rel_ab_dist <- geometric_distribution(decay)		# basic exponential distribution
	hS <- length(rel_ab_dist)
	#print(c(round(decay,5),hS))		# for debugging
	if (hS>=oS)	{
		raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
			# log likelihood
		lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
		}	else	{
		lnl <- oS*log(MINNO)
		}
	#print(c(round(decay,5),round(lnl,10)))		# for debugging
	return(lnl)
	}

# get best-fit geometric
optimize_geometric_abundance <- function(counts)	{
	oS <- length(counts)				# observed taxa
	observed <- fisher_plot(counts)
	decay <- (min(counts)/max(counts))^(1/(length(counts)-1))
	max_decay <- exp(log(decay)/2)
	min_decay <- exp(log(decay)*2)
	#
	cl <- list(fnscale=-1)
	nspec <- sum(counts)
	w <- optim(decay,fn=loglikelihood_geometric_rad,method="L-BFGS-B",nspec=nspec,oS=oS,observed=observed,lower=min_decay,upper=max_decay,control=cl)
	bH <- c(round(w$par,6),round(w$value,3),round(modified_AIC(w$value,1,nspec),3))
	names(bH) <- c("Geometric_decay","Geometric_log-likelihood","Geometric_AICc")
	return(bH)
	}

# calcualate lognormal for S entities that increase in magnitude by mag every standard deviation
lognormal_distribution <- function(mag, S)	{
	fi <- (S:1)/(S+1)
	prop <- mag^(qnorm(fi,0,1))/sum(mag^(qnorm(fi,0,1)))
	return(prop)
	}

# jackknife 5th-order richness estimator
jack5_Fisher<-function(observed)	{
	ntaxa<-sum(observed)
	ss<-0
	for (i in 1:length(observed))	ss<-ss+(i*observed[i])
	S <- round(ntaxa + (observed[1]*(((5*ss)-15)/ss)) - observed[2]*((10*(ss*ss)-(70*ss)+125)/(ss*(ss-1))) + observed[3]*(((10*(ss^3))-(120*(ss^2))+(485*ss)-660))/(ss*(ss-1)*(ss-2)) - observed[4]*((ss-4)^4)/(ss*(ss-1)*(ss-2)*(ss-3)) + observed[5]*((ss-5)^5)/(ss*(ss-1)*(ss-2)*(ss-3)*(ss-4)))
	return(S)
	}

# get minimum and maximum magnitude of increase for a lognormal given hS
accio_min_and_max_lognormal_mag_given_hS <- function(observed,counts,hS)	{
	oS <- sum(observed)				# observed taxa
	mxfinds <- length(observed)
	mnfinds <- min((1:mxfinds)[!observed %in% 0])
	rank_shifts <- sort(qnorm((1:hS)/(hS+1)),decreasing=TRUE)[1:oS]
	local_m <- mag_shifts <- vector(length=(oS-1))
	for (i in 1:(oS-1))		{
		mag_shifts[i] <- counts[i]/counts[i+1]
		local_m[i] <- exp(log(mag_shifts[i])/(rank_shifts[i]-rank_shifts[i+1]))
		}
	if (max(observed)==1)	{
		min_mag <- min(local_m)
		}	else {
		min_mag <- 1.01
		}
	max_mag <- prod(local_m)^(1/(length(local_m)))
	if (max_mag < min_mag)	{
		m <- min_mag
		min_mag <- max_mag
		max_mag <- m
		}
	return(c(min_mag,max_mag))
	}

# get log-likelihood for lognormal
loglikelihood_lognormal_rad_for_optim <- function(oS,nspec,rand_no,hS,observed,min_mag,max_mag)	{
	mag <- min_mag + rand_no*(max_mag - min_mag)
	#print(c(round(mag,9)))
	rel_ab_dist <- lognormal_distribution(mag=mag,S=hS)		# basic lognormal distribution
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
	return(round(lnl,2))
	}

# get loglikelihood of hypothesized lognormal
loglikelihood_lognormal_rad<-function(oS,nspec,mag,hS,observed)	{
	#get basic information
	#nspec<-0
	#for (i in 1:length(observed))	nspec<-nspec+(i*observed[i])
	#oS<-sum(observed)
	#optim(p0,fn=bestuniformpresrate,method="L-BFGS-B", lower=minr, 
	#	upper=maxr, control=cl, observed=observed, ncoll=ncoll)
	rad<-lognormal_distribution(mag,S=hS)		# basic lognormal distribution
	raw_expected<-expected_abundances(rad,nspec,S=hS)
	# I suspect the following function is meant to be replaced, as in loglikelihood_lognormal_rad_for_optim
	#lnl<-abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
	lnl<-distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
	return(lnl)
	}

# get most likely magnitude of increase for a lognormal distribution of hS entities
optimize_lognormal_abundance_given_hS <- function(observed,counts,hS)	{
	#print(hS)		# for debugging
	oS <- sum(observed)				# observed taxa
	nspec <- sum((1:length(observed))*observed)
	mxfinds <- length(observed)
	mnfinds <- min((1:mxfinds)[!observed %in% 0])
	mm <- accio_min_and_max_lognormal_mag_given_hS(observed,counts,hS)
	mag <- prod(mm)^0.5
	#inev <- mag <- exp(log(mxfinds/mnfinds)/(qnorm((hS-1)/(hS+1))-qnorm((hS-oS)/(hS+1))))
	cl <- list(fnscale=-1)
	rand_no <- (mag - mm[1])/(mm[2] - mm[1])
	w <- optim(rand_no,loglikelihood_lognormal_rad_for_optim,,method="L-BFGS-B",
		oS=oS,nspec=nspec,hS=hS,observed=observed,min_mag=mm[1],max_mag=mm[2],lower=0,upper=1,control=cl)
	w$par <- mm[1] + (w$par*(mm[2] - mm[1]))
	#w <- optim(mag,fn=loglikelihood_lognormal_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,
	#	observed=observed,lower=mm[1],upper=mm[2],control=cl)
	bH <- c(w$par,hS,w$value)
	names(bH) <- c("Lognormal_magnitude", "Lognormal_S","Lognormal_log-likelihood")
	return(bH)
	}

# get most likely lognormal distribution for a population of some sort
optimize_lognormal_abundance <- function(counts)	{
	stS <- oS <- length(counts)				# observed taxa
	observed <- fisher_plot(counts)
	nspec <- sum(counts)
	span <- 5
	enS <- stS+((span-1)*oS)
	incr <- floor(enS/span)
	cl <- list(fnscale=-1)
	peak <- 0
	while (incr>0)	{
		hS <- seq(stS,enS,by=incr)
		results <- sapply(hS,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
		mlnl <- max(results[3,])
		mlSc <- match(mlnl,results[3,])
		mlS <- hS[mlSc]
		bH <- results[,mlSc]
		spn <- length(results[3,])
		if (incr>1)	{
			# if runs so far are still producing higher likelihoods at higher richnesses:
			if (mlSc==spn)	{
				if (peak==0) {
					# if we have not yet found a peak, then keep increasing richness
					stS <- hS[spn]
					enS <- stS+((spn-1)*oS)
					}	else	{
					enS <- hS[spn]
					if (span<incr)	{
						incr <- floor((enS-hS[spn-1])/span)
						}	else {
						span <- enS-hS[spn-1]
						incr <- 1
						}
					stS <- hS[spn]-(span*incr)
					}	# end case where last number is best, but this is after finding a peak.
				} else if (mlSc==1)	{
				# if first richness is the best
				peak <- 1
				stS <- hS[mlSc]
				enS <- hS[mlSc+1]-1
				if (incr==1)	incr <- 0
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				} else {
				# if a richness in the middle is best & we still can find a better one
				peak <- 1
				hS2 <- c(mlS-1,mlS+1)
				results2 <- sapply(hS2,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
				if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
				# start just above 2nd best richness and go up to best
					stS <- hS[mlSc-1]+1
					enS <- hS[mlSc]
					if (span<incr)	{
						incr <- floor((enS-stS)/span)
						}	else {
						span <- enS-stS
						incr <- 1
						}
					}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
				# start at best richness and go just below 2nd best
					stS <- hS[mlSc]
					enS <- hS[mlSc+1]-1
					if (incr>span)	{
						incr <- floor((enS-stS)/span)
						}	else {
						span <- enS-stS
						incr <- 1
						}
					}	else {
					# we already had the best, so just end it
					incr <- 0
					}
				}
			# end case where we have a better richness in middle somewhere.  
			}	else	{
			incr <- 0
			}
		}
	bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,2,nspec))
	names(bH)[4] <- "Lognormal_AICc"
	return(bH)
	}

# generate Zipf distribution for S individuals with log-decay = zm
zipf_distribution <- function(zm, S)	{
	return(((1:S)^-zm)/sum((1:S)^-zm))
	#zpf <- ranks^-zm
	#prop <- zpf/sum(zpf)
	#return(prop)
	}

# get the likelihood of a Zipf model with log-log decay = zm & hS entities given oS observed with nspec individuals
loglikelihood_zipf_rad <- function(oS,nspec,zm,hS,observed)	{
	#print(zm)
	rel_ab_dist <- zipf_distribution(zm=zm,S=hS)		# basic zipf distribution
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
	return(round(lnl,3))
	}

# get the most likely Zipf model for a particular richness
optimize_zipf_abundance_given_hS <- function(counts,observed,max_zipf,hS)	{
	#print(hS)		# for debugging
	oS <- sum(observed)				# observed taxa
	nspec <- sum((1:length(observed))*observed)
	mxfinds <- length(observed)
	mnfinds <- min((1:mxfinds)[!observed %in% 0])
	inzm <- zm <- log(mxfinds/mnfinds)/log(oS-1)
	cl <- list(fnscale=-1)
	w <- optim(zm,fn=loglikelihood_zipf_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=0,upper=max_zipf,control=cl)
	bH <- c(w$par,hS,w$value)
	names(bH) <- c("Zipf_log-log_decay","Zipf_hypothesized_richness","Zipf_log-likelihood")
	return(bH)
	}

# get the maximum zipf decay implied by the distribution of data points
accio_max_zipf_decay <- function(counts,oS,nspec)	{
	xx <- log(1:oS)
	yy <- log(counts/nspec)-log(counts[1]/nspec)
	zz <- -1*(yy/xx)[2:length(xx)]
	return(max(zz))
	}

# get the most likely Zipf model
optimize_zipf_abundance <- function(counts,span=5)	{
	stS <- oS <- length(counts)				# observed taxa
	observed <- fisher_plot(counts)
	enS <- stS+((span-1)*oS)
	incr <- floor(enS/span)
	cl <- list(fnscale=-1)
	peak <- 0
	nspec <- sum(counts)
	max_zipf <- 2*accio_max_zipf_decay(counts,oS,nspec)
	pa <- 1
	pz <- span
	while (incr>0)	{
		hS <- seq(stS,enS,by=incr)
		results <- sapply(hS[pa:pz],optimize_zipf_abundance_given_hS,counts=counts,observed=observed,max_zipf=max_zipf)
		if (pa==2)	{
			results <- cbind(bH,results)
			} else if (pz<span)	{
			results <- cbind(results,bH)
			}
		mlnl <- max(results[3,])
		mlSc <- match(mlnl,results[3,])
		mlS <- hS[mlSc]
		bH <- results[,mlSc]
		spn <- length(results[3,])
		if (incr>1)	{
			# if runs so far are still producing higher likelihoods at higher richnesses:
			if (mlSc==spn)	{
				if (peak==0) {
					# if we have not yet found a peak, then keep increasing richness
					stS <- hS[spn]
					pa <- 2
					pz <- span
					enS <- stS+((spn-1)*oS)
					}	else	{
					enS <- hS[spn]
					if (span<incr)	{
						incr <- floor((enS-hS[spn-1])/span)
						}	else {
						span <- enS-hS[spn-1]	# check this! TDay
						incr <- 1
						}
					stS <- hS[spn]-(span*incr)
					pa <- 1
					pz <- span-1
					}	# end case where last number is best, but this is after finding a peak.
				} else if (mlSc==1)	{
				# if first richness is the best
				peak <- 1
				stS <- hS[mlSc]
				enS <- hS[mlSc+1]-1
				if (incr==1)	incr <- 0
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- span
				} else {
				# if a richness in the middle is best & we still can find a better one
				peak <- 1
				hS2 <- c(mlS-1,mlS+1)
				results2 <- sapply(hS2,optimize_zipf_abundance_given_hS,counts=counts,observed=observed,max_zipf=max_zipf)
				if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
				# start just above 2nd best richness and go up to best
					bH <- results2[,1]
					stS <- hS[mlSc-1]+1
					enS <- hS2[1]
					if (span<incr)	{
						incr <- floor((enS-stS)/span)
						}	else {
						span <- 1+enS-stS
						incr <- 1
						}
					pa <- 2
					pz <- span
					}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
				# start at best richness and go just below 2nd best
					stS <- hS2[2]
					enS <- hS[mlSc+1]-1
					bH <- results2[,2]
					if (incr>span)	{
						incr <- floor((enS-stS)/span)
						}	else {
						span <- 1+enS-stS
						incr <- 1
						}
					pa <- 1
					pz <- span-1
					}	else {
					# we already had the best, so just end it
					incr <- 0
					}
				}
			# end case where we have a better richness in middle somewhere.  
			}	else	{
			incr <- 0
			}
		}
	bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,k=2,n=nspec))
	names(bH)[4] <- "Lognormal_AICc"
	return(bH)
	}

# get gamma distribution
gamma_distribution <- function(a, b, S)	{
	p0 <- vector(length=S)
	p0[1] <- S/(S+1)
	for(i in 2:S) p0[i] <- p0[i-1]-(1/(S+1))
	prop <- qgamma(p0,a,b)/sum(qgamma(p0,a,b))
	return(prop)
	}

# get minimum alpha
accio_minimum_alpha_for_gamma_one <- function(hS, MINNO = 5e-324)	{
	min_alpha <- 1
	iS <- hS
	while (iS==hS)	{
		rel_ab_dist <- gamma_distribution(a=min_alpha,b=min_alpha,S=hS)		# basic lognormal distribution
		iS <- sum(rel_ab_dist>MINNO)
		if (iS == hS)	min_alpha <- min_alpha/2
		}
	return(2*min_alpha)
	}

# get log-likelihood for gamma
loglikelihood_gamma_one_rad <- function(oS,nspec,alpha,hS,observed, MINNO = 5e-324, MAXNO = 1.797693e+308)	{
	#print(c(log10(alpha),hS))
	rel_ab_dist <- gamma_distribution(a=alpha,b=alpha,S=hS)		# basic lognormal distribution
	#print(rel_ab_dist)
	if (sum(rel_ab_dist>MINNO) == hS)	{
		raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
		#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
		lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
		}	else	lnl <- -1*MAXNO
	return(round(lnl,2))
	}

# get most likely magnitude of increase for a lognormal distribution of hS entities
optimize_gamma_one_abundance_given_hS <- function(observed,counts,hS)	{
	#print(hS)		# for debugging
	oS <- sum(observed)				# observed taxa
	nspec <- sum((1:length(observed))*observed)
	cl <- list(fnscale=-1)
	alpha <- 2
	min_alpha <- accio_minimum_alpha_for_gamma_one(hS)
	w <- optim(alpha,loglikelihood_gamma_one_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=min_alpha,upper=10000,control=cl)
	#w <- optim(mag,fn=loglikelihood_lognormal_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=mm[1],upper=mm[2],control=cl)
	bH <- c(w$par,hS,w$value)
	names(bH) <- c("Gamma_alpha=beta", "Gamma_S","Gamma_log-likelihood")
	return(bH)
	}

# get most likely gamma distribution for a population of some sort
optimize_gamma_one_abundance <- function(counts)	{
	stS <- oS <- length(counts)				# observed taxa
	observed <- fisher_plot(counts)
	nspec <- sum(counts)
	span <- 5
	enS <- stS+((span-1)*oS)
	incr <- floor(enS/span)
	cl <- list(fnscale=-1)
	peak <- 0
	while (incr>0)	{
		hS <- seq(stS,enS,by=incr)
		results <- sapply(hS,optimize_gamma_one_abundance_given_hS,observed=observed,counts=counts)
		mlnl <- max(results[3,])
		mlSc <- match(mlnl,results[3,])
		mlS <- hS[mlSc]
		bH <- results[,mlSc]
		spn <- length(results[3,])
		if (incr>1)	{
			# if runs so far are still producing higher likelihoods at higher richnesses:
			if (mlSc==spn)	{
				if (peak==0) {
					# if we have not yet found a peak, then keep increasing richness
					stS <- hS[spn]
					enS <- stS+((spn-1)*oS)
					}	else	{
					enS <- hS[spn]
					if (span<incr)	{
						incr <- floor((enS-hS[spn-1])/span)
						}	else {
						span <- enS-hS[spn-1]
						incr <- 1
						}
					stS <- hS[spn]-(span*incr)
					}	# end case where last number is best, but this is after finding a peak.
				} else if (mlSc==1)	{
				# if first richness is the best
				peak <- 1
				stS <- hS[mlSc]
				enS <- hS[mlSc+1]-1
				if (incr==1)	incr <- 0
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				} else {
				# if a richness in the middle is best & we still can find a better one
				peak <- 1
				hS2 <- c(mlS-1,mlS+1)
				results2 <- sapply(hS2,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
				if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
				# start just above 2nd best richness and go up to best
					stS <- hS[mlSc-1]+1
					enS <- hS[mlSc]
					if (span<incr)	{
						incr <- floor((enS-stS)/span)
						}	else {
						span <- enS-stS
						incr <- 1
						}
					}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
				# start at best richness and go just below 2nd best
					stS <- hS[mlSc]
					enS <- hS[mlSc+1]-1
					if (incr>span)	{
						incr <- floor((enS-stS)/span)
						}	else {
						span <- enS-stS
						incr <- 1
						}
					}	else {
					# we already had the best, so just end it
					incr <- 0
					}
				}
			# end case where we have a better richness in middle somewhere.  
			}	else	{
			incr <- 0
			}
		}
	bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,2,nspec))
	names(bH)[4] <- "Gamma_AICc"
	return(bH)
	}



##### GRAPHICAL GOODNESS
draw_symmetric_axes <- function(mxx,xsize)	{
	if (mxx<=100)	{
		brks <- c(10,5,1)
		axis(1,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		axis(2,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		tcs <- seq(0,mxx,by=brks[1])
		axis(1,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,cex=xsize/3)
		axis(2,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,las=2,cex=xsize/3)
		tcs <- seq(0,mxx,by=brks[2])[!seq(0,mxx,by=brks[2]) %in% seq(0,mxx,by=brks[1])]
		prn <- TRUE
		if (mxx>50)	prn <- FALSE
		axis(1,at=tcs,tcl=-0.2,labels=prn,lwd=0.0,lwd.ticks=1.1,cex=xsize/3)
		axis(2,at=tcs,tcl=-0.2,labels=prn,lwd=0.0,lwd.ticks=1.1,las=2,cex=xsize/3)
		tcs <- seq(0,mxx,by=brks[3])[!seq(0,mxx,by=brks[3]) %in% seq(0,mxx,by=brks[2])]
		axis(1,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		}	else if (mxx<=200)	{
		brks <- c(20,10,2)
		axis(1,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		axis(2,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		tcs <- seq(0,mxx,by=brks[1])
		axis(1,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,cex=xsize/3)
		axis(2,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,las=2,cex=xsize/3)
		tcs <- seq(0,mxx,by=brks[2])[!seq(0,mxx,by=brks[2]) %in% seq(0,mxx,by=brks[1])]
		axis(1,at=tcs,tcl=-0.2,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.2,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		tcs <- seq(0,mxx,by=brks[3])[!seq(0,mxx,by=brks[3]) %in% seq(0,mxx,by=brks[2])]
		axis(1,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		}	else if (mxx<500)	{
		brks <- c(50,25,5)
		axis(1,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		axis(2,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		tcs <- seq(0,mxx,by=brks[1])
		axis(1,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,cex=xsize/3)
		axis(2,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,las=2,cex=xsize/3)
		tcs <- seq(0,mxx,by=brks[2])[!seq(0,mxx,by=brks[2]) %in% seq(0,mxx,by=brks[1])]
		axis(1,at=tcs,tcl=-0.2,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.2,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		tcs <- seq(0,mxx,by=brks[3])[!seq(0,mxx,by=brks[3]) %in% seq(0,mxx,by=brks[2])]
		axis(1,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		}	else	{
		brks <- c(100,50,10)
		axis(1,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		axis(2,at=seq(0,mxx),tcl=0.0,labels=FALSE,lwd=1.1,lwd.ticks=0)
		tcs <- seq(0,mxx,by=brks[1])
		axis(1,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,cex=xsize/3)
		axis(2,at=tcs,tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=1.1,las=2,cex=xsize/3)
		tcs <- seq(0,mxx,by=brks[2])[!seq(0,mxx,by=brks[2]) %in% seq(0,mxx,by=brks[1])]
		axis(1,at=tcs,tcl=-0.2,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.2,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		tcs <- seq(0,mxx,by=brks[3])[!seq(0,mxx,by=brks[3]) %in% seq(0,mxx,by=brks[2])]
		axis(1,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1)
		axis(2,at=tcs,tcl=-0.1,labels=FALSE,lwd=0.0,lwd.ticks=1.1,las=2)
		}
	}

accio_matrix_tree_from_ape_tree <- function(ape)	{
	### ape is the output of read.nexus
	#	ape$edge gives apes' version of the tree, which alters the numbers
	#	This changes them back to the numbers in the original (1,(2,3)) format
	#	The information for doing this is in ape$tip.label
	base <- ape$edge[1,1]
	notu <- base-1	# the lowest htu number is 1 more than maximum otu number
	yyy <- ape$edge[order(ape$edge[,2]),]
	yyy[1:notu,2] <- as.numeric(ape$tip.label)
	zzz <- yyy[order(yyy[,2]),]
	vector_tree <- zzz[order(zzz[,2]),1]
	vector_tree <- c(vector_tree[1:(base-1)],-1,vector_tree[base:length(vector_tree)])
	mtree <- accio_matrix_tree_from_vector_tree(vector_tree)
	}
	
##############################################################

# old stuff



#ZERO  1e-323
#MINEXPN <- 10^-10
#MINNO = 5e-324
#MAXNO = 1.797693e+308
#gap <- inapplicableValue <- "-"
#missingCharValue <- missing <- "?"
#missingCharValue <- "?"		# replace "?" with "?"
#inapplicableValue <- "-"			# replace "-" with "-"

#library(paleotree)
#source("D://dave/workspace/paleotree/R/exhaustionFunctions.R")
#nexustreefile <- "Dicranograptidae_Tree.txt"
#phyloTree <- read.tree(nexustreefile)
#strat_file <- "Dicranograptidae_Ranges.txt"
#strat_data <- read.table(strat_file,header=FALSE,sep="\t")
#setwd("D:\\dave\\research\\0 wagner exhaustion code")
#nexus_file_name <- "Dicranograptidae_Song_&_Zhang_2014.nex"
#char_data <- accio_data_from_nexus_file(nexus_file_name,polymorphs=FALSE)
# rewrite read.tree output to match what I use
#nTips <- Ntip(phyloTree)
#
#charMat <- char_data$Matrix
#########################################################
# modify data
# what is this outgroup thing??
#if (char_data$Outgroup>0)	{
#	outgroup <- char_data$Outgroup
#	}	else	{
#	outgroup <- 1	
#	}
#
# modify strat data
#
#if (length(simplify2array(strat_data))>notu)	{
#	FAs <- strat_data[,1]
#}else{
#	FAs <- strat_data
#	}

	
	
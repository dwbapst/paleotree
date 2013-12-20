#constrain 




constrainParPaleo<-function(f, ..., formulae=NULL, names=parnames(f),bounds=parbounds(f),extra=NULL) {
	#based on Rich FitzJohn's constrain function for diversitree 10-22-13
		# comment lines with double ## indicate Rich's original comments
		#I claim all responsibility for how ugly I can make Rich's code
		#...its actually kind of a challenge - DWB
	##
	## For the first case, everything is OK on the lhs and rhs
	## For subsequent cases:
	## lhs cannot contain things that are
	## - constrained things (already lhs anywhere)
	## - things constrained to (things on the rhs anywhere)
	## rhs cannot contain things that are
	## - constrained things (already lhs anywhere)
	## It is possibly worth pulling out all the numerical constants and
	## the "paired" parameters here to avoid using eval where
	## unnecessary. However, this makes the function substantially uglier
	## for a very minor speedup.
	#
	#let's make some example data
		#f=function(pqr=c(p,q,r)){pqr[1]^2+2*pqr[2]+pqr[3]} 
		#f<-make_paleotreeFunc(f,c("p","q","r"),list(c(0,0,0),rep(Inf,3)))
		#formulae=c(p~q,list());names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
		#f=function(pqr=c(p,q,r)){p^2+2*q+r}; 
		#f<-make_paleotreeFunc(f,c("p","q","r"),list(c(0,0,0),rep(Inf,3)))
		#formulae=c(p~q,r~q,list());names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
		#f=function(pqr=c(p.1,p.2,q.1,q.2,r.1,r.2)){p.1^2+2*q.1+r.1/(p.2^2+2*q.2+r.2)}
		#f<-make_paleotreeFunc(f,c("p.1","q.1","r.1","p.2","q.2","r.2"),list(rep(0,6),rep(Inf,6)))
		#formulae=c(p.1~q.all,p.match~r.match,list());names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
		#FOR USE WITH known f
		#formulae = c(rRate~pRate);names=parnames(f);extra=NULL;bounds=parbounds(f)
	#
	if(!is(f,'paleotreeFunc')){
		stop("Error: Given function does not appear to be a paleotree likelihood function")}
	if ( inherits(f, "constrained") ) {	#this thing checks to see if its already constrained 
		formulae <- c(attr(f, "formulae"), formulae)
		f <- attr(f, "origfunction")
		}
	rels <- list()				#rels are the things we're gonna constrain to something else
	names.lhs <- names.rhs <- names	#lhs is untouched pars??, rhs is the final pars??
	formulae <- c(formulae, list(...))	#adding the ... to formulae
	#expand formulae in case they contain systematic constraints here! here's some examples:
		#names=c("p.1.1","q.1.1","p.2.1","q.2.1","p.1.2","q.1.2","p.2.2","q.2.2")
		#formulae = c(p.all.match~q.all.match,list())
		#formulae = c(p.1.1~p.all.all,list())
		#formulae = c(p.1.match~q.all.match,list())
		#formulae = c(p.1.match~q.all.match,p.1.1~p.all.all,list())
	breakTerms<-lapply(formulae,function(x) unlist(strsplit(all.vars(as.formula(x)),".",fixed=TRUE)))
	needExpand<-sapply(breakTerms,function(x) any("match"==x)|any("all"==x))
	if(any(needExpand)){
		breakNames<-t(rbind(sapply(names,function(x) unlist(strsplit(x,".",fixed=TRUE)))))
		nparcat<-ncol(breakNames)
		newFormulae<-list()
		for(i in which(needExpand)){
			newFormulae<-c(newFormulae,expandConstrainForm(formulae[[i]],breakNames,nparcat))
			}
		formulae<-c(formulae[-which(needExpand)],newFormulae) #formulae		
		formulae<-unique(formulae)	#any duplicates?
		}
	for( formula in formulae ) {
		res <- constrainParsePaleo(formula, names.lhs, names.rhs, extra)
		if ( attr(res, "lhs.is.target") ) {
			i <- try(which( sapply(rels,function(x) identical(x, res[[1]]))),silent=TRUE)
			if(inherits(i,"try-error")){
				stop(sprintf("Error parsing constraint with %s on lhs",as.character(res[[1]])))
				}
			rels[i] <- res[[2]]	#DWB: gives warning message that symbol cannot be coerced to list
			## This will not work with *expressions* involving the LHS; that
			## would require rewriting the expressions themselves (which
			## would not be too hard to do). But for now let's just cause
			## an error...
			lhs.txt <- as.character(res[[1]])
			if ( any(sapply(rels, function(x) lhs.txt %in% all.vars(x))) ){
				stop(sprintf("lhs (%s) is in an expression and can't be constrained",lhs.txt))
				}
			}
		names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
		names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
		rels <- c(rels, structure(res[2], names=as.character(res[[1]])))
	  	}
	#in a function (p,q,r), with constraint p~q, names.lhs="r", names.rhs="q""r", rels = list(p=r)
	#okay, now we know which ones will be lhs, rhs and rels
	#need to test that all the bounds for the equivalencies are the same
		#check each rels for consistency with bounds
	relsIsPar<-rels[sapply(rels,function(x) names==x)]	#test to see which are pars
	if(length(relsIsPar)>0){
		termBound<-t(sapply(names(relsIsPar),function(x) 
			sapply(bounds,function(y) y[which(names==x)])))
		relsBound<-t(sapply(relsIsPar,function(x) 
			sapply(bounds,function(y) y[which(names==x)])))
		colnames(relsBound)<-colnames(termBound)<-NULL
		if(!identical(relsBound,termBound)){
			noMatch<-which(!apply(termBound==relsBound,1,all))
			noMatch<-paste(i,"~",names(relsIsPar)[noMatch])
			stop(paste("Upper and Lower bounds do not match for",noMatch))
			}
		}
	#back to usual diversitree code for a moment
	i <- match(unique(sapply(rels, as.character)), extra)	#match
	final <- c(extra[sort(i[!is.na(i)])], names.rhs)
	npar <- length(final)	
	#need to update the bounds at the same time the pars get updated
	bounds<-lapply(bounds,function(x) x[sapply(final,function(x) which(x==names))])
	## "free" are the parameters that have nothing special on their RHS
	## and are therefore passed directly through
	free <- setdiff(names.rhs, names(rels))
	free.i <- match(free, names) # index in full variables
	free.j <- match(free, final) # index in given variables.
	## Targets are processed in the same order as given by formulae.
	target.i <- match(names(rels), names)
	pars.out <- rep(NA, length(names))
	names(pars.out) <- names
	g <- function(pars, ..., pars.only=FALSE) {
		if ( length(pars) != npar ){
			stop(sprintf("Incorrect parameter length: expected %d, got %d",npar, length(pars)))
			}
	    pars.out[free.i] <- pars[free.j]
	    e <- structure(as.list(pars), names=final)
	    pars.out[target.i] <- unlist(lapply(rels, eval, e))
		if(pars.only){
			res<-pars.out
		}else{
			res<-f(pars.out, ...)
			}
		return(res)
		}
	class(g) <- c("constrained", class(f))
	attr(g, "parnames") <- final
	attr(g, "parbounds") <- bounds
	attr(g, "formulae") <- formulae
	attr(g, "extra") <- extra
	attr(g, "origfunction") <- f
	return(g)
	}
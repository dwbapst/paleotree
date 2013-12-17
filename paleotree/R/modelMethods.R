#a large number of model modifier functions, mostly having to do with parameters

#parnames

parnames <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parnames")
	}

parnames.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parnames")
	}

parnames.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parnames")
	}

`parnames<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parnames<-")
	}

`parnames<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parnames on a constrained function")
	}

`parnames<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	parnames<-attr(x,"parnames")
	value<-as.character(value)
	if(length(value)!=np){stop("length of new parnames not equal to number of parameters")}
	if(any(is.na(value))){stop("NA values in parnames replacement")}
	if(any(duplicated(value))){stop("Duplicated names in parnames replacement")}
	attr(x,"parnames")<-value
	return(x)
	}

#parbounds

parbounds <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parbounds")
	}

parbounds.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parbounds")
	}

parbounds.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parbounds")
	}

`parbounds<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parbound<-")
	}

`parbounds<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parbounds on a constrained function")
	}

`parbounds<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	if(!is.list(values) | !length(value)==2){stop("parbounds needs to be a list composed of two vectors")}
	lower<-as.numeric(value[[1]])
	if(length(lower)!=np){stop("length of new lower parbounds not equal to number of parameters")}
	if(any(is.na(lower))){stop("NA values in lower parbounds replacement")}
	upper<-as.numeric(value[[2]])
	if(length(upper)!=np){stop("length of new upper parbounds not equal to number of parameters")}
	if(any(is.na(upper))){stop("NA values in upper parbounds replacement")}
	attr(x,"parbounds")<-value
	return(x)
	}

#get upper and lower bounds for parameters

parLower <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parLower")
	}

parLower.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parbounds")[[1]]
	}

parLower.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parbounds")[[1]]
	}

`parLower<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parLower<-")
	}

`parLower<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parLower on a constrained function")
	}

`parLower<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	lower<-as.numeric(value)
	if(length(lower)!=np){stop("length of new lower parbounds not equal to number of parameters")}
	if(any(is.na(lower))){stop("NA values in lower parbounds replacement")}
	attr(x,"parbounds")[[1]]<-value
	return(x)
	}

parUpper <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parUpper")
	}

parUpper.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	attr(x, "parbounds")[[2]]
	}

parUpper.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
  	attr(x,"parbounds")[[2]]
	}

`parUpper<-` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parUpper<-")
	}

`parUpper<-.constrained` <- function(x, value){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	stop("Cannot set parUpper on a constrained function")
	}

`parUpper<-.paleotreeFunc` <- function(x, value) {
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	np<-attr(x,"np")	#number of parameters 
	#original uses base, not current number of params? I'm not following why...
	upper<-as.numeric(value)
	if(length(upper)!=np){stop("length of new upper parbounds not equal to number of parameters")}
	if(any(is.na(upper))){stop("NA values in upper parbounds replacement")}
	attr(x,"parbounds")[[2]]<-value
	return(x)
	}

#initial parameter values

parInit <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	UseMethod("parInit")
	}

parInit.constrained <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	res<-(attr(x, "parbounds")[[2]]-attr(x, "parbounds")[[1]])/2
	#infinite bounds probably too far from actual param value; use lower bound instead
	if(any(is.infinite(res))){res<-attr(x, "parbounds")[[1]]}
	return(res)
	}

parInit.paleotreeFunc <- function(x, ...){
	#based on Rich FitzJohn's argnames function for diversitree 10-22-13
	res<-(attr(x, "parbounds")[[2]]-attr(x, "parbounds")[[1]])/2
	#infinite bounds probably too far from actual param value; use lower bound instead
	if(any(is.infinite(res))){res<-attr(x, "parbounds")[[1]]}
	return(res)
	}
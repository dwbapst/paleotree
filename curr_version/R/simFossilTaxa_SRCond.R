simFossilTaxa_SRCond<-function(r,avgtaxa,p,q,anag.rate=0,prop.bifurc=0,prop.cryptic=0,nruns=1,mintime=1,maxtime=1000,
	minExtant=0,maxExtant=NULL,count.cryptic=FALSE,print.runs=FALSE,plot=FALSE){
	#wrapper for simulating clades large enough for 
		#getting some avg number of taxa under a given sampling parameter
	#for sampling: proper conditioning on number of taxa
		#how many taxa are needed in true clade to sample on average X taxa??
		#avgtaxa= average number of taxa you want to recover
		#mintaxa and maxtaxa will be set within 20% +/- of values necc for given avgtaxa
			#These should work well if avgtaxa is large (~50 or so)
	#simFossilTaxa_SRCond(r=0.1,p=0.1,q=0.1,nruns=10,avgtaxa=50,maxExtant=0)
	N<-avgtaxa/(r/(r+(q+anag.rate+(prop.bifurc*p))))
	results<-simFossilTaxa(p=p,q=q,anag.rate=anag.rate,prop.bifurc=prop.bifurc,prop.cryptic=prop.cryptic,nruns=nruns,mintaxa=N,
			maxtaxa=2*N,mintime=mintime,maxtime=maxtime,minExtant=minExtant,maxExtant=maxExtant,
			count.cryptic=FALSE,print.runs=print.runs,plot=plot)
	#if(nruns==1){results<-results[[1]]}
	return(results)
	}
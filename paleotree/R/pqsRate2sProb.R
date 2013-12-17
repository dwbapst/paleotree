pqsRate2sProb<-function(r,p,q,int.length=1){
	#A more accurate estimat of R given r, p and q
	#assuming p,q,r are constant and the timespan is infinte
		#dt is interval length for R
	#USES equations 26-29 from appendix to Foote (2000)
	#prob of samp for lineages that cross both boundaries
		#note typo in Foote (2000), eq 26, corrected version below
	dt<-int.length
	PDbt<-function(r,dt){1-exp(-r*dt)}
	#prob of samp for lineages that only cross bottom boundary
	PDbL<-function(q,r,dt){
		(((r+(q*exp(-(q+r)*dt)))/(q+r))-exp(-q*dt))/(1-exp(-q*dt))
		}
	#prob of samp for lineages that only cross upper boundary
	PDFt<-function(p,r,dt){
		(((r+(p*exp(-(p+r)*dt)))/(p+r))-exp(-p*dt))/(1-exp(-p*dt))
		}
	#prob of samp for lineages that cross neither boundary
		#29b corrected with addition sign!
	PDFL<-function(p,q,r,dt){
		if(p==q){
			NbNFL<-1/(exp(-q*dt)+(p*dt)-1)		#N(b)/N(FL) based on eq 1b and 6b
			term1<-(r*dt)/(p+r)				#first term in square brackets in eq 29b
			term2<-(1-exp(-p*dt))/p				#second term
			term3<-(p*(1-exp(-(p+r)*dt)))/((p+r)^2)	#third term
			terms<-term1-term2+term3			#full terms in square brackets
			res<-(NbNFL)*p*terms					#P(D|FL)
		}else{
			NbNFL<-1/(((q*exp((p-q)*dt))+((p-q)*exp(-q*dt))-p)/(p-q))
			term1<-(p*r*(exp((p-q)*dt)-1))/((q+r)*(p-q))
			term2<-(p*q*exp(-(q+r)*dt)*(exp((p+r)*dt)-1))/((p+r)*(q+r))
			term3<-exp(-q*dt)*(exp(p*dt)-1)
			terms<-term1+term2-term3
			res<-(NbNFL)*terms
			}
		res
		}
	#need to weight the PDs by the P of those taxon classes
		#use N equations from Foote (2000), relative to Nb to be probs
	Pbt<-exp(-q*dt)	
	PbL<-(1-exp(-q*dt))
	PFt<-exp((p-q)*dt)*(1-exp(-p*dt))
	if(p==q){PFL<-exp(-q*dt)+(p*dt)-1
		}else{PFL<-((q*exp((p-q)*dt))+((p-q)*exp(-q*dt))-p)/(p-q)}
	res<-sum(PDbt(r,dt)*Pbt,PDbL(q,r,dt)*PbL,
		PDFt(p,r,dt)*PFt,PDFL(p,q,r,dt)*PFL)
	names(res)<-NULL
	return(res)
	}
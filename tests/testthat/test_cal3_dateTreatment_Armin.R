# testing consistency issues with date treatment under cal3 (05-31-19)
	# making sure date treatment issues as reported by Armin Elsler in October 2018
	# are avoided in current code base 

test_that("cal3 datetreatment is consistent with arguments", {

#library(paleotree)

###########################################################

test_identical_tip_ages <- function(trees){
	# a function to test if tip ages are different
	#
	# checks
	if(!inherits(trees,"multiPhylo")){
		stop("trees is not class multiphylo")
		}
	if(length(trees)<2){
		stop("fewer than two trees given")
		}
	###########
	# get tip ages
	node_dates <- lapply(trees,dateNodes)
	tip_dates <- lapply(1:length(node_dates),function(i){
		x <- node_dates[[i]]
		x <- x[1:Ntip(trees[[i]])]
		names(x) <- trees[[i]]$tip.label
		x <- x[order(names(x))]
		return(x)
		})
	# test if tip ages are nearly the same
	all_nearly_identical <- all(sapply(tip_dates,
		function(x) isTRUE(all.equal(x,tip_dates[[1]]))
		))
	return(all_nearly_identical)
	}

################################

#Simulate some fossil ranges with simFossilRecord
set.seed(444)
record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
                        nTotalTaxa=c(30,40), nExtant=0)
taxa<-fossilRecord2fossilTaxa(record)
#simulate a fossil record with imperfect sampling with sampleRanges
rangesCont <- sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram <- taxa2cladogram(taxa,plot=TRUE)

# set up a set of ranges with identical FADs and LADs
	# just repeat FADs
rangesCont_pointocc <- rangesCont
rangesCont_pointocc[,2] <- rangesCont[,1]

#this library allows one to use
	# rate calibrated type time-scaling methods (Bapst, 2014.)
#to use these, we need an estimate of the sampling rate
# (we set it to 0.5 above)
likFun<-make_durationFreqCont(rangesCont)
srRes<-optim(parInit(likFun),likFun,
	lower=parLower(likFun),
	upper=parUpper(likFun),
    method="L-BFGS-B",
	control=list(maxit=1000000)
	)
	
sRate <- srRes[[1]][2]

# we also need extinction rate and branching rate
# we can get extRate from getSampRateCont too
#we'll assume extRate=brRate (ala Foote et al., 1999);
# may not always be a good assumption
divRate<-srRes[[1]][1]


##################################################################################
# 1) What is the behavior of `cal3TimePaleoPhy` when each taxon is supposed
	# to be a single occurrence with a precise age? 
# This would be the default usage with `firstLast`, with precise
	# first and last appearance times being identical (`rangesCont_pointocc`),
	# and `FAD.only = FALSE`.
# - Expectation: Tip ages shouldn't change from one tree to another,
	# as there is no tip age uncertainty expressed.

# #cal3TimePaleoPhy method using "firstLast" 
	# ancestors excluded
	# FAD.only = FALSE
	# FADs and LADs are identical (point occurrences)
ttrees_1 <- cal3TimePaleoPhy(
	cladogram, rangesCont_pointocc,
	brRate=divRate, extRate=divRate,
	sampRate=sRate, dateTreatment="firstLast",
	FAD.only=FALSE, ntrees=2,
	anc.wt=0,plot=FALSE)

expect_true(	
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_1)
	)

# Result
# Expectation met, treatment 'firstLast' has no impact on resulting tip ages.

##################################################################################
# 2) Is this behavior identical to `cal3TimePaleoPhy` with `minMax`,
	# using identical minimum and maximum ages?
# This would be the default usage with `firstLast`, with precise
	# first and last appearance times being identical(`rangesCont_pointocc`),
	# and `FAD.only = FALSE`.
# - Expectation: Tip ages shouldn't change from one tree to another,
	# as (again) there is no tip age uncertainty expressed.

# #cal3TimePaleoPhy method using "minMax" 
	# ancestors excluded
	# FAD.only = FALSE
	# FADs and LADs are identical (point occurrences)
ttrees_2<- cal3TimePaleoPhy(
	cladogram, rangesCont_pointocc,
	brRate=divRate, extRate=divRate,
	sampRate=sRate, dateTreatment="minMax",
	FAD.only=FALSE, ntrees=2,
	anc.wt=0,plot=FALSE)

expect_true(		
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_2)
	)

# RESULT
# Expectation met, treatment 'minMax' has no impact on resulting tip ages.

#################################################################################
# 3a) What is the difference in behavior of 
	# `cal3TimePaleoPhy` with `firstLast` versus `minMax`?
	# when the oldest and youngest ages given for each taxon are non-identical?
# - Expectation: Tip ages for `firstLast` should always be equal
	# to the given last appearance times, with no variation across analyses.
	
# #cal3TimePaleoPhy method using "firstLast" 
	# ancestors excluded
	# FAD.only = FALSE
	# FADs and LADs can differ
ttrees_3a <- cal3TimePaleoPhy(
	cladogram, rangesCont,
	brRate=divRate, extRate=divRate,
	sampRate=sRate, dateTreatment="firstLast",
	FAD.only=FALSE, ntrees=2,
	anc.wt=0,plot=FALSE)

expect_true(		
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_3a)
	)

# Result
# Expectation met, treatment 'firstLast' has no impact on resulting tip ages.

#################################################################################
# 3b) What is the difference in behavior of 
	# `cal3TimePaleoPhy` with `firstLast` versus `minMax`?
	# when the oldest and youngest ages given for each taxon are non-identical?
# - Expectation: Tip ages for `minMax` should be somewhere between the given
	# minimum and maximum bounds, with variation across runs.

# #cal3TimePaleoPhy method using "minMax" 
	# ancestors excluded
	# FAD.only = FALSE
	# FADs and LADs can differ
ttrees_3b <- cal3TimePaleoPhy(
	cladogram, rangesCont,
	brRate=divRate, extRate=divRate,
	sampRate=sRate, dateTreatment="minMax",
	FAD.only=FALSE, ntrees=2,
	anc.wt=0,plot=FALSE)

expect_false(		
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_3b)
	)

# RESULT
# Tip ages DO differ among the trees, as expected. Huzzah!

################################################################################
# 4a) When `FAD.only = TRUE`, what is the difference in behavior
	#  of `cal3TimePaleoPhy` with `firstLast` versus `minMax`?
	# With FADs and LADs *not* being identical.
# - Expectation: Tip ages for `firstLast` should always be equal
	# to the given first appearance times, with no variation across analyses.

# #cal3TimePaleoPhy method using "firstLast" 
	# ancestors excluded
	# FAD.only = TRUE
	# FADs and LADs can differ
ttrees_4a <- cal3TimePaleoPhy(
	cladogram, rangesCont,
	brRate=divRate, extRate=divRate,
	sampRate=sRate, dateTreatment="firstLast",
	FAD.only=TRUE, ntrees=2,
	anc.wt=0,plot=FALSE)

expect_true(	
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_4a)
	)

# Result
# Expectation met, treatment 'firstLast' has no impact on resulting tip ages.

################################################################################
# 4b) When `FAD.only = TRUE`, what is the difference in behavior
	#  of `cal3TimePaleoPhy` with `firstLast` versus `minMax`?
	# With FADs and LADs *not* being identical.
# - Expectation: Tip ages for `minMax` should be somewhere between the given
	# minimum and maximum bounds, with variation across runs.

expect_error(	
	# #cal3TimePaleoPhy method using "minMax" 
		# ancestors excluded
		# FAD.only = TRUE
		# FADs and LADs can differ
	ttrees_4b<- cal3TimePaleoPhy(
		cladogram, rangesCont,
		brRate=divRate, extRate=divRate,
		sampRate=sRate, dateTreatment="minMax",
		FAD.only=TRUE, ntrees=2,
		anc.wt=0,plot=FALSE)
	)

expect_error(
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_4b)
	)

# RESULT
# Returns an error, as it finds minMax incongruent with FAD.only = TRUE.

###########################################################################
# 5a) What is the difference in behavior with `randObs` between 
	# `FAD.only = FALSE` and `FAD.only = TRUE` ?
	# With FADs and LADs *not* being identical.
# - Expectation: When `FAD.only` is `FALSE`, tip ages from `randObs`
	# should be somewhere between the given minimum and maximum bounds,
	# with variation across runs, much like `minMax`. 

# #cal3TimePaleoPhy method using "randObs" 
	# ancestors excluded
	# FAD.only = FALSE
	# FADs and LADs can differ
ttrees_5a<- cal3TimePaleoPhy(
	cladogram, rangesCont,
	brRate=divRate, extRate=divRate,
	sampRate=sRate, dateTreatment="randObs",
	FAD.only=FALSE, ntrees=2,
	anc.wt=0,
	plot=FALSE
	)

expect_false(	
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_5a)
	)

# RESULT
# Expectation met, there is variation as expected.

###########################################################################
# 5b) What is the difference in behavior with `randObs` between 
	# `FAD.only = FALSE` and `FAD.only = TRUE` ? 
	# With FADs and LADs *not* being identical.
# - Expectation: When `FAD.only` is `TRUE`, tip agesfrom `randObs`
	# should be equal to the given first appearance times, with no variation
	# across runs, as with `firstLast`.

# #cal3TimePaleoPhy method using "randObs" 
	# ancestors excluded
	# FAD.only = TRUE
	# FADs and LADs can differ
expect_error(
	ttrees_5b<- cal3TimePaleoPhy(
		cladogram, rangesCont,
		brRate=divRate, extRate=divRate,
		sampRate=sRate, dateTreatment="randObs",
		FAD.only=TRUE, ntrees=2,
		anc.wt=0,
		plot=FALSE
		)
	)

expect_error(	
	# test if tip ages have changed
	test_identical_tip_ages(trees = ttrees_5b)
	)

# RESULT
# Returns an error, as it finds randObs incongruent with FAD.only = TRUE.

})

\name{footeValues}
\alias{footeValues}
\title{Calculates Values for Foote's Inverse Survivorship Analyses}
\usage{
footeValues(p, q, r, PA_n = 0, PB_1 = 0, p_cont = TRUE, q_cont = TRUE,
  Nb = 1)
}
\arguments{
  \item{p}{Instantaneous origination/branching rate of taxa
  under a continuous model; under a pulsed mode
  (p_cont=FALSE), a per-interval probability instead. Given
  as a vector with length equal to the number of intervals,
  so a different value may be given for each separate
  interval. Must be the same length as q and r.}

  \item{q}{Instantaneous extinction rate of taxa under a
  continuous model; under a pulsed mode (q_cont=FALSE), a
  per-interval 'rate' instead, which cannot exceed 1. Given
  as a vector with length equal to the number of intervals,
  so a different value may be given for each separate
  interval. Must be the same length as p and r.}

  \item{r}{Instantaneous sampling rate of taxa. Given as a
  vector with length equal to the number of intervals, so a
  different value may be given for each separate interval.
  Must be the same length as p and q.}

  \item{PA_n}{The probability of sampling a taxon before
  the first interval included in a survivorship study.
  Usually zero.}

  \item{PB_1}{The probability of sampling a taxon after the
  last interval included in a survivorship study. Usually
  zero for extinct groups, although more logically has the
  value of 1 when there are still extant taxa (i.e., if the
  last interval is the Holocene and the group is still
  alive, the probability of sampling them later is probably
  1...).}

  \item{p_cont}{If TRUE (the default), then origination is
  assumed to be a continuous time process with an
  instantaneous rate. If FALSE, the origination is treated
  as a pulsed discrete-time process with a probability.}

  \item{q_cont}{If TRUE (the default), then extinction is
  assumed to be a continuous time process with an
  instantaneous rate. If FALSE, the extinction is treated
  as a pulsed discrete-time process with a probability.}

  \item{Nb}{The number of taxa that enter an interval (b is
  for 'bottom'). This is an arbitrary constant used to
  scale other values in these calculations and can be
  safely set to 1.}
}
\value{
Returns a matrix with number of rows equal to the number of
intervals (i.e. the length of p, q and r) and named columns
representing the different values calculated by the
function: "Nb", "Nbt", "NbL", "NFt", "NFL", "PD_bt",
"PD_bL", "PD_Ft", "PD_FL", "PA", "PB", "Xbt", "XbL", "XFt"
and "XFL".
}
\description{
This function calculates the intermediary values needed for
fitting Foote's inverse survivorship analyses, as listed in
the table of equations in Foote (2003), with the analyses
themselves described further in Foote (2001) and Foote
(2005).
}
\details{
Although most calculations in this function agree with the
errata for Foote's 2003 table (see references), there were
some additional corrections for Prob(D|FL) made as part of
a personal communication in 2013 between the package author
and Michael Foote.
}
\examples{
#very simple example with three intervals, same value for all parameters
rate<-rep(0.1,3) 			#example rates (for the most part)
footeValues(rate,rate,rate)	#all continuous
footeValues(rate,rate,rate,p_cont=FALSE)	#origination pulsed
footeValues(rate,rate,rate,q_cont=FALSE) #extinction pulsed
footeValues(rate,rate,rate,p_cont=FALSE,q_cont=FALSE) #all pulsed
export
}
\author{
David W. Bapst, with advice from Michael Foote.
}
\references{
Foote, M. 2001. Inferring temporal patterns of
preservation, origination, and extinction from taxonomic
survivorship analysis. \emph{Paleobiology} 27(4):602-630.

Foote, M. 2003. Origination and Extinction through the
Phanerozoic: A New Approach. \emph{The Journal of Geology}
111(2):125-148.

Foote, M. 2003. Erratum: Origination and Extinction through
the Phanerozoic: a New Approach. \emph{The Journal of
Geology} 111(6):752-753.

Foote, M. 2005. Pulsed origination and extinction in the
marine realm. \emph{Paleobiology} 31(1):6-20.
}


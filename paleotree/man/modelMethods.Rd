\name{modelMethods}
\alias{modelMethods}
\alias{parInit}
\alias{parInit.constrained}
\alias{parInit.paleotreeFunc}
\alias{parLower}
\alias{parLower.constrained}
\alias{parLower.paleotreeFunc}
\alias{parLower<-}
\alias{parLower<-.constrained}
\alias{parLower<-.paleotreeFunc}
\alias{parUpper}
\alias{parUpper.constrained}
\alias{parUpper.paleotreeFunc}
\alias{parUpper<-}
\alias{parUpper<-.constrained}
\alias{parUpper<-.paleotreeFunc}
\alias{parbounds}
\alias{parbounds.constrained}
\alias{parbounds.paleotreeFunc}
\alias{parbounds<-}
\alias{parbounds<-.constrained}
\alias{parbounds<-.paleotreeFunc}
\alias{parnames}
\alias{parnames.constrained}
\alias{parnames.paleotreeFunc}
\alias{parnames<-}
\alias{parnames<-.constrained}
\alias{parnames<-.paleotreeFunc}
\title{Model Function Methods: Parameter Names, Bounds and Initial Values}
\usage{
parnames(x, ...)

\method{parnames}{paleotreeFunc}(x, ...)

\method{parnames}{constrained}(x, ...)

parnames(x) <- value

\method{parnames}{constrained}(x) <- value

\method{parnames}{paleotreeFunc}(x) <- value

parbounds(x, ...)

\method{parbounds}{paleotreeFunc}(x, ...)

\method{parbounds}{constrained}(x, ...)

parbounds(x) <- value

\method{parbounds}{constrained}(x) <- value

\method{parbounds}{paleotreeFunc}(x) <- value

parLower(x, ...)

\method{parLower}{constrained}(x, ...)

\method{parLower}{paleotreeFunc}(x, ...)

parLower(x) <- value

\method{parLower}{constrained}(x) <- value

\method{parLower}{paleotreeFunc}(x) <- value

parUpper(x, ...)

\method{parUpper}{constrained}(x, ...)

\method{parUpper}{paleotreeFunc}(x, ...)

parUpper(x) <- value

\method{parUpper}{constrained}(x) <- value

\method{parUpper}{paleotreeFunc}(x) <- value

parInit(x, ...)

\method{parInit}{constrained}(x, ...)

\method{parInit}{paleotreeFunc}(x, ...)
}
\arguments{
  \item{x}{A function of S3 class 'paleotreeFunc' with all
  necessary attributes expected of that class, which
  include parameter names and upper and lower bounds. As I
  have deliberately not exported the function which creates
  this class, it should be impossible for regular users to
  obtain such objects easily without using one of the
  'make' functions, which automatically output a function
  of the appropriate class and attributes.}

  \item{...}{'Ignored arguments to future methods' (i.e.
  for diversitree). Kept here only so constrainParPaleo is
  kept as close to the parent method in diversitree as
  possible.}

  \item{value}{The new value with which to replace the
  parameter names or bounds. Must be a vector of the same
  length as the number of parameters. For parbounds, must
  be a list composed of two vectors.}
}
\value{
Returns the sought parameter names, bounds or initial
values or (for the replacement methods) returns a modified
function with the respective attributes altered.
}
\description{
A large number of functions for obtaining and modifying the
parameters of likelihood models made in paleotree. These
functions allow users to obtain or set parameter names, or
obtain and set parameter bounds, both of which are treated
as an attribute of the function class used by paleotree. In
practice, this allows users to quickly obtain parameter
names and upper and lower values for use in bounded
optimizers, including reasonable starting values.
}
\details{
Parameter names cannot be changed for a constrained
function.

The parInit function calls the bounds for each parameter
and gives the midway value for each, giving a users a
shorthand to quickly generate initial parameter values
which are within the set bounds, for use in functions such
as \code{\link{optim}}. These values may be too high (due
to the liberal upper bounds I set for parameters in many of
the likelihood functions) and thus users should always try
slightly different values to see if the resulting maximum
likelihood parameter values change.

As parInit depends on the upper and lower bounds attribute,
no function is offered to allow it to be replaced (as there
is nothing to replace!).
}
\examples{
#example with make_durationFreqCont
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
rangesCont <- sampleRanges(taxa,r=0.5)
likFun <- make_durationFreqCont(rangesCont)

#get parameter names
parnames(likFun)

#get the bounds for those parameters
parbounds(likFun)

#can also get these seperately
parLower(likFun)
parUpper(likFun)

#initial parameter values
parInit(likFun)   #arbitrary midway value between par bounds

#can then use these in optimizers, such as optim with L-BFGS-B
#see the example for make_durationFreqCont

#renaming parameter names
likFun2 <- likFun
parnames(likFun2) <- c("extRate","sampRate")
parnames(likFun2)
#test if reset correctly
parnames(likFun2)==c("extRate","sampRate")
#also works for constrained functions
constrainFun<-constrainParPaleo(likFun,q.1~r.1)
parnames(constrainFun)
#also modified the parameter bounds, see!
parbounds(constrainFun)
parInit(constrainFun)
#but cannot rename parameter for constrained function!
}
\author{
These functions are strongly based on or inspired by the
\code{argnames} functions provided for handling models in
Rich Fitzjohn's library \code{diversitree}, but the
functions presented here are derviations written by David
Bapst.
}
\seealso{
These model methods were introduced to interact with the
new model framework introduced in paleotree v1.9, in
particular to interface with
\code{\link{constrainParPaleo}}.
}


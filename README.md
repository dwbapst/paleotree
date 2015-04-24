paleotree: Paleontological and Phylogenetic Analyses of Evolution

paleotree is an R library for analysing, time-scaling and simulating phylogenies of extinct/fossil lineages. In particular, most users are interested in the functions timePaleoPhy, bin_timePaleoPhy, cal3TimePaleoPhy and bin_cal3TimePaleoPhy, which time-scales cladograms of fossil taxa into dated phylogenies.

This package also contains a large number of functions associated with simulating diversification in incompletely sampled fossil records, under various models of morphological differentiation, and translating such simulated fossil records into their respective phylogenies. There are also likelihood functions offered for estimating sampling rates from different types of observables available from fossil record data.

This code is authored by David Bapst and offered under CC0. The most recent public release of the code is on CRAN at:

	http://cran.r-project.org/web/packages/paleotree/index.html

You can install the most recent public release version of paleotree in R using:

	install.packages("paleotree")

The latest pre-release version of paleotree can be found at github:

	https://github.com/dwbapst/paleotree
	
The most recent commit is currently: [![Build Status](https://travis-ci.org/dwbapst/paleotree.svg?branch=master)](https://travis-ci.org/dwbapst/paleotree) (Travis CI)
	
You can install this latest development version using the R function install_github in the package 'devtools':

	library(devtools)
	install_github("dwbapst/paleotree")

Once installed, you can check the version number of your paleotree install using the R function packageVersion:

	packageVersion("paleotree")

As of 01/13/15, the paleotree repository was restructured so the package no longer sat within a subdirectory. An old directory composed of deprecated, previously public functions (more than a year old) was removed.

This code is authored by David Bapst and offered under GPL >2.0.

The current total number of downloads of the paleotree package from the RStudio CRAN mirror is: [![Number of Downloads](http://cranlogs.r-pkg.org/badges/paleotree)]

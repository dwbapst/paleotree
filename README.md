paleotree: Paleontological and Phylogenetic Analyses of Evolution

paleotree is an R package for transforming, time-scaling, and modifying phylogenies containing extinct (i.e. fossil) lineages. In particular, most users are interested in the functions timePaleoPhy, bin_timePaleoPhy, cal3TimePaleoPhy and bin_cal3TimePaleoPhy, which time-scale cladograms of fossil taxa into dated phylogenies. This package also contains a large number of likelihood functions for estimating sampling and diversification rates from different types of data available from the fossil record (e.g. range data, occurrence data, etc). paleotree users can also simulate diversification and sampling in the fossil record using the function simFossilRecord, which is a detailed simulator for branching birth-death-sampling processes composed of discrete taxonomic units arranged in ancestor-descendant relationships. Users can use simFossilRecord to simulate diversification in incompletely sampled fossil records, under various models of morphological differentiation (i.e. the various patterns by which morphotaxa originate from one another), and with time-dependent, longevity-dependent and/or diversity-dependent rates of diversification, extinction and sampling. Additional functions allow users to translate simulated ancestor-descendant data from simFossilRecord into standard time-scaled phylogenies or unscaled cladograms that reflect the relationships among taxon units.

The most recent public release of the code is on CRAN at:

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

This code is authored by David Bapst and offered under CC0.

The current total number of downloads of the paleotree package from the RStudio CRAN mirror is: [![Number of Downloads](http://cranlogs.r-pkg.org/badges/grand-total/paleotree)](https://github.com/metacran/cranlogs.app)

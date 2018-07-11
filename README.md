## paleotree: Paleontological and Phylogenetic Analyses of Evolution

paleotree is an R package for transforming, 'a posteriori' time-scaling, and modifying phylogenies containing extinct (i.e. fossil) lineages. In particular, most users are interested in the functions timePaleoPhy, bin_timePaleoPhy, cal3TimePaleoPhy and bin_cal3TimePaleoPhy, which a posteriori time-scale cladograms of fossil taxa into dated phylogenies. This package also contains a large number of likelihood functions for estimating sampling and diversification rates from different types of data available from the fossil record (e.g. range data, occurrence data, etc). paleotree users can also simulate diversification and sampling in the fossil record using the function simFossilRecord, which is a detailed simulator for branching birth-death-sampling processes composed of discrete taxonomic units arranged in ancestor-descendant relationships. Users can use simFossilRecord to simulate diversification in incompletely sampled fossil records, under various models of morphological differentiation (i.e. the various patterns by which morphotaxa originate from one another), and with time-dependent, longevity-dependent and/or diversity-dependent rates of diversification, extinction and sampling. Additional functions allow users to translate simulated ancestor-descendant data from simFossilRecord into standard time-scaled phylogenies or unscaled cladograms that reflect the relationships among taxon units.

The most recent public release of the code is on CRAN at:

	[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/paleotree)](https://cran.r-project.org/package=paleotree)

You can install the most recent public release version of `paleotree` in R from CRAN using:

```
install.packages("paleotree")
```
	
The latest pre-release version of `paleotree` can be found at github:

	https://github.com/dwbapst/paleotree
	
The most recent commit is currently: [![Build Status](https://travis-ci.org/dwbapst/paleotree.svg?branch=master)](https://travis-ci.org/dwbapst/paleotree) (Travis CI)
	
You can install this latest development version using the R function `install_github` in the package `devtools`:

```
devtools::install_github("dwbapst/paleotree")
```
	
The above command will install the **master** branch (as the default option), which should be the most recent stable release (ideally the CRAN release, and thus should be identical to installing `paleotree` from CRAN). There will also generally be an in-development branch where the code may be in a state of change, with new features and bug-fixes; for example right now the development branch is named **developmentBranch** and can be specified with the `ref = ` argument:

```
devtools::install_github("dwbapst/paleotree", ref="developmentBranch")
```	
	
Once installed, you can check the version number of your paleotree install using the R function packageVersion:

```
packageVersion("paleotree")
```

If you use paleotree in your research, you can cite my paper describing paleotree in Methods in Ecology and Evolution:

	Bapst, D.W. 2012. paleotree: an R package for paleontological and phylogenetic analyses of evolution. Methods in Ecology and Evolution. 3: 803-807. doi: 10.1111/j.2041-210X.2012.00223.x
	
	http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2012.00223.x/abstract
	
You can also call the citation for paleotree from within R, using the citation function:
	
```
citation("paleotree")
```
	
This code is mainly authored by David Bapst, with some functions authored by Peter Wagner, and offered under CC0.

The current total number of downloads of the paleotree package from the RStudio CRAN mirror is: [![Number of Downloads](http://cranlogs.r-pkg.org/badges/grand-total/paleotree)](https://github.com/metacran/cranlogs.app)

[![Research software impact](http://depsy.org/api/package/cran/paleotree/badge.svg)](http://depsy.org/package/r/paleotree)

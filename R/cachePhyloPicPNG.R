## Cache phylopic images for a tree
## 
## This function takes a tree and saves phylopic thumbnails locally,
## in order to dramatically speed up later rendering.
## 

## @param tree A phylogeny of class \code{phylo} which will be
## plotted, with the terminal tip taxa replaced by silhouettes.
## The tree will be plotted with edge lengths.

## @param taxaDataPBDB  A \code{data.frame} of taxonomic data from
## the Paleobiology Database containing an \code{$image_no} variable,
## as returned when \code{show = "img"} is used. See \emph{Details}.

## @param cacheDir Where to save the output PNGs

## @export

cachePhyloPicPNG <- function(
		tree, taxaDataPBDB=tree$taxaDataPBDB,
		cacheDir
		){
	#################
	ids <- getPhyloPicIDNum(taxaData=taxaDataPBDB, tree=tree)
	for(i in seq_along(ids)) {
		picPNG <- getPhyloPicPNG(ids[i])
		png::writePNG(picPNG,
			target=file.path(cacheDir,
				paste0(ids[i], ".png")
				)
			)
	}
}


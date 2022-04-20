#' Creating a Taxon-Tree from Taxonomic Data Downloaded from the Paleobiology Database
#' 
#' The function \code{makePBDBtaxonTree} creates phylogeny-like 
#' object of class \code{phylo} from the taxonomic information
#' recorded in a taxonomy download from the PBDB for
#' a given group. Two different algorithms are provided,
#' the default being based on parent-child taxon relationships,
#' the other based on the nested Linnean hierarchy. The function
#' \code{plotTaxaTreePBDB} is also provided as a minor helper
#' function for optimally plotting the labeled topologies that are
#' output by \code{makePBDBtaxonTree}.
#' 
 
#' @details
#' This function should not be taken too seriously.
#' Many groups in the Paleobiology Database have
#' out-of-date or very incomplete taxonomic information.
#' This function is meant to help visualize
#' what information is present, and by use of time-scaling
#' functions, allow us to visualize the intersection
#' of temporal and phylogenetic, mainly to look for incongruence
#' due to either incorrect taxonomic placements,
#' erroneous occurrence data or both. 
#' 
#' Note however that, contrary to common opinion among some
#' paleontologists, taxon-trees may be just as useful for 
#' macroevolutionary studies as reconstructed phylogenies
#' (Soul and Friedman, 2015).

#' @param taxaDataPBDB A table of taxonomic data collected from
#' the Paleobiology Database, using the taxa list option
#' with \code{show = class}. Should work with versions 1.1-1.2 of
#' the API, with either the \code{pbdb} or \code{com} vocab. However,
#' as \code{accepted_name} is not available in API v1.1, the resulting
#' tree will have a taxon's *original* name and not
#' any formally updated name.

#' @param rankTaxon The selected taxon rank; must be one of \code{'species'},
#' \code{'genus'}, \code{'family'}, \code{'order'}, \code{'class'} or \code{'phylum'}.

#' @param cleanTree When \code{TRUE} (the default), the tree is run through a series of
#' post-processing, including having singles collapsed,
#' nodes reordered and being written out as a Newick string and read
#' back in, to ensure functionality with ape functions
#' and ape-derived functions. 
#' If \code{FALSE}, none of this post-processing is done and
#' users should beware, as such trees can lead to hard-crashes of R.

#' @param method Controls which algorithm is used for calculating
#' the taxon-tree. The default option is \code{method  = "parentChild"}
#' which converts the listed binary parent-child taxon relationships in
#' the Paleobiology Database- these parent-child relationships (if missing
#' from the input dataset) are autofilled using API calls to the
#' Paleobiology Database. Alternatively, users may use
#' \code{method = "Linnean"}, which converts the table of Linnean taxonomic
#' assignments (family, order, etc as provided by \code{show = class} in
#' PBDB API calls) into a taxon-tree. Two methods formerly both implemented
#' under \code{method  = "parentChild"} are also available as
#' \code{method = "parentChildOldMergeRoot"} and \code{method = "parentChildOldQueryPBDB"}
#' respectively. Both of these use similar algorithms as the current
#' \code{method  = "parentChild"} but differ in how they treat taxa with
#' parents missing from the input taxonomic dataset.
#' \code{method = "parentChildOldQueryPBDB"} behaves most similar
#' to \code{method = "parentChild"}  in that it queries the Paleobiology
#' Database via the API , but repeatedly does so for information on parent
#' taxa of the 'floating' parents, and continues within a \code{while}
#' loop until only one such unassigned parent taxon remains. This latter
#' option may talk a long time or never finish, depending on the
#' linearity and taxonomic structures encountered in the PBDB taxonomic
#' data; i.e. if someone a taxon was ultimately its own indirect child
#' in some grand loop by mistake, then under this option
#' \code{makePBDBtaxonTree} might never finish. In cases where taxonomy
#' is bad due to weird and erroneous taxonomic assignments reported by
#' the PBDB, this routine may search all the way back to a very ancient
#' and deep taxon, such as the \emph{Eukaryota} taxon.
#' \code{method = "parentChildOldMergeRoot"} will combine these disparate
#' potential roots and link them to an artificially-constructed
#' pseudo-root, which at least allows for visualization of the taxonomic
#' structure in a limited dataset. This latter option will be fully
#' offline, as it does not do any additional API calls
#' of the Paleobiology Database, unlike other options.



#' @param tipSet This argument only impacts analyses where 
#' \code{method  = "parentChild"} is used. This \code{tipSet} argument controls
#' which taxa are selected as tip taxa for the output tree. 
#' \code{tipSet  = "nonParents"} selects all child taxa which
#' are not listed as parents in \code{parentChild}.
#' Alternatively, \code{tipSet = "all"} will add a tip to every
#' internal node with the parent-taxon name encapsulated in parentheses.
#' The default is \code{NULL} - if \code{tipSet = NULL} and \code{method  = "parentChild"},
#' then \code{tipSet} will be set so \code{tipSet = "nonParents"}.

#' @param APIversion Version of the Paleobiology Database API used by
#' \code{makePBDBtaxonTree} when \code{method  = "parentChild"} or
#' \code{method  = "parentChildOldQueryPBDB"} is used. The current default
#' is \code{APIversion = "1.2"}, the most recent API version as of 12/11/2018.

#' @param annotatedDuplicateNames A logical determining whether duplicate taxon names,
#' when found in the Paleobiology Database for taxa (presumably reflecting an issue with
#' taxa being obsolete but with incomplete seniority data), should be annotated to include
#' sequential numbers so to modify them, via function\code{base}'s
#' \code{\link[base]{make.unique}}. This only applies to
#' \code{method = "parentChild"}, with the default option being
#' \code{annotatedDuplicateNames = TRUE}. If more than 26 duplicates are found, an error
#' is issued. If this argument is \code{FALSE}, an error is issued if duplicate taxon
#' names are found.

#' @param taxaTree A phylogeny of class \code{phylo}, presumably a taxon tree as output from
#' \code{makePBDBtaxonTree} with higher-taxon names as node labels.

#' @param edgeLength The edge length that the plotted tree should be plotted
#' with (\code{plotTaxaTreePBDB} plots phylogenies as non-ultrametric,
#' not as a cladogram with aligned tips).

#' @return
#' A phylogeny of class \code{phylo}, where each tip is a taxon of the given \code{rankTaxon}. See additional details
#' regarding branch lengths can be found in the sub-algorithms used to create the taxon-tree by this function:
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}.
#' 
#' Depending on the \code{method}
#' used, either the element \code{$parentChild} or \code{$taxonTable} is added to the list structure of
#' the output phylogeny object, which was used as input for one of the two algorithms mentioned above.
#' 
#' Please note that when applied to output from the taxa option of the API version 1.1, the taxon names
#' returned are the \emph{original} taxon names as 'accepted_name' is not available in API v1.1, while
#' under API v1.2, the returned taxon names should be the most up-to-date formal names for those taxa.
#' Similar issues also effect the identification of parent taxa, as the accepted name of the
#' parent ID number is only provided in version 1.2 of the API.

#' @seealso
#' Two other functions in paleotree are used as sub-algorithms by \code{makePBDBtaxonTree}
#' to create the taxon-tree within this function,
#' and users should consult their manual pages for additional details:
#' 
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}
#' 
#' Closely related functions for 
#' 
#' Other functions for manipulating PBDB data can be found at \code{\link{taxonSortPBDBocc}},
#' \code{\link{occData2timeList}}, and the example data at \code{\link{graptPBDB}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.
#' 
#' Soul, L. C., and M. Friedman. 2015. Taxonomy and Phylogeny Can Yield
#' Comparable Results in Comparative Palaeontological Analyses. \emph{Systematic Biology} 
#' (\href{http://sysbio.oxfordjournals.org/content/early/2015/03/23/sysbio.syv015.abstract}{Link})


#' @examples
#' # Note that most examples here use argument 
#'     # failIfNoInternet = FALSE so that functions do
#'     # not error out but simply return NULL if internet
#'     # connection is not available, and thus
#'     # fail gracefully rather than error out (required by CRAN).
#' # Remove this argument or set to TRUE so functions DO fail
#'     # when internet resources (paleobiodb) is not available.
#' 
#' set.seed(1)
#' \donttest{
#' 
#' #get some example occurrence and taxonomic data
#' data(graptPBDB)
#' 
#' #get the taxon tree: Linnean method
#' graptTreeLinnean <- makePBDBtaxonTree(
#'     taxaDataPBDB = graptTaxaPBDB,
#'     rankTaxon = "genus",
#'     method = "Linnean", 
#'     failIfNoInternet = FALSE)
#' 
#' #get the taxon tree: parentChild method
#' graptTreeParentChild <- makePBDBtaxonTree(
#'     taxaDataPBDB = graptTaxaPBDB,
#'     rankTaxon = "genus",
#'     method = "parentChild", 
#'     failIfNoInternet = FALSE)
#'     
#' if(!is.null(graptTreeParentChild) & 
#'         !is.null(graptTreeLinnean)){
#'     # if those functions worked...
#'     # let's plot these and compare them! 
#'     plotTaxaTreePBDB(graptTreeParentChild)
#'     plotTaxaTreePBDB(graptTreeLinnean)
#'     }
#' 
#' 
#' # pause 3 seconds so we don't spam the API
#' Sys.sleep(3)
#' 
#' ####################################################
#' # let's try some other groups
#' 
#' ###################################
#' #conodonts
#' 
#' conoData <- getCladeTaxaPBDB("Conodonta")
#' 
#' conoTree <- makePBDBtaxonTree(
#'     taxaDataPBDB = conoData,
#'     rankTaxon = "genus",
#'     method = "parentChild", 
#'     failIfNoInternet = FALSE)
#' 
#' if(!is.null(conoTree)){    
#'     # if it worked, plot it!
#'     plotTaxaTreePBDB(conoTree)
#'     }
#' 
#' # pause 3 seconds so we don't spam the API
#' Sys.sleep(3)
#' 
#' #############################
#' #asaphid trilobites
#' 
#' asaData <- getCladeTaxaPBDB("Asaphida")
#' 
#' asaTree <- makePBDBtaxonTree(
#'     taxaDataPBDB = asaData,
#'     rankTaxon = "genus",
#'     method = "parentChild", 
#'     failIfNoInternet = FALSE)
#' 
#' if(!is.null(asaTree)){    
#'     # if it worked, plot it!
#'     plotTaxaTreePBDB(asaTree)
#'     }
#' 
#' # pause 3 seconds so we don't spam the API
#' Sys.sleep(3)
#' 
#' ###############################
#' #Ornithischia
#' 
#' ornithData <- getCladeTaxaPBDB("Ornithischia")
#' 
#' ornithTree <- makePBDBtaxonTree(
#'     taxaDataPBDB = ornithData,
#'     rankTaxon = "genus",
#'     method = "parentChild", 
#'     failIfNoInternet = FALSE)
#' 
#' if(!is.null(ornithTree)){    
#'     # if it worked, plot it!
#'     plotTaxaTreePBDB(ornithTree)
#'     }
#' 
#' # pause 3 seconds so we don't spam the API
#' Sys.sleep(3)
#' 
#' #try Linnean!
#' 
#' #but first... need to drop repeated taxon first: Hylaeosaurus
#'     # actually this taxon seems to have been repaired 
#'     # as of September 2019 !
#' # findHylaeo <- ornithData$taxon_name == "Hylaeosaurus"
#' # there's actually only one accepted ID number
#' # HylaeoIDnum <- unique(ornithData[findHylaeo,"taxon_no"])
#' # HylaeoIDnum 
#' # so, take which one has occurrences listed
#' # dropThis <- which((ornithData$n_occs < 1) & findHylaeo)
#' # ornithDataCleaned <- ornithData[-dropThis,]
#' 
#' ornithTree <- makePBDBtaxonTree(
#'     ornithData,
#'     rankTaxon = "genus",
#'     method = "Linnean", 
#'     failIfNoInternet = FALSE)
#' 
#' if(!is.null(ornithTree)){    
#'     # if it worked, plot it!
#'     plotTaxaTreePBDB(ornithTree)
#'     }
#' 
#' # pause 3 seconds so we don't spam the API
#' Sys.sleep(3)
#' 
#' #########################
#' # Rhynchonellida
#' 
#' rhynchData <- getCladeTaxaPBDB("Rhynchonellida")
#' 
#' rhynchTree <- makePBDBtaxonTree(
#'     taxaDataPBDB = rhynchData,
#'     rankTaxon = "genus",
#'     method = "parentChild", 
#'     failIfNoInternet = FALSE)
#' 
#' if(!is.null(rhynchTree)){    
#'     # if it worked, plot it!
#'     plotTaxaTreePBDB(rhynchTree)
#'     }
#' 
#' #some of these look pretty messy!
#' 
#' }
#' 


# # plot it!
# # let's make a simple helper function
#    # for plotting these taxon trees
# plotPBDBtaxonTree <- function(tree){
#    plot(tree,show.tip.label = FALSE,
#        no.margin = TRUE,edge.width = 0.35)
#    nodelabels(tree$node.label,adj = c(0,1/2))
#       }
# 


#' @name makePBDBtaxonTree
#' @rdname makePBDBtaxonTree
#' @export
makePBDBtaxonTree <- function(
                      taxaDataPBDB, 
                      rankTaxon,
                      method = "parentChild", 
                      #solveMissing = NULL,
                      tipSet = NULL, 
                      cleanTree = TRUE,
                      annotatedDuplicateNames = TRUE,
                      APIversion = "1.2",
                      failIfNoInternet = TRUE
                      ){        
    ############################################################
    ############################################################
    # library(paleotree);data(graptPBDB);
    # taxaDataPBDB <- graptTaxaPBDB; rankTaxon = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
    # taxaDataPBDB <- graptTaxaPBDB; rankTaxon = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
    # taxaDataPBDB <- graptTaxaPBDB; rankTaxon = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
    # taxaDataPBDB <- graptTaxaPBDB; rankTaxon = "genus"; method = "Linnean"; 
    #
    #CHECKS
    if(length(method) != 1 | !is.character(method)){
        stop("method must be a single character value")
        }
    if(!any(method == c("Linnean", "parentChild"))){
        stop("method must be one of either 'Linnean' or 'parentChild'")
        }
    #if(!is.null(solveMissing)){
    #    if(length(solveMissing)>1 | !is.character(solveMissing)){
    #        stop("solveMissing must be either NULL or a single character value")
    #        }
    #    if(is.na(match(solveMissing,c("queryPBDB","mergeRoots")))){
    #        stop('solveMissing but be either NULL or "queryPBDB" or "mergeRoots"')
    #        }
    #    }
    if(!is(taxaDataPBDB,"data.frame")){
        stop("taxaDataPBDB isn't a data.frame")
        }
    if(length(rankTaxon) != 1 | !is.character(rankTaxon)){
        stop("rankTaxon must be a single character value")
        }
    if(!any(sapply(c("species","genus","family","order","class","phylum"),
            function(x) x == rankTaxon))){
        stop("rankTaxon must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")
        }
    #
    #########################################
    # CLEAN DATA
    #
    #translate to a common vocabulary
    dataTransform <- translatePBDBtaxa(taxaDataPBDB)
    dataTransform <- unique(dataTransform)
    #
    if(method == "parentChild"){
        # need two things: a table of parent-child relationships as IDs
            #and a look-up table of IDs and taxon names
        # 
        if(is.null(tipSet)){
            tipSet <- "nonParents"
            }
        ##############################
        # FIND ALL PARENTS FIRST
            # three column matrix with taxon name, taxon ID, parent ID
            # (in that order)
        parData <- getAllParents(
            dataTransform,
            status="all", 
            annotatedDuplicateNames = annotatedDuplicateNames,
            convertAccepted = FALSE,
            stopIfSelfParent = FALSE,
            failIfNoInternet = failIfNoInternet
            )
        if(is.null(parData) & !failIfNoInternet){return(NULL)}
        #print(parData)
        #
        #######################################
        # NOW FILTER OUT TIP TAXA WE WANT
        #tipIDs <- getTaxaIDsDesiredRank(data=dataTransform, rank=rankTaxon)
        tipIDs <- dataTransform$taxon_no[dataTransform$taxon_rank==rankTaxon]
        # check that all of these are unique values
        tipIDs <- unique(tipIDs)
        # figure out which taxon numbers match tip IDs
        whichTip <- match(tipIDs, parData$taxon_no)        
        #
        ###############################
        # BUILD PARENT-CHILD matrix
        # get parent-child matrix for just desired OTUs 
            # start matrix with those parent-child relationships
            # subset these from parData using the taxon ID numbers
        pcMat <- subsetParDataPBDB(
            subsetNum = tipIDs,
            parData = parData
            )            
        # starting from desired tip OTUs, work backwards to a common ancestor from the full parData
        pcMat<-constructParentChildMatrixPBDB(
            initPCmat=pcMat,
            parData = parData
            )    
        #################################
        # convert parent-child matrix to accepted taxon names
        pcMat <- apply(pcMat, 1:2, 
            convertParentChildMatNames,
            parData = parData
            )
        # 03-24-19 
        # why are some taxa gettign mysterious ".1" additions?
        #print(sort(unique(pcMat[,1])))
        # made it so that original names get annotated, not newly added parents
        # this issue is now avoided
        #
        # remove "NODE" root-stem 
        pcMat <- pcMat[pcMat[,1] != "NODE", ]
        ############################
        # Calculate the taxon tree!
        tree <- parentChild2taxonTree(
            parentChild = pcMat,
            tipSet = tipSet,
            cleanTree = cleanTree)
        #
        tree$parentChild <- pcMat
        }
    ##############
    #
    if(method == "Linnean"){
        tree <- getLinneanTaxonTreePBDB(
            dataTransform = dataTransform, 
            tipSet = tipSet,
            cleanTree = cleanTree,
            rankTaxon = rankTaxon
            )
        }
    #########
    #
    if(method == "parentChildOldMergeRoot" | method == "parentChildOldQueryPBDB"){
        tree <- parentChildPBDBOld(
            dataTransform = dataTransform, 
            tipSet = tipSet,
            cleanTree = cleanTree,
            method = method, 
            APIversion = APIversion,
            failIfNoInternet = failIfNoInternet
            )
        if(is.null(tree)){return(NULL)}
        }
    ####################
    tree$taxaDataPBDB <- taxaDataPBDB
    return(tree)
    }        
    
    
#' @rdname makePBDBtaxonTree
#' @export    
plotTaxaTreePBDB<-function(taxaTree, edgeLength = 1){
    taxaTree$edge.length <- rep(edgeLength ,Nedge(taxaTree))
    taxaTree$edge.length [taxaTree$edge[,2] <= Ntip(taxaTree)] <- edgeLength/5
    taxaTree$root.edge <- edgeLength*1.5
    plot(taxaTree,
        show.tip.label = FALSE,
            no.margin = TRUE,
        root.edge=TRUE,
        edge.width = 0.2)
    nodelabels(taxaTree$node.label,
        cex=0.5, 
        adj = c(1.1,0.5))
    }


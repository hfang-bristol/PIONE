#' Function to do prioritisation through random walk techniques
#'
#' \code{oPier} is supposed to prioritise nodes given an input graph and a list of seed nodes. It implements Random Walk with Restart (RWR) and calculates the affinity score of all nodes in the graph to the seeds. The priority score is the affinity score. Parallel computing is also supported for Linux-like or Windows operating systems. It returns an object of class "pNode". 
#'
#' @param seeds a named input vector containing a list of seed nodes. For this named vector, the element names are seed/node names (e.g. gene symbols), the element (non-zero) values used to weight the relative importance of seeds. Alternatively, it can be a matrix or data frame with two columns: 1st column for seed/node names, 2nd column for the weight values
#' @param g an object of class "igraph" to represent network. It can be a weighted graph with the node attribute 'weight'
#' @param seeds.inclusive logical to indicate whether non-network seed genes are included for prioritisation. If TRUE (by default), these genes will be added to the netowrk 
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' an object of class "pNode", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 5 containing node priority information, where nNode is the number of nodes in the input graph, and the 5 columns are "name" (node names), "node" (1 for network genes, 0 for non-network seed genes), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores)}
#'  \item{\code{g}: an input "igraph" object}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{oPier}}
#' @include oPier.r
#' @examples
#' # a) provide the input nodes/genes with the significance info
#' sig <- rbeta(500, shape1=0.5, shape2=1)
#' \dontrun{
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' data <- data.frame(symbols=org.Hs.eg$gene_info$Symbol[1:500], sig)
#' 
#' # b) provide the network
#' g <- xRDataLoader(RData.customised='org.Hs.PCommons_UN', RData.location=RData.location)
#'
#' # c) perform priority analysis
#' pNode <- oPier(seeds=data, g=g, restart=0.75)
#' }

oPier <- function(seeds, g, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    data <- seeds
    if(is.null(data)){
        setSeeds <- data
        stop("The input seeds must be not NULL.\n")
    }else{
		if (is.vector(data)){
			if(length(data)>1){
				# assume a vector
				if(is.null(names(data))){
					warning("The input seeds do not have node names (assuming equal weights).\n")
					tmp <- rep(1/length(data), length(data))
					names(tmp) <- data
					data <- tmp
				}
			}else{
				# assume a file
				data <- utils::read.delim(file=data, header=FALSE, row.names=NULL, stringsAsFactors=FALSE)
			}
		}
		if (is.vector(data)){
			scores <- data
		}else if(is.matrix(data) | is.data.frame(data)){
			data <- as.matrix(data)
			data_list <- split(x=data[,2], f=as.character(data[,1]))
			res_list <- lapply(data_list, function(x){
				x <- as.numeric(x)
				x <- x[!is.na(x)]
				## maximum weight if multiple per gene
				if(length(x)>0){
					max(x)
				}else{
					NULL
				}
			})
			scores <- unlist(res_list)
		}
		
		scores[scores<0] <- 0
		setSeeds <- data.frame(scores)
	}
    
    if(is(g,"igraph")){
    	ig <- g
    
		if(verbose){
			message(sprintf("The input graph has %d nodes and %d edges (%s) ...", vcount(ig),ecount(ig),as.character(Sys.time())), appendLF=TRUE)
		}
		
		####################
		if(seeds.inclusive){
			ind <- match(rownames(setSeeds), V(ig)$name)
			if(sum(is.na(ind))>0){
				newnodes <- rownames(setSeeds)[is.na(ind)]
				ig <- igraph::add_vertices(ig, nv=length(newnodes), name=newnodes)
				
				if(verbose){
					message(sprintf("In the resulting graph, %d non-network seed genes are added as isolated nodes (%s) ...", length(newnodes),as.character(Sys.time())), appendLF=TRUE)
				}
			}
		}
		####################
		
	}else{
		stop("The function must apply to the 'igraph' object.\n")
	}

    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################"))
        message(sprintf("'xRWR' is being called (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
    PTmatrix <- suppressWarnings(oRWR(g=ig, normalise=normalise, setSeeds=setSeeds, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose))
	rownames(PTmatrix) <- V(ig)$name
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################"))
        message(sprintf("'xRWR' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    if(1){
    	seeds <- rep(0, nrow(PTmatrix))
    	flag <- match(rownames(setSeeds), rownames(PTmatrix))
    	seeds[flag[!is.na(flag)]] <- 1
    	weights <- rep(0, nrow(PTmatrix))
    	weights[flag[!is.na(flag)]] <- setSeeds[!is.na(flag),1]
    	
    	########
    	nodes <- rep(0, nrow(PTmatrix))
    	flag <- match(V(g)$name, rownames(PTmatrix))
    	nodes[flag[!is.na(flag)]] <- 1
    	########
    	
    	df <- data.frame(name=rownames(PTmatrix), node=nodes, seed=seeds, weight=weights, priority=as.matrix(PTmatrix), stringsAsFactors=FALSE)
    	df <- df[order(-df$priority,-df$seed,-df$node), ]
    	df <- cbind(df, rank=rank(-df$priority,ties.method='min'))
    }
    
    ####################################################################################

    pNode <- list(g=ig,
                  priority = df
                 )
    class(pNode) <- "pNode"   
    
    invisible(pNode)
}

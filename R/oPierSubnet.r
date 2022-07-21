#' Function to identify a gene network from top prioritised genes
#'
#' \code{oPierSubnet} is supposed to identify maximum-scoring gene subnetwork from a graph with the node information on priority scores, both are part of an object of class "pNode". It returns an object of class "igraph". 
#'
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget")
#' @param priority.quantile the quantile of the top priority genes. By default, 10% of top prioritised genes will be used for network analysis. If NULL or NA, all prioritised genes will be used
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param test.permutation logical to indicate whether the permutation test is perform to estimate the significance of identified network with the same number of nodes. By default, it sets to false
#' @param num.permutation the number of permutations generating the null distribution of the identified network
#' @param respect how to respect nodes to be sampled. It can be one of 'none' (randomly sampling) and 'degree' (degree-preserving sampling)
#' @param aggregateBy the aggregate method used to aggregate edge confidence p-values. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RDS files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph". It has ndoe attributes: signficance, score, type, priority (part of the "pNode" object). If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence)
#' @note The priority score will be first scaled to the range x=[0 100] and then is converted to pvalue-like significant level: 10^(-x). Next, \code{\link{oSubneterGenes}} is used to identify a maximum-scoring gene subnetwork that contains as many highly prioritised genes as possible but a few lowly prioritised genes as linkers. An iterative procedure of scanning different priority thresholds is also used to identify the network with a desired number of nodes/genes. Notably, the preferential use of the same network as used in gene-level prioritisation is due to the fact that gene-level affinity/priority scores are smoothly distributed over the network after being walked. In other words, the chance of identifying such a gene network enriched with top prioritised genes is much higher.
#' @export
#' @seealso \code{\link{oPierSubnet}}
#' @include oPierSubnet.r
#' @examples
#' \dontrun{
#' # a) perform network analysis
#' # find maximum-scoring subnet with the desired node number=50
#' subnet <- oPierSubnet(pNode, priority.quantile=0.1, subnet.size=50, RData.location=RData.location)
#'
#' # b) save subnet results to the files called 'subnet_edges.txt' and 'subnet_nodes.txt'
#' output <- igraph::get.data.frame(subnet, what="edges")
#' utils::write.table(output, file="subnet_edges.txt", sep="\t", row.names=FALSE)
#' output <- igraph::get.data.frame(subnet, what="vertices")
#' utils::write.table(output, file="subnet_nodes.txt", sep="\t", row.names=FALSE)
#'
#' # c) visualise the identified subnet
#' ## do visualisation with nodes colored according to the priority
#' xVisNet(g=subnet, pattern=V(subnet)$priority, vertex.shape="sphere")
#' ## do visualisation with nodes colored according to pvalue-like signficance
#' xVisNet(g=subnet, pattern=-log10(as.numeric(V(subnet)$significance)), vertex.shape="sphere", colormap="wyr")
#' 
#' # d) visualise the identified subnet as a circos plot
#' library(RCircos)
#' xCircos(g=subnet, entity="Gene", RData.location=RData.location)
#' }

oPierSubnet <- function(pNode, priority.quantile=0.1, network=NA, STRING.only=NA, network.customised=NULL, subnet.significance=0.01, subnet.size=NULL, test.permutation=FALSE, num.permutation=100, respect=c("none","degree"), aggregateBy=c("Ztransform","fishers","logistic","orderStatistic"), verbose=TRUE, placeholder=NULL, guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    respect <- match.arg(respect)
    aggregateBy <- match.arg(aggregateBy)
    
    if(is(pNode,"pNode")){
        df_priority <- pNode$priority[, c("seed","weight","priority")]
        
		network <- network[1]
		if(!is.na(network)){
			network <- match.arg(network)
		}else{
			if(is.null(network.customised)){
				network.customised <- pNode$g
			}
		}
		
		priority <- df_priority$priority
		names(priority) <- rownames(df_priority)
		
		## scale to the range [0 100] and then convert to pvalue-like signficant level
		x <- priority
		y <- (x - min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))
		pval <- 10^(-100*y)
		
	}else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	df_priority <- pNode$priority[, c("name","rank","rating")]
    	df_priority$priority <- df_priority$rating
    	
    	network <- network[1]
		if(!is.na(network)){
			network <- match.arg(network)
		}else{
			if(is.null(network.customised)){
				network.customised <- pNode$metag
			}
		}
		
		priority <- df_priority$priority
		names(priority) <- rownames(df_priority)
		
		##############
		# convert into pvalue by 10^(-x*2)
		x <- priority
		pval <- 10^(-x*2)
		##############
				
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }

    
	if(verbose){
		now <- Sys.time()
		message(sprintf("The '%s' object contains %d prioritised genes", class(pNode), length(priority)), appendLF=TRUE)
	}
	
	## only keep the top priority (quantile)
	## priority quantile
	priority.quantile <- as.numeric(priority.quantile)
	if(length(priority.quantile>0 & priority.quantile<1) & !is.na(priority.quantile)){
		cf <- stats::quantile(pval, priority.quantile, na.rm=TRUE)
		ind <- which(pval<cf)
		pval <- pval[ind]
		
		priority <- priority[ind]
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Among prioritised genes, %d genes are used for network analysis", length(pval)), appendLF=TRUE)
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t maximum priority: %1.2e; minimum priority: %1.2e", max(priority), min(priority)), appendLF=TRUE)
		message(sprintf("\t minimum p-value: %1.2e; maximum p-value: %1.2e", min(pval), max(pval)), appendLF=TRUE)
	}
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################"))
        message(sprintf("oSubneterGenes is being called (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
    if(is.na(network)){
    	subg <- oSubneterGenes(data=pval, network.customised=network.customised, seed.genes=TRUE, subnet.significance=subnet.significance, subnet.size=subnet.size, test.permutation=test.permutation, num.permutation=num.permutation, respect=respect, aggregateBy=aggregateBy, verbose=verbose, placeholder=placeholder, guid=guid)
    }else{
    	subg <- oSubneterGenes(data=pval, network=network, STRING.only=STRING.only, network.customised=network.customised, seed.genes=TRUE, subnet.significance=subnet.significance, subnet.size=subnet.size, test.permutation=test.permutation, num.permutation=num.permutation, respect=respect, aggregateBy=aggregateBy, verbose=verbose, placeholder=placeholder, guid=guid)
	}
	
	# extract relevant info
	if(ecount(subg)>0 && is(subg,"igraph")){
		relations <- igraph::get.data.frame(subg, what="edges")[,c(1,2)]
		if(!is.null(subg$combinedP)){
			relations$edgeConfidence <- igraph::get.data.frame(subg, what="edges")[,"edgeConfidence"]
		}
		nodes <- igraph::get.data.frame(subg, what="vertices")
		nodes <- cbind(name=nodes$name, description=nodes$description, significance=nodes$significance, score=nodes$score, type=nodes$type, priority=priority[rownames(nodes)])
		if(igraph::is.directed(subg)){
			subnet <- igraph::graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
		}else{
			subnet <- igraph::graph.data.frame(d=relations, directed=FALSE, vertices=nodes)
		}
		if(!is.null(subg$combinedP)){
			subnet$combinedP <- subg$combinedP
		}
	}else{
		subnet <- NULL
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################"))
        message(sprintf("oSubneterGenes has finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    return(subnet)
}

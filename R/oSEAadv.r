#' Function to conduct enrichment analysis given a list of sets and a list of ontologies
#'
#' \code{oSEAadv} is supposed to conduct enrichment analysis given a list of sets and a list of ontologies. It is an advanced version of \code{xSEA}, returning an object of the class 'eSAD'.
#'
#' @param list_vec an input vector containing gene symbols. Alternatively it can be a list of vectors, representing multiple groups of genes
#' @param sets a tibble with two columns 'onto' and 'set' (each is an SET object)
#' @param background a background vector containing gene symbols as the test background. If NULL, by default all annotatable are used as background
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param which.distance which terms with the distance away from the ontology root (if any) is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the test statistic used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > fisher's exact test > binomial test. In other words, in terms of the calculated p-value, hypergeometric test < fisher's exact test < binomial test
#' @param background.annotatable.only logical to indicate whether the background is further restricted to the annotatable. By default, it is NULL: if ontology.algorithm is not 'none', it is always TRUE; otherwise, it depends on the background (if not provided, it will be TRUE; otherwise FALSE). Surely, it can be explicitly stated
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to account for the hierarchy of the ontology. It can be one of "none", "pc", "elim" and "lea". For details, please see 'Note' below
#' @param elim.pvalue the parameter only used when "ontology.algorithm" is "elim". It is used to control how to declare a signficantly enriched term (and subsequently all genes in this term are eliminated from all its ancestors)
#' @param lea.depth the parameter only used when "ontology.algorithm" is "lea". It is used to control how many maximum depth is used to consider the children of a term (and subsequently all genes in these children term are eliminated from the use for the recalculation of the signifance at this term)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @return 
#' an object of class "eSAD", a tibble with 5 columns; they are "group" (the input group names), "onto" (input ontologies), "info" [a list-column each containing 14 columns including including "id" (set ID), "name" (set name), "nA" (number of annotations), "nO" (number of overlaps), "fc" (fold change), "zscore" (z-score), "pvalue" (p-value), "adjp" (adjusted p-value), "or" (odds ratio), "CIl" (lower bound confidence interval for the odds ratio), "CIu" (upper bound confidence interval for the odds ratio), "distance" (distance to the root), "namespace", "overlap" (the detailed members for overlaps)], "data" (a list-column each containing 1 column), and "background" (a list-column each containing 1 column).
#' @note none
#' @export
#' @seealso \code{\link{oSEA}}
#' @include oSEAadv.r
#' @examples
#' \dontrun{
#' sets <- tibble(onto=c('GOMF','GOBP','KEGG','Bucket','PSG')[c(1)]) %>% mutate(set=map(onto,~oRDS(str_c("org.Hs.eg",.x),placeholder=placeholder)))
#' list_vec <- BioGRID_HCoV %>% nest(data=-from_tax) %>% mutate(gene=map(data,~pull(.x,to))) %>% select(-data) %>% deframe()
#' 
#' # basic analysis
#' esad <- oSEAadv(list_vec, sets)
#' esad %>% oSEAextract() %>% write_delim('results_esad.txt','\t')
#' gp <- oSEAballoon(esad, top=10, adjp.cutoff=0.05, zlim=c(0,10), slim=NULL)
#' gp + theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=6),strip.text.y=element_text(size=6,angle=0))
#'
#' # advanced analysis
#' ig <- oRDS("ig.GOMF", placeholder=placeholder)
#' V(ig)$name <- V(ig)$term_name
#' # advanced visual
#' gp <- oSEAggraph(esad, ig, fixed=F, leave=F, layout='dendrogram', node.label.size=0)
#' gp$gp_template + ggraph::geom_node_text(aes(filter=leaf,label=name,angle=node_angle(x,y)), repel=F,hjust=0,size=1.5) + ggraph::geom_node_text(aes(filter=!leaf,label=name), repel=T,size=1.5,color='red',alpha=0.5,check_overlap=T)+ expand_limits(x=c(-2,2), y=c(-2,2)) -> gp_template
#' gp + geom_edge_diagonal2(aes(color=node.term_namespace, alpha=stat(index)),width=0.2) + scale_edge_alpha('direction', guide='edge_direction') + theme(legend.box='vertical') -> gp1
#' gp1 + guides(edge_alpha=F, edge_colour=guide_legend('namespace',direction="vertical"), size=guide_legend('-log10(FDR)','top',direction="horizontal",ncol=3), color=guide_colorbar('Z-score','top',direction="horizontal",barheight=0.5)) + theme(legend.position='right')
#' }

oSEAadv <- function(list_vec, sets, background=NULL, size.range=c(10,2000), min.overlap=5, which.distance=NULL, test=c("fisher","hypergeo","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=FALSE, verbose=TRUE, silent=FALSE)
{
    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    test <- match.arg(test)
    p.tail <- match.arg(p.tail)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    p.tail <- match.arg(p.tail)
    
    ################################################
    ## support enrichment analysis for modular genes from the cModule object
    if(is(list_vec,'cModule')){
    	list_vec <- split(x=list_vec$mem$nodes, f=list_vec$mem$modules)
    }
    ################################################
    
    ############
    if(length(list_vec)==0){
    	return(NULL)
    }
    ############
    if(is.vector(list_vec) & !is(list_vec,"list")){
    	list_vec <- list(list_vec)
	}else if(is(list_vec,"list")){
		## Remove null elements in a list
		list_vec <- base::Filter(base::Negate(is.null), list_vec)
		if(length(list_vec)==0){
			return(NULL)
		}
    }else{
        stop("The input data must be a vector or a list of vectors.\n")
    }
    
	list_names <- names(list_vec)
	if(is.null(list_names)){
		list_names <- paste0('G', 1:length(list_vec))
		names(list_vec) <- list_names
	}
    
    onto <- info <- NULL
    
    pb <- dplyr::progress_estimated(length(list_vec))
    ls_df <- lapply(seq_len(length(list_vec)), function(i){
		if(verbose){
			message(sprintf("\nAnalysing group %d ('%s') (%s) ...", i, names(list_vec)[i], as.character(Sys.time())), appendLF=T)
		}
		
		eset <- set <- data <- NULL

    	sets %>% dplyr::mutate(eset=purrr::map2(onto, set, function(x,y){
			if(verbose){
				message(sprintf("\tontology '%s' (%s) ...", x, as.character(Sys.time())), appendLF=T)
			}
    		eSet <- oSEA(data=list_vec[[i]], set=y, ig=NULL, background=background, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=FALSE)
    	})) %>% dplyr::filter(purrr::map_lgl(eset,~!is.null(.x))) %>% dplyr::mutate(info=purrr::map(eset,~.x$info)) %>% dplyr::mutate(data=purrr::map(eset,~.x$data)) %>% dplyr::mutate(background=purrr::map(eset,~.x$background)) %>% dplyr::select(onto, info, data, background) -> df_res
    	
    	pb$tick()$print()
    	
    	tibble::tibble(group=names(list_vec)[i], df_res)
	})
    df_res_all <- do.call(rbind, ls_df)
    
    eSAD <- df_res_all
    class(eSAD) <- c("eSAD",class(eSAD))
     ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (oSEAadv): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    return(eSAD)
}

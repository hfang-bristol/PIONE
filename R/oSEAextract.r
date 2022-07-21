#' Function to extract enrichment results
#'
#' \code{oSEAextract} is supposed to extract results of enrichment analysis. 
#'
#' @param obj an object of class "eSET" or "eSAD"
#' @param sortBy which statistics will be used for sorting and viewing gene sets (terms). It can be "adjp" for adjusted p value (FDR), "pvalue" for p value, "zscore" for enrichment z-score, "fc" for enrichment fold change, "nA" for the number of sets (terms), "nO" for the number in overlaps, "or" for the odds ratio, "name" for set name, "id" for set ID and "distance" for set distance
#' @return
#' a tibble with 14 columns:
#' \itemize{
#'  \item{\code{name}: set name}
#'  \item{\code{id}: set ID}
#'  \item{\code{nA}: number of annotations}
#'  \item{\code{nO}: number of overlaps}
#'  \item{\code{fc}: enrichment fold change}
#'  \item{\code{zscore}: enrichment z-score}
#'  \item{\code{pvalue}: nominal p-value}
#'  \item{\code{adjp}: adjusted p-value (FDR)}
#'  \item{\code{or}: odds ratio}
#'  \item{\code{CIl}: lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: upper bound confidence interval for the odds ratio}
#'  \item{\code{distance}: set distance to the root}
#'  \item{\code{namespace}: set namespace}
#'  \item{\code{overlap}: semi-colon separated members for overlaps}
#' }
#' @note none
#' @export
#' @seealso \code{\link{oSEAextract}}
#' @include oSEAextract.r
#' @examples
#' \dontrun{
#' oSEAextract(obj)
#' }

oSEAextract <- function(obj, sortBy=c("adjp","pvalue","zscore","fc","nA","nO","or","id","name","distance")) 
{
    
    sortBy <- match.arg(sortBy)
    
    if(is(obj,"eSET") | is(obj,"eSAD")){
		
		overlap <- NULL
		group <- onto <- info <- NULL
		pvalue <- adjp <- zscore <- fc <- or <- nA <- nO <- id <- name <- distance <- NULL
		
		if(is(obj,"eSET")){
			obj$info %>% mutate(overlap=purrr::map_chr(overlap,~str_c(.x,collapse=", "))) -> tab
			
    	}else if(is(obj,"eSAD")){
    		obj %>% dplyr::select(group,onto,info) %>% tidyr::unnest(info) %>% mutate(overlap=purrr::map_chr(overlap,~str_c(.x,collapse=", "))) -> tab
    		
    	}
    
		switch(sortBy, 
			adjp={res <- tab %>% dplyr::arrange(pvalue,-zscore)},
			pvalue={res <- tab %>% dplyr::arrange(pvalue,-zscore)},
			zscore={res <- tab %>% dplyr::arrange(-zscore,pvalue)},
			fc={res <- tab %>% dplyr::arrange(-fc,pvalue)},
			or={res <- tab %>% dplyr::arrange(-or,pvalue)},
			nA={res <- tab %>% dplyr::arrange(-nA,pvalue,-zscore)},
			nO={res <- tab %>% dplyr::arrange(-nO,pvalue,-zscore)},
			name={res <- tab %>% dplyr::arrange(name)},
			id={res <- tab %>% dplyr::arrange(id)},
			distance={res <- tab %>% dplyr::arrange(distance,pvalue,-zscore)}
		)
		
		if(is(obj,"eSAD")){
			res <- res %>% dplyr::arrange(group,onto)
		}
    	
    }else{
        warnings("There is no enrichment in the input object.\n")
        res <- NULL
    }
    
    res
}

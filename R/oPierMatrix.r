#' Function to extract priority or evidence matrix from a list of pNode objects
#'
#' \code{oPierMatrix} is supposed to extract priority or evidence matrix from a list of pNode objects. Also supported is the aggregation of priority matrix (similar to the meta-analysis) generating the priority results; we view this functionality as the discovery mode of the prioritisation.
#'
#' @param list_pNode a list of "pNode" objects or a "pNode" object
#' @param displayBy which priority will be extracted. It can be "score" for priority score/rating (by default), "rank" for priority rank, "weight" for seed weight, "pvalue" for priority p-value, "evidence" for the evidence (seed info)
#' @param combineBy how to resolve nodes/targets from a list of "pNode" objects. It can be "intersect" for intersecting nodes (by default), "union" for unionising nodes
#' @param aggregateBy the aggregate method used. It can be either "none" for no aggregation, or "orderStatistic" for the method based on the order statistics of p-values, "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'fishers' and 'logistic' are preferred here. Also supported are methods summing up evidence 'sum', taking the maximum of evidence ('max') or sequentially weighting evidence 'harmonic'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' If displayBy is 'evidence', an object of the class "eTarget", a list with following components:
#' \itemize{
#'  \item{\code{evidence}: a data frame of nGene X 6 containing gene evidence information, where nGene is the number of genes, and the 7 columns are seed info including "Overall" for the total number of different types of seeds, followed by details on individual type of seeds (that is, "dGene", "pGene", "fGene", "nGene", "eGene", "cGene")}
#'  \item{\code{metag}: an "igraph" object}
#' }
#' Otherwise (if displayBy is not 'evidence'), if aggregateBy is 'none' (by default), a data frame containing priority matrix, with each column/predictor for either priority score, or priorty rank or priority p-value. If aggregateBy is not 'none', an object of the class "dTarget", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a data frame of n X 4+7 containing gene priority (aggregated) information, where n is the number of genes, and the 4 columns are "name" (gene names), "rank" (ranks of the priority scores), "rating" (the 5-star score/rating), "description" (gene description), and 7 seed info columns including "seed" (whether or not seed genes), "nGene" (nearby genes), "cGene" (conformation genes), "eGene" (eQTL gens), "dGene" (disease genes), "pGene" (phenotype genes), and "fGene" (function genes)}
#'  \item{\code{predictor}: a data frame containing predictor matrix, with each column/predictor for either priority score/rating, or priorty rank or priority p-value}
#'  \item{\code{metag}: an "igraph" object}
#'  \item{\code{list_pNode}: a list of "pNode" objects}
#' }
#' @note none
#' @export
#' @seealso \code{\link{oPierMatrix}}
#' @include oPierMatrix.r
#' @examples
#' \dontrun{
#' # get predictor matrix for targets
#' df_score <- oPierMatrix(ls_pNode)
#' # get evidence for targets
#' eTarget <- oPierMatrix(ls_pNode, displayBy="evidence")
#' # get target priority in a discovery mode
#' dTarget <- oPierMatrix(ls_pNode, displayBy="pvalue", aggregateBy="fishers")
#' }

oPierMatrix <- function(list_pNode, displayBy=c("score","rank","weight","pvalue","evidence"), combineBy=c('union','intersect'), aggregateBy=c("none","fishers","logistic","Ztransform","orderStatistic","harmonic","max","sum"), verbose=TRUE, placeholder=NULL, guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################

    displayBy <- match.arg(displayBy)
    combineBy <- match.arg(combineBy)
    aggregateBy <- match.arg(aggregateBy) 
    
   	if(is(list_pNode,"pNode")){
		list_pNode <- list(list_pNode)
	}else if(is(list_pNode,"list")){
		## Remove null elements in a list
		list_pNode <- base::Filter(base::Negate(is.null), list_pNode)
		if(length(list_pNode)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'pNode' objects or a 'pNode' object.\n")
	}
	
	## check list_names
	list_names <- names(list_pNode)
	if(is.null(list_names)){
		list_names <- paste('Predictor', 1:length(list_pNode), sep=' ')
		names(list_pNode) <- list_names
	}
	
	## get nodes involved
	ls_nodes <- lapply(list_pNode, function(x){
		V(x$g)$name
	})
	if(combineBy=='intersect'){
		nodes <- base::Reduce(intersect, ls_nodes)
	}else if(combineBy=='union'){
		nodes <- base::Reduce(union, ls_nodes)
	}
	nodes <- sort(nodes)
	
	if(displayBy=='evidence' | (displayBy=='pvalue' & aggregateBy!="none") | (displayBy=='score' & (aggregateBy %in% c("harmonic","max","sum")))){
		############
		## seed info
		predictor_names <- names(list_pNode)
		predictor_names <- gsub('_.*', '', predictor_names)
		ls_df <- lapply(1:length(list_pNode), function(i){
			pNode <- list_pNode[[i]]
			genes <- rownames(pNode$priority)[pNode$priority$seed==1]
			df <- data.frame(Gene=genes, Predictor=rep(predictor_names[i], length(genes)), stringsAsFactors=FALSE)
		})
		df <- do.call(rbind, ls_df)
		mat_evidence <- as.matrix(oSparseMatrix(df, rows=nodes, columns=NULL, verbose=FALSE))
		###
		# sorted by the number of seed gene types, followed by the total number of seed genes
		if(ncol(mat_evidence)>1){
        	mat_evidence <- mat_evidence[order(mat_evidence[,1],apply(mat_evidence,1,sum),decreasing=TRUE),]
        }else{
			####
			# deal with only 1 predictor
			####
			tmp <- mat_evidence[order(mat_evidence[,1],apply(mat_evidence,1,sum),decreasing=TRUE),]
			tmp_matrix <- matrix(tmp, ncol=ncol(mat_evidence))
			colnames(tmp_matrix) <- colnames(mat_evidence)
			rownames(tmp_matrix) <- names(tmp)
			mat_evidence <- tmp_matrix
			###
		}
		
		############
		if(0){
			## get edges involved
			ls_edges <- pbapply::pblapply(list_pNode, function(x){
				relations <- igraph::get.data.frame(x$g, what="edges")
			})
			edges <- do.call(rbind, ls_edges) %>% unique()
			## get metag
			metag <- igraph::graph.data.frame(d=edges, directed=FALSE, vertices=nodes)
		}else{
			#############
			### very fast
			#############
			edges <- igraph::get.data.frame(list_pNode[[1]]$g, what="edges")
			## get metag
			metag <- igraph::graph.data.frame(d=edges, directed=FALSE, vertices=nodes)
		}
		##############
		
	}
	
	if(displayBy!='evidence'){
		## Combine into a data frame called 'df_predictor'
		ls_priority <- pbapply::pblapply(list_pNode, function(pNode){
			p <- pNode$priority
			ind <- match(nodes, rownames(p))
			#ind <- ind[!is.na(ind)]
			if(displayBy=='score' | displayBy=='pvalue'){
				res <- p[ind, c("priority")]
			}else if(displayBy=='rank'){
				res <- p[ind, c("rank")]
			}else if(displayBy=='weight'){
				res <- p[ind, c("weight")]
			}
		})
		df_predictor <- do.call(cbind, ls_priority)
		rownames(df_predictor) <- nodes
	
		## replace NA with worst value
		if(displayBy=='score' | displayBy=='weight' | displayBy=='pvalue'){
			df_predictor[is.na(df_predictor)] <- 0
		}else if(displayBy=='rank'){
			df_predictor[is.na(df_predictor)] <- length(nodes)
		}

	}
	
	## only when displayBy=='pvalue'
	## Convert into p-values by computing an empirical cumulative distribution function
	if(displayBy=='pvalue'){
		ls_pval <- lapply(1:ncol(df_predictor), function(j){
			x <- df_predictor[,j]
			my.CDF <- stats::ecdf(x)
			pval <- 1 - my.CDF(x)
		})
		df_pval <- do.call(cbind, ls_pval)
		rownames(df_pval) <- rownames(df_predictor)
		colnames(df_pval) <- colnames(df_predictor)
		df_predictor <- df_pval
		
		## aggregate p values
		if(aggregateBy != "none"){
			df_ap <- dnet::dPvalAggregate(pmatrix=df_predictor, method=aggregateBy)
			df_ap <- sort(df_ap, decreasing=FALSE)
			
			## get rank
			df_rank <- rank(df_ap, ties.method="min")
			######
			df_ap[df_ap==0] <- min(df_ap[df_ap!=0])
			######
			## adjp
			df_adjp <- stats::p.adjust(df_ap, method="BH")
			######
			## rating: first log10-transformed ap and then being rescaled into the [0,5] range
			rating <- -log10(df_ap)
			####
			rating <- sqrt(rating)
			####
			rating <- 5 * (rating - min(rating))/(max(rating) - min(rating))
			
			## df_priority
			df_priority <- data.frame(name=names(df_ap), rank=df_rank, pvalue=df_ap, fdr=df_adjp, rating=rating, stringsAsFactors=FALSE)
			### add description
			df_priority$description <- oSymbol2GeneID(df_priority$name, details=TRUE, verbose=verbose, placeholder=placeholder, guid=guid)$description
			###
			
			## df_predictor
			ind <- match(names(df_ap), rownames(df_predictor))
			if(ncol(df_predictor)>1){
				df_predictor <- df_predictor[ind,]
			}else{
				####
				# deal with only 1 predictor
				####
				tmp <- df_predictor[ind,]
				tmp_matrix <- matrix(tmp, ncol=ncol(df_predictor))
				colnames(tmp_matrix) <- colnames(df_predictor)
				rownames(tmp_matrix) <- names(tmp)
				df_predictor <- tmp_matrix
				####				
			}
			
			## reorder mat_evidence
			ind_row <- match(df_priority$name, rownames(mat_evidence))
			ind_col <- match(unique(predictor_names), colnames(mat_evidence))
			if(ncol(mat_evidence)>1){
				mat_evidence <- mat_evidence[ind_row,ind_col]
			}else{
				####
				# deal with only 1 predictor
				####
				tmp <- mat_evidence[ind_row,ind_col]
				tmp_matrix <- matrix(tmp, ncol=ncol(mat_evidence))
				colnames(tmp_matrix) <- colnames(mat_evidence)
				rownames(tmp_matrix) <- names(tmp)
				mat_evidence <- tmp_matrix
				###
			}
			overall <- apply(mat_evidence!=0, 1, sum)
			
			## return dTarget
			#priority <- cbind(df_priority,Overall=overall, mat_evidence)
			#priority <- data.frame(df_priority[,c("name","rank","rating","description")], seed=ifelse(overall!=0,'Y','N'), mat_evidence[,c("nGene","cGene","eGene","dGene","pGene","fGene")], stringsAsFactors=FALSE)
			
			if(0){
				## previous version
				ind <- match(c("nGene","cGene","eGene","dGene","pGene","fGene"), colnames(mat_evidence))
			}else{
				#################
				## 20211115
				## now generalised version
				#################
				ind <- match(unique(predictor_names), colnames(mat_evidence))
			}
			priority <- data.frame(df_priority[,c("name","rank","rating","description")], seed=ifelse(overall!=0,'Y','N'), mat_evidence[,ind[!is.na(ind)]], stringsAsFactors=FALSE)
			dTarget <- list(priority  = priority,
							predictor = df_predictor,
							metag	  = metag,
							list_pNode  = list_pNode
						 )
			class(dTarget) <- "dTarget"
			
			df_predictor <- dTarget
		}
		
	}else if(displayBy=='evidence'){
		## reorder mat_evidence
		ind_col <- match(unique(predictor_names), colnames(mat_evidence))
		mat_evidence <- mat_evidence[, ind_col]
		overall <- apply(mat_evidence!=0, 1, sum)

		## return eTarget
		eTarget <- list(evidence  = data.frame(Overall=overall, mat_evidence),
						metag	  = metag
						)
		class(eTarget) <- "eTarget"
		df_predictor <- eTarget
		
	}else if(displayBy=='score' & (aggregateBy %in% c("harmonic","max","sum"))){
		if(aggregateBy=="max"){
			summaryFun <- max
		}else if(aggregateBy=="sum"){
			summaryFun <- sum
		}else if(aggregateBy=="harmonic"){
			summaryFun <- function(x){
				base::sum(x / base::rank(-x,ties.method="min")^2)
			}
		}
		
		df_harmonic <- apply(df_predictor, 1, summaryFun)
		
		if(1){
			df_harmonic <- sort(df_harmonic, decreasing=T)
			
			## get rank
			df_rank <- rank(-df_harmonic, ties.method="min")
			## rating: first sqrt and then being rescaled into the [0,5] range
			####
			rating <- sqrt(df_harmonic)
			####
			rating <- 5 * (rating - min(rating))/(max(rating) - min(rating))
			
			## df_priority
			df_priority <- data.frame(name=names(df_harmonic), rank=df_rank, harmonic=df_harmonic, rating=rating, stringsAsFactors=FALSE)
			### add description
			df_priority$description <- oSymbol2GeneID(df_priority$name, details=TRUE, verbose=verbose, placeholder=placeholder, guid=guid)$description
			###
			
			## df_predictor
			ind <- match(names(df_harmonic), rownames(df_predictor))
			if(ncol(df_predictor)>1){
				df_predictor <- df_predictor[ind,]
			}else{
				####
				# deal with only 1 predictor
				####
				tmp <- df_predictor[ind,]
				tmp_matrix <- matrix(tmp, ncol=ncol(df_predictor))
				colnames(tmp_matrix) <- colnames(df_predictor)
				rownames(tmp_matrix) <- names(tmp)
				df_predictor <- tmp_matrix
				####				
			}
			
			## reorder mat_evidence
			ind_row <- match(df_priority$name, rownames(mat_evidence))
			ind_col <- match(unique(predictor_names), colnames(mat_evidence))
			if(ncol(mat_evidence)>1){
				mat_evidence <- mat_evidence[ind_row,ind_col]
			}else{
				####
				# deal with only 1 predictor
				####
				tmp <- mat_evidence[ind_row,ind_col]
				tmp_matrix <- matrix(tmp, ncol=ncol(mat_evidence))
				colnames(tmp_matrix) <- colnames(mat_evidence)
				rownames(tmp_matrix) <- names(tmp)
				mat_evidence <- tmp_matrix
				###
			}
			overall <- apply(mat_evidence!=0, 1, sum)
			
			## return dTarget
			#priority <- cbind(df_priority,Overall=overall, mat_evidence)
			#priority <- data.frame(df_priority[,c("name","rank","rating","description")], seed=ifelse(overall!=0,'Y','N'), mat_evidence[,c("nGene","cGene","eGene","dGene","pGene","fGene")], stringsAsFactors=FALSE)
			if(0){
				## previous version
				ind <- match(c("nGene","cGene","eGene","dGene","pGene","fGene"), colnames(mat_evidence))
			}else{
				#################
				## 20211115
				## now generalised version
				#################
				ind <- match(unique(predictor_names), colnames(mat_evidence))
			}
			priority <- data.frame(df_priority[,c("name","rank","rating","description")], seed=ifelse(overall!=0,'Y','N'), mat_evidence[,ind[!is.na(ind)]], stringsAsFactors=FALSE)
			dTarget <- list(priority  = priority,
							predictor = df_predictor,
							metag	  = metag,
							list_pNode  = list_pNode
						 )
			class(dTarget) <- "dTarget"
			
			df_predictor <- dTarget
		}
	}
	
	if(verbose){
		
		if(displayBy=="evidence"){
			message(sprintf("A matrix of %d genes x %d evidence are generated", nrow(mat_evidence), ncol(mat_evidence)), appendLF=TRUE)
			
		}else{
			if(displayBy=="pvalue" & aggregateBy!="none"){
				message(sprintf("A total of %d genes are prioritised, combined by '%s' and aggregated by '%s' from %d predictors", nrow(df_predictor$priority), combineBy, aggregateBy, length(list_pNode)), appendLF=TRUE)
			}else{
				message(sprintf("A matrix of %d genes x %d predictors are generated, displayed by '%s' and combined by '%s'", nrow(df_predictor), ncol(df_predictor), displayBy, combineBy), appendLF=TRUE)
			}

		}
	}
	
####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
	
    invisible(df_predictor)
}

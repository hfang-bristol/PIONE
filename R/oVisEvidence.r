#' Function to visualise evidence for prioritised genes in a gene network
#'
#' \code{oVisEvidence} is supposed to visualise evidence for prioritised genes in a gene network. It returns an object of class "igraph". 
#'
#' @param xTarget an object of class "dTarget", "sTarget" or "eTarget"
#' @param g an object of class "igraph". If NA, the 'metag' will be used, which is part of the input object "xTarget". If provided, it must have a node attribute called 'priority'
#' @param nodes which node genes are in query. If NULL, the top gene will be queried
#' @param node.info tells the additional information used to label nodes. It can be one of "none" (only gene labeling), "smart" for (by default) using three pieces of information (if any): genes, 5-star ratings, and associated ranks (marked by an @ icon)
#' @param neighbor.order an integer giving the order of the neighborhood. By default, it is 1-order neighborhood
#' @param neighbor.seed logical to indicate whether neighbors are seeds only. By default, it sets to true
#' @param neighbor.top the top number of the neighbors with the highest priority. By default, it sets to NULL to disable this parameter
#' @param largest.comp logical to indicate whether the largest component is only retained. By default, it sets to true for the largest component being left
#' @param show logical to indicate whether to show the graph
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param legend.position the legend position. If NA, the legend is not shown
#' @param legend.horiz logical specifying the legend horizon. If TRUE, set the legend horizontally rather than vertically
#' @param mtext.side the side of marginal text. If NA, it is not shown
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param edge.width the width of the edge. If NULL, the width edge is proportional to the 'weight' edge attribute (if existed)
#' @param vertex.size the size of each vertex. If null, each vertex has the size proportional to the degree of nodes
#' @param vertex.size.nonseed the size of each nonseed vertex. If null, each vertex has the size proportional to the degree of nodes
#' @param vertex.label.color the color of vertex labels
#' @param vertex.label.color.nonseed the color of nonseed vertex labels
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return
#' a subgraph, an object of class "igraph".
#' @export
#' @seealso \code{\link{oColormap}}
#' @include oVisEvidence.r
#' @examples
#' \dontrun{
#' ## TNFRSF1A
#' oVisEvidence(xTarget, nodes="TNFRSF1A", neighbor.order=1, neighbor.seed=TRUE, neighbor.top=NULL, vertex.label.color="black", vertex.label.cex=0.7, vertex.label.dist=0.6, vertex.label.font=1, vertex.label.family="Arial", legend.position="bottomleft", legend.horiz=TRUE, newpage=FALSE)
#' ## UBA52
#' oVisEvidence(xTarget, nodes="UBA52", neighbor.order=1, neighbor.seed=TRUE, neighbor.top=20, vertex.label.color="black", vertex.label.cex=0.7, vertex.label.dist=0.6, vertex.label.font=1, legend.position="bottomleft", legend.horiz=TRUE, newpage=FALSE)
#' }

oVisEvidence <- function(xTarget, g=NA, nodes=NULL, node.info=c("smart","simple","none"), neighbor.order=1, neighbor.seed=TRUE, neighbor.top=NULL, largest.comp=TRUE, show=TRUE, colormap="ggplot2", legend.position="topleft", legend.horiz=FALSE, mtext.side=3, verbose=TRUE, edge.width=NULL, vertex.size=NULL, vertex.size.nonseed=NULL, vertex.label.color="steelblue", vertex.label.color.nonseed=NULL, ...)
{

    node.info <- match.arg(node.info)

    if(is(xTarget,"dTarget")){
    	if(is.null(xTarget$pPerf)){
    		df_evidence <- xTarget$priority[, 5:ncol(xTarget$priority)]
    	}else{
    		df_evidence <- xTarget$priority[, 6:ncol(xTarget$priority)]
    	}
        df_priority <- xTarget$priority[, c("rank","rating")]
        df_priority$priority <- df_priority$rating
		
    }else if(is(xTarget,"sTarget")){
        df_evidence <- as.data.frame(xTarget$evidence$evidence)
        df_priority <- xTarget$priority[, c("rank","rating")]
        df_priority$priority <- df_priority$rating
		
    }else if(is(xTarget,"eTarget")){
        df_evidence <- as.data.frame(xTarget$evidence)
        # here, sorted by the number of seed gene types
        df_priority <- df_evidence[order(df_evidence[,1],decreasing=TRUE),]
		
		#neighbor.top <- NULL
		
    }else{
    	stop("The function must apply to a 'dTarget' or 'sTarget' or 'eTarget' object.\n")
    }
	
	if(is(g,"igraph")){
		g <- g
	}else{
		if(is(xTarget,"dTarget") | is(xTarget,"eTarget")){
			g <- xTarget$metag
		}else if(is(xTarget,"sTarget")){
			g <- xTarget$evidence$metag
		}
		
		V(g)$priority <- df_priority[V(g)$name, 'priority']
	}
	if(!is(g,'igraph')){
		stop("The input 'g' must be provided!\n")
	}
	
	if(is.null(nodes)){
		nodes <- rownames(df_evidence)[1]
	}else{
		ind <- match(nodes, rownames(df_evidence))
		ind <- which(!is.na(ind))
		if(length(ind)>=1){
			nodes <- nodes[ind]
		}else{
			#nodes <- rownames(df_evidence)[1]
			warning(sprintf("\tNo found for queried %s", nodes), appendLF=TRUE)
			return(NULL)
		}
	}
    
    neighs.out <- igraph::neighborhood(g, order=neighbor.order, nodes=nodes, mode="all")
	neighbors <- names(unlist(neighs.out))
	if(neighbor.seed){
		# restrict to seeds
		ind <- neighbors %in% rownames(df_evidence)[df_evidence$seed=='Y']
		neighbors <- neighbors[ind]
	}
	if(!is.null(neighbor.top)){
		neighbor.top <- as.integer(neighbor.top)
		if(neighbor.top > length(neighbors)){
			neighbor.top <- length(neighbors)
		}
		
		ind <- match(rownames(df_priority), neighbors)
		df_neighbors <- df_priority[!is.na(ind),]
		neighbors <- rownames(df_neighbors)[1:neighbor.top]
	}
	vids <- union(neighbors, nodes)
    subg <- dnet::dNetInduce(g, nodes_query=vids, knn=0, remove.loops=TRUE, largest.comp=largest.comp)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("The %d-order graph induced by %d input gene(s): %d nodes and %d edges", neighbor.order, length(nodes), vcount(subg), ecount(subg)), appendLF=TRUE)
	}
	
	ind <- match(V(subg)$name, rownames(df_evidence))
	df_val <- df_evidence[ind,]
	ls_val <- lapply(1:nrow(df_val), function(i){
		if(sum(df_val[i,-1])==0){
			NULL
		}else{
			as.numeric(df_val[i,-1]!=0)
		}
	})
	names(ls_val) <- rownames(df_val)
	
	## vertex.label
	vertex.label <- V(subg)$name
	if(!is(xTarget,"eTarget")){
		ind <- match(vertex.label, rownames(df_priority))
		df_nodes <- df_priority[ind, ]
		if(node.info=="smart"){
			vertex.label <- paste0(vertex.label, "\n[", signif(df_nodes$priority,digits=3), "@", df_nodes$rank, "]")
		}else if(node.info=="simple"){
			vertex.label <- paste0(vertex.label, "\n(", df_nodes$rank, ")")
		}
	}
	
	## nodes NULL are drawn as circles
	vertex.shape <- rep("pie", length(ls_val))
	vertex.shape[sapply(ls_val, is.null)] <- "circle"
	## pie color
	pie.color <- oColormap(colormap)(ncol(df_val)-1)
	## legend text
	legend.text <- colnames(df_val)[-1]
	if(0){
		legend.text[grep('OMIM|disease',legend.text,ignore.case=TRUE)] <- "dGene"
		legend.text[grep('Phenotype',legend.text,ignore.case=TRUE)] <- "pGene"
		legend.text[grep('Function',legend.text,ignore.case=TRUE)] <- "fGene"
		legend.text[grep('nearbyGenes',legend.text,ignore.case=TRUE)] <- "nGene"
		legend.text[grep('eQTL',legend.text,ignore.case=TRUE)] <- "eGene"
		legend.text[grep('HiC|Hi-C',legend.text,ignore.case=TRUE)] <- "cGene"
	}
	## vertex size
	if(is.null(vertex.size)){
		vertex.size <- igraph::degree(subg)
		if(min(vertex.size) == max(vertex.size)){
			vertex.size <- 8
		}else{
			vertex.size <- 5 * (vertex.size - min(vertex.size))/(max(vertex.size) - min(vertex.size)) + 8
		}
		
		if(!is.null(vertex.size.nonseed)){
			vertex.size[sapply(ls_val, is.null)] <- vertex.size.nonseed
		}
		
	}else{
		if(!is.null(vertex.size.nonseed)){
			vertex.size <- rep(vertex.size, length(ls_val))
			vertex.size[sapply(ls_val, is.null)] <- vertex.size.nonseed
		}
	}
	
	## vertex.label.color
	if(!is.null(vertex.label.color)){
		if(!is.null(vertex.label.color.nonseed)){
			vertex.label.color <- rep(vertex.label.color, length(ls_val))
			vertex.label.color[sapply(ls_val, is.null)] <- vertex.label.color.nonseed
		}
	}
	
	## edge.width
	if(is.null(edge.width) & !is.null(E(subg)$weight)){
		## extract edge weight
		x <- as.numeric(E(subg)$weight)
		if(length(x)>0){
			if(max(x)-min(x)>0){
				## rescale into an interval [1,4] as edge width
				edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
			}
		}
	}
	
	#############################################################
	if(show){
		## draw graph
		oVisNet(subg, vertex.shape=vertex.shape, vertex.pie=ls_val, vertex.pie.color=list(pie.color), vertex.pie.border="grey", vertex.label=vertex.label, vertex.color="grey", vertex.size=vertex.size, signature=FALSE, edge.width=edge.width, vertex.label.color=vertex.label.color, ...)
		if(!is.na(legend.position)){
			graphics::legend(legend.position, legend=legend.text, col=pie.color, pch=13, bty="n", pt.cex=1.2, cex=1, text.col="darkgrey", text.font=4, horiz=legend.horiz)
		}
		if(!is.na(mtext.side)){
			graphics::mtext(paste0("Interacting partners for ", paste0(nodes,collapse=',')), side=mtext.side, adj=0, cex=0.8, font=4, family="sans")
		}
	}
	#############################################################
	
	## append evidence node attributes
	df_tmp <- df_val[,-1]
	colnames(df_tmp) <- legend.text
	for(i in 1:ncol(df_tmp)){
		igraph::vertex_attr(subg, colnames(df_tmp)[i]) <- df_tmp[,i]
	}
	V(subg)$vertex.label <- vertex.label
	
    return(subg)
}

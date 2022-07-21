#' Function to visualise evidence and priority scores for prioritised genes in a gene network
#'
#' \code{oVisEvidenceAdv} is supposed to visualise evidence and priority scores for prioritised genes in a gene network. It returns an object of class "ggplot". 
#'
#' @param xTarget an object of class "dTarget", "sTarget" or "eTarget"
#' @param g an object of class "igraph". If NA, the 'metag' will be used, which is part of the input object "xTarget". If provided, it must have node attributes ('priority','xcoord','ycoord')
#' @param nodes which node genes are in query. If NULL, the top gene will be queried
#' @param node.info tells the additional information used to label nodes. It can be one of "none" (only gene labeling), "smart" for (by default) using three pieces of information (if any): genes, 5-star ratings, and associated ranks (marked by an @ icon)
#' @param neighbor.order an integer giving the order of the neighborhood. By default, it is 1-order neighborhood
#' @param neighbor.seed logical to indicate whether neighbors are seeds only. By default, it sets to true
#' @param neighbor.top the top number of the neighbors with the highest priority. By default, it sets to NULL to disable this parameter
#' @param largest.comp logical to indicate whether the largest component is only retained. By default, it sets to true for the largest component being left
#' @param node.label.size a vector specifying node size or a character specifying which node attribute used for node label size
#' @param node.label.color the node label color
#' @param node.label.alpha the 0-1 value specifying transparency of node labelling
#' @param node.label.padding the padding around the labeled node
#' @param node.label.arrow the arrow pointing to the labeled node
#' @param node.label.force the repelling force between overlapping labels
#' @param node.shape an integer specifying node shape
#' @param node.color.title a character specifying the title for node coloring
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param node.size.range the range of actual node size
#' @param title a character specifying the title for the plot
#' @param edge.color a character specifying the edge colors
#' @param edge.color.alpha the 0-1 value specifying transparency of edge colors
#' @param edge.curve a numeric value specifying the edge curve. 0 for the straight line
#' @param edge.arrow.gap a gap between the arrow and the node
#' @param pie.radius the radius of a pie. If NULL, it equals roughly 1/75
#' @param pie.color the border color of a pie
#' @param pie.color.alpha the 0-1 value specifying transparency of pie border colors
#' @param pie.thick the pie border thickness
#' @param ... additional graphic parameters for oGGnetwork
#' @return
#' a ggplot object.
#' @export
#' @seealso \code{\link{oVisEvidence}}
#' @include oVisEvidenceAdv.r
#' @examples
#' \dontrun{
#' ## TNFRSF1A
#' oVisEvidenceAdv(xTarget, nodes="TNFRSF1A", neighbor.order=1, neighbor.seed=TRUE, neighbor.top=NULL)
#' }

oVisEvidenceAdv <- function(xTarget, g=NA, nodes=NULL, node.info=c("smart","simple","none"), neighbor.order=1, neighbor.seed=TRUE, neighbor.top=NULL, largest.comp=TRUE, node.label.size=2, node.label.color='steelblue4', node.label.alpha=0.9, node.label.padding=0.5, node.label.arrow=0, node.label.force=0.1, node.shape=19, node.color.title='Pi\nrating', colormap='white-yellow-red', ncolors=64, zlim=c(0,5), node.size.range=5, title='', edge.color="steelblue", edge.color.alpha=0.2, edge.curve=0, edge.arrow.gap=0.025, pie.radius=NULL, pie.color='black', pie.color.alpha=1, pie.thick=0.1,...)
{

    node.info <- match.arg(node.info)
	
	if(is(g,"igraph")){
		# If provided, it must have node attributes ('priority','xcoord','ycoord')
		if(any(is.null(V(g)$xcoord), is.null(V(g)$ycoord), is.null(V(g)$priority))){
			return(NULL)
		}
	}
	
	subg <- oVisEvidence(xTarget=xTarget, g=g, nodes=nodes, node.info=node.info, neighbor.order=neighbor.order, neighbor.seed=neighbor.seed, neighbor.top=neighbor.top, largest.comp=largest.comp, show=FALSE)
	
	if(!is(subg,'igraph')){
		return(NULL)
	}
	if(ecount(subg)<=1){
		return(NULL)
	}
	
	## layout
	if(any(is.null(V(subg)$xcoord), is.null(V(subg)$ycoord))){
		#glayout <- igraph::layout_as_tree(subg,root=dnet::dDAGroot(subg),circular=TRUE,flip.y=TRUE)
		if(0){
			glayout <- igraph::layout_with_kk(subg)
			V(subg)$xcoord <- glayout[,1]
			V(subg)$ycoord <- glayout[,2]
		}else{
			subg <- subg %>% oLayout("graphlayouts.layout_with_stress")
		}
	}
	
	#################
	gp <- oGGnetwork(g=subg, node.label='vertex.label', node.label.size=node.label.size, node.label.color=node.label.color, node.label.alpha=node.label.alpha, node.label.padding=node.label.padding, node.label.arrow=node.label.arrow, node.label.force=node.label.force, node.shape=node.shape, node.xcoord='xcoord', node.ycoord='ycoord', node.color='priority', node.color.title=node.color.title, colormap=colormap, ncolors=ncolors, zlim=zlim, node.size.range=node.size.range, title=title, edge.color=edge.color, edge.color.alpha=edge.color.alpha, edge.curve=edge.curve, edge.arrow.gap=edge.arrow.gap,...)
    #################

	#df <- gp$data_nodes
	df <- gp$data
	
	if(0){
		## previously
		columns <- c('dGene','pGene','fGene','nGene','eGene','cGene')
		df_sub <- as.data.frame(matrix(0, nrow=nrow(df), ncol=length(columns)+2))
		colnames(df_sub) <- c('x','y',columns)
		ind <- match(colnames(df_sub), colnames(df))
		df_sub[,!is.na(ind)] <- df[,ind[!is.na(ind)]]
		ind <- which(apply(df_sub[,c(-1,-2)],1,sum)!=0)
		if(length(ind)>0){
			df_sub <- df_sub[ind, ]
			df_sub_1 <- df_sub[,c(1,2)]
			df_sub_2 <- df_sub[,c(-1,-2)]
			df_sub_2[df_sub_2>=1] <- 1
			df_sub <- cbind(df_sub_1, df_sub_2)

			gp <- oPieplot(df_sub, columns, colormap='ggplot2', pie.radius=pie.radius, pie.color=pie.color, pie.color.alpha=pie.color.alpha, pie.thick=pie.thick, legend.title='Seed gene', gp=gp)
		
		}
		
	}else{
		#################
		## now (20211115)
		#################
		x <- y <- ycoord <- vertex.label <- . <- NULL
		
		df_sub_1 <- df %>% dplyr::select(x,y)
		df_sub_2 <- df %>% dplyr::select(ycoord:vertex.label) %>% dplyr::select(-1,-length(.))
		df_sub_2[df_sub_2>=1] <- 1
		ind <- which(apply(df_sub_2,1,sum)!=0)
		df_sub <- cbind(df_sub_1, df_sub_2) %>% dplyr::slice(ind)
		
		columns <- colnames(df_sub_2)
		gp <- oPieplot(df_sub, columns, colormap='ggplot2', pie.radius=pie.radius, pie.color=pie.color, pie.color.alpha=pie.color.alpha, pie.thick=pie.thick, legend.title='Seed gene', gp=gp)
	}
	
    invisible(gp)
}

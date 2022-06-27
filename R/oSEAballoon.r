#' Function to visualise enrichment results using a balloon plot
#'
#' \code{oSEAballoon} is supposed to visualise enrichment results using a balloon plot. It returns an object of class "ggplot".
#'
#' @param obj an object of class "eSET" or "eSAD". Alterntively, it can be a tibble having all these columns (named as 'name','adjp'); optionally ('group', 'namespace', 'onto', 'nO', 'distance')
#' @param top the number of the top terms (sorted according to adjp). For the "eSET" object, if it is 'auto' (for eSET), only the significant terms (see below adjp.cutoff) will be displayed
#' @param adjp.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05. Only works when top is 'auto' above
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the -log10(adjp)
#' @param color.title the legend title for colorbar
#' @param shape the point shape
#' @param size.range the point size range
#' @param slim the minimum and maximum point size values, defaulting to the range of the nO
#' @param size.title the size legend title
#' @param wrap.width a positive integer specifying wrap width of name
#' @param legend.direction the legend guide direction. It can be "horizontal", "vertical" and "auto" ("vertical" when having >15 rows; otherwise "horizontal")
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oSEAballoon}}
#' @include oSEAballoon.r
#' @examples
#' \dontrun{
#' obj %>% oSEAballoon()
#' #obj %>% oSEAextract() %>% filter(onto=='KEGG',group==694009) %>% ggdotchart('name',y='zscore',color='nO', facet='group',add='segments',rotate=T,label='or',font.label=list(color='white',vjust=0.5,size=9),dot.size=6) + theme_pubr(base_size=6) + scale_y_continuous(position="top")
#' }

oSEAballoon <- function(obj, top=10, adjp.cutoff=0.05, colormap="brewer.Reds", zlim=NULL, color.title="-log10(adjP)", shape='Q', size.range=c(0.5,4), slim=NULL, size.title="Overlap", wrap.width=NULL, legend.direction=c("auto","horizontal","vertical"))
{
    
	## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    legend.direction <- match.arg(legend.direction)
    
    if(is.null(obj)){
        warnings("There is no enrichment in the input object.\n")
        return(NULL)
    }
    
    group <- onto <- name <- nO <- adjp <- distance <- namespace <- NULL
    
    df <- NULL
    if(is(obj,'eSAD')){
    	obj %>% oSEAextract() -> df
	}else if(is(obj,'eSET')){
    	df <- obj$info
	}else{
		df <- obj
	}
	
	if(is(df,'data.frame')){
		
		if(all(c('name','adjp') %in% colnames(df))){
		
			if(!('nO' %in% colnames(df))){
				df %>% dplyr::mutate(nO=1) -> df
			}
			if(!('group' %in% colnames(df))){
				df %>% dplyr::mutate(group='group') -> df
			}
			if(!('onto' %in% colnames(df))){
				df %>% dplyr::mutate(onto='onto') -> df
			}
			if(!('distance' %in% colnames(df))){
				df %>% dplyr::mutate(distance=0) -> df
			}
			if(!('namespace' %in% colnames(df))){
				df %>% dplyr::mutate(namespace='namespace') -> df
			}
		}
	}
	
	if(is(top,'numeric')){
		top <- as.integer(top)
		df %>% dplyr::group_by(group,namespace) %>% dplyr::top_n(top,-adjp) %>% dplyr::ungroup() %>% dplyr::select(name) -> df_name
		df %>% dplyr::semi_join(df_name,by="name") -> df
	}
	
	df %>% dplyr::filter(adjp<adjp.cutoff) %>% dplyr::arrange(onto,namespace,dplyr::desc(distance),name) %>% dplyr::mutate(name=forcats::fct_inorder(name)) -> df
	
	n <- NULL
	if(df %>% dplyr::count(namespace) %>% dplyr::pull(n) %>% min() >1){
		df %>% dplyr::mutate(namespace=stringr::str_replace(namespace,' ','\n')) -> df
	}
	
	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################

	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(df$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste0(y[1], '...')
			}else{
				y
			}
		})
		df$name <- unlist(res_list)
	}
	
	name <- color <- or <- CIl <- CIu <- size <- NULL
	group <- ontology <- NULL
	
	df$color <- -log10(df$adjp)
	if(is.null(zlim)){
		tmp <- df$color
		zlim <- c(0, ceiling(max(tmp[!is.infinite(tmp)])))
	}
	df$color[df$color<=zlim[1]] <- zlim[1]
	df$color[df$color>=zlim[2]] <- zlim[2]
	
	df$size <- df$nO
	if(is.null(slim)){
		tmp <- df$size
		slim <- c(0, ceiling(max(tmp[!is.infinite(tmp)])))
	}
	df$size[df$size<=slim[1]] <- slim[1]
	df$size[df$size>=slim[2]] <- slim[2]
	
	###########################################
	df %>% ggplot(aes(x=group, y=name)) + geom_point(aes(color=color, size=size), shape=shape) -> gp
	
	gp <- gp + theme_classic() + theme(legend.position="right", legend.title=element_text(size=7), legend.text=element_text(size=6), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=6,angle=0), axis.text.x=element_text(size=7,angle=45,hjust=0,vjust=1))
	gp <- gp + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
	
	if(legend.direction=="auto"){
		if(df %>% dplyr::count(name) %>% nrow() <= 15){
			legend.direction <- "horizontal"
		}else{
			legend.direction <- "vertical"
		}
	}
	if(legend.direction=="vertical"){
		## color
		gp <- gp + scale_colour_gradientn(colors=oColormap(colormap)(64), limits=zlim, guide=guide_colorbar(title=color.title, title.position="top", barwidth=0.4,order=2))
		## size
		gp <- gp + scale_size_continuous(limits=slim, range=size.range, guide=guide_legend(size.title,title.position="top",keywidth=0.4,ncol=1,order=1))
		if(length(size.range)==1 | df %>% dplyr::count(nO) %>% nrow() == 1){
			gp <- gp + guides(size=FALSE)
		}
	}else if(legend.direction=="horizontal"){
		## color
		gp <- gp + scale_colour_gradientn(colors=oColormap(colormap)(64), limits=zlim, guide=guide_colorbar(title=color.title, title.position="top", barheight=0.4,direction="horizontal",order=2))
		## size
		gp <- gp + scale_size_continuous(limits=slim, range=size.range, guide=guide_legend(size.title,title.position="top",keyheight=0.4,ncol=3,byrow=TRUE,order=1))
		if(length(size.range)==1 | df %>% dplyr::count(nO) %>% nrow() == 1){
			gp <- gp + guides(size=FALSE)
		}
	}
	
	## strip
	if(df %>% dplyr::count(namespace) %>% nrow() > 1){
		gp <- gp + facet_grid(namespace~.,scales="free_y",space="free_y")
		gp <- gp + theme(strip.placement='outside', strip.background=element_rect(fill="transparent",color="grey80"), strip.text.y=element_text(size=7,face="italic",angle=0,hjust=0))
	}
	
	## x-axis position
	gp <- gp + scale_x_discrete(position="top")
	
	gp
}


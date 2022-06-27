#' Function to visualise enrichment results using a forest plot
#'
#' \code{oSEAforest} is supposed to visualise enrichment results using a forest plot. A point is colored by the significance level, and a horizontal line for the 95% confidence interval (CI) of odds ratio (OR; the wider the CI, the less reliable). It returns an object of class "ggplot".
#'
#' @param obj an object of class "eSET" or "eSAD". Alterntively, it can be a tibble having all these columns (named as 'name','adjp','or','CIl','CIu'); optionally ('group', 'namespace', 'onto', 'nO')
#' @param top the number of the top terms (sorted according to OR). For the "eSET" object, if it is 'auto' (for eSET), only the significant terms (see below adjp.cutoff) will be displayed
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
#' @param sortBy It can be "or" or "none"
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oSEAforest}}
#' @include oSEAforest.r
#' @examples
#' \dontrun{
#' obj %>% oSEAforest()
#' obj %>% oSEAextract() %>% filter(onto=='PSG') %>% oSEAforest()
#' }

oSEAforest <- function(obj, top=10, adjp.cutoff=0.05,  colormap="brewer.Reds", zlim=NULL, color.title="-log10(adjP)", shape=18, size.range=c(0,1), slim=NULL, size.title="Overlap", wrap.width=NULL, legend.direction=c("auto","horizontal","vertical"), sortBy=c("or","none"))
{
    
	## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    legend.direction <- match.arg(legend.direction)
    sortBy <- match.arg(sortBy)
    
    if(is.null(obj)){
        warnings("There is no enrichment in the 'obj' object.\n")
        return(NULL)
    }
    
    group <- onto <- name <- nO <- adjp <- or <- CIl <- CIu <- namespace <- NULL
    
    df <- NULL
    if(is(obj,'eSAD')){
    	obj %>% oSEAextract() -> df
	}else if(is(obj,'eSET')){
    	df <- obj$info
	}else{
		df <- obj
	}
	
	if(is(df,'data.frame')){
		
		if(all(c('name','adjp','or','CIl','CIu') %in% colnames(df))){
		
			if(!('nO' %in% colnames(df))){
				df %>% dplyr::mutate(nO=1) -> df
			}
			if(!('group' %in% colnames(df))){
				df %>% dplyr::mutate(group='group') -> df
			}
			if(!('onto' %in% colnames(df))){
				df %>% dplyr::mutate(onto='onto') -> df
			}
			if(!('namespace' %in% colnames(df))){
				df %>% dplyr::mutate(namespace='namespace') -> df
			}
		}
	}
	
	if(is(top,'numeric')){
		top <- as.integer(top)
		df %>% dplyr::group_by(group,namespace) %>% dplyr::top_n(top,or) %>% dplyr::ungroup() %>% dplyr::select(name) -> df_name
		df %>% dplyr::semi_join(df_name,by="name") -> df
	}
	
	if(sortBy=='or'){
		df %>% dplyr::filter(adjp<adjp.cutoff) %>% dplyr::arrange(group,onto,namespace,or,name) %>% dplyr::mutate(name=forcats::fct_inorder(name)) -> df
		
	}else{
		df %>% dplyr::filter(adjp<adjp.cutoff) %>% dplyr::arrange(group,onto,namespace,rev(name)) %>% dplyr::mutate(name=forcats::fct_inorder(name)) -> df
		
	}
	
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
	
	color <- size <- NULL
	
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
	
	df %>% ggplot(aes(y=name, x=log2(or), xmin=log2(CIl), xmax=log2(CIu), color=color)) + geom_pointrange(aes(size=size),shape=shape) -> gp
	#gp + geom_vline(xintercept=0, color='grey80') -> gp
	
	gp <- gp + theme_classic() + theme(legend.position="right", legend.title=element_text(size=7), legend.text=element_text(size=6), axis.title.y=element_blank(), axis.title.x=element_text(size=7), axis.text.y=element_text(size=6,angle=0), axis.text.x=element_text(size=7))
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
	
	## facet_grid + strip
	n_group <- df %>% dplyr::count(group) %>% nrow()
	n_namespace <- df %>% dplyr::count(namespace) %>% nrow()
	if(n_group!=1 | n_namespace!=1){
		if(n_group==1){
			gp <- gp + facet_grid(namespace~.,scales="free_y",space="free_y")	
		}else if(n_namespace==1){
			gp <- gp + facet_grid(.~group,scales="free_y",space="free_y")	
		}else{
			gp <- gp + facet_grid(namespace~group,scales="free_y",space="free_y")	
		}
		gp <- gp + theme(strip.placement='outside', strip.background=element_rect(fill="transparent",color="grey80"), strip.text.y=element_text(size=7,face="italic",angle=0,hjust=0), strip.text.x=element_text(size=7,face="bold",angle=0,hjust=0.5))
	}
	
	## x-axis position
	#gp <- gp + scale_x_continuous(position="top")
	#gp <- gp + ylab("log2(odds ratio)")
	
	gp
}


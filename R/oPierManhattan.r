#' Function to visualise prioritised genes using manhattan plot
#'
#' \code{oPierManhattan} is supposed to visualise prioritised genes using manhattan plot. Genes with the top priority are highlighed. It returns an object of class "ggplot".
#'
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget")
#' @param color a character vector for colors to alternate chromosome colorings. If NULL, ggplot2 default colors will be used. If a single character is provided, it can be "jet" (jet colormap) or "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta)
#' @param point.size the point size
#' @param top the number of the top targets to be labelled/highlighted
#' @param top.label.type how to label the top targets. It can be "box" drawing a box around the labels , and "text" for the text only
#' @param top.label.size the highlight label size
#' @param top.label.col the highlight label color
#' @param top.label.query which top genes in query will be labelled. By default, it sets to NULL meaning all top genes will be displayed. If labels in query can not be found, then all will be displayed
#' @param label.query.only logical to indicate whether only genes in query will be displayed. By default, it sets to FALSE. It only works when labels in query are enabled/found
#' @param chromosome.only logical to indicate whether only genes from input data will be displayed. By default, it sets to TRUE
#' @param y.scale how to transform the y scale. It can be "normal" for no transformation, "sqrt" for square root transformation, and "log" for log-based transformation
#' @param y.lab the y labelling. If NULL (by default), it shows the column of input data
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "placeholder"
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @param ... additional paramters associated with ggrepel::geom_text_repel
#' @return an object of class "ggplot", appended by an GR object called 'gr'
#' @note none
#' @export
#' @seealso \code{\link{oRDS}}, \code{\link{oColormap}}
#' @include oPierManhattan.r
#' @examples
#' \dontrun{
#' mp <- oPierManhattan(pNode, placeholder=placeholder)
#' mp$gr
#' ## control visuals
#' mp <- oPierManhattan(pNode, color='ggplot2', top=50, top.label.col="black", y.scale="sqrt", placeholder=placeholder)
#' mp
#' ## control labels
#' # only IL genes will be labelled
#' ind <- grep('^IL', rownames(pNode$priority))
#' top.label.query <- rownames(pNode$priority)[ind]
#' mp <- oPierManhattan(pNode, top.label.query=top.label.query, placeholder=placeholder)
#' mp
#' # only IL genes will be displayed
#' mp <- oPierManhattan(pNode, top.label.query=top.label.query, label.query.only=TRUE, placeholder=placeholder)
#' mp
#' }

oPierManhattan <- function(pNode, color=c("darkred","darkgreen"), point.size=0.5, top=50, top.label.type=c("box","text"), top.label.size=2, top.label.col="darkblue", top.label.query=NULL, label.query.only=FALSE, chromosome.only=TRUE, y.scale=c("normal","sqrt","log"), y.lab=NULL, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), font.family="sans", signature=TRUE, verbose=TRUE, placeholder=NULL, guid=NULL, ...)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    top.label.type <- match.arg(top.label.type)
    y.scale <- match.arg(y.scale)

    if(is(pNode,"pNode")){
        df_priority <- pNode$priority[, c("name","weight","priority")]
    }else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	df_priority <- pNode$priority[, c("name","rank","rating")]
    	df_priority$priority <- df_priority$rating
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=TRUE)
	}
    gr_Gene <- oRDS(GR.Gene[1], verbose=verbose, placeholder=placeholder, guid=guid)
    if(is.null(gr_Gene)){
    	GR.Gene <- "UCSC_knownGene"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.Gene), appendLF=TRUE)
		}
    	gr_Gene <- oRDS(GR.Gene, verbose=verbose, placeholder=placeholder, guid=guid)
    }
    
    ## ONLY restricted to genes with genomic locations
	#ind <- match(rownames(df_priority), names(gr_Gene))
	ind <- match(df_priority$name, names(gr_Gene))
	p_gr <- gr_Gene[ind[!is.na(ind)],]
	p_matrix <- df_priority[!is.na(ind),]
	
	## append genomic locations to GR object
	gr <- p_gr
	GenomicRanges::mcols(gr) <- cbind(GenomicRanges::mcols(gr), p_matrix)
	
	########################################
	if(label.query.only){
		if(!is.null(top.label.query)){
			top.label.query <- as.vector(t(top.label.query)) # just in case converting data.frame to vector
			ind <- match(names(gr), top.label.query)
			if(sum(!is.na(ind)) >= 1){
				gr <- gr[!is.na(ind)]
			}
		}
	}
	########################################
	
	
	## for sorting
	chrlabs <- paste('chr', as.character(c(1:22,'X','Y')), sep='')
	#######
	if(chromosome.only){
		ind <- chrlabs %in% unique(as.character(gr@seqnames@values))
		chrlabs <- chrlabs[ind]
	}
	#######	
	#eval(parse(text="seqlevels(gr) <- chrlabs"))
	GenomeInfoDb::seqlevels(gr) <- chrlabs
	
	## highlight points
	if(!is.null(top)){
		top <- as.integer(top)
		if(top > length(gr)){
			top <- length(gr)
		}
	}
	
	priority <- seqnames <- priority <- NULL
	###############################
	## calling ggbio::autoplot
	suppressMessages(ggp <- ggbio::autoplot(object=gr, aes(y=priority,color=seqnames,alpha=priority), coord="genome", geom='point', space.skip=0.01, size=point.size))
	
	## extract ggplot
	bp <- ggp@ggplot
	df <- bp$data
	
	## alternative colors
	if(!is.null(color)){
		if(length(color)>=2){
			alternative_colors <- color
			chrs <- levels(df[,1])
			N <- length(chrs)
			cols <- rep(alternative_colors, round(N/length(alternative_colors)) + 1)[1:N]
			names(cols) <- chrs
			bp <- bp + scale_color_manual(values=cols) + theme(legend.position="none")
		}else if(length(color)==1){
			chrs <- levels(df[,1])
			N <- length(chrs)
			cols <- oColormap(color)(N)
			names(cols) <- chrs
			bp <- bp + scale_color_manual(values=cols) + theme(legend.position="none")
		}
	}else{
		bp <- bp + theme(legend.position="none")
	}
	
	## vline
  	if(TRUE){
		vline.df <- df
		vline.df <- do.call(rbind, by(vline.df, vline.df$seqnames, function(dd){
			data.frame(start=min(dd$start), end=max(dd$end))
		}))
		## compute gap
		gap <- (vline.df$start[-1] + vline.df$end[-nrow(vline.df)])/2
		bp <- bp + geom_vline(xintercept=gap, alpha=0.5, color='gray70') + theme(panel.grid=element_blank())
  	}
	
	#bp <- bp + ggforce::facet_zoom(x=(seqnames=="chr2"))
	
	############
	## highlight top label
	############
	if(!is.null(top)){
		df_highlight <- bp$data[1:top,]
		
		#############################		
		## restrict to top in query for labels
		if(!is.null(top.label.query)){
			ind <- match(df_highlight$Symbol, top.label.query)
			if(sum(!is.na(ind)) >= 1){
				df_highlight <- df_highlight[!is.na(ind), ]
			}else{
				df_highlight <- NULL
			}
		}
		#############################		
		
		###########
		## potentially controlling only labels those in specific chromosome
		if(FALSE){
			ind <- match(df_highlight$seqnames, "chr1")
			if(sum(!is.na(ind)) >= 1){
				df_highlight <- df_highlight[!is.na(ind), ]
			}else{
				df_highlight <- NULL
			}
		}
		###########
		
		midpoint <- priority <- Symbol <- NULL
		if(!is.null(df_highlight)){
			if(top.label.type=="text"){
				bp <- bp + ggrepel::geom_text_repel(data=df_highlight, aes(x=midpoint,y=priority,label=Symbol), size=top.label.size, color=top.label.col, fontface='bold.italic', point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')), max.overlaps=Inf, ...)
			}else if(top.label.type=="box"){
				bp <- bp + ggrepel::geom_label_repel(data=df_highlight, aes(x=midpoint,y=priority,label=Symbol), size=top.label.size, color=top.label.col, fontface='bold.italic', box.padding=unit(0.2,"lines"), point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')), max.overlaps=Inf, ...)
			}
		}
	}
	
	## y scale
    if(y.scale=="sqrt"){
    	x <- NULL
    	bp <- bp + scale_y_continuous(trans=scales::sqrt_trans(), breaks=scales::trans_breaks("log10", function(x) 10^x, n=2))
    }else if(y.scale=="log"){
    	x <- NULL
    	bp <- bp + scale_y_continuous(trans=scales::log_trans(), breaks=scales::trans_breaks("log10", function(x) 10^x, n=2)) + annotation_logticks(sides='l')
    }
	
	if(!is.null(y.lab)){
		bp <- bp + ylab(y.lab)
	}
	
	bp <- bp + theme(axis.title.y=element_text(size=12), axis.text.y=element_text(color="black",size=8), axis.text.x=element_text(angle=45, hjust=1,color="black",size=10), panel.background=element_rect(fill=grDevices::rgb(0.98,0.98,0.98,1)))
	
	## caption
    if(signature){
    	caption <- paste("Created by oPierManhattan from Pi version", utils::packageVersion("Pi"))
    	bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	## put arrows on y-axis and x-axis
	bp <- bp + theme(axis.line.y=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")), axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))

    mp <- bp
    mp$gr <- gr
    
    invisible(mp)
}



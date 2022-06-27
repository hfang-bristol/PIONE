#' Function to visualise enrichment results using a upset plot
#'
#' \code{oUpsetAdv} is supposed to visualise enrichment results using a upset plot. The top is the kite plot, and visualised below is the combination matrix for overlapped genes. It returns an object of class "ggplot".
#'
#' @param data a data frame. It contains all these columns including  'value' and 'member' (separated by ',' or ', ')
#' @param sortBy how to sort bars/combinations. It can be 'degree' (number of non-zeros for a combination) or 'value' (value for a combination) or 'auto' (the same order as provided; by default)
#' @param color the point color
#' @param shape the point shape
#' @param size the point size
#' @param label.height.unit the unit specifying the per-row height of the combination matrix. By default, it is NULL. If specified (such as 8), it will be used to decide 'combmatrix.label.height' and 'combmatrix.label.text'
#' @param member.levels how to define the levels of members in the combination matri. It can be 'num' (the number of non-zeros for a member) or 'customised' (see 'levels.customised' below) or 'auto' (alphabeta order; by default)
#' @param levels.customised the customised levels of memebers in the combination matri. Only works when member.levels is 'customised' above
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oUpsetAdv}}
#' @include oUpsetAdv.r
#' @examples
#' \dontrun{
#' obj %>% transmute(value=nO, member=overlap) -> data
#' #data %>% arrange(-value) %>% mutate(member=fct_inorder(member)) -> data
#' data %>% oUpsetAdv() -> gp
#' gp + geom_point(aes(fill=as.factor(i)), shape=22, size=2, color="transparent") + geom_line(aes(group=i,color=as.factor(i)))
#' 
#' # customised levels for members
#' data %>% separate_rows(member,sep=',') %>% arrange(member) %>% pull(member) %>% unique() -> levels
#' data %>% oUpsetAdv(member.levels='customised',levels.customised=levels)
#' }

oUpsetAdv <- function(data, sortBy=c("auto","value","degree"), color="cyan4", shape=18, size=2, label.height.unit=NULL, member.levels=c("auto","num","customised"), levels.customised=NULL)
{
    
    sortBy <- match.arg(sortBy)
    member.levels <- match.arg(member.levels)
    
	member <- value <- i <- n <- NULL
    
    df <- NULL
	if(is(data,'data.frame')){
		if(all(c('value','member') %in% colnames(data))){
			df <- data %>% mutate(i=purrr::map_int(member,~stringr::str_split(.x, ',|, ', simplify=TRUE) %>% length()))
		}
	}

	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################
	
    if(sortBy=='degree'){
    	df %>% dplyr::arrange(i,value) -> df
    }else if(sortBy=='value'){
    	df %>% dplyr::arrange(value,i) -> df
    }
	
	df %>% mutate(member=fct_inorder(member)) %>% ggplot(aes(x=member,y=value)) + geom_col(fill="steelblue", color='transparent', width=0.1, alpha=0.1) + geom_point(shape=shape, color=color, size=size, fill='white', alpha=1) + theme_classic() + theme(axis.title.x=element_blank()) -> gp
	
	########################
	levels <- NULL
	if(member.levels=='num'){
		df %>% tidyr::separate_rows(member, sep=', |,') %>% dplyr::count(member) %>% dplyr::arrange(-n,member) %>% dplyr::pull(member) -> levels
	}else if(member.levels=='customised'){
		levels <- levels.customised
	}
	########################
	gp + ggupset::axis_combmatrix(sep=", |,", levels=levels) -> gp
	
	if(is.null(label.height.unit)){
		gp <- gp + ggupset::theme_combmatrix(combmatrix.panel.point.color.fill="steelblue", combmatrix.panel.point.color.empty="grey90", combmatrix.panel.point.size=1.5, combmatrix.panel.line.size=0.1)
	}else{
	
		if(member.levels=='auto'){
			df %>% tidyr::separate_rows(member, sep=', |,') %>% dplyr::count(member) %>% dplyr::arrange(member) %>% dplyr::pull(member) -> levels
		}
	
		#label.height.unit=8
		gp <- gp + ggupset::theme_combmatrix(combmatrix.panel.point.color.fill="steelblue", combmatrix.panel.point.color.empty="grey90", combmatrix.panel.point.size=1.5, combmatrix.panel.line.size=0.1, combmatrix.label.height=unit(label.height.unit*(length(levels)+1),"pt"), combmatrix.label.text=element_text(size=label.height.unit))
	}
	
	gp
}


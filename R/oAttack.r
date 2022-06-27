#' Function to perform targeted or combinatory attack for an igraph object
#'
#' \code{oAttack} is supposed to perform targeted or combinatory attack for an igraph object, where the effect of node removal is defined as the fraction of network nodes disconnected from the giant component (the largest connected component) remained. There are two types of attack analysis: removing a single node (accordingly we define the attackness per a node), removing nodes sequentially (targeted attack: Nodes ranked by for example betweenness centrality, and then the top nodes successively removed to calculate the effect of node removal), and removing nodes in combinations (combinatorial attack, often used based on the concept of combinatorial optimisation). Also attackness for a node is defined as the effect estimated when a single node removed.
#'
#' @param ig an object of class "igraph" with node attribute 'name'
#' @param measure how to perform attack analysis. This can be targeted attack (nodes ranked by 'degree', 'betweenness' or preranked ('prerank')) or combinatorial attack 'comb'
#' @param nodes.prerank a vector containing nodes preranked for targeted attack. Only works when measure is 'prerank'
#' @param nodes.combine a list containing nodes combined for combinatorial attack. Only works when measure is 'combine'
#' @return 
#' a tibble with 4 columns, including 'measure', 'frac.disconnected' (fraction of network nodes disconnected from the giant component), 'frac.removed' (fraction of network nodes removed), 'nodes.removed' (nodes removed, provided separated by ',').
#' @note none
#' @export
#' @seealso \code{\link{oAttack}}
#' @include oAttack.r
#' @examples
#' \dontrun{
#' set.seed(825)
#' ig <- sample_pa(20)
#' V(ig)$name <- paste0('n',seq(1,vcount(ig)))
#' 
#' ig %>% oAttack(measure='betweenness')
#' ig %>% oAttack(measure='degree')
#' 
#' # attackness for single nodes
#' nodes.combine <- utils::combn(V(ig)$name, 1,simplify=F) -> nodes.combine
#' ig %>% oAttack(measure='combine', nodes.combine=nodes.combine) -> df_attackness
#' # targeted attack (sequential removal of nodes preranked by attackness)
#' nodes.prerank <- df_attackness %>% arrange(-frac.disconnected) %>% pull(nodes.removed)
#' ig %>% oAttack(measure='prerank', nodes.prerank=nodes.prerank) -> df_res
#' 
#' # combinatorial attack
#' # 1) given combinations
#' nodes.combine <- list(c('n1'),c('n1','n3','n5'))
#' ig %>% oAttack(measure='combine', nodes.combine=nodes.combine)
#' # 2) combinations for any twos
#' V(ig)$name %>% combn(2, simplify=F) -> nodes.combine
#' oAttack(ig, measure="combine", nodes.combine=nodes.combine) %>% transmute(value=frac.disconnected, member=nodes.removed) -> data
#' data %>% arrange(-value) %>% top_n(5, value) %>% oUpsetAdv(member.levels="num")
#' # 3) a series of optimal combinations
#' levels <- V(ig)$name
#' tibble(i=seq(3)) %>% mutate(combine=map(i,~levels %>% combn(.x, simplify=F))) %>% mutate(res=map(combine, ~oAttack(ig, measure="combine", nodes.combine=.x) %>% select(frac.disconnected,nodes.removed) %>% top_n(1,frac.disconnected))) %>% select(i,res) %>% unnest(res) -> data
#' data %>% transmute(value=frac.disconnected, member=nodes.removed) %>% oUpsetAdv(member.levels='customised',levels.customised=levels, label.height.unit=6.5) + scale_y_continuous(limits=c(0,1)) + geom_line(aes(group=i,color=as.factor(i))) + geom_point(aes(fill=as.factor(i)), shape=22, size=2, color="transparent")
#' # 4) optimal combinations in the context of predefined nodes
#' nodes_predefined <- c('n2','n5')
#' nodes_rest <- setdiff(V(ig)$name, nodes_predefined)
#' levels.customised <- c(nodes_predefined, nodes_rest)
#' tibble(i=seq(3)) %>% mutate(combine=map(i,~nodes_rest %>% combn(.x,simplify=F) %>% lapply(function(y) union(y,nodes_predefined)))) %>% mutate(res=map(combine, ~oAttack(ig, measure="combine", nodes.combine=.x) %>% select(frac.disconnected,nodes.removed) %>% top_n(1,frac.disconnected))) %>% select(i,res) %>% unnest(res) -> data
#' data %>% transmute(value=frac.disconnected, member=nodes.removed) %>% oUpsetAdv(member.levels='customised',levels.customised=levels.customised, label.height.unit=6.5) + scale_y_continuous(limits=c(0,1)) + geom_line(aes(group=i,color=as.factor(i))) + geom_point(aes(fill=as.factor(i)), shape=22, size=2, color="transparent")
#' }

oAttack <- function(ig, measure=c('degree','betweenness','prerank','combine'), nodes.prerank=NULL, nodes.combine=NULL)
{

    measure <- match.arg(measure)
	
	max.comp.orig <- max(igraph::components(ig)$csize)
	n <- igraph::vcount(ig)
    
    if(measure=="combine"){
    	m <- length(nodes.combine)
    	max.comp.removed <- rep(max.comp.orig, m)
    	nodes.removed <- rep(max.comp.orig, m)
    	removed.pct <- rep(max.comp.orig, m)
    	pb <- dplyr::progress_estimated(m)
   		for(i in seq_len(m)){
   			pb$tick()$print()
   			ind <- match(V(ig)$name, nodes.combine[[i]])
   			v <- V(ig)$name[!is.na(ind)]
   			nodes.removed[i] <- paste(v,collapse=',')
			g.manual <- igraph::delete_vertices(ig, v)
            max.comp.removed[i] <- max(igraph::components(g.manual)$csize)
            removed.pct[i] <- 1 - igraph::vcount(g.manual) / n
   		}
		comp.pct <- max.comp.removed/max.comp.orig
   		
    }else{
    	removed.pct <- seq(0, 1, length=n + 1)
    	
    	if(measure=="prerank"){
    		ind <- match(nodes.prerank, V(ig)$name)
    		nodes_to_attack_inorder <- nodes.prerank[!is.na(ind)]
    		if(length(nodes_to_attack_inorder) < n){
    			return(NULL)
    		}
    	}else{
			if(measure=="betweenness"){
				#val <- igraph::centr_betw(ig)$res
				val <- igraph::betweenness(ig)
			}else if(measure=="degree"){
				val <- igraph::degree(ig)
			}
			nodes_to_attack_inorder <- sort(val, decreasing=TRUE) %>% names()
		}
		
        removed.pct <- seq(0, 1, length=n+1)
        max.comp.removed <- rep(max.comp.orig, n+1)
        nodes.removed <- rep(NA, n+1)
        g <- ig
        for(i in seq_len(n-1)){
        	g <- igraph::delete_vertices(g, nodes_to_attack_inorder[i])
            max.comp.removed[i+1] <- max(igraph::components(g)$csize)
            nodes.removed[i+1] <- paste(nodes_to_attack_inorder[1:i],collapse=',')
        }
		comp.pct <- max.comp.removed/max.comp.orig
		comp.pct[length(comp.pct)] <- 0
		nodes.removed[length(comp.pct)] <- paste(nodes_to_attack_inorder,collapse=',')
		
    }
	
	res <- tibble::tibble(measure=measure, frac.disconnected=1-comp.pct, frac.removed=removed.pct, nodes.removed=nodes.removed)

    return(res)
}
######################################################################
# DR
######################################################################
#' @title Definition for S3 class \code{DR}
#' @description \code{DR} has 3 components: df, index, gp.
#' @param df a data frame
#' @param index a data frame
#' @param gp a ggplot object
#' @return an object of S3 class \code{DR}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' DR(df, index, gp)
#' }
DR <- function(df, index, gp){
	## integrity checks
	if(!is(df,'data.frame') | !is(index,'data.frame') | !is(gp,'ggplot')){
		stop("The S3 class 'DR' object failed to pass integrity checks!\n")
	}
	value <- list(df=df, index=index, gp=gp)
	class(value) <- "DR"
	return(value)
}
#' @param x an object of class \code{DR}
#' @param ... other parameters
#' @rdname DR
#' @export
#' @method print DR
print.DR <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $df: a data frame of %d rows X %d columns", dim(x$df)[1],dim(x$df)[2]), "\n", sep="")
	cat(sprintf("  $index: a data frame of %d rows X %d columns", dim(x$index)[1],dim(x$index)[2]), "\n", sep="")
	cat(sprintf("  $gp: a ggplot object"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$df:\n")
	print(x$df[1:2,])
	cat("......\n")
	cat("$index:\n")
	print(x$index[1:2,])
	cat("......\n")
}


######################################################################
# dTarget
######################################################################
#' @title Definition for S3 class \code{dTarget}
#' @description \code{dTarget} has 3 components: priority, predictor and metag.
#' @param priority a data frame
#' @param predictor a data frame
#' @param metag an 'igraph' object
#' @return an object of S3 class \code{dTarget}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' dTarget(priority, predictor, metag)
#' }
dTarget <- function(priority, predictor, metag){
	## integrity checks
	if(!is(priority,'data.frame') | !is(predictor,'data.frame') | !is(metag,'igraph')){
		stop("The S3 class 'dTarget' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, predictor=predictor, metag=metag)
	class(value) <- "dTarget"
	return(value)
}
#' @param x an object of class \code{dTarget}
#' @param ... other parameters
#' @rdname dTarget
#' @export
print.dTarget <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	if(!is.null(x$metag)){
		cat(sprintf("  $metag: an igraph object with %d nodes and %d edges", vcount(x$metag), ecount(x$metag)), "\n", sep="")
	}
	cat(sprintf("  $predictor: a data frame of %d rows X %d columns", dim(x$predictor)[1],dim(x$predictor)[2]), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	cat(sprintf("  $list_pNode: a list of %d 'pNode' objects", length(x$list_pNode)), "\n", sep="")
	if(!is.null(x$pPerf)){
		cat(sprintf("  $pPerf: an object of the class 'pPerf'"), "\n", sep="")
	}
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print(x$priority[1:2,], row.names=FALSE)
	cat("......\n")
	if(!is.null(x$pPerf)){
		cat("$pPerf:\n")
		print(x$pPerf)
		cat("......\n")
	}
}


######################################################################
# pNode
######################################################################
#' @title Definition for S3 class \code{pNode}
#' @description \code{pNode} has 7 components: priority, g, SNP, Gene2SNP, nGenes, eGenes and cGenes.
#' @param priority a data frame
#' @param g an 'igraph' object
#' @return an object of S3 class \code{pNode}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' pNode(evidence, metag)
#' }
pNode <- function(priority, g){
	## integrity checks
	if(!is(priority,'data.frame') | !is(g,'igraph')){
		stop("The S3 class 'pNode' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, g=g)
	class(value) <- "pNode"
	return(value)
}
#' @param x an object of class \code{pNode}
#' @param ... other parameters
#' @rdname pNode
#' @export
print.pNode <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $g: an igraph object with %d nodes and %d edges", vcount(x$g), ecount(x$g)), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print(x$priority[1:2,], row.names=FALSE)
}
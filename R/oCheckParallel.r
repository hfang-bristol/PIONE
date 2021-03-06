#' Function to check whether parallel computing should be used and how
#'
#' \code{oCheckParallel} is used to check whether parallel computing should be used and how
#'
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return TRUE for using parallel computing; FALSE otherwise
#' @note
#' Whether parallel computation with multicores is used is system-specific. Also, it will depend on whether these two packages "foreach" and "doParallel" have been installed.
#' @export
#' @seealso \code{\link{oCheckParallel}}
#' @include oCheckParallel.r
#' @examples
#' \dontrun{
#' oCheckParallel()
#' }

oCheckParallel <- function(multicores=NULL, verbose=TRUE)
{
    
    # @import doParallel
    # @import foreach
    
    flag_parallel <- FALSE
    pkgs <- c("doParallel","foreach")
    if(all(pkgs %in% rownames(utils::installed.packages()))){
        tmp <- sapply(pkgs, function(pkg) {
            #suppressPackageStartupMessages(require(pkg, character.only=TRUE))
            requireNamespace(pkg, quietly=TRUE)
        })
        
        if(all(tmp)){
            doParallel::registerDoParallel()
            cores <- foreach::getDoParWorkers()
            if(is.null(multicores)){
                multicores <- max(1, ceiling(cores))
            }else if(is.na(multicores)){
                multicores <- max(1, ceiling(cores))
            }else if(multicores < 1 | multicores > 2*cores){
                multicores <- max(1, ceiling(cores))
            }else{
                multicores <- as.integer(multicores)
            }
            doParallel::registerDoParallel(cores=multicores) # register the multicore parallel backend with the 'foreach' package
            
            if(verbose){
                message(sprintf("\tdo parallel computation using %d cores ...", multicores, as.character(Sys.time())), appendLF=TRUE)
            }
            flag_parallel <- TRUE
        }
        
    }
    
    return(flag_parallel)
}

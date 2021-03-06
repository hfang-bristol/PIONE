#' Function to create a sparse matrix for an input file with three columns
#'
#' \code{oSparseMatrix} is supposed to create a sparse matrix for an input file with three columns.
#'
#' @param input.file an input file containing three columns: 1st column for rows, 2nd for columns, and 3rd for numeric values. Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param rows a vector specifying row names. By default, it is NULL
#' @param columns a vector specifying column names. By default, it is NULL
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' an object of the dgCMatrix class (a sparse matrix)
#' @note If rows (or columns) are not NULL, the rows (or columns) of resulting sparse matrix will be union of those from input.file and those from rows (or columns).
#' None
#' @export
#' @seealso \code{\link{oSparseMatrix}}
#' @include oSparseMatrix.r
#' @examples
#' # create a sparse matrix of 4 X 2
#' input.file <- rbind(c('R1','C1',1), c('R2','C1',1), c('R2','C2',1), c('R3','C2',1), c('R4','C1',1))
#' res <- oSparseMatrix(input.file)
#' res
#' # get a full matrix
#' as.matrix(res)
#'
#' res <- oSparseMatrix(input.file, columns=c('C1','C2','C3'))
#' res

oSparseMatrix <- function(input.file, rows=NULL, columns=NULL, verbose=TRUE)
{
    
    if(is.matrix(input.file) | is.data.frame(input.file)){
        if(verbose){
            message(sprintf("Load the input data matrix of %d X %d ...", dim(input.file)[1], dim(input.file)[2]), appendLF=TRUE)
        }
        
        if(ncol(input.file)==2){
            input.file <- cbind(input.file, rep(1,nrow(input.file)))
        }
        
        if(is.data.frame(input.file)){
            input <- cbind(as.character(input.file[,1]), as.character(input.file[,2]), as.character(input.file[,3]))
        }else{
            input <- input.file
        }
    }else if(is.character(input.file) & sum(input.file!='')){
        if(verbose){
            message(sprintf("Read the input file '%s' ...", input.file), appendLF=TRUE)
        }
        input <- utils::read.delim(input.file, header=TRUE, sep="\t", colClasses="character")
        if(ncol(input)==2){
            input <- cbind(input, rep(1,nrow(input)))
        }
    }else{
        return(NULL)
    }
    
    x <- input
    if(!is.null(x)){
    	
    	## row names
    	if(is.null(rows)){
        	x_row <- sort(unique(x[,1]))
        }else{
        	#x_row <- sort(union(rows, unique(x[,1])))
        	x_row <- rows
        }
        
    	## column names
    	if(is.null(columns)){
        	x_col <- sort(unique(x[,2]))
        }else{
        	#x_col <- sort(union(columns, unique(x[,2])))
        	x_col <- columns
        }
        
        #############################
        ind_row <- match(x[,1], x_row)
        ind_col <- match(x[,2], x_col)
        x.sparse <- Matrix::sparseMatrix(i=ind_row, j=ind_col, x=as.numeric(x[,3]), dims=c(length(x_row),length(x_col)))
        rownames(x.sparse) <- x_row
        colnames(x.sparse) <- x_col
        
        if(verbose){
            message(sprintf("There are %d entries, converted into a sparse matrix of %d X %d.", nrow(x), dim(x.sparse)[1], dim(x.sparse)[2]), appendLF=TRUE)
        }
    }else{
        return(NULL)
    }
    
    invisible(x.sparse)
}

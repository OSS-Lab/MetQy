#' Creates labels to add to a character vector to make all entries unique
#'
#' Creates a label vector of length N, such that they are unique and help distinguish them
#'
#' @param N - numeric. Indicates the length of the label vector to be return.
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @return A character vector of length N with unique set of alpha-numeric labels.
#'
#' @export

############################################################################################################################################

misc_create_labels <- function(N,...){
  stopifnot(is.numeric(N))
  stopifnot(length(N)==1)

  if(N>length(LETTERS)){
    how_much_larger <- ceiling(N/length(LETTERS)) #N - length(LETTERS)

    use_Labels <- NULL
    for(h in 1:how_much_larger)  use_Labels <- c(use_Labels, paste(LETTERS,h,sep=""))
    use_Labels <- use_Labels[1:N]
  } else {
    use_Labels <- paste(LETTERS,1,sep="")
    use_Labels <- use_Labels[1:N]
  }
  return(use_Labels)
}

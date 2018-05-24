#' Check for duplicated strings and append letter code to distinguish them.
#'
#' @param NAME_VECTOR - character vector.
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @export
#'

############################################################################################################################################

misc_check_duplicate_names <- function(NAME_VECTOR,...){
  duplicated_rows <- unique(NAME_VECTOR[which(duplicated(NAME_VECTOR))])
  if(length(duplicated_rows)>0){
    # warning("Duplicated dataset names provided. Appending \\\\s**[A-Z]\\\\d to distinguish them ")

    for(D in duplicated_rows){
      index <- which(NAME_VECTOR==D)
      NAME_VECTOR[index] <- paste(NAME_VECTOR[index]," **",misc_create_labels(length(index)),sep="")
    }
  }
  return(NAME_VECTOR)
}

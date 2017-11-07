#' Scatter plot with overlapping factors/groups.
#'
#' Represent different groups as defined by FACTOR(S) in a scatter plot.
#' This function is used by \link{analysis_genomes_module_output()} to plot the first two Principal Components (PCs) from the PCA analysis,
#' overlapping different factors or groups (as determined by \code{FACTOR}).
#' It is also used by \link{analysis_pca_mean_distance_grouping()} to highlight a single group on PC plot.
#'
#' @param MATRIX      - two column matrix to plot. Only the first two columns will be used.
#'
#' @param FACTOR      - character vector or list of character vectors used to split the data into groups.
#'
#' @param factor_labs - optional. Character vector to distinguish \code{FACTOR} groups. Default (\code{NULL}). See Details.
#'
#' @param Filename    - optional. Character vector containing the file path(s) to save image(s) WITH extension.
#' The device to save is determined by file extension. Default  ("plot_scatter_byFactors.pdf").
#' If set to \code{""}, no file will be written.
#'
#' @details \code{factor_labs} is used as an extension for the filename for the plot files.
#' The names of the \code{FACTOR} object OR "factor" followed by a number
#'  will be used if \code{factor_labs} is not specified (i.e. \code{factor_labs = NULL}).
#'
#' @return
#' A list with as many entires as factors is generated (one for each factor) using \code{factor_labs} for the names.
#' For every \code{FACTOR}, a list will contain:
#' \preformatted{
#'  $nGroups - numeric. The number of groups found for that FACTOR (NOTE: That entries with "" will be excluded.)
#'  $table   - data frame with the number of entries for each group found.
#'  $file    - character vector of file name(s).
#' }
#'
#' @export
#'
#' @seealso \link{analysis_genomes_module_output},\link{analysis_pca_mean_distance_grouping}

############################################################################################################################################

plot_scatter_byFactors <- function(MATRIX,FACTOR,factor_labs = NULL,Filename = "plot_scatter_byFactors.pdf",Width = 9, Height = 7){

  #### MANAGE INPUT ----
  stopifnot(is.matrix(MATRIX)||ncol(MATRIX)<2)
  stopifnot(is.character(FACTOR)||is.list(FACTOR))
  if(ncol(MATRIX)>2) MATRIX <- MATRIX[,1:2]
  # FORMAT DECIMAL PLACES
  MATRIX <- matrix(as.numeric(format(MATRIX, digits = 3)),dim(MATRIX))

  # FILENAME
  if(Filename!=""&&is.character(Filename)){
    extension <- gsubfn::strapplyc(Filename,"([.][a-z]+)$",simplify = T)
    Filename_temp <- Filename
    for(E in extension) Filename_temp <- gsub(E,"",Filename_temp)
  }

  # DEFINE LISTS OF FACTORS
  if(is.character(FACTOR)){
    group_by_list <- list(FACTOR = FACTOR)
  }
  if(is.list(FACTOR)){
    group_by_list <- FACTOR
  }
  if(!exists("group_by_list")) stop("FACTOR not recognized. Most be a character vector or a list.")

  # ADD NAMES TO LISTS
  nFactors <- length(group_by_list)

  if(is.null(factor_labs)||nFactors!=length(factor_labs)){
    if(!is.null(names(group_by_list))){
      factor_labs <- names(group_by_list)
    } else {
      factor_labs <- paste("factor",1:nFactors,sep="")
    }
  }

  # OUTPUT
  output        <- vector(mode = "list",nFactors)
  names(output) <- factor_labs

  # if(!exists("width"))  width <- 9
  # if(!exists("height")) height <- 7

  for(G in 1:nFactors){
    this_factor   <- group_by_list[[G]]
    Groups        <- sort(unique(this_factor))
    if(length(which(Groups==""))>0) Groups <- Groups[-which(Groups=="")]
    nGroups       <- length(Groups)
    group_counts  <- numeric(nGroups)
    # OUTPUT
    output[[G]]$nGroups <- nGroups

    # GROUP COUNTS
    for(n in 1:nGroups) group_counts[n] <- length(which(this_factor==Groups[n]))

    # SAVE DATA
    output[[G]]$table   <- data.frame(Groups = Groups,  group_counts = group_counts,stringsAsFactors = F)

    group_legend  <- paste(Groups," (",group_counts,")",sep="")
    maxChar       <- max(nchar(group_legend))

    # PLOT ONLY IF LESS THAN 46 GROUPS
    if(nGroups<46){
      if(is.character(this_factor)&&length(this_factor)==nrow(MATRIX)){
        # FILE
        if(Filename!=""&&is.character(Filename)){
          # extension <- gsubfn::strapplyc(Filename,"[.][a-z]+",simplify = T)
          # Filename_temp <- Filename
          # for(E in extension) Filename_temp <- gsub(E,"",Filename_temp)
          File <-  paste(Filename_temp,"_",factor_labs[G],extension,sep="")
          output[[G]]$file <- File


          for(FF in File){
            pdf(file = FF,width = Width,height = Height)
            par(xpd=TRUE)
            # par(mar = c(5.1, 4.1, 1, 2.1+maxChar/3.5))
            if(nGroups > 1) {
              # room for legend
              par(mar = c(5.1, 5.1, 1, 11.25)) #c(bottom, left, top, right)
            }else{
              par(mar = c(5.1, 5.1, 1.1, 1.1)) 
            }
            plot(MATRIX,xlab = "PC1",ylab = "PC2",pch = 16, col = "grey95") # PC1 v PC2

            ### FACTOR OF ALL DATASETS ----
            if(nGroups < 3) {
              COLOURS <- c("red","blue")
            } else {
              if(nGroups < 13){
                COLOURS     <- RColorBrewer::brewer.pal(nGroups,"Paired")
              } else {
                COLOURS     <- rainbow(nGroups,end = 0.9)
              }
            }

            for(n in 1:nGroups){
              group_index <- which(this_factor==Groups[n])
              group_counts[n] <- length(group_index)
              points(MATRIX[group_index,1],MATRIX[group_index,2],col = COLOURS[n],pch = 16)
            }

            if(nGroups > 1) legend(max(MATRIX[,1])*1.075,max(MATRIX[,2])*1.1,legend = group_legend,pch = 16,col = COLOURS, horiz = F, cex = 0.7, bty = "n")
            par(xpd=FALSE)
            dev.off()
            par(mar = c(5.1, 4.1, 4.1, 2.1))
          }
        }
      }
    } else {
      # output[[G]]$table  <- NULL
      output[[G]]$file   <- NULL
      warning(paste("FACTOR '",factor_labs[G], "' skipped. Too many groups for plotting: ",nGroups))
      next
    }
  }

  return(output)
}

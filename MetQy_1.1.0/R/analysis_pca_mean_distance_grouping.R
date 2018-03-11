#' Given a set of coordinates and a grouping variable, the mean distance is calculated for each group.
#'
#' @param MATRIX           - matrix with the point coordinates (two column minimum).
#'
#' @param FACTOR           - used to split the data into groups. Character vector (or list of vetors) with the same length as rows
#'                           in \code{MATRIX}.
#'
#' @param nDim             - optional. Numeric vector of dimensions to use for the point coordinates. Default (\code{2}).
#'
#' @param factor_labs      - optional. Character vector of the length of \code{FACTOR} if \code{FACTOR} is a list or
#'                           length 1 if \code{FACTOR} is a character vector. Default (\code{NULL}; \code{FACTOR} names or letters used).
#'
#' @param plot_mean_dist   - logical. Should a plot of the mean distances be generated? Default (\code{TRUE}).
#'
#' @param Filename         - optional. Filename with path and extension. String added to distinguish \code{FACTOR} groups.
#'                            Default ("plot_mean_dist.pdf" saved to working directory).
#'
#' @param plot_top_percent - optional. Numeric vector between 0 and 1 or integer. Default (NULL, i.e. not plotted).
#'
#' @param ...              - further arguments for \code{plot_scatter()}.
#'
#' @details
#' The mean distance of groups with one member (derived from the \code{FACTOR} labels) cannot be calculated.
#'
#' \code{nDim} must be larger than 2.
#'
#' \code{plot_top_percent} refers to the fraction or number of groups with the highest calculated mean distance
#' that should be plotted on a scatter plot. Items are recycled if not as many as the number of grouping factors is provided.
#'
#' @return List containing a table of the mean distances (columns for group name, mean distance and size of group (N))
#' and a ggplot object of the data (NAs removed) if \code{plot_mean_dist} is TRUE.
#'
#' @seealso \link{analysis_pca_mean_distance_calculation}
#'
#' @examples
#' data(data_example_moduleIDs)
#' data(data_example_genomeIDs)
#'
#' # Calculate the module completion fraction (mcf) for the genomes and modules contained in the data objects above.
#' OUT         <- query_genomes_to_modules(data_example_genomeIDs,MODULE_ID = data_example_moduleIDs)
#'
#' pca <- prcomp(OUT$MATRIX)
#'
#' # Group data
#' this_FACTOR      <- rep(LETTERS[1:5],length(data_example_genomeIDs)/5)
#' mean_dist_output <- analysis_pca_mean_distance_grouping(pca$x,this_FACTOR,xLabs_angle = F,Width = 2, Height = 1.5,Filename = "../GitHub/MetQy/fig/plot_pca_scatter.png")
#'
#' @export

############################################################################################################################################

analysis_pca_mean_distance_grouping <- function(MATRIX,FACTOR,factor_labs = NULL,nDim = 2, plot_mean_dist = T, Filename = "plot_mean_dist.pdf",plot_top_percent = NULL,...){

  #### MANAGE INPUT ----
  stopifnot(is.matrix(MATRIX),ncol(MATRIX)>1)
  stopifnot(is.character(FACTOR)||is.list(FACTOR))
  stopifnot(nDim>1)

  # FILENAME
  if(Filename!=""&&is.character(Filename)){
    extension <- gsubfn::strapplyc(Filename,"[.][a-z]+",simplify = T)
    Filename_temp <- Filename
    for(E in extension) Filename_temp <- gsub(E,"",Filename_temp)
  }

  # MATCH LENGTH plot_top_percent with length FACTOR
  if(length(plot_top_percent)<length(FACTOR))
    plot_top_percent <- rep(plot_top_percent,1+ceiling((length(FACTOR)-length(plot_top_percent))/length(plot_top_percent)))

  # DEFINE LISTS OF FACTORS
  if(is.character(FACTOR)){
    group_by_list <- list(FACTOR = FACTOR)
  }
  if(is.list(FACTOR)){
    group_by_list <- FACTOR
  }

  output <- vector(mode = "list",length(group_by_list))

  # ADD NAMES TO LISTS
  nFactors <- length(group_by_list)
  if(is.null(factor_labs)||nFactors!=length(factor_labs)){
    if(!is.null(names(group_by_list))){
      factor_labs <- names(group_by_list)
    } else {
      factor_labs <- LETTERS[1:length(group_by_list)]
    }
  }
  names(output) <- factor_labs

  for(G in 1:length(group_by_list)){
    this_factor <- group_by_list[[G]]
    # CHECK DIMENSIONS
    stopifnot(length(this_factor)==nrow(MATRIX))

    Groups  <- sort(unique(this_factor))
    nGroups <- length(Groups)
    mean_dist_table <- data.frame(GROUPS = Groups,
                                  N = NA,
                                  # mean_dist = NA,
                                  stringsAsFactors = F)
    nDim_col_index <- NULL
    for(D in 1:length(nDim)){
      nDim_col_index <- c(nDim_col_index,ncol(mean_dist_table)+1)
      mean_dist_table[[ncol(mean_dist_table)+1]] <- NA
      names(mean_dist_table)[ncol(mean_dist_table)] <- paste("mean_dist_nDim_",nDim[D],sep="")
    }

    # GET GROUPS AND CALCULATE MEAN DISTANCE
    for(n in 1:nGroups){
      group_index <- which(this_factor==Groups[n])
      mean_dist_table$N[n] <- length(group_index)
      if(length(group_index)>1){
        for(D in 1:length(nDim)){
          mean_dist_table[n,nDim_col_index[D]] <- analysis_pca_mean_distance_calculation(MATRIX[group_index,1:nDim[D]])
        }
      }
    }

    # NORMILIZE BY NUMBER OF DATASETS IN GROUP (N)
    # mean_dist_table$mean_dist_norm <- mean_dist_table$mean_dist/mean_dist_table$N

    # ORDER BY INCREASING MEAN DISTANCE
    mean_dist_table <- mean_dist_table[order(mean_dist_table[,nDim_col_index[1]]),]

    output[[G]]$mean_dist_table <- mean_dist_table

    # PLOT MEAN DISTANCE
    if(plot_mean_dist){
      # MAIN FILE NAME
      File <-  paste(Filename_temp,"_",factor_labs[G],extension,sep="")
      output[[G]]$mean_dist_file <- File

      # DATA WITHOUT NAs
      plot_DATA <- mean_dist_table[which(mean_dist_table$N>1),] # data; remove NAs
      plot_DATA <- plot_DATA[,-2] # remove column with number of group members

      if(ncol(plot_DATA)>2){
        order_index <- order(apply(plot_DATA[,2:ncol(plot_DATA)],1,mean))
      }else{
        order_index <- order(plot_DATA[,2])
      }

      plot_DATA <- plot_DATA[order_index,]
      # head(plot_DATA)
      plot_data <- reshape2::melt(plot_DATA,id = "GROUPS")
      plot_data[,2] <- gsub("mean_dist_","",plot_data[,2])

      ## PLOT
      # y label
      if(length(nDim)>1){
        yLab = "Mean distance"
      }else{
        yLab = paste("Mean distance (nDim = ",nDim,")",sep="")
      }
      plot_scatter(plot_data = plot_data,colBy = 2,xLab = paste("Groups (by ",factor_labs[G],")",sep=""),yLab = yLab,
                   File = output[[G]]$mean_dist_file,...)
      # ZOOM IN ON TAIL
      if(nrow(plot_DATA)>50){
        plot_data <- reshape2::melt(plot_DATA[(nrow(plot_DATA)-19):nrow(plot_DATA),],id = "GROUPS")
        plot_data[,2] <- gsub("mean_dist_","",plot_data[,2])

        output[[G]]$mean_dist_file_tail <- paste(Filename_temp,"_",factor_labs[G],"_tail",extension,sep="")
        plot_scatter(plot_data = plot_data,colBy = 2,xLab = paste("Groups (by ",factor_labs[G],", top 20)",sep=""),
                     yLab = yLab,File = output[[G]]$mean_dist_file_tail,...)
      }

      # TOP PERCENT
      if(is.numeric(plot_top_percent[G])){
          remove_index <- c(which(mean_dist_table$N<2)) # ,which(mean_dist_table$GROUPS=="UNKNOWN")

          if(length(remove_index)>0){
            top_data <- mean_dist_table[-remove_index,] # data; remove NAs
          }else{
            top_data <- mean_dist_table
          }
          # OPTION TO SPECIFY % PER FACTOR

          if(plot_top_percent[G]>0&&plot_top_percent[G]<1){
            N          <- ceiling(nrow(top_data)*plot_top_percent[G])
            top_groups <- as.character(top_data$GROUPS[(nrow(top_data)-N+1):nrow(top_data)]) # ordered
          } else {
            if(plot_top_percent[G]<nrow(top_data)&&plot_top_percent[G]!=0){
              N           <- ceiling(plot_top_percent[G])
              top_groups  <- as.character(top_data$GROUPS[(nrow(top_data)-N+1):nrow(top_data)])

            } else {
             top_groups <- ""
             if(!plot_top_percent[G]<nrow(top_data)) warning(paste("plot_top_percent ",G," larger than number of datasets with 'mean dist'. Scatter plot top percent not done for factor",factor_labs[G],sep=""))
           }
          }
         # cat("top_groups", top_groups,fill = T)
         if(top_groups[1] != "") {
             output[[G]]$top_files <- paste(Filename_temp,"_",factor_labs[G],"_scatter_top_",plot_top_percent[G],"_",length(top_groups):1,"_",top_groups,extension,sep="")
           for(FF in 1:length(top_groups)){
             # FF = 1
             temp_factor <- as.character(this_factor)
             temp_factor[which(this_factor!=top_groups[FF])] <- "" # PLOT ONLY ONE GROUP ON TOP OF SCATTER OF ALL DATA

             plot_scatter_byFactors(MATRIX,FACTOR = temp_factor,factor_labs = top_groups[FF],
                                    Filename = paste(Filename_temp,"_",factor_labs[G],"_scatter_top_",plot_top_percent[G],"_",c(length(top_groups):1)[FF],extension,sep=""))
           }
         }
      }

    }else{
      output[[G]]$mean_dist_file <- output[[G]]$mean_dist_file_end <- NULL
    }

  }
  return(output)
}

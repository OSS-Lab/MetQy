#' Scatter plot with overlapping factors/groups.
#'
#' Represent different groups as defined by FACTOR(S) in a scatter plot.
#'
#' @details
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
#' @param Filename    - optional. Character vector containing the file path to save the image. Must have an extension (see Details).
#'                          Default ("plot_scatter_byFactors.pdf"). If set to \code{""}, no file will be written.
#'
#' @param Width     - Width size for file. Default (7 in).
#'
#' @param Height    - Height size for file. Default (5 in).
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @details \code{factor_labs} is used as an extension for the filename for the plot files.
#' The names of the \code{FACTOR} object OR "factor" followed by a number
#' will be used if \code{factor_labs} is not specified (i.e. \code{factor_labs = NULL}).
#'
#' A plot is only generated for FACTORS with less than 60 groups.
#' All the data is plotted in grey in the background, with groups being overlayed for each Factor.
#'
#'  Only the following file types are allowed: pdf, png, svg and jpeg.
#'
#' @return
#' A list with as many entires as factors is generated (one for each factor) using \code{factor_labs} for the names.
#' For every \code{FACTOR}, a list will contain:
#' \preformatted{
#'  $nGroups - numeric. The number of groups found for that FACTOR
#'                      NOTE: That entries with "" will be excluded.
#'  $table   - data frame with the number of entries for each group found.
#'  $file    - character vector of file name(s).
#' }
#'
#' @export
#'
#' @examples
#' data(data_example_moduleIDs)
#' data(data_example_genomeIDs)
#' # length(data_example_genomeIDs) # [1] 25
#'
#' # Calculate the module completion fraction (mcf) for the genomes
#' #       and modules contained in the data objects above.
#' OUT         <- query_genomes_to_modules(data_example_genomeIDs,
#'                                         MODULE_ID = data_example_moduleIDs)
#'
#' pca <- prcomp(OUT$MATRIX)
#'
#' # Make boxplots of the mcf output from query_genomes_to_modules
#' this_FACTOR <- rep(LETTERS[1:5],length(data_example_genomeIDs)/5)
#' plot_output <- plot_scatter_byFactors(pca$x[,1:2],FACTOR = this_FACTOR,
#'                                       factor_labs = "random",
#'                                       Filename = "plot_scatter_byFactors.png")
#'
#' # NAs are ommitted, so a single group can be contrasted with overall data
#' this_FACTOR <- c(rep(NA,20),rep(LETTERS[1],5))
#' plot_output <- plot_scatter_byFactors(pca$x[,1:2],FACTOR = this_FACTOR,
#'                                       factor_labs = "group_A",
#'                                       Filename = "plot_scatter_byFactors_single.png")
#'
#' @seealso \link{analysis_genomes_module_output},\link{analysis_pca_mean_distance_grouping}

############################################################################################################################################

plot_scatter_byFactors <- function(MATRIX,FACTOR,factor_labs = NULL,xLab = NULL,yLab = NULL,Filename = "plot_scatter_byFactors.pdf",Width = 7, Height = 5,...){

  #### MANAGE INPUT ----
  stopifnot(is.matrix(MATRIX)||ncol(MATRIX)<2)
  stopifnot(is.character(FACTOR)||is.list(FACTOR))

  # DATA AND AXIS LABLES ----
  if(is.null(xLab)) xLab <- colnames(MATRIX)[1]
  if(is.null(yLab)) yLab <- colnames(MATRIX)[2]
  if(ncol(MATRIX)>2) MATRIX <- MATRIX[,1:2]
  # FORMAT DECIMAL PLACES
  MATRIX <- matrix(as.numeric(format(MATRIX, digits = 3)),dim(MATRIX))


  # FILENAME ----
  MAKE_PLOT <- FALSE
  if(Filename!=""&&is.character(Filename)){
    if(length(Filename)>1){
      Filename <- Filename[1]
      warning("Only first Filename being used.")
    }

    extension <- gsubfn::strapplyc(Filename,"([.][a-z]+)$",simplify = T)
    if(length(extension)<1){
      warning("Extension not identified. Plot will not be generated.")
      MAKE_PLOT <- FALSE
    }else{
      MAKE_PLOT <- TRUE
    }

    if(length(match(extension,c(".pdf",".png",".svg",".jpeg")))<1){
      extension <- ".pdf"
      warning(paste(extension,"is not a valid file type. Generating pdf file."))
    }

    if(MAKE_PLOT){
      # TRUNCATE FILENAME TO ADD FACTOR GROUP INFO
      Filename_temp <- Filename
      Filename_temp <- gsub(extension,"",Filename_temp)
    }
  }

  # DEFINE LISTS OF FACTORS ----
  if(is.character(FACTOR)){
    if(length(FACTOR)>nrow(MATRIX)){
      FACTOR <- FACTOR[1:nrow(MATRIX)]
      warning("More")
    }
    group_by_list <- list(FACTOR = FACTOR)
  }
  if(is.list(FACTOR)){
    group_by_list <- FACTOR
  }
  if(!exists("group_by_list")) stop("FACTOR not recognized. Most be a character vector or a list.")

  # ADD NAMES TO LISTS ----
  nFactors <- length(group_by_list)

  if(is.null(factor_labs)||nFactors!=length(factor_labs)){
    if(!is.null(names(group_by_list))){
      factor_labs <- names(group_by_list)
    } else {
      factor_labs <- paste("factor",1:nFactors,sep="")
    }
  }

  # OUTPUT ----
  output        <- vector(mode = "list",nFactors)
  names(output) <- factor_labs

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
    if(nGroups<60){
      if(is.character(this_factor)&&length(this_factor)==nrow(MATRIX)){
        # FILE
        if(MAKE_PLOT){
          File <-  paste(Filename_temp,"_",factor_labs[G],extension,sep="")

          output[[G]]$file <- File

          if(length(grep("pdf",extension)>0)){
            pdf(file = File,width = Width,height = Height)
          }else if(length(grep("png",extension)>0)){
            png(filename = File,width = Width,height = Height,units = "in",res = 250)
          }else if(length(grep("svg",extension)>0)){
            svg(filename = File,width = Width,height = Height)
          }else if(length(grep("jpeg",extension)>0)){
            jpeg(filename = File,width = Width,height = Height,units = "in",res = 250)
          }

          par(xpd=TRUE) # all plotting is clipped to the figure region

          # MARGIN ----
          if(nGroups > 1) {
            nChars <- max(nchar(group_legend))
            if(nChars<0) nChars <- 0

            # room for legend
            par(mar = c(4.5, 4.5, 1.1, 1.5+nChars/3)) #c(bottom, left, top, right)
          }else{
            par(mar = c(4.5, 4.5, 1.1, 1.1))
          }

          ### COLOUR FACTOR OF ALL DATASETS ----
          if(nGroups < 3){
            COLOURS <- c('#e31a1c','#1f78b4')# red and blue
          }else if(nGroups < 8){
              # COLOURS     <- RColorBrewer::brewer.pal(nGroups,"Paired")
              COLOURS <- c('#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a',"#B15928","#000000")
              COLOURS <- COLOURS[1:nGroups]
          }else if(nGroups >=8){
              COLOURS     <- rainbow(nGroups,end = 0.9)
          }

          # BACKGROUND PLOT --> ALL DATA ------ (PC1 v PC2)
          plot(MATRIX,xlab = xLab,ylab = yLab,pch = 16, col = "grey80",
               xlim = c(min(MATRIX[,1])*0.99,max(MATRIX[,1])*1.01),
               ylim = c(min(MATRIX[,2])*0.99,max(MATRIX[,2])*1.01))


          for(n in 1:nGroups){
            group_index <- which(this_factor==Groups[n])
            group_counts[n] <- length(group_index)
            points(MATRIX[group_index,1],MATRIX[group_index,2],col = COLOURS[n],pch = 16)
          }

          if(nGroups > 1)
            legend(x = max(MATRIX[,1])*1.15,y = max(MATRIX[,2])*1.05,legend = group_legend,pch = 16,col = COLOURS,
                   horiz = F, cex = ifelse(nGroups<46,0.7,0.55), bty = "n")

          par(xpd=FALSE)

          while(dev.cur()>1) dev.off()
          # RESTORE DEFAULT
          par(mar = c(5.1, 4.1, 4.1, 2.1))
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

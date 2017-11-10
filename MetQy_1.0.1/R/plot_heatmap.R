#' Plot a heatmap
#' 
#' Plot a heatmap of the the module fraction present across genomes or of the variance across groups.
#'
#' @param data_in     - numeric matrix or 3-column data frame containing the module fraction match of the genomes and modules. See Details.
#'
#' @param row_lab     - optional. String to specify the axis label corresponding to the row values. Default ("Genomes").
#'
#' @param col_lab     - optional. String to specify the axis label corresponding to the column values. Default ("Modules").
#'
#' @param ORDER_MATRIX - logical. Should rows and columns be reorder according to the dendogram (hierarchical clustering). Default (\code{FALSE}).
#'
#' @param legend_name - optional. String to specify the legend title. Default ("Fraction\\n matched\\n").
#'
#' @param set_yLim    - optional. Logical or numeric (length 2) indicating whether to set the y limit of the heatmap scale. Default (\code{FALSE}). See Details.
#'
#' @param Filename    - optional. A character vector containing the file path to save image(s). The device to save is determined by file extension. Default  ("", i.e. file not written).
#'
#' @param Width       - optional. Width  size for file. Default (24 in).
#'
#' @param Height      - optional. Height size for file. Default (18 in).
#'
#' @details 
#' If \code{data_in} is a data frame, the heatmap will be made using the first column as the row entry labels, the second column as the column entry labels 
#' and the third as the actual value to be plotted. 
#' 
#' If a numeric matrix is provided, it will be 'melted' into this type of data frame by using the \code{reshape2::melt} function.
#' 
#' \code{set_yLim} - #' the default (\code{FALSE}), sets the y limit to \code{c(0, 1)}, while \code{TRUE} to the \code{c(min, max)} of the data.
#' If numeric, the values given are used for the lower and upper limits respectively.
#'
#' @return  ggplot object. Saves figures if Filename is specified.
#'
#' @export

############################################################################################################################################

plot_heatmap <- function(data_in, row_lab = "Genomes", col_lab = "Modules",ORDER_MATRIX = FALSE,legend_name = "Module\ncompleteness\nfraction\n", set_yLim = FALSE, Filename ="", Width = 24, Height = 18){

  # MANAGE INPUT ----
  if(!is.character(row_lab)||!length(row_lab)==1){
    warning("Incorrect 'row_lab' format. Using default values")
    row_lab <- "Genomes"
  }
  if(!is.character(col_lab)||!length(col_lab)==1){
    warning("Incorrect 'col_lab' format. Using default values")
    col_lab <- "Modules"
  }
  # set_yLim
  stopifnot((is.logical(set_yLim)&&length(set_yLim)==1)||(is.numeric(set_yLim)&&length(set_yLim)==2))

  # data_in - a 3-column DATA FRAME
  if(is.data.frame(data_in)){
    if(ncol(data_in)!=3){
      stop("'data_in' must be either a 3-column data frame or a numeric matrix")
    } else {
      if(is.numeric(data_in[,3])){
        DATA <- data_in
      } else {
        stop("if is.data.frame(data_in), the thrid column must be numeric ")
      }
    }
  }

  # data_in - a numeric MATRIX   - TRANSFORM to data frame (for ggplot2)
  # if(is.matrix(data_in)&&is.numeric(data_in)) {
  if(is.matrix(data_in)) {
    if(ORDER_MATRIX){
      # Use a hierarchical cluster to order the genomes
      if(nrow(data_in)>2){ reorder_rows <- order.dendrogram(as.dendrogram(hclust(dist(data_in))))   }else{reorder_rows <- 1:nrow(data_in)}
      # Use a hierarchical cluster to order the modules
      if(ncol(data_in)>2){ reorder_cols <- order.dendrogram(as.dendrogram(hclust(dist(t(data_in)))))}else{reorder_cols <- 1:ncol(data_in)}

      data_in <- data_in[reorder_rows,reorder_cols]
    }
    DATA <- reshape2::melt(data_in)
  }

  if(!exists("DATA")) stop(paste("'data_in' must be either a 3 column data frame or a matrix. 'data_in' is",data.class(data_in)))

  # CONSERVE ORDER ----
  DATA[,1] <- as.character(DATA[,1])
  DATA[,1] <- factor(DATA[,1], levels = unique(DATA[,1]))

  DATA[,2] <- as.character(DATA[,2])
  DATA[,2] <- factor(DATA[,2], levels = unique(DATA[,2]))

  if(row_lab != "") names(DATA)[1] <- row_lab
  if(col_lab != ""&&is.character(row_lab)&&length(row_lab)==1) names(DATA)[2] <- col_lab

  ## y lim
  if(is.logical(set_yLim)){
    yLim_min <- ifelse(set_yLim,min(DATA[,3]),0)
    yLim_max <- ifelse(set_yLim,max(DATA[,3]),1)
  }
  if(is.numeric(set_yLim)){
    yLim_min <- set_yLim[1]
    yLim_max <- set_yLim[2]
  }
  if(!exists("yLim_max")){
    yLim_min <- 0
    yLim_max <- 1
  }

  #### PLOT ----
  p <-
    ggplot2::ggplot(DATA, ggplot2::aes(DATA[,2],DATA[,1])) +
    ggplot2::geom_tile(data = DATA,ggplot2::aes(DATA[,2],DATA[,1],fill = DATA[,3])) +
    # ggplot2::scale_fill_gradient(low = "white", high = "black",guide = "colourbar", name = legend_name, limits = c(0,1)) +
    ggplot2::scale_fill_distiller(palette ="Spectral",guide = "colourbar", name = legend_name, limits = c(yLim_min,yLim_max)) +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill = 'white', colour = 'white')) +
    # AXES
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                  axis.title.y	= ggplot2::element_text(size = Width*1.3,margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)),
                  axis.title.x	= ggplot2::element_text(size = Width*1.3,margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0))) +
    ggplot2::xlab(names(DATA)[2]) +
    ggplot2::ylab(names(DATA)[1]) +
    # LEGEND
    ggplot2::theme(legend.position = "right",  # Legend
                legend.key = ggplot2::element_rect(colour = 'black'),
                legend.key.size = ggplot2::unit(Width/14, "cm"),
                legend.text  = ggplot2::element_text(size=Width*0.95),
                legend.title = ggplot2::element_text(size=Width*1))

  # MANAGE X-AXIS LABELS ----
  if(length(unique(DATA[,2]))>50){
    p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_blank())
  }else{
    nChars <- max(nchar(as.character(unique(DATA[,2])))[1:ifelse(length(unique(DATA[,2]))>2,3,length(unique(DATA[,2])))])-max(nchar(as.character(DATA[,1])))
    if(!exists("nChar")) nChar <- -1
    if(nChars<0) nChars <- 0
    
    p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_text(hjust = 1,vjust=1,size=Width*1.2, colour = "black",angle = 45)) +
              ggplot2::theme(plot.margin  = ggplot2::unit(c(0.3,0.5,1.5,0.5+nChars/3.5),"cm"))   # top, right, bottom, and left
  }

  # MANAGE Y-AXIS LABELS ----
  if(length(unique(DATA[,1]))>50){
    p <- p + ggplot2::theme(axis.text.y   = ggplot2::element_blank())
  }else{
    p <- p + ggplot2::theme(axis.text.y   = ggplot2::element_text(size=Width*1.2, colour = "black"))
  }

  # Save Figure ----
  if(Filename!=""&&is.character(Filename)){
    for(File in Filename){
      ggplot2::ggsave(filename =File,plot = p,width = Width, height = Height)
    }
  }

  # RETURN PLOT
  return(p)
}

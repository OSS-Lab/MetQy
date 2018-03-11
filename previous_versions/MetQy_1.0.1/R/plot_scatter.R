#' Scatter plot
#'
#' Scatter plot for pairs of categorical data with a numeric value.
#' 
#' @details 
#' This function is used by \link{analysis_genomes_module_output()} to plot the module variance accross genomes and
#' 
#' by \link{analysis_pca_mean_distance_grouping()} to plot the mean PCA distance of the groups analysed.
#'
#' @param plot_data - data frame. The columns to be plotted are indicated with \code{X} and \code{Y}.Order given is conserved.
#'
#' @param X,Y       - optional. Number indicating column to use for x axis and y axis, respectively. Default (1 and 3).
#'
#' @param colBy     - optional. Number indicating column to use for coloring grouping. Default (\code{NULL}).
#'
#' @param xLab,yLab - optional. String to use for x label and y label, respectively. Default 
#'                      (first and second column names, respectively).
#'
#' @param Filename  - optional. String(s) to use for file name. Default ("plot_scatter.pdf"). If set to \code{""} a file is not written.
#'
#' @param Width     - width for figure file. Default (12in).
#'
#' @param Height    - height for figure file. Default (8in).
#'
#' @details If \code{colBy} is set to \code{NULL}, no grouping will be done for colouring (a scale warning will be issued that can be safely ignored).
#'
#' @return
#' ggplot2 plot object
#'
#' @export
#' @examples
#' plot_data_example <- data.frame("Groups"=LETTERS[1:12],
#'                                 "Factor"=c(rep(1,4),rep(2,4),rep(3,4)),
#'                                 "Value"=runif(12,-10,10),stringsAsFactors = FALSE)
#' 
#' plot_scatter(plot_data_example,Filename = "")
#' 
#' # Change plot order to be according to the numeric value
#' plot_data_example <- plot_data_example[order(plot_data_example$Value),]
#' plot_scatter(plot_data_example,Filename = "")
#'
#' @seealso \link{analysis_genomes_module_output},\link{analysis_pca_mean_distance_grouping}

############################################################################################################################################

plot_scatter <- function(plot_data,X = 1, Y = 3,colBy = NULL,xLab = "",yLab = "",Filename = "plot_scatter.pdf",xLabs_angle = T,Width = 24, Height = 18){

  # INPUT
  stopifnot(is.data.frame(plot_data),is.character(xLab),is.character(yLab),is.numeric(X),is.numeric(Y))

  if(length(xLab)>1) xLab <- xLab[1]
  if(length(xLab)>1) yLab <- yLab[1]
  if(xLab == "")     xLab <- names(plot_data)[X]
  if(yLab == "")     yLab <- names(plot_data)[Y]

  # CONSERVE DATA ORDER
  # for(C in 1:ncol(plot_data)) plot_data[,C] <- factor(plot_data[,C], levels = unique(plot_data[,C]))
  
  if(is.null(colBy)){ 
    DATA        <- plot_data[,c(X,Y)]
    DATA$colBy <- 1
  }else{
    DATA        <- plot_data[,c(X,Y,colBy)]
  }
  names(DATA) <- c("X","Y","colBy")
  # DATA$colBy  <- as.factor(DATA$colBy)
  
  # CONSERVE ORDER
  DATA$X <- as.character(DATA$X)
  DATA$X <- factor(DATA$X,levels = unique(DATA$X))
  DATA$colBy <- as.character(DATA$colBy)
  DATA$colBy <- factor(DATA$colBy,levels = unique(DATA$colBy))
  
  if(length(unique(DATA$colBy))<=10){
    colours <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
    colours <- colours[1:length(unique(DATA$colBy))]
  }else{
    colours <- colors()[5:(length(unique(DATA$colBy))+5)]
  }
  # PLOT
  if(length(unique(DATA$colBy))>1){
    p <- ggplot2::ggplot(DATA,mapping = ggplot2::aes(x = DATA$X, y = DATA$Y,group = DATA$colBy,col = DATA$colBy)) + 
            ggplot2::geom_line() + ggplot2::scale_color_manual(values = colours)
  }else{
    p <- ggplot2::ggplot(DATA) + 
          ggplot2::geom_point(mapping = ggplot2::aes(x = DATA$X, y = DATA$Y,colour = DATA$colBy)) + ggplot2::scale_color_manual(values = colours)
  }
  p <- p + #ggplot2::scale_color_manual(values = colours) +
              ggplot2::xlab(xLab) + ggplot2::ylab(yLab)  + ggplot2::theme_bw() +
              # ggplot2::theme(axis.ticks.x   = ggplot2::element_blank()) +
              ggplot2::theme(axis.title.x   = ggplot2::element_text(size = 16)) +
              ggplot2::theme(axis.title.y   = ggplot2::element_text(size = 16)) +
              ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),panel.grid.minor.y = ggplot2::element_blank()) +
              ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),panel.grid.major.y = ggplot2::element_blank()) +
              ggplot2::theme(axis.title.y	= ggplot2::element_text(size = Width*1.3,margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)),
                             axis.title.x	= ggplot2::element_text(size = Width*1.3,margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0))) +
              ggplot2::theme(plot.margin= ggplot2::unit(c(0.5,0.5,0.5,0.5),"cm"))   # top, right, bottom, and left

  # IF MORE THAN ONE FACTOR IS PRESENT, ADD LEGEND
  if(length(unique(DATA$colBy))>1){
    p <- p + ggplot2::theme(legend.key = ggplot2::element_blank(),legend.title = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # X axis labels - if less than 50 X values, add labels
  if(length(DATA$X)>20){
    p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_blank())
  }else{
    nChars <- max(nchar(as.character(unique(DATA$X)))[1:ifelse(length(unique(DATA$X))>2,3,length(unique(DATA$X)))])-max(nchar(as.character(DATA$Y)))
    if(!exists("nChar")) nChar <- -1
    if(nChars<0) nChars <- 0
    
    p <- p + ggplot2::theme(plot.margin  = ggplot2::unit(c(0.3,0.5,0.5,0.25+nChars/5),"cm"))   # top, right, bottom, and left
    
    if(xLabs_angle){ 
      p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_text(hjust = 1,vjust=1,colour = "black",size = Width*1.1,angle = 45))
    }else{
      p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_text(hjust = 1,vjust=1,colour = "black",size = Width*1.1))
    }
  }
  
  # Y axis labels
  p <- p + ggplot2::theme(axis.text.y   = ggplot2::element_text(hjust = 1,vjust=1,colour = "black",size = Width*1.2))

  # Save Figure ----
  if(Filename!=""&&is.character(Filename)){
    for(File in Filename){
      ggplot2::ggsave(filename =File,plot = p,width = Width, height = Height)
    }
  }
  return(p)
}

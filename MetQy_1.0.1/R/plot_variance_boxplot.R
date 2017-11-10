#' Make a boxplot
#'
#' @param MATRIX_IN - numeric matrix containing the module fraction match of the genomes (rows) and modules (columns).
#'                        #The columns are used as factor for the box plot grouping.
#'
#' @param x_lab     - optional. String to specify the x-axis label corresponding to the row values. Default ("Modules").
#'
#' @param y_lab     - optional. String to specify the legend title. Default ("Fraction\\n matched").
#'
#' @param Filename  - optional. A character vector containing the file path to save image(s). The device to save is determined by file extension. Default  ("", i.e. file not written).
#'
#' @param Width     - Width size for file. Default (24in).
#'
#' @param Height    - Height size for file. Default (18in).
#'
#'
#' @return  ggplot object of plot
#' @export

############################################################################################################################################

plot_variance_boxplot <- function(MATRIX_IN, x_lab = "Modules", y_lab = "Fraction matched", Filename ="", Width = 24, Height = 18){
  # GIVE MATRIX AND IT CALCULATES DE VARIATION OF THE ROWS

  # MANAGE INPUT ----
  if(!is.matrix(MATRIX_IN)||!is.numeric(MATRIX_IN)) stop("'MATRIX_IN' must be a numeric matrix")

  if(!is.character(x_lab)||!length(x_lab)==1){
    warning("Incorrect 'x_lab' format. Using default values")
    x_lab <- ""
  }
  if(!is.character(y_lab)||!length(y_lab)==1){
    warning("Incorrect 'y_lab' format. Using default values")
    y_lab <- ""
  }

  # TRANSFORM MATRIX TO DATA FRAME
  DATA <- reshape2::melt(MATRIX_IN)

  # GROUP BY COLUMNS OF ORIGINAL DATA
  DATA[,2]       <- as.factor(DATA[,2])
  names(DATA)[3] <- y_lab

  #### PLOT ----

  # BOXPLOT FOR VARIANCE
  p <-
    ggplot2::ggplot(data = DATA,ggplot2::aes(x = DATA[,2],y = DATA[,3])) + ggplot2::geom_boxplot() + ggplot2::theme_bw() + 
    # AXIS LABELS
    ggplot2::xlab(x_lab) + ggplot2::ylab(y_lab) + ggplot2::ylim(0,1) +
    # AXES FORMATTING
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.title.y	= ggplot2::element_text(size = Width*1.2,margin = margin(t = 0, r = 20, b = 0, l = 0)),
                   axis.title.x	= ggplot2::element_text(size = Width*1.2,,margin = margin(t = 20, r = 0, b = 0, l = 0))) +
    ggplot2::theme(axis.text.x   = ggplot2::element_text(hjust = 1,vjust=1,size=Width*0.95, colour = "black",angle = 45)) +
    ggplot2::theme(axis.text.y   = ggplot2::element_text(size=Width*0.95, colour = "black")) 
  # MANAGE X-AXIS LABELS ----
  if(length(unique(DATA[,2]))<50){

    nChars <- max(nchar(as.character(unique(DATA[,2])))[1:ifelse(length(unique(DATA[,2]))>2,3,length(unique(DATA[,2])))])-max(nchar(as.character(DATA[,1])))
    if(exists("nChar")==F) nChar <- -1
    if(nChars<0) nChars <- 2

    p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_text(hjust = 1,vjust=1,size=Width*1, colour = "black",angle = 45)) +
      ggplot2::theme(plot.margin= ggplot2::unit(c(0.3,0.4,0.5,nChars*0.3),"cm"))   # top, right, bottom, and left
  }else{
    p <- p + ggplot2::theme(axis.text.x   = ggplot2::element_blank())
  }

  # MANAGE Y-AXIS LABELS ----
  # if(length(unique(DATA[,1]))<50){
    # p <- p + ggplot2::theme(axis.text.y   = ggplot2::element_text(size=Width*1, colour = "black"))
  # }else{
  #   p <- p + ggplot2::theme(axis.text.y   = ggplot2::element_blank())
  # }
  
  p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank(),plot.margin = ggplot2::unit(rep(0.5,4),units = "cm" )) # top, right, bottom, and left
  

  # Save Figure ----
  if(Filename!=""&&is.character(Filename)){
    for(file in Filename){
      ggplot2::ggsave(filename =file,plot = p,width = Width, height = Height)
    }
  }

  # RETURN PLOT
  return(p)
}

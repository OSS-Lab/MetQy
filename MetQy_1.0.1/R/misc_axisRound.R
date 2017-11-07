#' Find best value to round axis to.
#'
#' @param Vector - values that are being ploted.
#'
#' @param roundBy - optional. The value to use for rounding. Default (\code{NULL}).
#'
#' @param Min - optional. The value to use as the starting value for the axis. Default (\code{NULL}).
#'
#' @return List containing the minimum and maximum axis values (\code{$min} and \code{$max}, respectivey) and the array that can be used to indicate the axis breaks (\code{$array})
#'
#' @export

misc_axisRound <- function(Vector, roundBy = NULL, Min = NULL){

  # roundBy = NULL; Min = NULL

  # MANAGE INPUT
  stopifnot(is.numeric(Vector),is.null(roundBy)||is.numeric(roundBy),is.null(Min)||is.numeric(Min))

  if(is.null(roundBy)){
    roundBy_values   <- c(1:10*10^-5,2:10*10^-4,2:10*10^-3,2:10*10^-2,2:10*10^-1, 2:10*10^0,2:10*10^1,2:10*10^2,2:10*10^3,2:10*10^4,2:10*10^5,0.0025,0.025,0.25)
                          # c(1e-05,2e-05,5e-05,1e-04,2e-04,5e-04,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1:10,1:10*10,1:10*100,1:10*1000,1:10*10000,1:10*100000)
                          # ,50,100,200,300,500,1000,2000,3000,4000,5000,10000,20000,50000,100000)

    if(max(Vector)<1){
      # roundBy <- sort(roundBy_values[which(round(max(Vector)/5,digits = 5)%/%roundBy_values==1)],decreasing = T)[1]
      roundBy <- roundBy_values[which(round(max(Vector)/5,digits = 5)%/%roundBy_values==1)]
      # if(sum(roundBy==0.0025)>0) roundBy <- 0.0025
      if(sum(roundBy==0.025)>0) roundBy <- 0.025
      if(sum(roundBy==0.05)>0)  roundBy <- 0.05
      if(sum(roundBy==0.05)>0)  roundBy <- 0.0005
      if(sum(roundBy==0.25)>0)  roundBy <- 0.25
      if(sum(roundBy==0.5)>0)   roundBy <- 0.5
      if(length(roundBy)>1)     roundBy <- sort(roundBy,decreasing = T)[1]
    } else {
      roundBy <- roundBy_values[which(ceiling(max(Vector)/5)%/%roundBy_values==1)]

      # IF A GIVEN roundBy NUMBER ENDS WITH THE MAX, USE THAT
      max_index <- which((max(Vector)/roundBy)/round(max(Vector)/roundBy)==1)
      if(length(max_index)>0){
        roundBy <- sort(roundBy[max_index],decreasing = T)[1]
        Max     <- max(Vector)
      }else{
        if(sum(roundBy==0.05)>0)roundBy <- 0.05
        if(sum(roundBy==0.5)>0) roundBy <- 0.5
        if(sum(roundBy==5)>0)   roundBy <- 5
        if(sum(roundBy==50)>0)  roundBy <- 50
        if(sum(roundBy==100)>0) roundBy <- 100
        if(sum(roundBy==200)>0) roundBy <- 200
        if(sum(roundBy==500)>0) roundBy <- 500
        if(length(roundBy)>1)   roundBy <- sort(roundBy,decreasing = T)[1]
      }
    }

    if(roundBy==0){
      roundBy <- sort(roundBy_values[which(ceiling(max(Vector)/5)%/%roundBy_values==2)],decreasing = T)[1]
    }
    if(roundBy==0){
      stop("No appropriate rounding value found")
    }
  }

  if(is.null(Min)) Min <- (min(Vector)%/%roundBy)*roundBy

  if(!exists("Max")) Max   <- (max(Vector)%/%roundBy+1)*roundBy
  Array <- seq(Min,Max,by = roundBy)

  # OUTPUT
  output <- list(min = Min,
                 max = Max,
                 array = Array)
  return(output)
  # print(output)
}

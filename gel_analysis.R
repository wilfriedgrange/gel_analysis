rm(list=ls())# clean memory
library(baseline) 
image_file <- '/Users/wgrange/desktop/im.txt' # 8 bits file
################################
################################
#### functions #################
################################
################################

################################
# rotate matrix 
################################
rotate_matrix<- function(x) {t(apply(x, 2, rev))}

################################
# Local minima (maxima) function
# https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r
################################
find_peaks <- function (x, m = 5){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

################################
# Peak Fitting Function
################################
peak_fit<-function(Peak,getc,counter,num_peak,all) {
  half_size<-18 # half fitting range
  f1 <- function(par) # Gaussian fit
  {
    m <- par[1]
    sd <- par[2]
    k <- par[3]
    rhat <- abs(k) * exp(-0.5 * ((x - m)/sd)^2) 
    sum((r - rhat)^2)
  }
  # This is the section where we fit
  r<-getc[(Peak-half_size):(Peak+half_size)]
  x <- seq_along(r)
  max<-which(r==max(r))  
  num<-optim(c(max, 1, 1), f1, method="BFGS", control=list(reltol=1e-9))   # return loss (error)
  # plot
  x_int<-seq(1,length(r),0.01)
  predicted_values<-num$par
  g1<- abs(predicted_values[3]) *exp(-(x_int- predicted_values[1])**2/(2 *  predicted_values[2]**2))
  lines(x_int+Peak-half_size-1,g1,col = 'blue' ,lwd =2)
  A0<-round(abs(predicted_values[2]*predicted_values[3]*sqrt(2*pi)),0)
  text(predicted_values[1]+Peak-0.3*half_size, 0.9*abs(predicted_values[3]), A0, col = "blue",cex =1.5)
}


################################
################################
########## Main ################
################################
################################

# Read .txt file (8 bits image)
image0<-read.csv(image_file,header=F,sep='\t')
# convert to a matrix  
image0<-data.matrix(image0)
# Average slices
avslicex<-apply(image0,2,mean)
t<-find_peaks(-avslicex,m=18)


for (z in 1:(t-1))
{
  print("Lane : ")
  print(z)
  avslicex<-apply(image0[,t[z]:t[z+1]],1,mean)
  # Baseline correction
  bc.rollingBall <- baseline(t(matrix(avslicex)),wm=10, ws=10,method='rollingBall')
  getb<-as.vector(getBaseline(bc.rollingBall))
  getc<-as.vector(getCorrected(bc.rollingBall))->getc_duplicated
  # Plots
  layout(matrix(c(1,2), 1, 2, byrow = TRUE),
         widths =c(1,4))
  lane<-image0[,t[z]:t[z+1]]
  image(rotate_matrix(lane), axes = FALSE, col = grey(seq(0, 1, length = 255)),asp=nrow(lane)/ncol(lane))
  plot(getc,xlab='', ylab='Intensity (baseline corrected)')
  lines(getc)
  x_ticks <- axis(1, seq(0,length(getc),10), labels = FALSE)
  grid(lwd = 1, ny = NULL, nx = NA)
  abline(v = x_ticks, lwd = 2, lty = 3, col = "lightgray")
  getc_duplicated[getc_duplicated<5]<-0 # discard values smaller than 5
  t2<-find_peaks(getc_duplicated,18) # get peaks location
  t2<-t2[which(t2>18)] # and just keep those which are high enough
  # and fit
  values<-sapply(t2, function(w){
    peak_fit(w,getc,counter,length(t2),1)
  }
  )
}
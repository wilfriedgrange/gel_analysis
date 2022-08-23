rrm(list=ls())# clean memory
if (!require(baseline)) install.packages('baseline') # baseline package
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set wd
image_file <- 'IM1.txt' # 8 bits file

################################
################################
#### parameters ################
################################
################################
################################
ws<-10 # smooth baseline (change !)
wm_lowMW<-10# wm for low Molecular Weight (change !)
wm_highMW<-15 #  wm for high Molecular Weight (change !)
half_size<-18 # half fitting range (change !)
m_cut<- 18 # used to separate lanes (change !)

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
peak_fit<-function(Peak,getc,num_peak,half_size = 10) {
  f1 <- function(par) # Gaussian fit
  {
    m <- par[1]
    sd <- par[2]
    k <- par[3]
    rhat <- abs(k) * exp(-0.5 * ((x - m)/sd)^2) 
    sum((r - rhat)^2) # as in Least Squares
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
# Average slices (this determines lanes)
avslicex<-apply(image0,2,mean)
t<-find_peaks(-avslicex,m_cut)

for (z in 1:(length(t)-1))
{
  print('Analyzing new lane ....')
  avslicex<-apply(image0[,t[z]:t[z+1]],1,mean)
  # Baseline correction
  cut<-length(avslicex)/2
  bc.rollingBall_highMW <- baseline(t(matrix(avslicex[1:cut])),wm_highMW, ws,method='rollingBall')
  getb_highMW<-as.vector(getBaseline( bc.rollingBall_highMW ))
  getc_highMW<-as.vector(getCorrected( bc.rollingBall_highMW ))->getc_duplicated_highMW
  bc.rollingBall_lowMW <- baseline(t(matrix(avslicex[cut:length(avslicex)])),wm_lowMW, ws,method='rollingBall')
  getb_lowMW <-as.vector(getBaseline( bc.rollingBall_lowMW  ))
  getc_lowMW <-as.vector(getCorrected( bc.rollingBall_lowMW  ))->getc_duplicated_lowMW 
  getc<-c(getc_highMW,getc_lowMW)->getc_duplicated

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
    peak_fit(w,getc,length(t2),half_size)
  }
  )
}
print('done ....')

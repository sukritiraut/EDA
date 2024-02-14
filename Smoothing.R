## Part #1 deterministic equidistant design
## Generate n=101 equidistant points in [-2\pi, 2\pi]

##EQUIDISTANT DATA FOR PART 1
x <- 2*pi*seq(-1, 1, length=n)

##NON-EQUIDISTANT DATA FOR PART 2
x2 <- round(2*pi*sort(c(0.5, -1 + rbeta(50,2,2), rbeta(50,2,2))), 8)

##EDA OF MEXICAN HAT FUNCTION
#y <- (1-x^2)*exp(-0.5*x^2) 
y <- ((1-x^2)*exp(-0.5*x^2) + rnorm(length(x), sd=0.2))
#y2 <- (1-x2^2)*exp(-0.5*x2^2)
y2 <- ((1-x2^2)*exp(-0.5*x2^2) + rnorm(length(x), sd=0.2))

dmin= min(y,y2)
dmax = max(y,y2)
  
matplot(x, y, "l", ylim=c(dmin, dmax), ylab="Response", col="black",lwd=2, lty=1, main="Mexican hat function f(x) for equidistant and non-equidistant data points (with error)")
matlines(x, y2, col="red",lwd=2, lty=1)
legend("topleft", legend=c("f(x) for equidistant points", "f(x) for non-equidistant points"),
       col=c("black", "red"), lty=1:2, cex=0.8)

########################################################################
#PART 1

m <- 1000
n <- 101

## Initialize the matrix of fitted values for three methods
fvlp <- fvnw <- fvss <- matrix(0, nrow= n, ncol= m)


##Generate data, fit the data and store the fitted values

for (j in 1:m){
  ## simulate y-values

  y <- ((1-x^2)*exp(-0.5*x^2) + rnorm(length(x), sd=0.2));   #Mexican hat function
  
  ## Get the estimates and store them
  fvlp[,j] <- predict(loess(y ~ x, span = 0.75), newdata = x); #LOESS with span = 0.75
  fvnw[,j] <- ksmooth(x, y, kernel="normal", bandwidth= 0.2, x.points=x)$y; #nadaraya-watson (NW) kernel smoothing with gaussian kernel and bandwith=0.2
  fvss[,j] <- predict(smooth.spline(y ~ x), x=x)$y #spline smoothing with default tuning parameter
}

##true y value from mexican hat function 
y <- ((1-x^2)*exp(-0.5*x^2))


par(mfrow=c(2,2))

## Below is the sample R code to plot the mean of three estimators in a single plot
meanlp = apply(fvlp,1,mean); 
meannw = apply(fvnw,1,mean);
meanss = apply(fvss,1,mean);

dmin = min(meanlp, meannw, meanss);
dmax = max( meanlp, meannw, meanss);

matplot(x, meanlp, "l", ylim=c(dmin, dmax), ylab="Response", col="black",lwd=2, lty=1, main="Average estimated values vs true values")
matlines(x, meannw, col="red",lwd=2, lty=1)
matlines(x, meanss, col="blue",lwd=2, lty=1)
matlines (x, y, col="green",lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS","true"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)



## plot the empirical bias/variance/MSE



##bias

Biaslp <-  meanlp - y 
Biasnw <- meannw - y
Biasss <- meanss - y

Biasmin = min( Biaslp, Biasnw, Biasss);
Biasmax = max( Biaslp, Biasnw, Biasss);

matplot(x, Biaslp, "l", ylim=c(Biasmin, Biasmax), ylab="Response",col="black",lwd=2, lty=1,main="Bias")
matlines(x, Biasnw, col="red",lwd=2, lty=1)
matlines(x, Biasss, col="blue",lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)


##variance

varlp <- 1/1000*(apply((fvlp- meanlp)^2,1,sum))
varnw <- 1/1000*(apply((fvnw- meannw)^2,1,sum))
varss <- 1/1000*(apply((fvss- meanss)^2,1,sum))


varmin = min( varlp, varnw, varss);
varmax = max( varlp, varnw, varss);

matplot(x, varlp, "l", ylim=c(varmin, varmax), ylab="Response",col="black",lwd=2, lty=1,main="Variance")
matlines(x, varnw, col="red",lwd=2, lty=1)
matlines(x, varss, col="blue",lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)



##MSE

MSElp <- 1/1000*(apply((fvlp - y)^2,1,sum))
MSEnw <- 1/1000*(apply((fvnw - y)^2,1,sum))
MSEss <- 1/1000*(apply((fvss - y)^2,1,sum))

MSEmin = min(MSElp, MSEnw, MSEss);
MSEmax = max(MSElp, MSEnw, MSEss);

matplot(x, MSElp, "l", ylim=c(MSEmin, MSEmax), ylab="Response",col="black",lwd=2, lty=1,main="MSE")
matlines(x, MSEnw, col="red",lwd=2, lty=1)
matlines(x, MSEss, col="blue",lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)




###############################################################################


## Part #2 non-equidistant design


set.seed(79)
x2 <- round(2*pi*sort(c(0.5, -1 + rbeta(50,2,2), rbeta(50,2,2))), 8)


## assume you save the file "HW04part2.x.csv" in the local folder "C:/temp",
#x2 <- read.table(file= "C:/temp/HW04part2.x.csv", header=TRUE);




##Generate data, fit the data and store the fitted values
m <- 1000
n <- 101
#x2 <- read.table(file= "HW04part2-1.x.csv", sep = ",", header=TRUE);

## Initialize the matrix2 of fitted values for three methods
fvlp <- fvnw <- fvss <- matrix(0, nrow= n, ncol= m)


##Generate data, fit the data and store the fitted values

for (j in 1:m){
  ## simulate y-values
  
  y <- ((1-x2^2)*exp(-0.5*x2^2) + rnorm(length(x2), sd=0.2));   #Mex2ican hat function
  
  ## Get the estimates and store them
  fvlp[,j] <- predict(loess(y ~ x2, span = 0.3365), newdata = x2); #LOESS with span = 0.75
  fvnw[,j] <- ksmooth(x2, y, kernel="normal", bandwidth= 0.2, x.points=x2)$y; #nadaraya-watson (NW) kernel smoothing with gaussian kernel and bandwith=0.2
  fvss[,j] <- predict(smooth.spline(y ~ x2, spar= 0.7163), x=x2)$y #spline smoothing with default tuning parameter
}

##true y value from mex2ican hat function 
y <- ((1-x2^2)*exp(-0.5*x2^2))


par(mfrow=c(2,2))

## Below is the sample R code to plot the mean of three estimators in a single plot
meanlp = apply(fvlp,1,mean); 
meannw = apply(fvnw,1,mean);
meanss = apply(fvss,1,mean);

dmin = min(meanlp, meannw, meanss);
dmax = max( meanlp, meannw, meanss);

matplot(x2, meanlp, "l", ylim=c(dmin, dmax), ylab="Response",col="black", lwd=2, lty=1, main="Average estimated values vs true values")
matlines(x2, meannw, col="red", lwd=2, lty=1)
matlines(x2, meanss, col="blue", lwd=2, lty=1)
matlines (x2,y, col="green", lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS","true"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)




## plot the empirical bias/variance/MSE

##bias

Biaslp <-  meanlp - y 
Biasnw <- meannw - y
Biasss <- meanss - y

Biasmin = min( Biaslp, Biasnw, Biasss);
Biasmax = max( Biaslp, Biasnw, Biasss);

matplot(x2, Biaslp, "l", ylim=c(Biasmin, Biasmax), ylab="Response", col="black", lwd=2, lty=1,main="Bias")
matlines(x2, Biasnw, col="red", lwd=2, lty=1)
matlines(x2, Biasss, col="blue", lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)



##variance

varlp <- 1/1000*(apply((fvlp- meanlp)^2,1,sum))
varnw <- 1/1000*(apply((fvnw- meannw)^2,1,sum))
varss <- 1/1000*(apply((fvss- meanss)^2,1,sum))


varmin = min( varlp, varnw, varss);
varmax = max( varlp, varnw, varss);

matplot(x2, varlp, "l", ylim=c(varmin, varmax), ylab="Response", col="black", lwd=2, lty=1,main="Variance")
matlines(x2, varnw, col="red", lwd=2, lty=1)
matlines(x2, varss, col="blue", lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)




##MSE

MSElp <- 1/1000*(apply((fvlp - y)^2,1,sum))
MSEnw <- 1/1000*(apply((fvnw - y)^2,1,sum))
MSEss <- 1/1000*(apply((fvss - y)^2,1,sum))

MSEmin = min(MSElp, MSEnw, MSEss);
MSEmax = max(MSElp, MSEnw, MSEss);

matplot(x2, MSElp, "l", ylim=c(MSEmin, MSEmax), ylab="Response",col="black", lwd=2, lty=1,main="MSE")
matlines(x2, MSEnw, col="red", lwd=2, lty=1)
matlines(x2, MSEss, col="blue", lwd=2, lty=1)
legend("topleft", legend=c("LOWESS", "NW","SS"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)







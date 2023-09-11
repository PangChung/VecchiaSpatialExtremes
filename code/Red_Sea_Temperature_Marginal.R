#################
### LIBRARIES ###
#################

library(splines)
library(ismev)
library(fields)
library(ggplot2)
library(ggmap)
library(cowplot)

### Libraries that needed for the functions ###
library(evd)
library(parallel)
library(methods)
library(mvtnorm)
library(Rfast) 
library(partitions)
library(rgdal)

###########################################
### SET WORKING DIRECTORY AND LOAD DATA ###
###########################################

rm(list=ls()) # clear the environment
wd <- "."
setwd(wd)

load("data/data_original.RData")

x <- sort(unique(loc.sub[,1]))
y <- sort(unique(loc.sub[,2]))

nx <- length(x)
ny <- length(y)

###########################
### EXTRACT SUB-DATASET ###
###########################

n.times <- nrow(data.sub)
n.sites <- ncol(data.sub)
n.years <- length(unique(year))

#######################
### MARGINAL TRENDS ###
#######################

#BS <- bs(c(1:365),degree=3,df=14,intercept=FALSE)
BS <- cSplineDes(c(1:365),knots=quantile(c(1:365),c(0:12)/12),ord=4) ### This defines cyclic cubic B-splines defined on 12 knots (one for each month)

## This is to visualize the 12 cyclic splines
par(mfrow=c(3,4))
for(i in 1:ncol(BS)){
  plot(c(1:365),BS[,i],type="l")
}

## to see the chosen knots: attr(BS,"knots")

## We create the design matrix X
X <- cbind(1,c(1:n.times)/(365*100)) ## here, we add the intercept and "time" as a covariate
for(i in 1:ncol(BS)){ ## here we add all the 12 splines as covariates
  X <- cbind(X,rep(BS[,i],n.years))
}

n.par <- ncol(X) ## this gives the total number of parameters

est.par.mean <- matrix(nrow=n.par,ncol=n.sites) ## matrix with estimated parameters to describe the mean temperature
est.mean <- matrix(nrow=n.times,ncol=n.sites) ## fitted mean temperatures
residuals1 <- matrix(nrow=n.times,ncol=n.sites) ## residuals = temperature - fitted mean
est.par.sd <- matrix(nrow=2,ncol=n.sites) ## matrix with estimated parameters to describe the temperature standard deviation (on log-scale). We only estimate 2 parameters here (intercept+time trend) but do not consider seasonality.
est.sd <- matrix(nrow=n.times,ncol=n.sites) ## fitted temperature standard deviation
residuals2 <- matrix(nrow=n.times,ncol=n.sites) ## residuals = (temperature - fitted mean)/(fitted standard deviation)

distmat <- rdist.earth(loc.sub,miles=FALSE)
radius <- 50 ##in km (defines the local neighborhoods for estimating trends and the GEV marginal distribution for maxima)

par(mfrow=c(3,2))

### Estimation of marginal parameters: (1) Intercept + linear time trend + seasonality in mean temperature and parameters; and (2) Intercept + linear time trend in temperature standard deviation.
### CAREFUL! The following for-loop is quite INTENSIVE! (takes 5-10min to complete)
### Alternatively: can load the results by running the line below starting with "load(...)" 
for(j in 1:n.sites){
  print(j)
  ind.j <- which(distmat[,j]<=radius) ## select all sites that are within 50km of the j-th location 
  data.j <- as.vector(data.sub[,ind.j]) ## extract all data at all these sites and stack them in a big vector (here, the data are assumed to be independent across sites)
  X.j <- X ## we define the design matrix for the j-th location to be fitted
  for(i in 1:(length(ind.j)-1)){
    X.j <- rbind(X.j,X)
  }
  X.j.mean <- as.matrix(X.j) ## Design matrix to fit the mean tempature, which includes the intercept, linear time trend, and 12 cyclic splines for seasonality
  X.j.sd <- as.matrix(X.j[,1:2]) ## Design matrix to fit the tempature standard deviation, which only includes the intercept and linear time trend
  
  fit.mean.j <- lm(data.j~X.j.mean-1) ## linear model fit (by least squares) to estimate mean temperature
  est.par.mean[,j] <- fit.mean.j$coefficients ## estimated parameters
  est.mean.j <- fit.mean.j$fitted.values ## fitted mean temperatures
  residuals1.j <- data.j-est.mean.j ## residuals = temperature - fitted mean
  est.mean[,j] <- est.mean.j[nrow(X)*(which(ind.j==j)-1) + c(1:nrow(X))] ## we here extract the fitted mean that corresponds to the target location and add it to the result matrix
  residuals1[,j] <- data.sub[,j]-est.mean[,j] ## same for the residuals
  
  ## we here define the likelihood function for the j-th site based on a N(mu_j,sigma_j^2) assumption with mu_j fixed to what was estimated above. We here simply estimate sigma_j^2.
  nllk.j <- function(param,residuals1.j,X.j.sd){
    sds <- exp(X.j.sd%*%param) ## standard deviation has a log-linear specification 
    return(sum(log(sds)+0.5*(residuals1.j/sds)^2)) ## this returns -log(likelihood), without some additive constants
  }
  init <- c(sd(residuals1[,j]),0) ## these are the initial values
  ## we here minimize the function nllk.j (i.e., maximize the likelihood for the j-th site) using the Nelder-Mead method
  fit.sd.j <- optim(par=init,fn=nllk.j,residuals1.j=residuals1.j,X.j.sd=X.j.sd,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)   
  
  est.par.sd[,j] <- fit.sd.j$par ## estimated parameters for the standard deviation
  est.sd.j <- exp(X.j.sd%*%est.par.sd[,j]) ## fitted standard deviations
  est.sd[,j] <- exp(X[,1:2]%*%est.par.sd[,j]) ## we here extract the fitted standard deviation that corresponds to the target location and add it to the result matrix
  residuals2.j <- residuals1.j/est.sd.j ## we compute the residuals (temperature - mean temperature)/(temperature standard deviation)
  residuals2[,j] <- residuals1[,j]/est.sd[,j] ## same as above
  
  #plot(time,data.sub[,j],type="l",xlab="Time",ylab="Temperature [ºC]",main=paste("Mean at the ",j,"-th location",sep=""))
  #lines(time,est.mean[,j],col="red")
  #plot(time,est.sd[,j],type="l",,xlab="Time",ylab="SD of Temperature [ºC]",main=paste("SD at the ",j,"-th location",sep=""))
}

#save(list=c("BS","X","distmat","radius","est.par.mean","est.mean","residuals1","est.par.sd","est.sd","residuals2"),file="data/Trends_fits.RData")
#load(file="data/Trends_fits.RData")

#####################
### YEARLY MAXIMA ###
#####################

## We compute yearly maxima for each site, to which the GEV distribution will be fitted
maxima <- matrix(nrow=n.years,ncol=n.sites) 
for(k in 1:n.years){
  maxima[k,] <- apply(residuals2[c(((k-1)*365+1):(k*365)),],2,max)
}

################
### GEV FITS ###
################

est.par.GEV <- matrix(nrow=3,ncol=n.sites) ## matrix of estimated GEV parameters (fitted to yearly maxima)
sd.par.GEV <- matrix(nrow=3,ncol=n.sites) ## matrix of standard deviations estimated for the GEV parameters (fitted to yearly maxima)
maxima.frechet <- matrix(nrow=n.years,ncol=n.sites) ## matrix of yearly maxima transformed to the unit Fréchet distribution

## Estimation of GEV parameters
## This is quite fast!
for(j in 1:n.sites){
  print(j)
  ind.j <- which(distmat[,j]<=radius) ## Again, we pool stations together within a 50km distance of the j-th site to increase the sample size at each location.
  maxima.j <- as.vector(maxima[,ind.j]) ## we put all the data at all these sites together
  fit.j <- gev.fit(maxima.j,show=FALSE) ## we fit a stationary GEV(mu,sigma,xi) distribution to the pooled sample at the j-th site
  est.par.GEV[,j] <- fit.j$mle ## these are the GEV estimated parameters
  sd.par.GEV[,j] <- fit.j$se ## GEV estimated standard deviations
  maxima.frechet[,j] <- -1/log(pgev(maxima[,j],loc=fit.j$mle[1],scale=fit.j$mle[2],shape=fit.j$mle[3])) ## maxima transformed to the unit Fréchet scale
}

#save(list=c("maxima","maxima.frechet","est.par.GEV","sd.par.GEV"),file="data/GEV_fits.RData")
#load(file="data/GEV_fits.RData")

#save(list=ls(),file="data/Marginal_fits_ALL.RData")
#load(file="data/Marginal_fits_ALL.RData")


################
### QQ-PLOTS ###
################

## We check that the marginal fits are good at a set of selected locations

## QQ-plots for data on unit Fréchet scale at 9 FIXED sites
pdf("figure/qqplot_9_sites.pdf",width = 4*3,height=4*3)
par(mfrow=c(3,3))
for(j in c(1,100,250,400,500,600,750,800,1000)){
  qqplot(maxima.frechet[,j],qfrechet(c(1:n.years)/(n.years+1)),log="xy",xlab="Empirical quantiles",ylab="Theoretical Fréchet(1) quantiles",main=paste(j,"-th location",sep=""),pch=20)
  abline(0,1,col="red")
}
dev.off()

## QQ-plots for data on unit Fréchet scale at 16 RANDOM sites
pdf("figure/qqplot_16_sites.pdf",width = 4*4,height=4*4)
par(mfrow=c(4,4))
for(j in sort(sample(c(1:ncol(maxima.frechet)),16,replace=FALSE))){
  qqplot(maxima.frechet[,j],qfrechet(c(1:n.years)/(n.years+1)),log="xy",xlab="Empirical quantiles",ylab="Theoretical Fréchet(1) quantiles",main=paste(j,"-th location",sep=""),pch=20)
  abline(0,1,col="red")
}
dev.off()
######################
### MAPS AND PLOTS ###
######################

### TEMPORAL plot of the temperature time series at a given location, with estimated mean in red, estimated standard deviation, and residuals (data-fitted.mean)/fitted.SD
pdf("figure/marginal_diagnostic_3_sites.pdf",width = 4*3,height = 4*3)
my.loc <- c(100,400,700) ## North to South
par(mfrow=c(3,3),cex.main=2,cex.lab=2,mar=c(4,5,3,1))
for(j in 1:3){
  plot(time,data.sub[,my.loc[j]],type="l",xlab="Time",ylab="Temperature [ºC]",main=paste("Mean at the ",my.loc[j],"-th location",sep=""))
  lines(time,est.mean[,my.loc[j]],col="red")
}

for(j in 1:3){
  plot(time,est.sd[,my.loc[j]],type="l",,xlab="Time",ylab="SD of Temperature [ºC]",main=paste("SD at the ",my.loc[j],"-th location",sep=""))
}

for(j in 1:3){
  plot(time,residuals2[,my.loc[j]],type="l",,xlab="Time",ylab="Residuals",main=paste("Residuals at the ",my.loc[j],"-th location",sep=""))
  abline(h=0,col="red")
}
dev.off()

### SPATIAL map of temperature field corresponding to the "most extreme" day overall (based on sum of raw data)
my.day <- which.max(rowSums(data.sub))

data.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = data.sub[my.day,]) ## Create data frame for plotting SST data
data.mean.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.mean[my.day,]) ## Create data frame for plotting mean SST
data.sd.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.sd[my.day,]) ## Create data frame for plotting SD SST
data.residuals.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = residuals2[my.day,]) ## Create data frame for plotting residual SST data
data.mean.trend.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.par.mean[2,]) ## Create data frame for plotting mean trend SST
data.sd.trend.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.par.sd[2,]) ## Create data frame for plotting SD trend SST

my.point <- data.frame(lon = loc.sub[my.loc,1], lat = loc.sub[my.loc,2], temp = data.sub[my.day,my.loc]) ## If wants to add selected points on map

#key = 'AIzaSyABn_DX4l8KocC6QCV_lzP-1M4E3qR0Ki0'
#register_google(key = key)
xlim <- c(32,44)
ylim <- c(12,31)
## the following line retrieves the map from Google Map (CAREFUL: don't run this line too often...!)
#map <- get_map(location=c(xlim[1],ylim[1],xlim[2],ylim[2]),zoom=5,source="google",color="color",maptype="satellite",force=TRUE)
#save(map,file="data/GoogleMap.RData")
load(file="data/GoogleMap.RData")
ewbreaks <- c(32, 36, 40, 44)
nsbreaks <- c(12, 18, 24, 30)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))

### 1st: plot the raw data
gg1 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.plot, alpha = 1, na.rm = T, shape = 15, size = 1)
gg1 <- gg1 + scale_color_gradientn(colours = tim.colors())
gg1 <- gg1 + ggtitle("SST raw data (Most extreme day)")
gg1 <- gg1 + xlab("Longitude")
gg1 <- gg1 + ylab("Latitude")
gg1 <- gg1 + labs(colour = expression(paste("[ºC]")))
gg1 <- gg1 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg1 <- gg1 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg1 <- gg1 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg1 <- gg1 + theme(axis.ticks = element_line(size = 1))
gg1 <- gg1 + theme(axis.ticks.length = unit(.1, "cm"))
gg1 <- gg1 + theme(legend.text = element_text(size = 12))
gg1 <- gg1 + theme(legend.key.width = unit(0.5, "cm"))
gg1 <- gg1 + theme(legend.title = element_text(size = 12))
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg1

### 2nd: plot the estimated temperature mean
gg2 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.mean.plot, alpha = 1, na.rm = T, shape = 15, size = 1)
gg2 <- gg2 + scale_color_gradientn(colours = tim.colors())
gg2 <- gg2 + ggtitle("SST estimated mean (Most extreme day)")
gg2 <- gg2 + xlab("Longitude")
gg2 <- gg2 + ylab("Latitude")
gg2 <- gg2 + labs(colour = expression(paste("[ºC]")))
gg2 <- gg2 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg2 <- gg2 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg2 <- gg2 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg2 <- gg2 + theme(axis.ticks = element_line(size = 1))
gg2 <- gg2 + theme(axis.ticks.length = unit(.1, "cm"))
gg2 <- gg2 + theme(legend.text = element_text(size = 12))
gg2 <- gg2 + theme(legend.key.width = unit(0.5, "cm"))
gg2 <- gg2 + theme(legend.title = element_text(size = 12))
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg2

### 3rd: plot the estimated temperature SD
gg3 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.sd.plot, alpha = 1, na.rm = T, shape = 15, size = 1)
gg3 <- gg3 + scale_color_gradientn(colours = tim.colors())
gg3 <- gg3 + ggtitle("SST estimated SD (Most extreme day)")
gg3 <- gg3 + xlab("Longitude")
gg3 <- gg3 + ylab("Latitude")
gg3 <- gg3 + labs(colour = expression(paste("[ºC]")))
gg3 <- gg3 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg3 <- gg3 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg3 <- gg3 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg3 <- gg3 + theme(axis.ticks = element_line(size = 1))
gg3 <- gg3 + theme(axis.ticks.length = unit(.1, "cm"))
gg3 <- gg3 + theme(legend.text = element_text(size = 12))
gg3 <- gg3 + theme(legend.key.width = unit(0.5, "cm"))
gg3 <- gg3 + theme(legend.title = element_text(size = 12))
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg3

### 4th: plot the temperature residuals
gg4 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg4 <- gg4 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.residuals.plot, alpha = 1, na.rm = T, shape = 15, size = 1)
gg4 <- gg4 + scale_color_gradientn(colours = tim.colors())
gg4 <- gg4 + ggtitle("SST residuals (Most extreme day)")
gg4 <- gg4 + xlab("Longitude")
gg4 <- gg4 + ylab("Latitude")
gg4 <- gg4 + labs(colour = expression(paste("")))
gg4 <- gg4 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg4 <- gg4 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg4 <- gg4 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg4 <- gg4 + theme(axis.ticks = element_line(size = 1))
gg4 <- gg4 + theme(axis.ticks.length = unit(.1, "cm"))
gg4 <- gg4 + theme(legend.text = element_text(size = 12))
gg4 <- gg4 + theme(legend.key.width = unit(0.5, "cm"))
gg4 <- gg4 + theme(legend.title = element_text(size = 12))
gg4 <- gg4 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg4


### 5th: plot the temperature mean trend
gg5 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg5 <- gg5 + geom_point(aes(x = lon, y = lat, colour = est.par.mean[1,]), data = data.mean.trend.plot, alpha = 1, na.rm = T, shape = 15, size = 1)
gg5 <- gg5 + scale_color_gradientn(colours = tim.colors())
gg5 <- gg5 + ggtitle(expression(hat(beta)['0,i']~'in'~mean))
gg5 <- gg5 + xlab("Longitude")
gg5 <- gg5 + ylab("Latitude")
gg5 <- gg5 + labs(colour = expression(paste("[ºC/cent.]")))
gg5 <- gg5 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg5 <- gg5 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg5 <- gg5 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg5 <- gg5 + theme(axis.ticks = element_line(size = 1))
gg5 <- gg5 + theme(axis.ticks.length = unit(.1, "cm"))
gg5 <- gg5 + theme(legend.text = element_text(size = 12))
gg5 <- gg5 + theme(legend.key.width = unit(0.5, "cm"))
gg5 <- gg5 + theme(legend.title = element_text(size = 12))
gg51 <- gg5 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")

gg5 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg5 <- gg5 + geom_point(aes(x = lon, y = lat, colour = est.par.mean[2,]), data = data.mean.trend.plot, alpha = 1, na.rm = T, shape = 15, size = 1)
gg5 <- gg5 + scale_color_gradientn(colours = tim.colors())
gg5 <- gg5 + ggtitle(expression(hat(beta)['1,i']~'in'~mean))
gg5 <- gg5 + xlab("Longitude")
gg5 <- gg5 + ylab("Latitude")
gg5 <- gg5 + labs(colour = expression(paste("[ºC/cent.]")))
gg5 <- gg5 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg5 <- gg5 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg5 <- gg5 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg5 <- gg5 + theme(axis.ticks = element_line(size = 1))
gg5 <- gg5 + theme(axis.ticks.length = unit(.1, "cm"))
gg5 <- gg5 + theme(legend.text = element_text(size = 12))
gg5 <- gg5 + theme(legend.key.width = unit(0.5, "cm"))
gg5 <- gg5 + theme(legend.title = element_text(size = 12))
gg52 <- gg5 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")

### 6th: plot the temperature SD trend
gg61 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg61 <- gg61 + geom_point(aes(x = lon, y = lat, colour = est.par.sd[1,]), data = data.sd.trend.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg61 <- gg61 + scale_color_gradientn(colours = tim.colors())
gg61 <- gg61 + ggtitle(expression(hat(gamma)['0,i']~'in'~SD))
gg61 <- gg61 + xlab("Longitude")
gg61 <- gg61 + ylab("Latitude")
gg61 <- gg61 + labs(colour = expression(paste("[ºC/cent.]")))
gg61 <- gg61 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg61 <- gg61 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg61 <- gg61 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg61 <- gg61 + theme(axis.ticks = element_line(size = 1))
gg61 <- gg61 + theme(axis.ticks.length = unit(.1, "cm"))
gg61 <- gg61 + theme(legend.text = element_text(size = 12))
gg61 <- gg61 + theme(legend.key.width = unit(0.5, "cm"))
gg61 <- gg61 + theme(legend.title = element_text(size = 12))
gg61 <- gg61 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg61 

gg6 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg6 <- gg6 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.sd.trend.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg6 <- gg6 + scale_color_gradientn(colours = tim.colors())
gg6 <- gg6 + ggtitle(expression(hat(gamma)['1,i']~'in'~SD))
gg6 <- gg6 + xlab("Longitude")
gg6 <- gg6 + ylab("Latitude")
gg6 <- gg6 + labs(colour = expression(paste("[ºC/cent.]")))
gg6 <- gg6 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg6 <- gg6 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg6 <- gg6 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg6 <- gg6 + theme(axis.ticks = element_line(size = 1))
gg6 <- gg6 + theme(axis.ticks.length = unit(.1, "cm"))
gg6 <- gg6 + theme(legend.text = element_text(size = 12))
gg6 <- gg6 + theme(legend.key.width = unit(0.5, "cm"))
gg6 <- gg6 + theme(legend.title = element_text(size = 12))
gg62 <- gg6 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg62

pdf("figure/marginal_est_coef.pdf",width = 8,height = 8,onefile = TRUE)
plot_grid(gg51,gg52,gg61,gg62,ncol=2)
dev.off()

pdf(file="figure/MapData_MostExtreme.pdf",width=8,height=8,onefile=TRUE)
plot_grid(gg1,gg2,gg3,gg4,ncol=2)
dev.off()

### SPATIAL map of temperature field corresponding to the "most unusual" day overall (based on sum of residuals)
my.day <- which.max(rowSums(residuals2))
data.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = data.sub[my.day,]) ## Create data frame for plotting SST data
data.mean.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.mean[my.day,]) ## Create data frame for plotting mean SST
data.sd.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.sd[my.day,]) ## Create data frame for plotting SD SST
data.residuals.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = residuals2[my.day,]) ## Create data frame for plotting residual SST data
data.mean.trend.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.par.mean[2,]) ## Create data frame for plotting mean trend SST
data.sd.trend.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = est.par.sd[2,]) ## Create data frame for plotting SD trend SST

my.point <- data.frame(lon = loc.sub[my.loc,1], lat = loc.sub[my.loc,2], temp = data.sub[my.day,my.loc]) ## If wants to add selected points on map

#key = 'AIzaSyABn_DX4l8KocC6QCV_lzP-1M4E3qR0Ki0'
#register_google(key = key)
xlim <- c(32,44)
ylim <- c(12,31)
#map <- get_map(location=c(xlim[1],ylim[1],xlim[2],ylim[2]),zoom=5,source="google",color="color",maptype="satellite",force=TRUE)
#save(map,file="data/GoogleMap.RData")
#load(file="data/GoogleMap.RData")
ewbreaks <- c(32, 36, 40, 44)
nsbreaks <- c(12, 18, 24, 30)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))

### 1st: plot the raw data
gg1 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg1 <- gg1 + scale_color_gradientn(colours = tim.colors())
gg1 <- gg1 + ggtitle("SST raw data (Most 'unusual' day)")
gg1 <- gg1 + xlab("Longitude")
gg1 <- gg1 + ylab("Latitude")
gg1 <- gg1 + labs(colour = expression(paste("[ºC]")))
gg1 <- gg1 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg1 <- gg1 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg1 <- gg1 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg1 <- gg1 + theme(axis.ticks = element_line(size = 1))
gg1 <- gg1 + theme(axis.ticks.length = unit(.1, "cm"))
gg1 <- gg1 + theme(legend.text = element_text(size = 12))
gg1 <- gg1 + theme(legend.key.width = unit(0.5, "cm"))
gg1 <- gg1 + theme(legend.title = element_text(size = 12))
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg1

### 2nd: plot the estimated temperature mean
gg2 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.mean.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg2 <- gg2 + scale_color_gradientn(colours = tim.colors())
gg2 <- gg2 + ggtitle("SST estimated mean (Most 'unusual' day)")
gg2 <- gg2 + xlab("Longitude")
gg2 <- gg2 + ylab("Latitude")
gg2 <- gg2 + labs(colour = expression(paste("[ºC]")))
gg2 <- gg2 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg2 <- gg2 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg2 <- gg2 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg2 <- gg2 + theme(axis.ticks = element_line(size = 1))
gg2 <- gg2 + theme(axis.ticks.length = unit(.1, "cm"))
gg2 <- gg2 + theme(legend.text = element_text(size = 12))
gg2 <- gg2 + theme(legend.key.width = unit(0.5, "cm"))
gg2 <- gg2 + theme(legend.title = element_text(size = 12))
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg2

### 3rd: plot the estimated temperature SD
gg3 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.sd.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg3 <- gg3 + scale_color_gradientn(colours = tim.colors())
gg3 <- gg3 + ggtitle("SST estimated SD (Most 'unusual' day)")
gg3 <- gg3 + xlab("Longitude")
gg3 <- gg3 + ylab("Latitude")
gg3 <- gg3 + labs(colour = expression(paste("[ºC]")))
gg3 <- gg3 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg3 <- gg3 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg3 <- gg3 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg3 <- gg3 + theme(axis.ticks = element_line(size = 1))
gg3 <- gg3 + theme(axis.ticks.length = unit(.1, "cm"))
gg3 <- gg3 + theme(legend.text = element_text(size = 12))
gg3 <- gg3 + theme(legend.key.width = unit(0.5, "cm"))
gg3 <- gg3 + theme(legend.title = element_text(size = 12))
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg3

### 4th: plot the temperature residuals
gg4 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg4 <- gg4 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.residuals.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg4 <- gg4 + scale_color_gradientn(colours = tim.colors())
gg4 <- gg4 + ggtitle("SST residuals (Most 'unusual' day)")
gg4 <- gg4 + xlab("Longitude")
gg4 <- gg4 + ylab("Latitude")
gg4 <- gg4 + labs(colour = expression(paste("")))
gg4 <- gg4 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg4 <- gg4 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg4 <- gg4 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg4 <- gg4 + theme(axis.ticks = element_line(size = 1))
gg4 <- gg4 + theme(axis.ticks.length = unit(.1, "cm"))
gg4 <- gg4 + theme(legend.text = element_text(size = 12))
gg4 <- gg4 + theme(legend.key.width = unit(0.5, "cm"))
gg4 <- gg4 + theme(legend.title = element_text(size = 12))
gg4 <- gg4 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg4

pdf(file="figure/MapData_MostUsual.pdf",width=8,height=8,onefile=TRUE)
plot_grid(gg1,gg2,gg3,gg4,ncol=2)
dev.off()


### SPATIAL map of annual temperature maxima corresponding to the "most extreme" years (in terms of max/sum/min)
data.maxima.max.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = log(maxima.frechet[which.max(apply(maxima.frechet,1,max)),])) ## Create data frame for plotting SST maxima
data.maxima.sum.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = log(maxima.frechet[which.max(apply(maxima.frechet,1,sum)),])) ## Create data frame for plotting SST maxima
data.maxima.min.plot <- data.frame(lon = loc.sub[,1], lat = loc.sub[,2], temp = log(maxima.frechet[which.max(apply(maxima.frechet,1,min)),])) ## Create data frame for plotting SST maxima

my.point <- data.frame(lon = loc.sub[my.loc,1], lat = loc.sub[my.loc,2], temp = data.sub[my.day,my.loc]) ## If wants to add selected points on map

#key = 'AIzaSyABn_DX4l8KocC6QCV_lzP-1M4E3qR0Ki0'
#register_google(key = key)
xlim <- c(32,44)
ylim <- c(12,31)
#map <- get_map(location=c(xlim[1],ylim[1],xlim[2],ylim[2]),zoom=5,source="google",color="color",maptype="satellite",force=TRUE)
#save(map,file="data/GoogleMap.RData")
#load(file="data/GoogleMap.RData")
ewbreaks <- c(32, 36, 40, 44)
nsbreaks <- c(12, 18, 24, 30)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))

### 1st: plot annual maxima for 2003 (sum of maxima is larger)
gg1 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.maxima.sum.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg1 <- gg1 + scale_color_gradientn(colours = tim.colors())
gg1 <- gg1 + ggtitle("2003 SST annual maxima on Gumbel scale")
gg1 <- gg1 + xlab("Longitude")
gg1 <- gg1 + ylab("Latitude")
gg1 <- gg1 + labs(colour = expression(paste("[ºC]")))
gg1 <- gg1 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg1 <- gg1 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg1 <- gg1 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg1 <- gg1 + theme(axis.ticks = element_line(size = 1))
gg1 <- gg1 + theme(axis.ticks.length = unit(.1, "cm"))
gg1 <- gg1 + theme(legend.text = element_text(size = 12))
gg1 <- gg1 + theme(legend.key.width = unit(0.5, "cm"))
gg1 <- gg1 + theme(legend.title = element_text(size = 12))
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg1

### 2nd: plot annual maxima for 2009 (max of maxima is larger)
gg2 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.maxima.max.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg2 <- gg2 + scale_color_gradientn(colours = tim.colors())
gg2 <- gg2 + ggtitle("2009 SST annual maxima on Gumbel scale")
gg2 <- gg2 + xlab("Longitude")
gg2 <- gg2 + ylab("Latitude")
gg2 <- gg2 + labs(colour = expression(paste("[ºC]")))
gg2 <- gg2 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg2 <- gg2 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg2 <- gg2 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg2 <- gg2 + theme(axis.ticks = element_line(size = 1))
gg2 <- gg2 + theme(axis.ticks.length = unit(.1, "cm"))
gg2 <- gg2 + theme(legend.text = element_text(size = 12))
gg2 <- gg2 + theme(legend.key.width = unit(0.5, "cm"))
gg2 <- gg2 + theme(legend.title = element_text(size = 12))
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg2

### 3rd: plot annual maxima for 2010 (min of maxima is larger)
gg3 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.maxima.min.plot, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg3 <- gg3 + scale_color_gradientn(colours = tim.colors())
gg3 <- gg3 + ggtitle("2010 SST annual maxima on Gumbel scale")
gg3 <- gg3 + xlab("Longitude")
gg3 <- gg3 + ylab("Latitude")
gg3 <- gg3 + labs(colour = expression(paste("[ºC]")))
gg3 <- gg3 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg3 <- gg3 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg3 <- gg3 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg3 <- gg3 + theme(axis.ticks = element_line(size = 1))
gg3 <- gg3 + theme(axis.ticks.length = unit(.1, "cm"))
gg3 <- gg3 + theme(legend.text = element_text(size = 12))
gg3 <- gg3 + theme(legend.key.width = unit(0.5, "cm"))
gg3 <- gg3 + theme(legend.title = element_text(size = 12))
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = "white"), data = my.point, alpha = 1, shape = 16, size = 3, col = "black")
gg3

pdf(file="figure/MapAnnualMaxima.pdf",width=4*3,height=4,onefile=TRUE)
plot_grid(gg1,gg2,gg3,ncol=3)
dev.off()

## plot the fitted GEV parameters
#key = 'AIzaSyABn_DX4l8KocC6QCV_lzP-1M4E3qR0Ki0'
#register_google(key = key)
xlim <- c(32,44)
ylim <- c(12,31)
#map <- get_map(location=c(xlim[1],ylim[1],xlim[2],ylim[2]),zoom=5,source="google",color="color",maptype="satellite",force=TRUE)
#save(map,file="data/GoogleMap.RData")
#load(file="data/GoogleMap.RData")
ewbreaks <- c(32, 36, 40, 44)
nsbreaks <- c(12, 18, 24, 30)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))

### 1st: plot location parameter
data.df = data.frame(lon=loc.sub[,1],lat=loc.sub[,2],temp=est.par.GEV[1,])
gg1 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg1 <- gg1 + geom_point(aes(x = lon, y = lat, colour = temp),data=data.df, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg1 <- gg1 + scale_color_gradientn(colours = tim.colors())
gg1 <- gg1 + ggtitle(expression(hat(mu)[i]~'in GEV distribution'))
gg1 <- gg1 + xlab("Longitude")
gg1 <- gg1 + ylab("Latitude")
gg1 <- gg1 + labs(colour = expression(paste("Value")))
gg1 <- gg1 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg1 <- gg1 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg1 <- gg1 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg1 <- gg1 + theme(axis.ticks = element_line(size = 1))
gg1 <- gg1 + theme(axis.ticks.length = unit(.1, "cm"))
gg1 <- gg1 + theme(legend.text = element_text(size = 12))
gg1 <- gg1 + theme(legend.key.width = unit(0.5, "cm"))
gg1 <- gg1 + theme(legend.title = element_text(size = 12))
gg1

### 2nd: plot scale parameter
data.df = data.frame(lon=loc.sub[,1],lat=loc.sub[,2],temp=est.par.GEV[2,])
gg2 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg2 <- gg2 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.df, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg2 <- gg2 + scale_color_gradientn(colours = tim.colors())
gg2 <- gg2 + ggtitle(expression(hat(sigma)[i]~'in GEV distribution'))
gg2 <- gg2 + xlab("Longitude")
gg2 <- gg2 + ylab("Latitude")
gg2 <- gg2 + labs(colour = expression(paste("Value")))
gg2 <- gg2 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg2 <- gg2 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg2 <- gg2 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg2 <- gg2 + theme(axis.ticks = element_line(size = 1))
gg2 <- gg2 + theme(axis.ticks.length = unit(.1, "cm"))
gg2 <- gg2 + theme(legend.text = element_text(size = 12))
gg2 <- gg2 + theme(legend.key.width = unit(0.5, "cm"))
gg2 <- gg2 + theme(legend.title = element_text(size = 12))
gg2

### 3rd: plot shape parameter
data.df = data.frame(lon=loc.sub[,1],lat=loc.sub[,2],temp=est.par.GEV[3,])
gg3 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg3 <- gg3 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.df, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg3 <- gg3 + scale_color_gradientn(colours = tim.colors())
gg3 <- gg3 + ggtitle(expression(hat(xi)[i]~'in GEV distribution'))
gg3 <- gg3 + xlab("Longitude")
gg3 <- gg3 + ylab("Latitude")
gg3 <- gg3 + labs(colour = expression(paste("Value")))
gg3 <- gg3 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg3 <- gg3 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg3 <- gg3 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg3 <- gg3 + theme(axis.ticks = element_line(size = 1))
gg3 <- gg3 + theme(axis.ticks.length = unit(.1, "cm"))
gg3 <- gg3 + theme(legend.text = element_text(size = 12))
gg3 <- gg3 + theme(legend.key.width = unit(0.5, "cm"))
gg3 <- gg3 + theme(legend.title = element_text(size = 12))
gg3

data.df = data.frame(lon=loc.sub[,1],lat=loc.sub[,2],temp=apply(maxima.frechet,2,function(x){ks.test(x,pfrechet,alternative = "two.sided")$p.value}))
gg4 <- ggmap(map, extent = "panel", legend = "bottomleft")
gg4 <- gg4 + geom_point(aes(x = lon, y = lat, colour = temp), data = data.df, alpha = 1, na.rm = T, shape = 15, size = 1.0)
gg4 <- gg4 + scale_color_gradientn(colours = tim.colors())
gg4 <- gg4 + ggtitle('P-values in KS test')
gg4 <- gg4 + xlab("Longitude")
gg4 <- gg4 + ylab("Latitude")
gg4 <- gg4 + labs(colour = expression(paste("Value")))
gg4 <- gg4 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = c(29,47))
gg4 <- gg4 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(9,34))
gg4 <- gg4 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg4 <- gg4 + theme(axis.ticks = element_line(size = 1))
gg4 <- gg4 + theme(axis.ticks.length = unit(.1, "cm"))
gg4 <- gg4 + theme(legend.text = element_text(size = 12))
gg4 <- gg4 + theme(legend.key.width = unit(0.5, "cm"))
gg4 <- gg4 + theme(legend.title = element_text(size = 12))
gg4

pdf(file="figure/Map_GEV_parameter.pdf",width=4*2,height=4*2,onefile=TRUE)
plot_grid(gg1,gg2,gg3,gg4,ncol=2)
dev.off()

########################################
### EMPIRICAL BIVARIATE COEFFICIENTS ###
########################################
D <- nrow(loc.sub)
all.pairs <- comb_n(D,2)
n.pairs <- choose(D,2)
M <- n.pairs ### number of selected pairs of sites
sub.pairs <- all.pairs[,sample(c(1:n.pairs),M,replace=FALSE)] ### selected pairs

### distance for selected pairs
sub.dists <- sapply(1:M,
                    FUN=function(i){
                      return(distmat[sub.pairs[1,i],sub.pairs[2,i]])
                    })
### empirical extremal coefficients for selected pairs
sub.extcoef <- sapply(1:M,
                    FUN=function(i){
                      return(min(2,max(1,1/mean(1/pmax(maxima.frechet[,sub.pairs[1,i]],maxima.frechet[,sub.pairs[2,i]])))))
                    })
### scatterplot of empirical bivariat extremal coefficient vs distance
par(mfrow=c(1,1))
plot(sub.dists,sub.extcoef,pch=20,ylim=c(1,2),xlab="Distance",ylab="Empirical bivariate extremal coefficients")

n.classes <- 50
width.classes <- max(distmat)/n.classes
classes.mid <- width.classes/2 + c(0:(n.classes-1))*width.classes
class.extcoef <- sapply(1:n.classes,
                         FUN=function(i){
                           return(sub.extcoef[which(sub.dists>classes.mid[i]-width.classes/2-0.5 & sub.dists<classes.mid[i]+width.classes/2+0.5)])   
                         })
grp.extcoef <- as.factor(sapply(1:M,
                         FUN=function(i){
                           return(ceiling(sub.dists[i]/width.classes))   
                         }))
plot(sub.dists,sub.extcoef,type="n",pch=20,ylim=c(1,2),xlab="Distance",ylab="Empirical bivariate extremal coefficients")
boxplot(sub.extcoef~grp.extcoef,add=TRUE,at=classes.mid,xaxt="n",boxwex=30)
abline(h=c(1,2),col="lightgrey")

#save(maxima.frechet,distmat,loc.sub,file="data/data.RData")

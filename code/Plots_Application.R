### load the data and plot ###
library(mvtnorm)
library(Rfast) 
library(ggplot2)
library(partitions)
library(evd)
library(cowplot)
library(ggmap)
library(fields)
library(rgdal)

source("code/MLE_BrownResnick.R"); 
load("data/data.RData");load("data/prediction_estimates4_1.Rdata")
tiylim <- c(11,31)
#map <- get_map(location=c(xlim[1],ylim[1],xlim[2],ylim[2]),zoom=5,source="google",
#	color="color",maptype="satellite",force=TRUE)
load("data/GoogleMap.RData")
ewbreaks <- c(32, 36, 40, 44)
nsbreaks <- c(12, 18, 24, 30)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))

pdf("figure/loc_ratio.pdf",width = 15,height = 20,onefile = TRUE)
ratio=c(1,1/2,1/4,1/8,1/16,9/10)
for(ratio.ind in 1:6){
ind.sub <- c()
ind.sub[1]<-c(which.min(distmat[which.min(colmeans(distmat))[1],])) ## Vecchia sequence based on max-min ordering
for(j in 2:nrow(loc.sub.trans)){
  ind.sub[j] <- c(1:nrow(loc.sub.trans))[-ind.sub[1:(j-1)] ][which.max(rowMins(matrix(distmat[-ind.sub[1:(j-1)],ind.sub[1:(j-1)]],ncol=(j-1)),value=TRUE))[1]]
}
if(ratio.ind==6){
ind.sub <- rev(ind.sub)[1:floor(nrow(loc.sub.trans)*ratio[ratio.ind])]
}else{
ind.sub <- ind.sub[1:floor(nrow(loc.sub.trans)*ratio[ratio.ind])]
}
neighbours <- function(ind,ind.sub,q){
  ind.neighbours <- ind.sub[order(distmat[ind,ind.sub])[c(1:q,sample((q+1):length(ind.sub),4))]]
  return(ind.neighbours)
}
neighbours2predict <- sapply((1:nrow(distmat))[-ind.sub],FUN=neighbours,ind.sub=ind.sub,q=4)
par(mfrow=c(1,1),cex=2,cex.axis=2,cex.lab=2,mar=c(5,5,1,1))
plot(loc.sub.trans[ind.sub,1],loc.sub.trans[ind.sub,2],xlab="Transformed x (km)",ylab="Transformed y (km)",col="black",pch=20)
if(ratio.ind!=1){
points(loc.sub.trans[-ind.sub,1],loc.sub.trans[-ind.sub,2],col="red",pch=20)
}
#ind.neigbours <- unique(c(neighbours2predict))
#points(loc.sub.trans[ind.neigbours,1],loc.sub.trans[ind.neigbours,2],col="blue",pch=20)
}
dev.off()


### plot the contour of the variance ###
vario.contour<-function(loc,par){
  D = nrow(loc)
  val <- sapply(1:D,
                FUN=function(i){
                  loc = rbind(loc[i,],c(0,0))
                  sigma <- vario.func(loc,par=par.val)
                  return(V.biv(matrix(c(1,1),1,2),sigma))
                })
  return(val)
}
x <- seq(xlim[1],xlim[2],length=500)
y <- seq(ylim[1],ylim[2],length=500)
loc.contour <-as.matrix(expand.grid(x,y))
loc.contour.trans <- project(loc.contour, "+proj=utm +zone=37 ellps=WGS84")/1000
m.trans = apply(loc.contour.trans,2,mean)
loc.contour.trans <- t(apply(loc.contour.trans,1,function(x){x-m.trans}))
par.val = estimates[which.min(nlog.pred.val),]
z <- vario.contour(loc.contour.trans,par.val)
ind.fit = rep(F,nrow(loc.sub))
ind.fit[ind.sub] = T
data.contour <- data.frame(lon=loc.contour[,1],lat=loc.contour[,2],z=z)
data.df = data.frame(lon=loc.sub[,1],lat=loc.sub[,2],ind.fit=ind.fit)

gg1 <- ggmap(map, extent = "device", legend = "right")
gg1 <- gg1 + geom_point(aes(x = lon, y = lat),colour="#34568B",data=data.df, alpha = 1, na.rm = T, shape = 16, size = 0.5)
gg1 <- gg1 + stat_contour(data=data.contour,aes(x=lon,y=lat,z=z,colour=after_stat(level)),binwidth = 0.1)
gg1 <- gg1 + scale_color_gradientn(colours = tim.colors(),limit=c(1,2),na.value ="white")
gg1 <- gg1 + ggtitle("Contours of fitted extremal coefficients")
gg1 <- gg1 + xlab("Longitude")
gg1 <- gg1 + ylab("Latitude")
gg1 <- gg1 + labs(colour = expression(hat(theta)))
gg1 <- gg1 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = xlim)
gg1 <- gg1 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = ylim)
gg1 <- gg1 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg1 <- gg1 + theme(axis.ticks = element_line(size = 1))
gg1 <- gg1 + theme(axis.ticks.length = unit(.1, "cm"))
gg1 <- gg1 + theme(legend.text = element_text(size = 12))
gg1 <- gg1 + theme(legend.key.width = unit(0.5, "cm"))
gg1 <- gg1 + theme(legend.title = element_text(size = 12))
print(gg1)

gg2 <- ggmap(map, extent = "device", legend = "right")
gg2 <- gg2 + geom_point(aes(x = lon, y = lat),colour="#34568B",data=subset(data.df,!ind.fit), alpha = 1, na.rm = T, shape = 16, size = 1.0)
gg2 <- gg2 + ggtitle("Validation sites")
gg2 <- gg2 + xlab("Longitude")
gg2 <- gg2 + ylab("Latitude")
gg2 <- gg2 + scale_x_continuous(breaks = ewbreaks, labels = ewlabels, expand = c(0, 0), limits = xlim)
gg2 <- gg2 + scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = ylim)
gg2 <- gg2 + theme(axis.text = element_text(size=12), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),plot.title = element_text(hjust = 0.5))
gg2 <- gg2 + theme(axis.ticks = element_line(size = 1))
gg2 <- gg2 + theme(axis.ticks.length = unit(.1, "cm"))
gg2 <- gg2 + theme(legend.text = element_text(size = 12))
gg2 <- gg2 + theme(legend.key.width = unit(0.5, "cm"))
gg2 <- gg2 + theme(legend.title = element_text(size = 12))
print(gg2)

pdf("figure/variogram_contour_3.pdf",width = 4.2*2,height = 5)
plot_grid(gg1,gg2,ncol=2,rel_widths = c(1.28,1))
dev.off()

## plot of extremal coeffficents ###
load("data/prediction_estimates4_1.Rdata")
#load("data/prediction_estimates_fixed.Rdata")
par.mat = estimates
D <- nrow(loc.sub.trans)
all.pairs <- comb_n(D,2)
n.pairs <- choose(D,2)
M <- n.pairs ### number of selected pairs of sites
sub.pairs <- all.pairs[,sample(c(1:n.pairs),M,replace=FALSE)] ### selected pairs

### empirical extremal coefficients for selected pairs
sub.extcoef <- sapply(1:M,
                      FUN=function(i){
                        return(min(2,max(1,1/mean(1/pmax(maxima.frechet[,sub.pairs[1,i]],maxima.frechet[,sub.pairs[2,i]])))))
                      })

ext.boxplot <- function(id,n.classes = 50){
  par.val = par.mat[id,]
  A = matrix(c(cos(par.val[4]),sin(par.val[4]),-sin(par.val[4]),cos(par.val[4])),2,2)
  Sigma <- A%*%diag(c(1,par.val[3]),2)%*%t(A)
  loc.sub.trans.back <- loc.sub.trans %*% chol(Sigma)
  distmat.trans <- as.matrix(dist(loc.sub.trans.back,method="euclidean"))
  sub.dists <- sapply(1:M,
                      FUN=function(i){
                        return(distmat.trans[sub.pairs[1,i],sub.pairs[2,i]])
                      })
  ### fitted extremal coefficients for selected pairs
  extcoef.Pair = matrix(NA,ncol=2,nrow=2000)
  extcoef.Pair[,2] <- dist.pairs <- seq(min(sub.dists),max(sub.dists),length.out = 2000)
  loc.pairs <- cbind(0,dist.pairs) %*% solve(chol(Sigma))
  extcoef.Pair[,1] <- sapply(1:nrow(loc.pairs),
                             FUN=function(i){
                               loc = rbind(0,loc.pairs[i,])
                               sigma <- vario.func(loc,par=par.val)
                               return(V.biv(matrix(c(1,1),1,2),sigma))
                             })
  ### distance for selected pairs
  width.classes <- diff(range(sub.dists))/n.classes
  min.dist = min(sub.dists)
  classes.mid <- width.classes/2 + min.dist + c(0:(n.classes-1))*width.classes
  grp.extcoef <- pmin(sapply(1:M,
                                       FUN=function(i){
                                         return(ceiling((sub.dists[i]- min.dist + 0.1)/width.classes))
                                       }),n.classes)
  classes.mid = sort(classes.mid[unique(grp.extcoef)])
  #pdf(paste0("figure/boxplot_ext_fixed/boxplot_distance_",scenarios[order.ind[id,1],1],"_",scenarios[order.ind[id,1],2],"_",scenarios[order.ind[id,1],3],"_",scenarios[order.ind[id,1],4],"_",order.ind[id,2],".pdf"))
  pdf(paste("figure/boxplot_ext_3/boxplot_distance",scenarios[id,1],scenarios[id,2]*2,scenarios[id,3],scenarios[id,4],scenarios[id,5],".pdf",sep = "_"))
  par(mfrow=c(1,1),mar=c(3.5,4,2,0),cex.lab=1.5,mgp=c(2,1,0))
  plot(sub.dists,sub.extcoef,type="n",pch=20,ylim=c(1,2),xlab="Distance [km]",ylab="Empirical bivariate extremal coefficients")
  boxplot(sub.extcoef~grp.extcoef,add=TRUE,at=classes.mid,xaxt="n",boxwex=15,outline=FALSE)
  lines(extcoef.Pair[,2],extcoef.Pair[,1],lty=1,lwd=3,col="red")
  abline(h=c(1,2),col="lightgrey")
  dev.off()
}
sapply(1:nrow(par.mat),ext.boxplot,n.classes=100)

### 
D = nrow(loc.sub.trans)
pairs <- t(combn(D,2))
angle <- function(id){
   loc.pair<- loc.sub.trans[pairs[id,1],] - loc.sub.trans[pairs[id,2],]
   val <- (atan2(loc.pair[1],loc.pair[2]) + pi) %% pi
   return(val)
}
angle.ip <- sapply(1:nrow(pairs),angle)
all.pairs <- comb_n(D,2)
sub.dists <- sapply(1:ncol(all.pairs),
                    FUN=function(i){
                      return(distmat[all.pairs[1,i],all.pairs[2,i]])
                    })
sub.extcoef <- sapply(1:ncol(all.pairs),
                      FUN=function(i){
                        return(min(2,max(1,1/mean(1/pmax(maxima.frechet[,all.pairs[1,i]],maxima.frechet[,all.pairs[2,i]])))))
                      })
directions = c(15,45,75,105,135,165)
dists <- seq(min(sub.dists),max(sub.dists),length.out=1000)
n.classes=20
#n.ratio = 6
#model.ind = sapply(2:n.ratio,function(i){which(order.ind[,2]==i)[which.min(nlog.pred.val[order.ind[,2]==i])]} )
fun <- function(id,dists,direction){
  par.val = estimates[id,]
  extcoef.Pair <- sapply(dists,
                         FUN=function(dist.i){
                           loc = matrix(c(0,0,dist.i*cos(direction/pi),dist.i*sin(direction/pi)),2,2,byrow = T)
                           sigma <- vario.func(loc,par=par.val)
                           return(V.biv(matrix(c(1,1),1,2),sigma))
                         })
  return(extcoef.Pair)
}

ind.angle <- sub.dists.angle <- width.classes <- classes.mid <-min.dist <- grp.extcoef <- list()
for(i in 1:length(directions)){
ind.angle[[i]]  <- ( angle.ip >= (directions[i]-1)/180*pi & angle.ip <= (directions[i]+1)/180*pi )
sub.dists.angle[[i]]  = sub.dists[ind.angle[[i]]]
width.classes[[i]]  <- diff(range(sub.dists))/n.classes
min.dist[[i]]  = min(sub.dists)
classes.mid[[i]]  <- width.classes[[i]] /2 + min.dist[[i]]  + c(0:(n.classes-1))*width.classes[[i]]
grp.extcoef[[i]] <- pmin(sapply(1:sum(ind.angle[[i]]),
                           FUN=function(id){
                             return(ceiling((sub.dists.angle[[i]][id] - min.dist[[i]] + 0.1)/width.classes[[i]]))
                           }),n.classes)
classes.mid[[i]] = sort(classes.mid[[i]][unique(grp.extcoef[[i]])])
}

index=unique(order.ind[,1])
#index = 1:nrow(estimates)
for(j in index){
model.ind=which(order.ind[,1] == j & order.ind[,2]!=6)
pdf(paste0("figure/boxplot_ext_fixed/boxplot_distance_angle_",scenarios[j,1],"_",scenarios[j,2]*2,"_",scenarios[j,3],"_",scenarios[j,4],".pdf"),width=4*3,height=4*2)
par(mfrow=c(2,3),mar=c(3.5,2,3.5,0),pty="s",cex.lab=1.5,mgp=c(2,1,0))
for(i in 1:length(directions)){
  plot(sub.dists.angle[[i]],sub.extcoef[ind.angle[[i]]],xlim=range(sub.dists),type="n",pch=20,ylim=c(1,2),xlab="Distance [km]",ylab="Empirical bivariate extremal coefficients",main=paste0(directions[i],"º"))
   
  boxplot(sub.extcoef[ind.angle[[i]]]~grp.extcoef[[i]],add=TRUE,at=classes.mid[[i]],xaxt="n",boxwex=25,outline=FALSE)
  extcoef.Pair = lapply(model.ind,fun,dists=dists,direction=directions[i])
  for(id in 1:length(extcoef.Pair)){
  lines(dists,extcoef.Pair[[id]],lty=1,col=id)
  }
  abline(h=c(1,2),col="lightgrey")
  }
dev.off()
}

index=unique(order.ind[,1])
library(RColorBrewer)
for(j in index){
  pdf(paste0("figure/boxplot_ext_fixed/boxplot_distance_angle_",scenarios[j,1],"_",scenarios[j,2]*2,"_",scenarios[j,3],"_",scenarios[j,4],".pdf"),width=4*3,height=4*2)
  par(mfrow=c(2,3),mar=c(3.5,2,3.5,0),pty="s",cex.lab=1.5,mgp=c(2,1,0),cex.main=1.5)
  model.ind=which(order.ind[,1] == j & order.ind[,2]!=6)
  for(i in 1:length(directions)){
    plot(sub.dists.angle[[i]],sub.extcoef[ind.angle[[i]]],xlim=range(sub.dists),type="n",pch=20,ylim=c(1,2),xlab="Distance [km]",ylab="Empirical bivariate extremal coefficients",main=paste0(directions[i],"º"))
    grp.table <- as.numeric(table(grp.extcoef[[i]]))
    grp.table = 1-grp.table/max(grp.table)
    box.colors <- grey(level=grp.table,alpha=0.5) 
    boxplot(sub.extcoef[ind.angle[[i]]]~grp.extcoef[[i]],add=TRUE,at=classes.mid[[i]],xaxt="n",boxwex=35,outline=FALSE,col=box.colors)
    extcoef.Pair = lapply(model.ind,fun,dists=dists,direction=directions[i])
    colors <- rev(brewer.pal(n=length(extcoef.Pair)+1,name="Reds")[-1])
    for(id in 1:length(extcoef.Pair)){
      lines(dists,extcoef.Pair[[id]],lty=1,lwd=1.5,col=colors[id])
    }
    abline(h=c(1,2),col="lightgrey")
  }
  dev.off()
}

library(RColorBrewer)
pdf(paste0("figure/boxplot_ext_fixed/boxplot_distance_angle.pdf"),width=4*3,height=4*2)
par(mfrow=c(2,3),mar=c(3.5,2,3.5,0),pty="s",cex.lab=1.5,mgp=c(2,1,0),cex.main=1.5)
for(i in 1:length(directions)){
    plot(sub.dists.angle[[i]],sub.extcoef[ind.angle[[i]]],xlim=range(sub.dists),type="n",pch=20,ylim=c(1,2),xlab="Distance [km]",ylab="Empirical bivariate extremal coefficients",main=paste0(directions[i],"º"))
    #box.colors <- grey.colors(n=length(unique(grp.extcoef[[i]]))+1,start=0.1,end=1,alpha=0.5,rev=TRUE)[-1][rank(table(grp.extcoef[[i]]))]
    grp.table <- as.numeric(table(grp.extcoef[[i]]))
    grp.table = 1-grp.table/max(grp.table)
    box.colors <- grey(level=grp.table,alpha=0.5) 
    boxplot(sub.extcoef[ind.angle[[i]]]~grp.extcoef[[i]],add=TRUE,at=classes.mid[[i]],xaxt="n",boxwex=35,outline=FALSE,col=box.colors)
    model.ind=which(order.ind[,1] == 1 & order.ind[,2]!=6)
    extcoef.Pair = lapply(model.ind,fun,dists=dists,direction=directions[i])
    colors <- rev(brewer.pal(n=length(extcoef.Pair)+1,name="Reds")[-1])
    for(id in 1:length(extcoef.Pair)){
      lines(dists,extcoef.Pair[[id]],lty=1,lwd=1.5,col=colors[id])
    }
    model.ind=which(order.ind[,1] == 14 & order.ind[,2]!=6)
    extcoef.Pair = lapply(model.ind,fun,dists=dists,direction=directions[i])
    colors <- rev(brewer.pal(n=length(extcoef.Pair)+1,name="Blues")[-1])
    for(id in 1:length(extcoef.Pair)){
      lines(dists,extcoef.Pair[[id]],lty=1,lwd=1.5,col=colors[id])
    }
    abline(h=c(1,2),col="lightgrey")
}
dev.off()


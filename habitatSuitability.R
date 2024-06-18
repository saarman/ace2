# load library
library(raster)
library(maps)
library(mapdata)
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library(fields)
library(mapplots)
library(sp)
library(classInt)

# import netCDF files from https://github.com/lanl/culexmaxentmodels
pip <- brick("culex_pipiens_meansuitability.nc")
quinq.full <- brick("culex_quinquefasciatus_meansuitability.nc")

# set extent to match
quinq <- setExtent(quinq.full,extent(pip),keepres = TRUE)

# plot both distributions
plot(pip)
plot(quinq)

#############################################
# Import site sampling points, make spatial for later
#############################################
# add sampled points
sites <- read.csv("ace2_by_site-id.csv")
sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))

# plot both distributions with sampled points
plot(pip)
points(sitesSp, pch = 1)

plot(quinq)
points(sitesSp, pch = 1)


#############################################
# Plot maps with complementary colors:
#############################################
colfuncX <- colorRampPalette(c("white", "#EFCA17", "red"))
#colfuncX <- colorRampPalette(c("#FFE50E","#821250"))
colfuncY <- colorRampPalette(c("white","#0096EC","#701363"))
#colfuncY <- colorRampPalette(c("#75AFCE","#612576"))

#pdf(file="X.pdf", width = 8, height = 9)
plot(pip, main="Cx. pipiens", col=colfuncY(25),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

#pdf(file="Y.pdf", width = 8, height = 9)
plot(quinq,main="Cx. quinquefasciatus", col=colfuncX(25),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

#############################################
#create a color matrix for Bivariate map:
#############################################
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=9, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")

#############################################
#Build bivariate map function:
#############################################
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(temp,na.rm=T, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}

#############################################
# Use the bivariate.map function with quinq and pip
#############################################

bivmap<-bivariate.map(quinq,pip,colormatrix=col.matrix, nquantiles=9)
# Note: go back to the function that generates color matrices and try out new color schemes.

pdf(file="pipiens_vs_quinqs.pdf", width = 8, height = 9)
plot(bivmap,main="Cx. pipiens vs Cx. quinquefasciatus", col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,add=F,legend=F)
points(sitesSp, pch = 1)
dev.off()

pdf(file="legend.pdf", width = 8, height = 9)
col.matrix<-colmat(nquantiles=9, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")
dev.off()

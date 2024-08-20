# load library
library(raster)
library(maps)
library(mapdata)
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library(fields)
library(mapplots)
library(sp)
library(classInt)
library(ggplot2)

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
sites <- read.csv("ace2_by_site_date.csv")
sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))

# plot both distributions with sampled points
plot(pip)
points(sitesSp, pch = 1)

plot(quinq)
points(sitesSp, pch = 1)


#############################################
# Extract suitability at each location, create data frame
#############################################

# start data frame and name fields
sitesDf <- as.data.frame(sites[,1:7])
names(sitesDf) <- c("locality","site","pp","pq","qq","latitude","longitude")

# extract habitat suitability value for pip and quinq layers at each point
sitesDf$pip_hs <- extract(pip, sitesSp)[,1]
sitesDf$quinq_hs <- extract(quinq, sitesSp)[,1]

# calculate + add: hybrid index,suitability index
sitesDf$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
sitesDf$h_index <- ((2*sitesDf$pp)+(sitesDf$pq))/sitesDf$alleles # p/p+q
sitesDf$s_index <- sitesDf$pip_hs/(sitesDf$pip_hs+sitesDf$quinq_hs) # pip suitability / total suitability
#plot(sitesDf$h_index,sitesDf$s_index, xlim = c(0,1), ylim = c(0,1))
#abline(0,1)

# name rows
rownames(sitesDf)<- sitesDf$site


#############################################
# Barplots of basic data and results
#############################################

#pdf(file="barplots.pdf", width = 8, height = 9)
par(mfrow = c(1, 4), mar=c(10,10,2,2))
barplot(sitesDf$h_index,names.arg = rownames(sitesDf),horiz = T,las=2,xlab="H-index")
barplot(sitesDf$s_index,names.arg = rownames(sitesDf),horiz = T,las = 2, xlab="S-index")
barplot(sitesDf$pip_hs,names.arg = rownames(sitesDf),horiz = T,las = 2, xlab = "HSM-pip")
barplot(sitesDf$quinq_hs,names.arg = rownames(sitesDf),horiz = T,las = 2, xlab="HSM-quinq")
#dev.off()

#with ggplot
data <- as.data.frame(sitesDf)
data$site <- rownames(sitesDf)

p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = h_index), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=rev) +
  ylab("H-index\nquinq --------> pip") +
  xlab("Collection Site") +
  coord_flip()

p2 <- ggplot(data=data, aes(x=site, y=s_index)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = h_index), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=rev) +
  ylab("S-index\nquinq ---------> pip") +
  xlab("Collection Site") +
  coord_flip()

p3.1 <- ggplot(data=data, aes(x=site, y=pip_hs)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = h_index), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=rev) +
  ylab("HSM-pip") +
  xlab(NULL) +
  coord_flip()

p3.2 <- ggplot(data=data, aes(x=site, y=quinq_hs)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = h_index), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=rev) +
  ylab("HSM-quinq") +
  xlab(NULL) +
  coord_flip()


library(patchwork)
#pdf(file="lineplots.pdf", width = 5, height = 5)
p1 + p2
#p1 + p3.1 + p2 + p3.2
dev.off()

############################################
# Linear model fit
############################################

s_index.lm <- lm(h_index~s_index, data = sitesDf)
summary(s_index.lm)

# Call:
#   lm(formula = h_index ~ s_index, data = sitesDf)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.62197 -0.10887  0.00215  0.12581  0.35882 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.2078     0.1236  -1.682    0.113    
# s_index       1.6911     0.2530   6.685 7.29e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.239 on 15 degrees of freedom
# Multiple R-squared:  0.7487,	Adjusted R-squared:  0.7319 
# F-statistic: 44.69 on 1 and 15 DF,  p-value: 7.29e-06

lat.lm <- lm(h_index~latitude, data = sitesDf)
summary(lat.lm)


pip.lm <- lm(h_index~pip_hs, data = sitesDf)
summary(pip.lm)

quinq.lm <- lm(h_index~quinq_hs, data = sitesDf)
summary(quinq.lm)


pip.quinq.lm <- lm(h_index~pip_hs*quinq_hs, data = sitesDf)
summary(pip.quinq.lm)


#############################################
# Scatter plot of hybrid index versus suitability
#############################################


p4 <- ggplot(data=data, aes(x=h_index,y=latitude)) + 
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  #scale_y_continuous(limits=c(0,1)) +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Latitude") +
  xlab("H-index\nquinq ---------------------> pip")

p5 <- ggplot(data=data, aes(x=h_index,y=s_index)) + 
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("S-index\nquinq ---------------------> pip") +
  xlab("H-index\nquinq ---------------------> pip")

p6.1 <- ggplot(data=data, aes(x=h_index, y=pip_hs)) + 
  geom_point(size=3) +
  #geom_text(aes(label = data$site), hjust = -.25) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("HSM-pip") +
  xlab("H-index\nquinq -------------------> pip")

p6.2 <- ggplot(data=data, aes(x=h_index,y=quinq_hs)) + 
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("HSM-quinq") +
  xlab("H-index\nquinq -------------------> pip")

#pdf(file="scatterplots.pdf", width = 10, height = 5)
p4 + p5
dev.off()

#############################################
# Plot maps with complementary colors:
#############################################

upperleft=rgb(0,150,235, maxColorValue=255)
upperright=rgb(130,0,80, maxColorValue=255)
bottomleft="white" #rgb(220,220,220,  maxColorValue=255)
bottomright=rgb(255,230,15, maxColorValue=255)

colfuncX <- colorRampPalette(c(bottomleft,upperleft))
#colfuncX <- colorRampPalette(c("white", "#EFCA17", "red"))
#colfuncX <- colorRampPalette(c("#FFE50E","#821250"))
colfuncY <- colorRampPalette(c(bottomleft,bottomright))
#colfuncY <- colorRampPalette(c("white","#0096EC","#701363"))
#colfuncY <- colorRampPalette(c("#75AFCE","#612576"))

#pdf(file="pipiens.pdf", width = 8, height = 9)
plot(pip, main="Cx. pipiens", col=colfuncX(4),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

#pdf(file="quinquefasciatus.pdf", width = 8, height = 9)
plot(quinq,main="Cx. quinquefasciatus", col=colfuncY(4),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

#############################################
# Truncate the pip and quinq anything below 0.25 --> NA
#############################################
pip.trunc <- pip
pip.trunc[pip.trunc < 0.25] <- NA

#pdf(file="pipiens.pdf", width = 8, height = 9)
plot(pip.trunc, main="Cx. pipiens", col=colfuncX(10),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

quinq.trunc <- quinq 
quinq.trunc[quinq.trunc < 0.25] <- NA

#pdf(file="quinquefasciatus.pdf", width = 8, height = 9)
plot(quinq.trunc,main="Cx. quinquefasciatus", col=colfuncY(10),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

#############################################
# Overlapping map with transparency
#############################################
library(scales)  #to use alpha() function

dev.off()
#pdf(file="pip_quinq_overlap_pies.pdf", width = 8, height = 9)
plot(pip, main="Cx. pipiens vs Cx. quinquefasciatus ", col=alpha(colfuncX(10),1),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(quinq,main="Cx. quinquefasciatus", col=alpha(colfuncY(10),0.75),frame.plot=F,axes=F,box=F,add=T,legend=F)

#add in pie charts source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
map("state", col = "grey", fill = FALSE, add=T)
# add from data table
for (i in 1:length(sitesDf$site)){
  add.pie(z = c(sitesDf[i,3],sitesDf[i,4],sitesDf[i,5]), x = sitesDf$longitude[i], y = sitesDf$latitude[i], clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 1, radius = 1 )
}
dev.off()

#pdf(file="pip_quinq_overlap.pdf", width = 8, height = 9)
plot(pip, main="Cx. pipiens vs Cx. quinquefasciatus ", col=alpha(colfuncX(10),1),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(quinq,main="Cx. quinquefasciatus", col=alpha(colfuncY(10),0.75),frame.plot=F,axes=F,box=F,add=T,legend=F)
dev.off()

#pdf(file="USA_ace2_pies.pdf", width = 8, height = 9)
#add in pie charts source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
map("state", col = "grey", fill = T, add=F)
# add from data table
for (i in 1:length(sitesDf$site)){
  add.pie(z = c(sitesDf[i,3],sitesDf[i,4],sitesDf[i,5]), x = sitesDf$longitude[i], y = sitesDf$latitude[i], clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 1, radius = 1 )
}
dev.off()









#    #############################################
#    # Bivariate map with new package bivariatemaps
#    #############################################
#    # https://www.youtube.com/watch?v=eJ40ffGQSCU
#    # https://rfunctions.blogspot.com/2015/03/bivariate-maps-bivariatemap-function.html
#    # now in a package https://rfunctions.blogspot.com/2022/05/new-r-package-on-cran-bivariatemaps.html
#    # The functions in bivariatemaps package are very similar to what we had from Giuseppe, maybe from the same source?
#    
#    # load dependencies
#    library(classInt)
#    library(raster)
#    library(rgdal)
#    library(dismo)
#    library(XML)
#    library(sp)
#    library(maps)
#    library(bivariatemaps)
#    
#    col.matrix<-colmat(nquantiles=10, xlab = "Cx. quinquefasciatus", ylab = "Cx. pipiens")
#    
#    bivmap<-bivariate.map(quinq, pip,colormatrix=col.matrix, nquantiles=10)
#    bivmap.trunc<-bivariate.map(quinq.trunc, pip.trunc, colormatrix=col.matrix, nquantiles=10)
#    
#    usa <- pip
#    usa[usa > 0] <- 1
#    plot(usa,col="lightgrey",main="Cx. pipiens vs Cx. quinquefasciatus", frame.plot=F,axes=F,box=F,legend=F)
#    plot(bivmap.trunc, col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,add=T,legend=F)
#    #plot(bivmap,col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,add=T,legend=F)
#    points(sitesSp, pch = 1)
#    
#    
#    # I think the quantiles of the col.matrix function
#    # is confusing me because the distribution of values is not even
#    # from 0-1, there are a lot of very low values
#    
#    hist(pip)
#    quantile(pip)
#    
#    hist(quinq)
#    quantile(quinq)
#    
#    # Yeah, the way this splits the distribution by quantiles makes it 
#    # difficult to understand for the bivmap. I would rather split
#    # by raw value instead of quantiles so that the scale was absolute
#    # rather than relative in the final product
#    
#    
#    #############################################
#    #create a color matrix for Bivariate map with Giuseppe's old code
#    #############################################
#    colmat<-function(nquantiles=4, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
#      my.data<-seq(0,1,.01)
#      my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
#      my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
#      my.pal.2<-findColours(my.class,c(upperright, bottomright))
#      col.matrix<-matrix(nrow = 101, ncol = 101, NA)
#      for(i in 1:101){
#        my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
#        col.matrix[102-i,]<-findColours(my.class,my.col)}
#      plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
#      for(i in 1:101){
#        col.temp<-col.matrix[i-1,]
#        points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
#      seqs<-seq(0,100,(100/nquantiles))
#      seqs[1]<-1
#      col.matrix<-col.matrix[c(seqs), c(seqs)]}
#    
#    col.matrix<-colmat(nquantiles=4, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")
#    
#    #############################################
#    #Build bivariate map function:
#    #############################################
#    bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=4){
#      
#      quanmean<-getValues(rasterx)
#      temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
#      brks<-with(temp, quantile(temp,na.rm=T, probs = c(seq(0,1,1/nquantiles))))
#      r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = FALSE))
#      quantr<-data.frame(r1[,2]) 
#      
#      quanvar<-getValues(rastery)
#      temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
#      brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
#      r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = FALSE))
#      quantr2<-data.frame(r2[,2])
#      
#      as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
#      col.matrix2<-colormatrix
#      cn<-unique(colormatrix)
#      for(i in 1:length(col.matrix2)){
#        ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
#      cols<-numeric(length(quantr[,1]))
#      for(i in 1:length(quantr[,1])){
#        a<-as.numeric.factor(quantr[i,1])
#        b<-as.numeric.factor(quantr2[i,1])
#        cols[i]<-as.numeric(col.matrix2[b,a])}
#      r<-rasterx
#      r[1:length(r)]<-cols
#      return(r)}
#    
#    #############################################
#    # Use the bivariate.map function with quinq and pip
#    #############################################
#    
#    bivmap<-bivariate.map(quinq,pip,colormatrix=col.matrix, nquantiles=4)
#    # Note: go back to the function that generates color matrices and try out new color schemes.
#    
#    #pdf(file="pipiens_vs_quinqs.pdf", width = 8, height = 9)
#    plot(bivmap,main="Cx. pipiens vs Cx. quinquefasciatus", col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,legend=F)
#    points(sitesSp, pch = 1)
#    #dev.off()
#    
#    #pdf(file="legend.pdf", width = 8, height = 9)
#    col.matrix<-colmat(nquantiles=4, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")
#    #dev.off()
#    
#    #############################################
#    # Use the bivariate.map function with quinq.trunc and pip.trunc
#    #############################################
#    bivmap.trunc<-bivariate.map(quinq.trunc,pip.trunc,colormatrix=col.matrix, nquantiles=4)
#    # Note: go back to the function that generates color matrices and try out new color schemes.
#    
#    #pdf(file="pipiens_vs_quinqs_truncated.pdf", width = 8, height = 9)
#    usa <- pip
#    usa[usa > 0] <- 1
#    plot(usa,col="grey",frame.plot=F,axes=F,box=F,legend=F)
#    plot(bivmap.trunc, main="Cx. pipiens vs Cx. quinquefasciatus", col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,add=T,legend=F)
#    points(sitesSp, pch = 1)
#    #dev.off()
#    
#    #pdf(file="legend_truncated.pdf", width = 8, height = 9)
#    col.matrix<-colmat(nquantiles=4, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")
#    #dev.off()
#    
#    
#    
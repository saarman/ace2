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
sites <- read.csv("ace2_by_site-id.csv")
sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))

# plot both distributions with sampled points
plot(pip)
points(sitesSp, pch = 1)

plot(quinq)
points(sitesSp, pch = 1)

#############################################
# Extract suitability at each location
#############################################
pip_points <- extract(pip, sitesSp)
quinq_points <- extract(quinq, sitesSp)

sites_suitability <- cbind(sites, pip_points,quinq_points)
names(sites_suitability) <- c("site","pp","pq","qq","latitude","longitude","pip_suitability","quinq_suitability")

# calculate habitat ratio and p frequency
sites_suitability$pip_rel.suitability <- sites_suitability$pip_suitability/(sites_suitability$pip_suitability+sites_suitability$quinq_suitability)
sites_suitability$pip_freq <- ((sites$pp*2)+(sites$pq))/((sites$pp*2)+(sites$pq)+(sites$qq*2))

# reorder and name rows
sites_suitability <- sites_suitability[rev(order(sites_suitability$latitude)),]
rownames(sites_suitability)<- c("01-BYR-WA",
                            "03-USU-UT",
                            "02-BAR-MA",
                            "04-SLC-UT",
                            "05-HUN-NJ",
                            "06-SOM-NJ",
                            "07-SUT-CA",
                            "09-ROC-MD",
                            "08-LIN-CA",
                            "10-MAR-AZ",
                            "11-COS-TX",
                            "12-ANA-FL",
                            "13-MIA-FL")

sites_suitability <- sites_suitability[rev(order(rownames(sites_suitability))),]

#pdf(file="barplots.pdf", width = 8, height = 9)
par(mfrow = c(1, 3), mar=c(10,10,2,2))
barplot(sites_suitability$pip_freq,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
dev.off()

#with ggplot, but i'm not sure I like how it looks
data <- as.data.frame(sites_suitability)
data$site <- rownames(sites_suitability)

p1 <- ggplot(data=data, aes(x=site, y=pip_freq)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = pip_freq), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=data$site) +
  ylab("Hybrid Index") +
  xlab("Collection Site") +
  coord_flip()

p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = pip_freq), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=data$site) +
  ylab("HSM Cx. pipiens") +
  xlab(NULL) +
  coord_flip()

p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = pip_freq), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=data$site) +
  ylab("HSM Cx. quinquefasciatus") +
  xlab(NULL) +
  coord_flip()


library(patchwork)
#pdf(file="lineplots.pdf", width = 8, height = 6)
p1 + p2 + p3
#dev.off()

#############################################
# Scatter plot of hybrid index versus suitability
#############################################

p4 <- ggplot(data=data, aes(y=pip_suitability,x=pip_freq)) + 
  geom_point() +
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  ylab("HSM Cx. pipiens") +
  xlab("Hybrid Index") +
  coord_flip()

p5 <- ggplot(data=data, aes(y=quinq_suitability,x=pip_freq)) + 
  geom_point() +
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  ylab("HSM Cx. quinquefasciatus") +
  xlab(NULL)  +
  coord_flip()

#pdf(file="scatterplots.pdf", width = 7.5, height = 4)
p4 + p5
#dev.off()

############################################
# Linear model fit
############################################

pip.lm <- lm(pip_freq~pip_suitability, data = sites_suitability)
summary(pip.lm)

#Call:
#  lm(formula = pip_freq ~ pip_suitability, data = sites_suitability)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.64121 -0.03795  0.04920  0.09788  0.19487 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -0.07059    0.12916  -0.546    0.596    
#pip_suitability  1.08937    0.17347   6.280 6.01e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.2263 on 11 degrees of freedom
#Multiple R-squared:  0.7819,	Adjusted R-squared:  0.7621 
#F-statistic: 39.43 on 1 and 11 DF,  p-value: 6.005e-05

quinq.lm <- lm(pip_freq~quinq_suitability, data = sites_suitability)
summary(quinq.lm)

#Call:
#  lm(formula = pip_freq ~ quinq_suitability, data = sites_suitability)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.67159 -0.34451 -0.01923  0.33020  0.52006 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)         1.2225     0.3476   3.517  0.00483 **
#  quinq_suitability  -0.8780     0.4913  -1.787  0.10150   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4265 on 11 degrees of freedom
#Multiple R-squared:  0.225,	Adjusted R-squared:  0.1545 
#F-statistic: 3.193 on 1 and 11 DF,  p-value: 0.1015


#############################################
# Plot maps with complementary colors:
#############################################
upperleft=rgb(0,150,235, maxColorValue=255)
upperright=rgb(130,0,80, maxColorValue=255)
bottomleft="grey"
bottomright=rgb(255,230,15, maxColorValue=255)

colfuncX <- colorRampPalette(c(bottomleft,upperleft))
#colfuncX <- colorRampPalette(c("white", "#EFCA17", "red"))
#colfuncX <- colorRampPalette(c("#FFE50E","#821250"))
colfuncY <- colorRampPalette(c(bottomleft,bottomright))
#colfuncY <- colorRampPalette(c("white","#0096EC","#701363"))
#colfuncY <- colorRampPalette(c("#75AFCE","#612576"))

pdf(file="pipiens.pdf", width = 8, height = 9)
plot(pip, main="Cx. pipiens", col=colfuncX(4),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
dev.off()

pdf(file="quinquefasciatus.pdf", width = 8, height = 9)
plot(quinq,main="Cx. quinquefasciatus", col=colfuncY(4),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
dev.off()


#############################################
# Trunkate the pip and quinq anything below 0.1 --> NA
#############################################
pip.trunk <- pip
pip.trunk[pip.trunk < 0.1] <- NA

#pdf(file="pipiens.pdf", width = 8, height = 9)
plot(pip.trunk, main="Cx. pipiens", col=colfuncX(4),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

quinq.trunk <- quinq 
quinq.trunk[quinq.trunk < 0.1] <- NA

#pdf(file="quinquefasciatus.pdf", width = 8, height = 9)
plot(quinq.trunk,main="Cx. quinquefasciatus", col=colfuncY(4),frame.plot=F,axes=F,box=F,add=F,legend=T)
points(sitesSp, pch = 1)
#dev.off()

#############################################
# Bivariate map with new example https://cran.r-project.org/web/packages/biscale/vignettes/biscale.html
#############################################





#############################################
#create a color matrix for Bivariate map with Giuseppe's old code
#############################################
colmat<-function(nquantiles=4, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
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

col.matrix<-colmat(nquantiles=4, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")

#############################################
#Build bivariate map function:
#############################################
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=4){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(temp,na.rm=T, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = FALSE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = FALSE))
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

bivmap<-bivariate.map(quinq,pip,colormatrix=col.matrix, nquantiles=4)
# Note: go back to the function that generates color matrices and try out new color schemes.

#pdf(file="pipiens_vs_quinqs.pdf", width = 8, height = 9)
plot(bivmap,main="Cx. pipiens vs Cx. quinquefasciatus", col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,legend=F)
points(sitesSp, pch = 1)
#dev.off()

#pdf(file="legend.pdf", width = 8, height = 9)
col.matrix<-colmat(nquantiles=4, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")
#dev.off()

#############################################
# Use the bivariate.map function with quinq.trunk and pip.trunk
#############################################
bivmap.trunk<-bivariate.map(quinq.trunk,pip.trunk,colormatrix=col.matrix, nquantiles=4)
# Note: go back to the function that generates color matrices and try out new color schemes.

#pdf(file="pipiens_vs_quinqs_trunkated.pdf", width = 8, height = 9)
usa <- pip
usa[usa > 0] <- 1
plot(usa,col="grey",frame.plot=F,axes=F,box=F,legend=F)
plot(bivmap.trunk, main="Cx. pipiens vs Cx. quinquefasciatus", col=as.vector(col.matrix), frame.plot=F,axes=F,box=F,add=T,legend=F)
points(sitesSp, pch = 1)
#dev.off()

#pdf(file="legend_trunkated.pdf", width = 8, height = 9)
col.matrix<-colmat(nquantiles=4, xlab="Cx. quinquefasciatus", ylab="Cx. pipiens")
#dev.off()



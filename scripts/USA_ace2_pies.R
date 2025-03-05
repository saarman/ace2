# This script starts after data wrangling, with a data input csv file.

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
library(terra)

#############################################
# Import site sampling points, create spatial
#############################################
# add sampled points
sites <- read.csv("ace2_by_site_date.csv")
sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))

#############################################
# Habitat suitability models as netCDF files
#############################################
# import netCDF files from https://github.com/lanl/culexmaxentmodels

# pip North America with chosen extent 
pip <- brick("culex_pipiens_meansuitability.nc")

# quinx, set extent to match
quinq.full <- brick("culex_quinquefasciatus_meansuitability.nc")
quinq <- setExtent(quinq.full,extent(pip),keepres = TRUE)
quinq <- resample(quinq, pip)
#identical(extent(pip),extent(quinq))

#############################################
# Plot both distributions with complementary colors:
#############################################
library(scales)  #to use alpha() function

upperleft=rgb(0,150,235, maxColorValue=255)
upperright=rgb(130,0,80, maxColorValue=255)
bottomleft="white" #rgb(220,220,220,  maxColorValue=255)
bottomright=rgb(255,230,15, maxColorValue=255)

colfuncX <- colorRampPalette(c(bottomleft,upperleft))
colfuncY <- colorRampPalette(c(bottomleft,bottomright))

# Set xpd to NA to allow for plotting in the margins
par(xpd = NA)

plot(pip, main="", col=alpha(colfuncX(10),1),frame.plot=F,axes=F,box=F,add=F,legend=F,mar = rep(10,4))
plot(quinq, col=alpha(colfuncY(10),0.75),frame.plot=F,axes=F,box=F,add=T,legend=F)
#points(sitesSp, pch = 20)

#############################################
# Extract suitability at each location, create data frame
#############################################

# start data frame and name fields
sitesDf <- as.data.frame(sites[,1:7])
names(sitesDf) <- c("locality","site","pp","pq","qq","latitude","longitude")

# extract habitat suitability value for pip and quinq layers at each point
sitesDf$pip_hs <- extract(pip, sitesSp)
sitesDf$quinq_hs <- extract(quinq, sitesSp)

# calculate + add: hybrid index,suitability index
sitesDf$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
sitesDf$h_index <- ((2*sitesDf$pp)+(sitesDf$pq))/sitesDf$alleles # p/p+q
sitesDf$s_index <- sitesDf$pip_hs/(sitesDf$pip_hs+sitesDf$quinq_hs) # pip suitability / total suitability
#plot(sitesDf$h_index,sitesDf$s_index, xlim = c(0,1), ylim = c(0,1))
#abline(0,1)

# name rows
rownames(sitesDf)<- sitesDf$site


#######################################
# Plot Pies onto habitat suitability with jitter
#######################################
# Start the PDF device driver
#pdf(file="Habitat_ace2_pies.pdf", width=11, height=8.5)

# Set xpd to NA to allow for plotting in the margins
par(xpd = NA)

plot(pip, main="", col=alpha(colfuncX(10),1),frame.plot=F,axes=F,box=F,add=F,legend=F,mar = rep(10,4))
plot(quinq, col=alpha(colfuncY(10),0.75),frame.plot=F,axes=F,box=F,add=T,legend=F)
#points(sitesSp, pch = 20)

# add from data table with x/y offset
x_offset <- c(0,0,0,5,5+(0.96788-0.82005),3,8.5,-6,0,-6-(.59845-.30633),0,-8,-2,0,0,5,5)
y_offset <- c(0,0,0,3,-.5,-6,-3,2.5,0,-2.5,0,-3,-4,0,0,1,.5)
for (i in 1:17){
  lines(x = c(sitesDf$longitude[i],sitesDf$longitude[i]+x_offset[i]), y = c(sitesDf$latitude[i],sitesDf$latitude[i]+y_offset[i]))
  add.pie(z = c(sitesDf$pp[i],sitesDf$pq[i],sitesDf$qq[i]), x = sitesDf$longitude[i]+x_offset[i], y = sitesDf$latitude[i]+y_offset[i], clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 2, radius = 2 )
}

dev.off()

#######################################
# Plot Pies onto Map with jitter
#######################################

# Start the PDF device driver
#pdf(file="USA_ace2_pies.pdf", width=11, height=8.5)

# Set xpd to NA to allow for plotting in the margins
par(xpd = NA)

# plot outline of states
map("state", col = "grey90", fill = TRUE, mar = rep(10, 4))

# add from data table with x/y offset
x_offset <- c(0,0,0,5,5+(0.96788-0.82005),3,8.5,-6,0,-6-(.59845-.30633),0,-8,-2,0,0,5,5)
y_offset <- c(0,0,0,3,-.5,-6,-3,2.5,0,-2.5,0,-3,-4,0,0,1,.5)
for (i in 1:17){
  lines(x = c(sitesDf$longitude[i],sitesDf$longitude[i]+x_offset[i]), y = c(sitesDf$latitude[i],sitesDf$latitude[i]+y_offset[i]))
  add.pie(z = c(sitesDf$pp[i],sitesDf$pq[i],sitesDf$qq[i]), x = sitesDf$longitude[i]+x_offset[i], y = sitesDf$latitude[i]+y_offset[i], clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 2, radius = 2 )
}

dev.off()

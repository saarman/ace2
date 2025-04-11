# This script starts after data wrangling, with a data input txt file.

# load libraries
library(maps)
library(mapdata)
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library(mapplots)
library(sp)
library(classInt)
library(ggplot2)
#library(fields)
#library(terra)
#library(raster)

#setwd("~/Dropbox/Mosquitoes/VecTech_DNA/OneHealth_MS/ace2_github")

#############################################
# Import site sampling points, create spatial
#############################################
# add sampled points
sites <- read.csv("./data/ThisStudy.txt", sep = "\t")
sitesSp <- SpatialPoints(coords = cbind(sites$long,sites$lat))


#############################################
# Create data frame
#############################################
colnames(sites)
# start data frame and name fields
sitesDf <- as.data.frame(sites[,c(2,1,4,5,6,9,10,7,8)])
names(sitesDf) <- c("locality","site","pp","pq","qq","latitude","longitude","year","h-index")

# name rows
rownames(sitesDf)<- sitesDf$site


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


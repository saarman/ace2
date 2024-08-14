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

#############################################
# Import site sampling points, make spatial for later
#############################################
# add sampled points
sites <- read.csv("ace2_by_site_date.csv")
sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))



#TODO: clean up!!! below this line


# plot both distributions with sampled points
plot(pip)
points(sitesSp, pch = 1)

plot(quinq)
points(sitesSp, pch = 1)


#######################################
# Plot Pies onto Map by named locality, with jittered data points, add in lines to locality
#######################################

# #create and plot coord = lat,long
length(rownames(df_id))

# plot coords onto map
map("state", col = "grey90", fill = TRUE, mar = rep(10, 4))

#points(df_id$longitude,df_id$latitude,col="black",pch = 19, cex=0.6)
# add in pie charts source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
# One at a time to allow for added lines:
#1
x_offset <- (5)
y_offset <- (2)
add.pie(z = c(df_id$pp[1],df_id$pq[1],df_id$qq[1]), x = df_id$longitude[1]+x_offset, y = df_id$latitude[1]+y_offset, clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 1, radius = 1 )
lines(x = c(df_id$longitude[1],df_id$longitude[1]+x_offset), y = c(df_id$latitude[1],df_id$latitude[1]+y_offset))
#2 
x_offset <- (5)
y_offset <- (2)
add.pie(z = c(df_id$pp[2],df_id$pq[2],df_id$qq[2]), x = df_id$longitude[2]+x_offset, y = df_id$latitude[2]+y_offset, clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 1, radius = 1 )
lines(x = c(df_id$longitude[2],df_id$longitude[2]+x_offset), y = c(df_id$latitude[2],df_id$latitude[2]+y_offset))

#etc.
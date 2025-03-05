

setwd("~/Dropbox/Caccone_Aksoy/Glossina-Spiroplasma/")
#install any packages you don't already have
# install.packages("LEA")

#load all the packages you might need:
library(LEA)
library(fields)
library(RColorBrewer)
library(mapplots)
library(maptools)
library(raster)
library(rgdal)
library(maptools)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

??pie

  #one at a time:
pie(c(2,3), main="OCA",col=rainbow(6),radius=.8,clockwise=TRUE,labels=c("Spiroplasma +", "Spiroplasma -"))


#put on the map:
dev.off()

#read in coordinate pop file which also contains the sample sizes for each pie chart
coord.pop = read.table("spiroplasma_updated.coord")
coord.pop

#read in the sample proportion spiroplasma positive and negative in two columns
qpop = read.table("spiroplasma-pies_updated.txt")
qpop = as.matrix(qpop)
qpop

#define the reaches of your map
coord = read.table("spiroplasma_updated.coord")
plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n")
plot(coord)

Uganda = read.csv("Uganda.csv")
Uganda
library(maps)
library(mapdata)
map("worldHires","Uganda")
points(Uganda$longitude, Uganda$latitude,col="red",cex=0.6)
points(coord)
#add background map of the world.
map(add = T, col = "grey90", fill = TRUE)


#add in pie charts
for (i in 1:33){
  add.pie(z = qpop[i,],coord.pop[i,3], x = coord.pop[i,1], y = coord.pop[i,2], clockwise=TRUE, labels = "", col = c("red","yellow"))
}


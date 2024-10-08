#"Ace2 PCR results, Culex pipiens complex"
library(maps)
library(mapdata)
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library(fields)
library(mapplots)


#######################################
## Tables of results:
#######################################

#ace_results <- read.csv("VecTechAce2Results_27Feb2024.csv")
ace_results <- read.csv("VecTechAce2Results_5Jun2024.csv")

#Allele frequencies for each locality
table(ace_results$Species,ace_results$Locality)

#Genotypes in p/q notation, pipiens = pp, quinquefasciatus = qq, 
# genotype df
df_geno <- data.frame(sampleID = ace_results$SampleID, 
                      site = ace_results$Locality,
                      locality = paste0(substr(ace_results$SampleID,1,3),
                      "-", substr(ace_results$SampleID,4,5)), 
                      genotype = ace_results$Species, 
                      latitude = ace_results$Latitude, 
                      longitude = ace_results$Longitude)
df_geno$genotype[df_geno$genotype == "Culex pipiens pipiens"] <- "pp"
df_geno$genotype[df_geno$genotype == "hybrid"] <- "pq"
df_geno$genotype[df_geno$genotype == "Culex quinquefasciatus"] <- "qq"

#View data frame
View(df_geno)

#create simple table by coordinates (sites)
table(df_geno[,c(2,4)])

#create simple table by named localities
table(df_geno[,c(3,4)])

#######################################
## Pie charts:
#######################################

## One pie chart at a time for coordinates (sites)
tab_sites <- table(df_geno[,c(2,4)])
dimnames(tab_sites)$site
for (i in dimnames(tab_sites)$site){
  pie(tab_sites[i,], main=i,col=c("black","grey","white"),radius=.8,clockwise=TRUE,labels = NA)
}

## One pie chart at a time for named localities
tab_id <- table(df_geno[,c(3,4)])
dimnames(tab_id)$locality
for (i in dimnames(tab_id)$locality){
  pie(tab_id[i,], main=i,col=c("black","grey","white"),radius=.8,clockwise=TRUE,labels = NA)
}

# # convert counts to frequencies
# frequencies <- tab
# 
# #read in the sample proportion of each genotype into three columns
# for (i in 1:length(tab[,1])){
#   tab[i,] <- tab[i,]/sum(tab[i,])
# }

#######################################
## Plot coordinates
#######################################

#create and plot coord = lat,long
coord <- unique(df_geno[,c(2,3,6,5)])
View(coord)

#plot on x-y
plot(coord$longitude, coord$latitude,xlab = "Longitude", ylab = "Latitude")

#plot coords onto map
map("state", col = "grey90", fill = TRUE)
points(coord$longitude, coord$latitude,col="red",cex=0.6, pch = 3)

#######################################
# Plot Pies onto Map by coordinates (sites)
#######################################

#get order straight in table
ordered <- tab_sites[coord$site,]

tab_sites[coord$site,]

#make sure coord is correct
ordered_coord <- coord[coord$site == coord$site,]

#add in pie charts source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
map("state", col = "grey90", fill = TRUE)
# add from data table
for (i in 1:length(ordered_coord$site)){
  add.pie(z = ordered[i,], x = ordered_coord$longitude[i], y = ordered_coord$latitude[i], clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 1, radius = 1 )
}
# add from Andrea's extra locality 
#add.pie(z = c(64/87,1/87,10/87), x = -78.42910058623431, y = 39.21137430019763, clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 3, radius = 3 )

#put into one dataframe
df_gps <- as.data.frame.matrix(ordered)
df_gps$latitude <- ordered_coord$latitude
df_gps$longitude <- ordered_coord$longitude

write.csv(df_gps,file = "ace2_by_gps.csv")

#######################################
# Plot Pies onto Map by named locality
#######################################
#get order straight in table
ordered <- tab_id[unique(coord$locality),]

#add gps locality for each name
ordered_coord <- as.data.frame(aggregate(latitude ~ locality, coord, mean))
ordered_coord$longitude <- as.data.frame(aggregate(longitude ~ locality, coord, mean))$longitude


#create and save plot as png
#png(file="map-pies.png", width=1200, height=700)

#add in pie charts source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
map("state", col = "grey90", fill = TRUE)
# add from data table
for (i in 1:length(ordered_coord$locality)){
  add.pie(z = ordered[i,1:3], x = ordered_coord$longitude[i], y = ordered_coord$latitude[i], clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 1, radius = 1 )
}
# add from Andrea's extra locality 
#add.pie(z = c(64/87,1/87,10/87), x = -78.42910058623431, y = 39.21137430019763, clockwise=TRUE, labels = "", col = c("black","grey","white"), cex = 3, radius = 3 )

dev.off()

#put into one dataframe
#put into one dataframe
df_id <- as.data.frame.matrix(ordered)
df_id$latitude <- ordered_coord$latitude
df_id$longitude <- ordered_coord$longitude

write.csv(df_id,file = "ace2_by_site-id.csv")

# Note, I added date to this csv file manually in excel after the fact. 
# This became "ace2_by_site_date.csv"


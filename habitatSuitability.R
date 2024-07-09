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
# Extract suitability at each location
#############################################
pip_points <- extract(pip, sitesSp)
quinq_points <- extract(quinq, sitesSp)

sites_suitability <- cbind(sites[,1:7], pip_points,quinq_points)
names(sites_suitability) <- c("locality","site","pp","pq","qq","latitude","longitude","pip_suitability","quinq_suitability")

# calculate h-index
sites_suitability$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
sites_suitability$h_index <- ((2*sites_suitability$pp)+(sites_suitability$pq))/sites_suitability$alleles

sites_suitability[c("pp","pq","qq","alleles","h_index")]

# name rows
rownames(sites_suitability)<- sites_suitability$site

#pdf(file="barplots.pdf", width = 8, height = 9)
par(mfrow = c(1, 3), mar=c(10,10,2,2))
barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
#dev.off()

#with ggplot, but i'm not sure I like how it looks
data <- as.data.frame(sites_suitability)
data$site <- rownames(sites_suitability)

p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = h_index), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=rev) +
  ylab("H-Index\nCxq ----------------> Cxp") +
  xlab("Collection Site") +
  coord_flip()

p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  geom_line() +
  geom_segment(aes(xend=site, yend=0), color="black") + 
  #geom_text(aes(label = h_index), hjust = -.5) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(limits=rev) +
  ylab("HSM-pip") +
  xlab(NULL) +
  coord_flip()

p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
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
#pdf(file="lineplots.pdf", width = 8, height = 6)
p1 + p2 + p3
#dev.off()

#############################################
# Scatter plot of hybrid index versus suitability
#############################################

p4 <- ggplot(data=data, aes(y=pip_suitability,x=h_index)) + 
  geom_point() +
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  ylab("HSM-pip") +
  xlab("H-Index\nCxq -------------------------------------------> Cxp")
  #coord_flip()

p5 <- ggplot(data=data, aes(y=quinq_suitability,x=h_index)) + 
  geom_point() +
  #geom_text(aes(label = data$site), hjust = -.25) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  ylab("HSM-quinq") +
  xlab("H-Index\nCxq -------------------------------------------> Cxp")
  #coord_flip()

#pdf(file="scatterplots.pdf", width = 7.5, height = 4)
p4 + p5
#dev.off()

data
############################################
# Linear model fit
############################################

# pip.lm <- lm(h_index~pip_suitability, data = sites_suitability)
# summary(pip.lm)
# Call:
#   lm(formula = h_index ~ pip_suitability, data = sites_suitability)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.54326 -0.05763  0.04165  0.06985  0.24012 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -0.1451     0.1027  -1.414    0.178    
# pip_suitability   1.0969     0.1450   7.567  1.7e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2172 on 15 degrees of freedom
# Multiple R-squared:  0.7924,	Adjusted R-squared:  0.7786 
# F-statistic: 57.26 on 1 and 15 DF,  p-value: 1.698e-06
# 
quinq.lm <- lm(h_index~quinq_suitability, data = sites_suitability)
summary(quinq.lm)
> # Import site sampling points, make spatial for later
  > #############################################
> # add sampled points
  > sites <- read.csv("ace2_by_site_date.csv")
> sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))
> # plot both distributions with sampled points
  > plot(pip)
> points(sitesSp, pch = 1)
> plot(quinq)
> points(sitesSp, pch = 1)
> #############################################
> # Extract suitability at each location
  > #############################################
> pip_points <- extract(pip, sitesSp)
> quinq_points <- extract(quinq, sitesSp)
> sites_suitability <- cbind(sites[,1:7], pip_points,quinq_points)
> names(sites_suitability) <- c("locality","site","pp","pq","qq","latitude","longitude","pip_suitability","quinq_suitability")
> # calculate habitat ratio and hindex
  > sites_suitability$pip_rel.suitability <- sites_suitability$pip_suitability/(sites_suitability$pip_suitability+sites_suitability$quinq_suitability)
> sites_suitability$h_index <- ((sites$pp*2)+(sites$pq))/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> sites_suitability$site
[1] "01-BYRWA" "02-COOIL" "03-BARMA" "04-USUUT" "05-SLCUT" "06-HUNNJ" "07-SOMNJ" "08-SUTCA" "09-ROCMD"
[10] "10-LINCA" "11-STGUT" "12-PHOAZ" "13MARAZ"  "14-COSTX" "15-SLILA" "16-ANAFL" "17-MIAFL"
> # name rows
  > rownames(sites_suitability)<- sites_suitability$site
> #pdf(file="barplots.pdf", width = 8, height = 9)
  > par(mfrow = c(1, 3), mar=c(10,10,2,2))
> barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
> barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
> barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
> #with ggplot, but i'm not sure I like how it looks
  > data <- as.data.frame(sites_suitability)
> data$site <- rownames(sites_suitability)
> p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("H-Index\nCxq --------------> Cxp") +
  +   xlab("Collection Site") +
  +   coord_flip()
> 
  > p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-pip") +
  +   xlab(NULL) +
  +   coord_flip()
> 
  > p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-quinq") +
  +   xlab(NULL) +
  +   coord_flip()
> 
  > 
  > library(patchwork)

Attaching package: ‘patchwork’

The following object is masked from ‘package:raster’:
  
  area

Warning message:
  package ‘patchwork’ was built under R version 4.0.5 
> #pdf(file="lineplots.pdf", width = 8, height = 6)
  > p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  > #dev.off()
  > p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("H-Index\nCxq --------------> Cxp") +
  +   xlab("Collection Site") +
  +   #coord_flip()
  + 
  + p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-pip") +
  +   xlab(NULL) +
  +   #coord_flip()
  + 
  + p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-quinq") +
  +   xlab(NULL) +
  +   #coord_flip()
  + 
  + 
  + library(patchwork)
Error in `ggplot_add()`:
  ! Can't add `library(patchwork)` to a <ggplot> object.
Run `rlang::last_trace()` to see where the error occurred.
> #pdf(file="lineplots.pdf", width = 8, height = 6)
> p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
> #dev.off()
> sites_suitability
                                            locality     site pp pq qq latitude  longitude pip_suitability
01-BYRWA                  Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900      0.96894002
02-COOIL                       Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087      0.98287416
03-BARMA                       Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427      0.98026067
04-USUUT                       Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005      0.82503319
05-SLCUT           Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788      0.96941930
06-HUNNJ                               Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462      0.76642287
07-SOMNJ                                Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610      0.89291251
08-SUTCA                             Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845      0.84323812
09-ROCMD                               Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055      0.93759370
10-LINCA                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633      0.80387235
11-STGUT St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000      0.53358889
12-PHOAZ                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414      0.22717942
13MARAZ                                 Maricopa, AZ  13MARAZ  0  0 16 33.07362 -111.97377      0.12688294
14-COSTX                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893      0.05475381
15-SLILA                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900      0.12959486
16-ANAFL                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074      0.09963039
17-MIAFL                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909      0.19097427
         quinq_suitability pip_rel.suitability    h_index
01-BYRWA         0.5034906          0.65805478 1.00000000
02-COOIL         0.9749584          0.50202157 1.00000000
03-BARMA         0.6294978          0.60894890 1.00000000
04-USUUT         0.2072268          0.79924938 1.00000000
05-SLCUT         0.9783818          0.49769933 0.51724138
06-HUNNJ         0.2833190          0.73010601 0.94444444
07-SOMNJ         0.4881530          0.64653883 1.00000000
08-SUTCA         0.8341666          0.50270402 0.92857143
09-ROCMD         0.8457388          0.52575371 1.00000000
10-LINCA         0.7629824          0.51304843 1.00000000
11-STGUT         0.5539052          0.49065912 0.00000000
12-PHOAZ         0.9860736          0.18724819 0.00000000
13MARAZ          0.8615471          0.12836816 0.00000000
14-COSTX         0.5350969          0.09282655 0.08108108
15-SLILA         0.6448851          0.16733147 0.02702703
16-ANAFL         0.7196956          0.12160043 0.00000000
17-MIAFL         1.0000000          0.16035129 0.00000000
> p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
+   geom_line() +
+   geom_segment(aes(xend=site, yend=0), color="black") + 
+   #geom_text(aes(label = h_index), hjust = -.5) +
+   geom_point(size=3) +
+   scale_y_continuous(limits=c(0,1)) +
+   scale_x_discrete(limits=data$site) +
+   ylab("H-Index\nCxq --------------> Cxp") +
+   xlab("Collection Site") +
+   coord_flip()
> #pdf(file="lineplots.pdf", width = 8, height = 6)
> p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
> p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
+   geom_line() +
+   geom_segment(aes(xend=site, yend=0), color="black") + 
+   #geom_text(aes(label = h_index), hjust = -.5) +
+   geom_point(size=3) +
+   scale_y_continuous(limits=c(0,1)) +
+   scale_x_discrete(limits=rev) +
+   ylab("H-Index\nCxq ----------------> Cxp") +
+   xlab("Collection Site") +
+   #coord_flip()
+ 
+ p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
+   geom_line() +
+   geom_segment(aes(xend=site, yend=0), color="black") + 
+   #geom_text(aes(label = h_index), hjust = -.5) +
+   geom_point(size=3) +
+   scale_y_continuous(limits=c(0,1)) +
+   scale_x_discrete(limits=data$site) +
+   ylab("HSM-pip") +
+   xlab(NULL) +
+   #coord_flip()
+ 
+ p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
+   geom_line() +
+   geom_segment(aes(xend=site, yend=0), color="black") + 
+   #geom_text(aes(label = h_index), hjust = -.5) +
+   geom_point(size=3) +
+   scale_y_continuous(limits=c(0,1)) +
+   scale_x_discrete(limits=data$site) +
+   ylab("HSM-quinq") +
+   xlab(NULL) +
+   #coord_flip()
+ 
+ 
+ library(patchwork)
Error in `ggplot_add()`:
! Can't add `library(patchwork)` to a <ggplot> object.
Run `rlang::last_trace()` to see where the error occurred.
> p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("H-Index\nCxq ----------------> Cxp") +
  +   xlab("Collection Site")
> p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-pip") +
  +   xlab(NULL) 
> p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-quinq") +
  +   xlab(NULL) 
> #pdf(file="lineplots.pdf", width = 8, height = 6)
  > p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  > p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("H-Index\nCxq ----------------> Cxp") +
  +   xlab("Collection Site") +
  +   coord_flip()
> p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-pip") +
  +   xlab(NULL) +
  +   coord_flip()
> p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=data$site) +
  +   ylab("HSM-quinq") +
  +   xlab(NULL) +
  +   coord_flip()
> #pdf(file="lineplots.pdf", width = 8, height = 6)
  > p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  > p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("H-Index\nCxq ----------------> Cxp") +
  +   xlab("Collection Site") +
  +   coord_flip()
> 
  > p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("HSM-pip") +
  +   xlab(NULL) +
  +   coord_flip()
> 
  > p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("HSM-quinq") +
  +   xlab(NULL) +
  +   coord_flip()
> #pdf(file="lineplots.pdf", width = 8, height = 6)
  > p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  > p4 <- ggplot(data=data, aes(y=pip_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM Cx. pipiens") +
  +   xlab("H-Index\nCxq ----------------> Cxp")
> p5 <- ggplot(data=data, aes(y=quinq_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM Cx. quinquefasciatus") +
  +   xlab("H-Index \n 0 = Cx. quinq.                                  1 = Cx. pip.")
> p4 <- ggplot(data=data, aes(y=pip_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM Cx. pipiens") +
  +   xlab("H-Index\nCxq -------------------------------------------> Cxp")
>   #coord_flip()
  > 
  > p5 <- ggplot(data=data, aes(y=quinq_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM Cx. quinquefasciatus") +
  +   xlab("H-Index\nCxq -------------------------------------------> Cxp")
>   #coord_flip()
  > 
  > #pdf(file="scatterplots.pdf", width = 7.5, height = 4)
  > p4 + p5
> p4 <- ggplot(data=data, aes(y=pip_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM-pip") +
  +   xlab("H-Index\nCxq -------------------------------------------> Cxp")
>   #coord_flip()
  > 
  > p5 <- ggplot(data=data, aes(y=quinq_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM-quinq") +
  +   xlab("H-Index\nCxq -------------------------------------------> Cxp")
>   #coord_flip()
  > 
  > #pdf(file="scatterplots.pdf", width = 7.5, height = 4)
  > p4 + p5
> #dev.off()
  > data
locality     site pp pq qq latitude  longitude pip_suitability
01-BYRWA                  Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900      0.96894002
02-COOIL                       Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087      0.98287416
03-BARMA                       Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427      0.98026067
04-USUUT                       Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005      0.82503319
05-SLCUT           Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788      0.96941930
06-HUNNJ                               Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462      0.76642287
07-SOMNJ                                Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610      0.89291251
08-SUTCA                             Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845      0.84323812
09-ROCMD                               Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055      0.93759370
10-LINCA                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633      0.80387235
11-STGUT St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000      0.53358889
12-PHOAZ                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414      0.22717942
13MARAZ                                 Maricopa, AZ  13MARAZ  0  0 16 33.07362 -111.97377      0.12688294
14-COSTX                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893      0.05475381
15-SLILA                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900      0.12959486
16-ANAFL                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074      0.09963039
17-MIAFL                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909      0.19097427
quinq_suitability pip_rel.suitability    h_index
01-BYRWA         0.5034906          0.65805478 1.00000000
02-COOIL         0.9749584          0.50202157 1.00000000
03-BARMA         0.6294978          0.60894890 1.00000000
04-USUUT         0.2072268          0.79924938 1.00000000
05-SLCUT         0.9783818          0.49769933 0.51724138
06-HUNNJ         0.2833190          0.73010601 0.94444444
07-SOMNJ         0.4881530          0.64653883 1.00000000
08-SUTCA         0.8341666          0.50270402 0.92857143
09-ROCMD         0.8457388          0.52575371 1.00000000
10-LINCA         0.7629824          0.51304843 1.00000000
11-STGUT         0.5539052          0.49065912 0.00000000
12-PHOAZ         0.9860736          0.18724819 0.00000000
13MARAZ          0.8615471          0.12836816 0.00000000
14-COSTX         0.5350969          0.09282655 0.08108108
15-SLILA         0.6448851          0.16733147 0.02702703
16-ANAFL         0.7196956          0.12160043 0.00000000
17-MIAFL         1.0000000          0.16035129 0.00000000
> sites$pq
[1]  2  0  0  0 11  2  3  8  3  2  0  0  0  3  1  0  0
> sites
locality     site pp pq qq latitude  longitude      collection.aprox
1                   Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900               9/26/23
2                        Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087        7/31/23-8/1/23
3                        Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427               8/18/23
4                        Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005               7/25/23
5            Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788               9/19/23
6                                Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462               9/12/23
7                                 Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610               9/20/23
8                              Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845              10/11/23
9                                Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055                6/6/23
10                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633               9/18/23
11 St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000              10/10/23
12                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414               2/24/24
13                                Maricopa, AZ  13MARAZ  0  0 16 33.07362 -111.97377              10/22/23
14                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893              11/11/23
15                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900               4/30/24
16                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074       10/2023-11/2023
17                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909 11/28/2023-11/29/2023
collection.year collection.month
1             2023              Sep
2             2023          Jul-Aug
3             2023              Aug
4             2023              Sep
5             2023              Sep
6             2023              Sep
7             2023              Sep
8             2023              Sep
9             2023              Jun
10            2023              Sep
11            2023              Oct
12            2024              Feb
13            2023              Oct
14            2023              Nov
15            2024              Apr
16            2023          Oct-Nov
17            2023              Nov
> #############################################
> # Import site sampling points, make spatial for later
  > #############################################
> # add sampled points
  > sites <- read.csv("ace2_by_site_date.csv")
> sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))
> # plot both distributions with sampled points
  > plot(pip)
> points(sitesSp, pch = 1)
> plot(quinq)
> points(sitesSp, pch = 1)
> #############################################
> # Extract suitability at each location
  > #############################################
> pip_points <- extract(pip, sitesSp)
> quinq_points <- extract(quinq, sitesSp)
> sites_suitability <- cbind(sites[,1:7], pip_points,quinq_points)
> names(sites_suitability) <- c("locality","site","pp","pq","qq","latitude","longitude","pip_suitability","quinq_suitability")
> #############################################
> # Import site sampling points, make spatial for later
  > #############################################
> # add sampled points
  > sites <- read.csv("ace2_by_site_date.csv")
> sites
locality     site pp pq qq latitude  longitude      collection.aprox
1                   Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900               9/26/23
2                        Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087        7/31/23-8/1/23
3                        Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427               8/18/23
4                        Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005               7/25/23
5            Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788               9/19/23
6                                Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462               9/12/23
7                                 Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610               9/20/23
8                              Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845              10/11/23
9                                Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055                6/6/23
10                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633               9/18/23
11 St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000              10/10/23
12                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414               2/24/24
13                                Maricopa, AZ 13-MARAZ  0  0 16 33.07362 -111.97377              10/22/23
14                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893              11/11/23
15                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900               4/30/24
16                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074       10/2023-11/2023
17                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909 11/28/2023-11/29/2023
collection.year collection.month
1             2023              Sep
2             2023          Jul-Aug
3             2023              Aug
4             2023              Sep
5             2023              Sep
6             2023              Sep
7             2023              Sep
8             2023              Sep
9             2023              Jun
10            2023              Sep
11            2023              Oct
12            2024              Feb
13            2023              Oct
14            2023              Nov
15            2024              Apr
16            2023          Oct-Nov
17            2023              Nov
> sites <- read.csv("ace2_by_site_date.csv")
> sitesSp <- SpatialPoints(coords = cbind(sites$longitude,sites$latitude))
> 
  > # plot both distributions with sampled points
  > plot(pip)
> points(sitesSp, pch = 1)
> 
  > plot(quinq)
> points(sitesSp, pch = 1)
> 
  > #############################################
> # Extract suitability at each location
  > #############################################
> pip_points <- extract(pip, sitesSp)
> quinq_points <- extract(quinq, sitesSp)
> 
  > sites_suitability <- cbind(sites[,1:7], pip_points,quinq_points)
> names(sites_suitability) <- c("locality","site","pp","pq","qq","latitude","longitude","pip_suitability","quinq_suitability")
> 
  > # calculate h-index
  > sites_suitability$h_index <- ((sites$pp*2)+(sites$pq))/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> sites_suitability
locality     site pp pq qq latitude  longitude pip_suitability
1                   Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900      0.96894002
2                        Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087      0.98287416
3                        Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427      0.98026067
4                        Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005      0.82503319
5            Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788      0.96941930
6                                Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462      0.76642287
7                                 Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610      0.89291251
8                              Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845      0.84323812
9                                Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055      0.93759370
10                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633      0.80387235
11 St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000      0.53358889
12                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414      0.22717942
13                                Maricopa, AZ 13-MARAZ  0  0 16 33.07362 -111.97377      0.12688294
14                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893      0.05475381
15                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900      0.12959486
16                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074      0.09963039
17                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909      0.19097427
quinq_suitability    h_index
1          0.5034906 1.00000000
2          0.9749584 1.00000000
3          0.6294978 1.00000000
4          0.2072268 1.00000000
5          0.9783818 0.51724138
6          0.2833190 0.94444444
7          0.4881530 1.00000000
8          0.8341666 0.92857143
9          0.8457388 1.00000000
10         0.7629824 1.00000000
11         0.5539052 0.00000000
12         0.9860736 0.00000000
13         0.8615471 0.00000000
14         0.5350969 0.08108108
15         0.6448851 0.02702703
16         0.7196956 0.00000000
17         1.0000000 0.00000000
> sites_suitability
locality     site pp pq qq latitude  longitude pip_suitability quinq_suitability    h_index
1                   Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900      0.96894002         0.5034906 1.00000000
2                        Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087      0.98287416         0.9749584 1.00000000
3                        Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427      0.98026067         0.6294978 1.00000000
4                        Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005      0.82503319         0.2072268 1.00000000
5            Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788      0.96941930         0.9783818 0.51724138
6                                Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462      0.76642287         0.2833190 0.94444444
7                                 Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610      0.89291251         0.4881530 1.00000000
8                              Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845      0.84323812         0.8341666 0.92857143
9                                Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055      0.93759370         0.8457388 1.00000000
10                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633      0.80387235         0.7629824 1.00000000
11 St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000      0.53358889         0.5539052 0.00000000
12                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414      0.22717942         0.9860736 0.00000000
13                                Maricopa, AZ 13-MARAZ  0  0 16 33.07362 -111.97377      0.12688294         0.8615471 0.00000000
14                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893      0.05475381         0.5350969 0.08108108
15                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900      0.12959486         0.6448851 0.02702703
16                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074      0.09963039         0.7196956 0.00000000
17                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909      0.19097427         1.0000000 0.00000000
> # calculate h-index
  > sites_suitability$h_index <- ((sites$pp*2)-(sites$pq))/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> # name rows
  > rownames(sites_suitability)<- sites_suitability$site
> #pdf(file="barplots.pdf", width = 8, height = 9)
  > par(mfrow = c(1, 3), mar=c(10,10,2,2))
> barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
Error in plot.new() : figure margins too large
> barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
Error in plot.new() : figure margins too large
> barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
Error in plot.new() : figure margins too large
> par(mfrow = c(1, 3), mar=c(10,10,2,2))
Warning messages:
  1: In doTryCatch(return(expr), name, parentenv, handler) :
  display list redraw incomplete
2: In doTryCatch(return(expr), name, parentenv, handler) :
  display list redraw incomplete
> barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
> barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
> barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
> # calculate h-index
  > sites_suitability$h_index <- (absolute((sites$pp*2)-(sites$pq)-(sites$qq*2)))/((sites$pp*2)+(sites$pq)+(sites$qq*2))
Error in absolute((sites$pp * 2) - (sites$pq) - (sites$qq * 2)) : 
  could not find function "absolute"
> # calculate h-index
  > sites_suitability$h_index <- (abs((sites$pp*2)-(sites$pq)-(sites$qq*2)))/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> # name rows
  > rownames(sites_suitability)<- sites_suitability$site
> #pdf(file="barplots.pdf", width = 8, height = 9)
  > par(mfrow = c(1, 3), mar=c(10,10,2,2))
> barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
> barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
> # calculate h-index
  > sites_suitability$h_index <- (abs((sites$pp*2)-(sites$pq)))/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> # name rows
  > rownames(sites_suitability)<- sites_suitability$site
> #pdf(file="barplots.pdf", width = 8, height = 9)
  > par(mfrow = c(1, 3), mar=c(10,10,2,2))
> barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
> barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
> barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
> sites_suitability
locality     site pp pq qq latitude  longitude pip_suitability quinq_suitability
01-BYRWA                  Byron, WA (Byron Pond mix) 01-BYRWA 37  2  0 46.19300 -119.89900      0.96894002         0.5034906
02-COOIL                       Cook County, IL (mix) 02-COOIL 17  0  0 42.03176  -87.93087      0.98287416         0.9749584
03-BARMA                       Barnstable County, MA 03-BARMA 12  0  0 41.79362  -69.99427      0.98026067         0.6294978
04-USUUT                       Cache, UT (Hyde Park) 04-USUUT 20  0  0 41.79696 -111.82005      0.82503319         0.2072268
05-SLCUT           Salt Lake City, UT (mix of sites) 05-SLCUT  2 11  7 40.74805 -111.96788      0.96941930         0.9783818
06-HUNNJ                               Hunterdon, NJ 06-HUNNJ 16  2  1 40.53959  -74.83462      0.76642287         0.2833190
07-SOMNJ                                Somerset, NJ 07-SOMNJ 12  3  0 40.53358  -74.58610      0.89291251         0.4881530
08-SUTCA                             Sutter-Yuba, CA 08-SUTCA  9  8  1 39.16647 -121.59845      0.84323812         0.8341666
09-ROCMD                               Rockville, MD 09-ROCMD 17  3  0 39.05775  -77.13055      0.93759370         0.8457388
10-LINCA                       Lincoln, CA (1-2 mix) 10-LINCA 38  2  0 38.90447 -121.30633      0.80387235         0.7629824
11-STGUT St. George, Washington County, UT (1-3 mix) 11-STGUT  0  0 20 37.17900 -113.32000      0.53358889         0.5539052
12-PHOAZ                Phoenix, Maricopa County, AZ 12-PHOAZ  0  0 20 33.44844 -112.07414      0.22717942         0.9860736
13-MARAZ                                Maricopa, AZ 13-MARAZ  0  0 16 33.07362 -111.97377      0.12688294         0.8615471
14-COSTX                         College Station, TX 14-COSTX  0  3 17 30.60044  -96.26893      0.05475381         0.5350969
15-SLILA                                 Slidell, LA 15-SLILA  0  1 18 30.32700  -89.74900      0.12959486         0.6448851
16-ANAFL                               Anastasia, FL 16-ANAFL  0  0 19 29.90311  -81.40074      0.09963039         0.7196956
17-MIAFL                 Miami-Dade County, FL (mix) 17-MIAFL  0  0 19 25.80141  -80.19909      0.19097427         1.0000000
h_index
01-BYRWA 0.94736842
02-COOIL 1.00000000
03-BARMA 1.00000000
04-USUUT 1.00000000
05-SLCUT 0.24137931
06-HUNNJ 0.83333333
07-SOMNJ 0.77777778
08-SUTCA 0.35714286
09-ROCMD 0.83783784
10-LINCA 0.94871795
11-STGUT 0.00000000
12-PHOAZ 0.00000000
13-MARAZ 0.00000000
14-COSTX 0.08108108
15-SLILA 0.02702703
16-ANAFL 0.00000000
17-MIAFL 0.00000000
> sites_suitability[c("pp","pq","qq","h_index")]
pp pq qq    h_index
01-BYRWA 37  2  0 0.94736842
02-COOIL 17  0  0 1.00000000
03-BARMA 12  0  0 1.00000000
04-USUUT 20  0  0 1.00000000
05-SLCUT  2 11  7 0.24137931
06-HUNNJ 16  2  1 0.83333333
07-SOMNJ 12  3  0 0.77777778
08-SUTCA  9  8  1 0.35714286
09-ROCMD 17  3  0 0.83783784
10-LINCA 38  2  0 0.94871795
11-STGUT  0  0 20 0.00000000
12-PHOAZ  0  0 20 0.00000000
13-MARAZ  0  0 16 0.00000000
14-COSTX  0  3 17 0.08108108
15-SLILA  0  1 18 0.02702703
16-ANAFL  0  0 19 0.00000000
17-MIAFL  0  0 19 0.00000000
> # calculate h-index
  > sites_suitability$N <- (sites$pp*2)+(sites$pq)+(sites$qq*2)
> sites_suitability$h_index <- ((sites$pp*2)+(sites$pq)+(sites$qq*2))-(sites$pp*2)-sites$pq
> sites_suitability[c("pp","pq","qq","N","h_index")]
pp pq qq  N h_index
01-BYRWA 37  2  0 76       0
02-COOIL 17  0  0 34       0
03-BARMA 12  0  0 24       0
04-USUUT 20  0  0 40       0
05-SLCUT  2 11  7 29      14
06-HUNNJ 16  2  1 36       2
07-SOMNJ 12  3  0 27       0
08-SUTCA  9  8  1 28       2
09-ROCMD 17  3  0 37       0
10-LINCA 38  2  0 78       0
11-STGUT  0  0 20 40      40
12-PHOAZ  0  0 20 40      40
13-MARAZ  0  0 16 32      32
14-COSTX  0  3 17 37      34
15-SLILA  0  1 18 37      36
16-ANAFL  0  0 19 38      38
17-MIAFL  0  0 19 38      38
> # calculate h-index
  > sites_suitability$N <- (sites$pp*2)+(sites$pq)+(sites$qq*2)
> sites_suitability$h_index <- (((sites$pp*2)+(sites$pq)+(sites$qq*2))-(sites$pp*2)-sites$pq)/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> sites_suitability[c("pp","pq","qq","N","h_index")]
pp pq qq  N    h_index
01-BYRWA 37  2  0 76 0.00000000
02-COOIL 17  0  0 34 0.00000000
03-BARMA 12  0  0 24 0.00000000
04-USUUT 20  0  0 40 0.00000000
05-SLCUT  2 11  7 29 0.48275862
06-HUNNJ 16  2  1 36 0.05555556
07-SOMNJ 12  3  0 27 0.00000000
08-SUTCA  9  8  1 28 0.07142857
09-ROCMD 17  3  0 37 0.00000000
10-LINCA 38  2  0 78 0.00000000
11-STGUT  0  0 20 40 1.00000000
12-PHOAZ  0  0 20 40 1.00000000
13-MARAZ  0  0 16 32 1.00000000
14-COSTX  0  3 17 37 0.91891892
15-SLILA  0  1 18 37 0.97297297
16-ANAFL  0  0 19 38 1.00000000
17-MIAFL  0  0 19 38 1.00000000
> sites_suitability$h_index <- (((sites$pp*2)+(sites$pq)+(sites$qq*2))-(sites$qq*2)-sites$pq)/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> sites_suitability[c("pp","pq","qq","N","h_index")]
pp pq qq  N   h_index
01-BYRWA 37  2  0 76 0.9736842
02-COOIL 17  0  0 34 1.0000000
03-BARMA 12  0  0 24 1.0000000
04-USUUT 20  0  0 40 1.0000000
05-SLCUT  2 11  7 29 0.1379310
06-HUNNJ 16  2  1 36 0.8888889
07-SOMNJ 12  3  0 27 0.8888889
08-SUTCA  9  8  1 28 0.6428571
09-ROCMD 17  3  0 37 0.9189189
10-LINCA 38  2  0 78 0.9743590
11-STGUT  0  0 20 40 0.0000000
12-PHOAZ  0  0 20 40 0.0000000
13-MARAZ  0  0 16 32 0.0000000
14-COSTX  0  3 17 37 0.0000000
15-SLILA  0  1 18 37 0.0000000
16-ANAFL  0  0 19 38 0.0000000
17-MIAFL  0  0 19 38 0.0000000
> # calculate h-index
  > sites_suitability$alleles <- (sites$pp*2)+(sites$pq)+(sites$qq*2)
> sites_suitability$h_index <- (((sites$pp*2)+(sites$pq)+(sites$qq*2))-(sites$qq*2)-sites$pq)/((sites$pp*2)+(sites$pq)+(sites$qq*2))
> 
  > sites_suitability[c("pp","pq","qq","alleles","h_index")]
pp pq qq alleles   h_index
01-BYRWA 37  2  0      76 0.9736842
02-COOIL 17  0  0      34 1.0000000
03-BARMA 12  0  0      24 1.0000000
04-USUUT 20  0  0      40 1.0000000
05-SLCUT  2 11  7      29 0.1379310
06-HUNNJ 16  2  1      36 0.8888889
07-SOMNJ 12  3  0      27 0.8888889
08-SUTCA  9  8  1      28 0.6428571
09-ROCMD 17  3  0      37 0.9189189
10-LINCA 38  2  0      78 0.9743590
11-STGUT  0  0 20      40 0.0000000
12-PHOAZ  0  0 20      40 0.0000000
13-MARAZ  0  0 16      32 0.0000000
14-COSTX  0  3 17      37 0.0000000
15-SLILA  0  1 18      37 0.0000000
16-ANAFL  0  0 19      38 0.0000000
17-MIAFL  0  0 19      38 0.0000000
> 3/36
[1] 0.08333333
> 1-(3/36)
[1] 0.9166667
> sites_suitability$alleles <- sites$alleles
> # calculate h-index
  > sites$alleles <- (sites$pp*2)+(sites$pq)+(sites$qq*2)
> sites_suitability$alleles <- sites$alleles
> # calculate h-index
  > sites$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
> sites_suitability$N_alleles <- sites$alleles
> sites_suitability$h_index <- (sites$alleles-(2*sites$pp)-sites$pq)/sites$alleles
> sites_suitability[c("pp","pq","qq","N_alleles","h_index")]
pp pq qq N_alleles    h_index
01-BYRWA 37  2  0        78 0.02564103
02-COOIL 17  0  0        34 0.00000000
03-BARMA 12  0  0        24 0.00000000
04-USUUT 20  0  0        40 0.00000000
05-SLCUT  2 11  7        40 0.62500000
06-HUNNJ 16  2  1        38 0.10526316
07-SOMNJ 12  3  0        30 0.10000000
08-SUTCA  9  8  1        36 0.27777778
09-ROCMD 17  3  0        40 0.07500000
10-LINCA 38  2  0        80 0.02500000
11-STGUT  0  0 20        40 1.00000000
12-PHOAZ  0  0 20        40 1.00000000
13-MARAZ  0  0 16        32 1.00000000
14-COSTX  0  3 17        40 0.92500000
15-SLILA  0  1 18        38 0.97368421
16-ANAFL  0  0 19        38 1.00000000
17-MIAFL  0  0 19        38 1.00000000
> 25/40
[1] 0.625
> 15/40
[1] 0.375
> 3/40
[1] 0.075
> # calculate h-index
  > sites$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
> sites_suitability <- cbind(sites[,1:7], pip_points,quinq_points)
> names(sites_suitability) <- c("locality","site","pp","pq","qq","latitude","longitude","pip_suitability","quinq_suitability")
> # calculate h-index
  > sites$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
> sites_suitability$h_index <- ((2*sites$pp)+(sites$pq))/sites$alleles
> # calculate h-index
  > sites_suitability$alleles <- (sites$pp*2)+(sites$pq*2)+(sites$qq*2)
> sites_suitability$h_index <- ((2*sites_suitability$pp)+(sites_suitability$pq))/sites_suitability$alleles
> sites_suitability[c("pp","pq","qq","alleles","h_index")]
pp pq qq alleles    h_index
1  37  2  0      78 0.97435897
2  17  0  0      34 1.00000000
3  12  0  0      24 1.00000000
4  20  0  0      40 1.00000000
5   2 11  7      40 0.37500000
6  16  2  1      38 0.89473684
7  12  3  0      30 0.90000000
8   9  8  1      36 0.72222222
9  17  3  0      40 0.92500000
10 38  2  0      80 0.97500000
11  0  0 20      40 0.00000000
12  0  0 20      40 0.00000000
13  0  0 16      32 0.00000000
14  0  3 17      40 0.07500000
15  0  1 18      38 0.02631579
16  0  0 19      38 0.00000000
17  0  0 19      38 0.00000000
> 38/78
[1] 0.4871795
> # name rows
  > rownames(sites_suitability)<- sites_suitability$site
> #pdf(file="barplots.pdf", width = 8, height = 9)
  > par(mfrow = c(1, 3), mar=c(10,10,2,2))
> barplot(sites_suitability$h_index,names.arg = rownames(sites_suitability),horiz = T,las=2,xlab="H-Index")
> barplot(sites_suitability$pip_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab = "HSM-pip")
> barplot(sites_suitability$quinq_suitability,names.arg = rownames(sites_suitability),horiz = T,las = 2, xlab="HSM-quinq")
> #with ggplot, but i'm not sure I like how it looks
  > data <- as.data.frame(sites_suitability)
> data$site <- rownames(sites_suitability)
> p1 <- ggplot(data=data, aes(x=site, y=h_index)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("H-Index\nCxq ----------------> Cxp") +
  +   xlab("Collection Site") +
  +   coord_flip()
> p2 <- ggplot(data=data, aes(x=site, y=pip_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("HSM-pip") +
  +   xlab(NULL) +
  +   coord_flip()
> p3 <- ggplot(data=data, aes(x=site, y=quinq_suitability)) + 
  +   geom_line() +
  +   geom_segment(aes(xend=site, yend=0), color="black") + 
  +   #geom_text(aes(label = h_index), hjust = -.5) +
  +   geom_point(size=3) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   scale_x_discrete(limits=rev) +
  +   ylab("HSM-quinq") +
  +   xlab(NULL) +
  +   coord_flip()
> library(patchwork)
> #pdf(file="lineplots.pdf", width = 8, height = 6)
  > p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  > pdf(file="lineplots.pdf", width = 8, height = 6)
> p1 + p2 + p3
`geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  `geom_line()`: Each group consists of only one observation.
ℹ Do you need to adjust the group aesthetic?
  > dev.off()
RStudioGD 
2 
> p4 <- ggplot(data=data, aes(y=pip_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM-pip") +
  +   xlab("H-Index\nCxq -------------------------------------------> Cxp")
> p5 <- ggplot(data=data, aes(y=quinq_suitability,x=h_index)) + 
  +   geom_point() +
  +   #geom_text(aes(label = data$site), hjust = -.25) +
  +   geom_point(size=3) +
  +   scale_x_continuous(limits=c(0,1)) +
  +   scale_y_continuous(limits=c(0,1)) +
  +   ylab("HSM-quinq") +
  +   xlab("H-Index\nCxq -------------------------------------------> Cxp")
> p4 + p5
> pdf(file="scatterplots.pdf", width = 7.5, height = 4)
> p4 + p5
> dev.off()
RStudioGD 
2 
> pip.lm <- lm(h_index~pip_suitability, data = sites_suitability)
> summary(pip.lm)

Call:
  lm(formula = h_index ~ pip_suitability, data = sites_suitability)

Residuals:
  Min       1Q   Median       3Q      Max 
-0.54326 -0.05763  0.04165  0.06985  0.24012 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)      -0.1451     0.1027  -1.414    0.178    
pip_suitability   1.0969     0.1450   7.567  1.7e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2172 on 15 degrees of freedom
Multiple R-squared:  0.7924,	Adjusted R-squared:  0.7786 
F-statistic: 57.26 on 1 and 15 DF,  p-value: 1.698e-06

> #   Estimate Std. Error t value Pr(>|t|)    
  > # (Intercept)      -0.1451     0.1027  -1.414    0.178    
  > # pip_suitability   1.0969     0.1450   7.567  1.7e-06 ***
  > #   ---
  > #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  > # 
  > # Residual standard error: 0.2172 on 15 degrees of freedom
  > # Multiple R-squared:  0.7924,	Adjusted R-squared:  0.7786 
  > # F-statistic: 57.26 on 1 and 15 DF,  p-value: 1.698e-06

quinq.lm <- lm(h_index~quinq_suitability, data = sites_suitability)
summary(quinq.lm)

#Call:
#  lm(formula = h_index ~ quinq_suitability, data = sites_suitability)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.6123 -0.4142  0.1082  0.3296  0.6589 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)         0.9690     0.3374   2.872   0.0116 *
#  quinq_suitability  -0.6440     0.4598  -1.400   0.1817  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4483 on 15 degrees of freedom
#Multiple R-squared:  0.1156,	Adjusted R-squared:  0.05668 
#F-statistic: 1.961 on 1 and 15 DF,  p-value: 0.1817

pip.quinq.lm <- lm(h_index~pip_suitability*quinq_suitability, data = sites_suitability)
summary(pip.quinq.lm)

# Call:
#   lm(formula = h_index ~ pip_suitability * quinq_suitability, data = sites_suitability)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.48189 -0.02010  0.02957  0.04941  0.27676 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                       -0.06076    0.43667  -0.139   0.8915  
# pip_suitability                    1.35361    0.57058   2.372   0.0338 *
#   quinq_suitability                 -0.11288    0.56034  -0.201   0.8435  
# pip_suitability:quinq_suitability -0.39619    0.73972  -0.536   0.6013  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2065 on 13 degrees of freedom
# Multiple R-squared:  0.8374,	Adjusted R-squared:  0.7999 
# F-statistic: 22.32 on 3 and 13 DF,  p-value: 2.099e-05

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
#pdf(file="pip_quinq_overlap.pdf", width = 8, height = 9)
plot(pip, main="Cx. pipiens vs Cx. quinquefasciatus ", col=alpha(colfuncX(10),0.5),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(quinq,main="Cx. quinquefasciatus", col=alpha(colfuncY(10),0.5),frame.plot=F,axes=F,box=F,add=T,legend=F)
points(sitesSp, pch = 1)
#dev.off()






















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
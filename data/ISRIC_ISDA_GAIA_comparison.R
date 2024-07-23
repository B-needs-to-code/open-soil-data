wd <- "C:/Users/bolen/OneDrive/Desktop/Academia/Code_Github/open-soil-data/data" #BMO
setwd(wd)

library(terra)
library(sp)
library(raster)
library(ggplot2)
library(remotes)
library(geodata)     #Extract, add to gaia
Gaia <- read.csv("gaia_trials.csv")
afph <- soil_af_isda_vsi("ph") #path=tempdir(), quiet=TRUE) are not working in fx
afph
xy = cbind(-115, 40)         # section sets projection, don't understand syntax 
v = vect(xy, crs=crs(rast()))    #where/what crs?! coordinate reference system?
v
vv = project(v, afph)
vv
coords = Gaia[,c("lng","lat")] #works
e = extract(afph, coords, raw=FALSE)   #creates values, but doesnt include coords, are they right?
colnames(e) <- c("IDcheck","ph_h2o_0_20cm_ISDA")

Gaia= cbind.data.frame(Gaia, 
                       IDcheck=e['IDcheck'],
                       ph_h2o_0_20cm_ISDA=e['ph_h2o_0_20cm_ISDA'] )
head(e)
Gaia_fid_uni = remove.duplicates(Gaia)

library(tidyr)
                Gaia$ph_h2o_0_20cm_ISDA = Gaia$ph_h2o_0_20cm_ISDA/10
               G2 =Gaia %>% drop_na(p_h) #2653
               G2= G2%>% drop_na(ph_h2o_0_20cm_ISDA)
                remove.duplicates()
#plot(subset[Gaia[is.na!Gaia$p_H]],subset[Gaia$ph_h2o_0_20cm_ISDA] <- incomplete syntax for dropping p_h NA without tidyr
plot(G2$p_H,G2$ph_h2o_0_20cm_ISDA, ylim=c(4,7),xlim=c(4,7),las=1)
abline(0,1)
model_isda = lm(ph_h2o_0_20cm_ISDA ~ pH , data=Gaia)
abline(model_isda,col="red")
mean(Gaia$pH, na.rm=TRUE)
mean(Gaia$ph_h2o_0_20cm_ISDA, na.rm=TRUE)  #check for spatial clusters, map the errors, compare for other variables  (which are more and less mgmt related texture stays the same, OM changes, compare to spectroscopy M3)
plot(Gaia$pH,Gaia$sgrids_ph) #both are poorly correlated
#trendline(Gaia$pH, Gaia$sgrids_ph, model = "line2P", plot = TRUE, linecolor = "red",
#trendline(Gaia$pH, Gaia$sgrids_ph, model = "line2P", plot = TRUE, linecolor = "red", #dont have trendline() installed
#          lty = 1, lwd = 1, summary = TRUE, ePos = "topleft", eDigit = 5,
#         eSize = 1)       #if using trendline
# example formatting plot(bodymass, height, pch = 16, cex = 1.3, col = "blue", main = "HEIGHT PLOTTED AGAINST BODY MASS", xlab = "BODY MASS (kg)", ylab = "HEIGHT (cm)")

@@ -58,19 +63,20 @@
  summary(resid(model_isda))
summary((model_isda))
xlim = c(2,7)
abline (0,1)



#how to pull adjacent cells (ad subtract by LatLong equivalent of 30m???) 4 directly adjacent, 4 diagonally adjacent?
#how to pull adjacent cells (ad subtract by LatLong equivalent of 30m???) 4 directly adjacent, 4 diagonally adjacent? #while loop +- based on count number
#direct comparisons first call on corresponding attributes out from Gaia, then use lat and long to compare with isDA.  
#create ggplot % correlation, linear regression for each variable, group by country. /
#create ggplot % correlation, linear regression for each variable, group by country. 
# create new data frame from Gaia data, ad +- long lats to data frame with loops, comparing each adjacent cell value in isDa to Gaia value, save closest value cell. Repeat comparison.
#create new data frame from Gaia data, ad +- long lats to data frame with loops, comparing each adjacent cell value in isDa to Gaia value, save farthest value cell Repeat comparison.
#expected values high and low create interval. % of actual soil data values within high and low values. 
#resample isda raster at 90*90
#regression 
#uncertainty comparison
#by country
#by unique observation


#which package for ISRIC Values.
#downloading ISDA only does a part of the data set ~49500 cells, 30.3 million km2 in Africa, only 44km2 if 30mx30m


#compare shifted cell values for ph with while loop, depth to 30cm (from various depth values), 

ph <- geodata::soil_af_isda_vsi("ph")
back-transformation: x/10
ph
#class       : SpatRaster
#dimensions  : 289306, 328563, 4  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -3502583, 6354307, -4149107, 4530073  (xmin, xmax, ymin, ymax)
#coord. ref. : WGS 84 / Pseudo-Mercator (EPSG:3857)
#source      : ph.tif
#names       : ph_0-20cm, ph_20-50cm, ph-sd_0-20cm, ph-sd_20-50cm

#This has 30x30 m spatial resolution. Processing is similar to that of soil grids. 

#You can get the dev version with
#install.packages("remotes")
install
library(remotes)
remotes::install_github("rspatial/geodata")
install.packages("rcpp")
library(rcpp)
install.packages(rspatial)
install.packages("jsonlite")

library(rspatial)
library(geodata)
ph_isda <- soil_af_isda_vsi("ph")
back-transformation: x/10
ph_isda
back-transformation: x/10

ph2_isric= soil_af_vsi("phh2o", 30, stat="mean") #also stat=#<- soil_af_elements_vsi("ph")
ph2_isric
back-transformation: x/10

# install.packages("stars") #  FOR ISRIC WCS 250m need to convert to use terra extraction
#library(stars)
# bb=c(29002683,-11500000,4781127,29002683) # bounding box (homolosine) for ethiopia,keya, tanzania
# igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection
# 14.8n-11.5n 29.002683e 47.811276e
# sg_url="/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"
# gdal_translate(paste0(sg_url,'ocs/ocs_0-30cm_mean.vrt'),
#                "./crop_roi_igh_r.tif",
#                tr=c(250,250),
#                projwin=bb,
#                projwin_srs =igh)
# 
# 
# gdalwarp("./crop_roi_igh_r.vrt",
#          "./crop_roi_ll_r.vrt", 
#          s_srs=igh, 
#          t_srs="EPSG:4326", 
#          of="VRT",
#         overwrite = TRUE)
#https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/webdav_from_R.md need to use terra instead gdal is not active

#Open ISDA Data  (also can do with ISRIC)
# 
# this <- system('hostname', TRUE)
# if (this == "LAPTOP-IVSPBGCA") {
#   wd <- "G:/.shortcut-targets-by-id/1mfeEftF_LgRcxOT98CBIaBbYN4ZHkBr_/share/gaiatrials//"
# } else if (this == "DESKTOP-JS0UEJ1") {
#   wd <- "C:/Users/bolen/OneDrive/Desktop/Academia/Code_Github/open-soil-data/data" #added correct directory -BMO
# } else {
#   wd <- "~/Documents/GitHub/gaiatrials/"
# }
# setwd(wd)
#install.packages("rgdal", dependencies = TRUE)
#eth_gaia <- subset( Gaia, country == "Ethiopia") , rwa_gaia <- subset( Gaia, country == "Rwanda"), tan_gaia <- subset( Gaia, country == "Tanzania")   
#check with View(tan_gaia)


###########################################################################
# Monitoring spatiotemporal dynamics of large herbivores across an African
# rangeland using hierarchical multi-species distance sampling

# Ramiro D. Crego, Harry B.M. Wells, Grant Connette, Jared A. Stabach,
# Naitareu Soit & Stewart Thompson

# Load libraries

library(tidyverse)
library(ggplot2)
library(sf)
library(tmap)

# Load Conservancy shapefile
aoi <- st_read("./Data/Naboisho.shp")

# Transect dates
Time <- c('2017-10-01','2017-11-01', '2017-12-01', '2018-01-01', '2018-02-01', '2018-03-01','2018-04-01', '2018-05-01','2018-06-01', '2018-07-01', '2018-08-01', "2018-09-01", "2018-10-01", "2018-11-01", "2018-12-01", "2019-01-01", "2019-02-01", "2019-03-01", "2019-04-01", "2019-05-01", "2019-06-01", "2019-07-01", "2019-08-01", "2019-09-01", "2019-10-01", "2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01")
Time <- as.Date(Time, format='%Y-%m-%d')
Time <- as.data.frame(Time)

#aoi <- cbind(aoi, t(Time))

## Extract NDVI 
library(lubridate)
library(rgee)
ee_Initialize()
ee_check()

#Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#Join Image and Points based on a maxDifference Filter within a temporal window

#Set temporal window in days for filter. This will depend on the remote sensing data used.
tempwin <- 16 

#Set the filter
maxDiffFilter<-ee$Filter$maxDifference(
  difference=tempwin*24*60*60*1000, #days * hr * min * sec * milliseconds
  leftField= "date_millis", #Timestamp of the telemetry data
  rightField="system:time_start" #Image date 
)

# Define the join. We implement the saveBest function for the join, which finds the image that best matches the filter (i.e., the image closest in time to the particular GPS fix location). 
saveBestJoin<-ee$Join$saveBest(
  matchKey="bestImage",
  measureKey="timeDiff"
)

#Function to add property with raster pixel value from the matched image
add_value<-function(feature){
  #Get the image selected by the join
  img1<-ee$Image(feature$get("bestImage"))$select(band)
  #Get mean pixel value within feature
  ft_value<-img1$reduceRegion(ee$Reducer$mean(),feature$geometry())
  ft_value2<-img1$reduceRegion(ee$Reducer$stdDev(),feature$geometry())
  #Return the data containing mean value and image date.
  feature$setMulti(list(Mean = ft_value$get(band), SD = ft_value2$get(band), DateTimeImage = img1$get('system:index')))
}

# Function to remove image property from features
removeProperty<- function(feature) {
  #Get the properties of the data
  properties = feature$propertyNames()
  #Select all items except images
  selectProperties = properties$filter(ee$Filter$neq("item", "bestImage"))
  #Return selected features
  feature$select(selectProperties)
}


################
# Extract NDVI

# Load image collection and set days
start<-"2017-08-30"
end<-"2020-12-31"
imagecoll<-ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end)
band <- "NDVI" 


########
# Extract values

finalNDVIlist <- list()

for(i in 1:30){
  data1 <- mutate(aoi, Date = as.factor(Time[i,]))
  # Send sf to GEE
  data <- sf_as_ee(data1)
  # Transform day into milliseconds
  data<-data$map(add_date)
  # Apply the join
  Data_match<-saveBestJoin$apply(data, imagecoll, maxDiffFilter)
  # Add pixel value to the data
  DataFinal<-Data_match$map(add_value)
  # Remove image property from the data
  DataFinal<-DataFinal$map(removeProperty)
  # Move GEE object into R
  temp<- ee_as_sf(DataFinal, via = 'getInfo')
  # Append to list
  finalNDVIlist[[i]] <- temp
}  

finalNDVIlist[[1]]

AOI_NDVI <- matrix(NA, nrow = 30, ncol = 3)
for(i in 1:30){
  temp <- finalNDVIlist[[i]]
  temp1 <- st_drop_geometry(temp)
  AOI_NDVI[i,] <- cbind(temp1[,'Mean']/10000, temp1[,'SD']/10000, Time[i,])
}
AOI_NDVI <- as.data.frame(AOI_NDVI)
names(AOI_NDVI) <- c("Mean", "SD", "Date")
AOI_NDVI$Date <- as.Date(AOI_NDVI$Date)
AOI_NDVI$yLower <- AOI_NDVI$Mean - AOI_NDVI$SD
AOI_NDVI$yUpper <- AOI_NDVI$Mean + AOI_NDVI$SD

aoindvi <- ggplot(AOI_NDVI, aes(x = Date, y = Mean)) + geom_ribbon(aes(ymin = yLower, ymax = yUpper), fill = "grey", alpha = 0.5)+ geom_line(colour = "black") + scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") + stat_smooth(method = lm, se = T) + labs(y = "NDVI")

ggsave(file = "./Figures/NDVI.jpg", plot = aoindvi, width = 8, height = 6)

hist(AOI_NDVI$Mean)
summary(lm(AOI_NDVI$Mean ~ AOI_NDVI$Date))
plot(lm(AOI_NDVI$Mean ~ AOI_NDVI$Date))

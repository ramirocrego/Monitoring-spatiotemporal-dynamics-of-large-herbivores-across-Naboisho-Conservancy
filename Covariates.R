###########################################################################
# Monitoring spatiotemporal dynamics of large herbivores across an African
# rangeland using hierarchical multi-species distance sampling

# Ramiro D. Crego, Harry B.M. Wells, Grant Connette, Jared A. Stabach,
# Naitareu Soit & Stewart Thompson

#########################
## Create covariates

# Load libraries
library(dplyr)
library(tidyverse)
library(plyr)
library(ggplot2)
library(sf)
library(tmap)
library(stars)

# Load transects
transects <- st_read("./Data/NaboishoTransects.shp")
transects$ID <- tolower(transects$ID)
head(transects)

transectp <- transects %>% st_transform(crs = "+proj=utm +zone=36 +datum=WGS84 +units=m") %>% st_simplify(F,60) %>% st_buffer(dist = 300, endCapStyle = "FLAT", joinStyle= "ROUND")
plot(transectp[1])
trans <- st_transform(transectp, crs = "+proj=longlat +datum=WGS84")

tmap_mode("view")
tm_shape(transects) +
  tm_lines(col = "ID", palette = "Set1", lwd=2, alpha = 0.5, title.col = "Transects") +
tm_shape(trans) +
  tm_polygons(col = "ID", palette = "Set1", lwd=2, alpha = 0.5, title.col = "Transects") 


# Extract time for each survey and attach to transects
load("./Data/DistData")
str(data)

nametransects <- unique(data$TransectID)
nametransects

t1 <- data %>% filter(TransectID == nametransects[1]) %>% select(Date) %>% unique() %>% mutate(DateT = Date)

for (i in 2:8){

  t2 <- data %>% filter(TransectID == nametransects[i]) %>% select(Date) %>% unique() %>% mutate(DateT = Date)

  t1 <- full_join(t1,t2, by = 'Date')
}

t1 <- t1 %>% arrange(Date)
names(t1) <- c('Date',nametransects)
t1$YM <- paste0(format(t1$Date, "%Y"),"_",format(t1$Date, "%m"))
head(t1)


t <- t1 %>% group_by(YM) %>% summarise_all(min, na.rm = TRUE)
tfinal <- t(as.data.frame(t))
tfinal <- tfinal[c(-1,-2,-11),-36]
tfinal <- as.data.frame(tfinal)

## Combine date of transects to the sf
trans <- trans %>% add_column(tfinal)
head(trans)


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
  #Return the data containing mean value and image date.
  feature$setMulti(list(FtVal = ft_value$get(band), DateTimeImage = img1$get('system:index')))
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
start<-"2016-05-30"
end<-"2021-12-31"
imagecoll<-ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end)
band <- "NDVI" 


########
# Extract values
colnam <- colnames(tfinal)
finalNDVIlist <- list()

for(i in 1:ncol(tfinal)){
data1 <- mutate(transects, Date = tfinal[,i])
# Send sf to GEE
data <- sf_as_ee(data1)
# Transform day into milliseconds
data<-data$map(add_date)
# Apply the join
Data_match<-saveBestJoin$apply(data, imagecoll, maxDiffFilter)
# Add pixel value to the data
data<-Data_match$map(add_value)
# Remove image property from the data
data<-data$map(removeProperty)
# Move GEE object into R
temp<- ee_as_sf(data, via = 'getInfo')
# Append to list
finalNDVIlist[[i]] <- temp
}  

finalNDVIlist[[1]]

NDVIs <- matrix(NA, nrow = nrow(tfinal), ncol = ncol(tfinal))
for(i in 1:ncol(tfinal)){
  temp <- finalNDVIlist[[i]]
  temp1 <- st_drop_geometry(temp)
  NDVIs[,i] <- temp1[,'FtVal']/10000
}
NDVIs


### Prepare for model

# Order in covariates
colnamesforShape <- c("rekero woodlands","ilkisieusieu woodlands","sampu enkare woodlands","olare sampu woodlands","balanite plains","payia plains","enoolera plains","rekero plains")

# Order in data
colnamesCovData <- c("balanite plains","enoolera plains" ,"ilkisieusieu woodlands", "olare sampu woodlands", "payia plains","rekero plains","rekero woodlands","sampu enkare woodlands")

NDVIs <- t(NDVIs)
head(NDVIs)
colnames(NDVIs) <- colnamesforShape
NDVIs <- NDVIs[,colnamesCovData]

meanNDVI <- mean(NDVIs) 
sdNDVI <- sd(NDVIs)
NDVI <- (NDVIs - meanNDVI)/sdNDVI
save(NDVI, file = "./Data/NDVI") # used for model


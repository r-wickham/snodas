
library(raster)
library(rgdal)
library(RNRCS)
library(dplyr)

baseDir <- getwd()

archiveDir <- paste0(baseDir,"\\archive")
vecDir <- paste0(baseDir,"\\vector")
rasDir <- paste0(baseDir,"\\raster")

#Reading in as RGDAL object
elevNWWFileName <- paste0(rasDir,"\\elev\\usgs_elev")


#loading a SWE raster
sweFileNames <- dir(paste0(rasDir,"\\precip_liquid"),pattern = "*.tif$",full.names = T)



dwrPoly <- readOGR(dsn = paste0(vecDir,"\\dwr_wgs84.shp"))

### Functions #######################

#testing
# rasList <- list(elevNWW, swe)
rasFileNames <- c(elevNWWFileName, sweFileNames)
samplePoly <- dwrPoly

#extract paired points from raster sets
extractPairedDataFromRasters <- function(rasFileNames, samplePoly){
  #Inputs are:
  # rasFileNames     a character vector of raster names
  # samplePoly      a 'SpatialPolygonsDataFrame' object
  #
  #Ouput is a dataframe where each column is named after the name
  #  respective raster.  Data are paired data extracted from grid cells
  #  of the overlain rasters.  The raster with the smallest dimesionsal
  #  grids (by area) is used to establish the points at which to extract
  #  data (via 'rasterToPoints').  Data are then extracted from all rasters
  #  that are contained within the 'samplePoly' object (via a subset 
  #  of raster grid points using 'point.in.polygon' function)
  require(raster)
  require(rgdal)
  require(lubridate)
  
  sTime <- Sys.time()
  cat(sprintf("\nProcessing raster data:\t%s", sTime))
  
  #First, some function definitions:
  #square root of sum of squares (i.e., area compute from xy-dimensions)
  sss <- function(...) sqrt(sum(unlist(...)^2))
  
  rasList <- sapply(rasFileNames, raster) #Reading in raster datasets
  
  #Retrieving polygon coordinates
  if(class(samplePoly) != "SpatialPolygonsDataFrame")
    stop("\nextractPairedDataFromRasters: 'samplePoly' argument must of of class 'SpatialPolygonsDataFrame'")
  polyCoords <- samplePoly@polygons[[1]]@Polygons[[1]]@coords 
  
  #QC checks that function will be properly executed:
  #  - raster objects are of the correct class
  if( !all(sapply(rasList,class) %in% "RasterLayer") )
    stop("\nextractPairedDataFromRasters: All elements in 'rasList' must be of the class 'RasterLayer'")
  
  #  - polyon and rasters all overlap to some extent
  #    checks for at least 1 point of polygon indices covered by raster
  rasExtents <- sapply(rasList, extent, simplify=F)
  rasOverlapPoly <- sapply(rasExtents,
         function(x) sum(point.in.polygon(point.x = polyCoords[,1], point.y = polyCoords[,2],
                                      pol.x = x[c(1:2,2:1)],
                                      pol.y = x[c(3,3,4,4)]))>0, simplify=T)
  if( !all(rasOverlapPoly) ) 
    stop(sprintf("\nextractPairedDataFromRasters: Raster(s) do not not overlap with sample polygon:\n\t%s",
                paste0(rasFileNames[!rasOverlapPoly],collapse="\n\t")))
  
  #  - all the coordinate reference systems are the same
  allCRS <- c(sapply(rasList,crs,simplify = T),crs(samplePoly))
  if( length(unique(allCRS)) != 1 )
    stop(paste0("\nextractPairedDataFromRasters: All coordinate reference systems of",
                " 'samplePoly' and elements in 'rasList' must be the same."))
  
  #1. Check which raster has the smallest grid area;
  #  res returns a vector of two, x-y dimensions, fed to sum of square function
  minRes <- which.min(sapply(rasList, function(x) sss(res(x)) ,simplify=T))
  cat(sprintf("\nRaster to be use for points extraction:\t%s", basename(rasFileNames)[minRes]))
  
  #2. Whichever raster has the smallest resolution, use as the template
  #    for the extraction of data; reduce smallest resolution raster to 
  #    points via 'rasterToPoints'
  cat("\nReducing raster to points")
  rasCoords <- rasterToPoints(rasList[[minRes]])[,1:2] #only using x-y, third column is vertical
  
  #3. Determine which of the reduce points are in the sample polygon via 'point.in.polygon'
  #    
  cat("\nEstablishing x-y points from rasters inside sample polygon")
  ptsInPolyIndex <- point.in.polygon(point.x = rasCoords[,1], point.y = rasCoords[,2],
                   pol.x = polyCoords[,1],
                   pol.y = polyCoords[,2])
  
  #4. Extract data from each of the raster sources given the 'rasCoords' definition,
  #    applied as an index to the 'ptsInPolyIndex' object
  cat("\nExtracting data from rasters")
  ptsInPoly <- rasCoords[ptsInPolyIndex==1,] #subsetting vectorized raster pts to those in poly
  
  out <- data.frame(ptsInPoly,  #Merging x-y with extracted raster data
                sapply(rasList,extract,y=ptsInPoly[,],simplify=T))
  names(out) <- c("x","y",basename(names(rasList))) #renaming using simplified basename of rasters
  
  #return list of the vectorized raster data and rasters clipped to the extents
  #  of the sample polygon
  rasList <- sapply(rasList, crop, y=samplePoly)
  names(rasList) <- basename(names(rasList))
  
  eTime <- Sys.time()
  cat(sprintf("\nProcessing raster data:\t%s", eTime))
  cat(sprintf("\nElapsed Time:\t%s", as.period(round(sTime %--% eTime,0))))

  return(list(rasPts=out, ras=rasList))
}


### Main #################

ras <- extractPairedDataFromRasters(rasFileNames = rasFileNames,samplePoly = dwrPoly)

xList <- list("usgs_elev")
yList <- as.list(basename(sweFileNames))

library(manipulate)


snotelMeta <- bind_rows(grabNRCS.meta(ntwrks = c("SNTL","SNTLT"), cnvrt.elev = FALSE))



#determining which snotel sites in polygon
polyCoords <- dwrPoly@polygons[[1]]@Polygons[[1]]@coords
sntlInBasin <- point.in.polygon(point.x = snotelMeta$longitude,
                 point.y = snotelMeta$latitude,
                 pol.x = polyCoords[,1],
                 pol.y = polyCoords[,2])

snotelMeta <- snotelMeta[as.logical(sntlInBasin),]
snotelMeta[, c("latitude","longitude","elev_ft")] <- sapply(snotelMeta[, c("latitude","longitude","elev_ft")],as.numeric)

getNRCSSWE <- function(snotelMeta, date){
  #returns vector of SWE for a single day given metadata
  #out <- vector("numeric",nrow(snotelMeta))
  snotelMeta$swe <- NA
  for( k in 1:nrow(snotelMeta)){
    snotelMeta$swe[k] <- grabNRCS.data(network = snotelMeta$ntwk[k],
                  site_id = strsplit(snotelMeta$site_id[k],":")[[1]][2],
                  timescale = "daily",
                  DayBgn = date, DayEnd = date)$Snow.Water.Equivalent..in..Start.of.Day.Values[1]
  }
  return(snotelMeta)
}


manipulate(
  {
    layout(matrix(1:2, nrow=2))
    date <- as.Date(yName, "swe_%Y-%m-%d")
    snotelMeta <-getNRCSSWE(snotelMeta,date)
    plot(ras$rasterPts$usgs_elev*3.28, rasterPts[,yName], pch=3, cex=0.1, main=yName,
         ylim=c(0,max(c(snotelMeta$swe,rasterPts[,yName]))) )
    points(x = as.numeric(snotelMeta$elev_ft), snotelMeta$swe, pch=3, cex=2, col="red")
    text(x = as.numeric(snotelMeta$elev_ft), y=snotelMeta$swe, labels=snotelMeta$site_name,col="red")
    

    # plot(ras$ras$usgs_elev)
    plot(ras$ras[[yName]])
    plot(dwrPoly,add=T)
    points(x = as.numeric(snotelMeta$longitude), as.numeric(snotelMeta$latitude), pch=3, cex=1, col="red")
    text(x = as.numeric(snotelMeta$longitude), y=as.numeric(snotelMeta$latitude), labels=snotelMeta$site_name)
    
  },
  yName=picker(yList))

#Plots showing:
# x-y of each grid point with an overlain histogram of averages every 500' elevation
#  for the swe and accumulated precip (maybe separate out histogram into solid/liquid)
#spatial overlay of SWE and accumulated precip
#show HUC8, and rivers, with the sample polygon boldened


#Possible in plotly?
# plot_ly(x=ras$rasPts$usgs_elev,
#         y=ras$rasPts$`swe_2005-3-15.tif`,
#         text=paste0()) %>% add_markers()

library(leaflet)
library(d3scatter)
library(crosstalk)
library(plotly)

# maybe with some crosstalk for the Snotel data?
#https://rstudio.github.io/crosstalk/using.html


#Adding SWE to snotel metadata
yName <- 'swe_2003-12-1.tif'
date <- as.Date(yName, "swe_%Y-%m-%d")
snotelMeta <-getNRCSSWE(snotelMeta,date)

# sharedSwe <- SharedData$new(snotelMeta)

names(ras$rasPts) <- gsub("-","_",names(ras$rasPts))

#color pallette
pal <- colorBin(palette = rainbow(20,end = 0.8),domain = ras$rasPts$swe_2003_12_1.tif, bins=20, na.color = "#00000000")  

p <- bscols(

  plot_ly(data = snotelMeta,x= ~elev_ft,
          y= ~swe, name="SNOTEL",
          marker=list(size=12,color="black"),
          text=~paste0("Site: ",site_name,"<br>Elev (ft):  ", elev_ft),
          type = "scatter", mode="markers")  %>%
    add_text(x = ~elev_ft, y=~swe, data = snotelMeta,
             text = ~site_name, textposition="top center",
             showlegend=F) %>%
    add_markers(data = ras$rasPts,
                x=~usgs_elev*3.28, y=~swe_2003_12_1.tif,
                name="SNODAS",
                marker=list(size=1,color=hsv(0,0,0,0.5)),
                hoverinfo="none", text="") %>%
    layout(yaxis=list(title="SWE (in)"), xaxis=list(title="Elevation (ft)")),
  


leaflet(data = sharedSwe, width = 800,height = 800) %>%
  addTiles() %>%
  addCircleMarkers(label=~site_name,
                   labelOptions =
                     labelOptions(noHide = T, textOnly=T,
                                  style=list("color"="black",
                                             "font-size"="16px"))
                   ) %>%
  addPolygons(lng = polyCoords[,1], lat = polyCoords[,2], fill=F) %>%
  addRasterImage(x = ras$ras$`swe_2003-12-1.tif`,opacity = 0.3,
                 colors=rainbow(20,end = 0.8) ) %>% 
  addLegend(pal = pal, values = ras$rasPts$swe_2003_12_1.tif)
)


htmltools::save_html(html = p, file = "C:/Temp/DWR_SWE_TEST.html")


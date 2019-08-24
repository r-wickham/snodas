

source(dir("scripts/util",full.names = T))

startDate <- as.Date("2015-09-29")
endDate <- as.Date("2018-09-29")
endDate <- Sys.Date()
dates <- seq(startDate,endDate,by="d")
saveDir <- "raster"
parameters = c("swe","liquid precipitation","solid precipitation")
# 
refExtentObject <-
  rgdal::readOGR(dsn = "vector\\DistrictBoundary_100mi_Buffer.shp")
# 
rList <- downloadSNODASparallel(dates = dates,
                                saveDir = saveDir,
                                createArchive = T,
                                parameters = parameters,
                                refExtentObject = refExtentObject)

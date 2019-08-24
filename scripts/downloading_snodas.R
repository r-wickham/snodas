

source(dir("scripts/util",full.names = T))

startDate <- as.Date("2007-03-25")
endDate <- as.Date("2018-09-29")
endDate <- Sys.Date()
dates <- seq(startDate,endDate,by="d")
saveDir <- "raster"
parameters = c("swe","liquid precipitation","solid precipitation")
# 
refExtentObject <-
  rgdal::readOGR(dsn = "vector\\DistrictBoundary_100mi_Buffer.shp")

for( wy in 2009:2019){
  
  startDate <- as.Date( sprintf("%g-10-01",wy-1) )
  endDate <- min(as.Date( sprintf("%g-09-30",wy) ), Sys.Date())
  dates <- seq(startDate,endDate,by="d")
  
  rList <- downloadSNODASparallel(dates = dates,
                                  saveDir = saveDir,
                                  createArchive = T,
                                  parameters = parameters,
                                  refExtentObject = refExtentObject)
  rm(rList); gc()
}

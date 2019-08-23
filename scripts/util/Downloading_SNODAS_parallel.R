### Downloading Files from FTP Site ####################################
#This is a function that, given a date(s), an appropriate data type(s),
#  and an extent object will download the corresponding SNODAS data and 
#  save it to a designated file location.



downloadSNODASparallel <- function(dates, saveDir=NULL, parameters,
                                   createArchive=F, maxDlAttempts=5, refExtentObject=NULL){
  #Input arguments:
  #'dates'  is a vector containing either dates (e.g., Date[1:2], format: "2018-04-04" "2018-04-05")
  #         or character strings of dates (e.g., chr [1:2] "2018-04-04" "2018-04-05")
  #'saveDir' is the fully qualified directory path to the location where data will be saved
  #           NOTE: data will be saved within a subfolder of 'saveDir' indicating the type of 
  #          data (e.g., 'swe', 'precip')
  #'parameters' is a vector of strings corresponding to parameters to be archived
  #            All parameter strings should correspond to those in snodasMeta$par
  #'createArchive' is a logical specifying whether or not an archive of geotiffs should be created
  #'maxDlAttempts' integer specifying the maximum number of tries to download from the ftp site
  #'refExtentObject' is any object that has an extent that can be extracted.  This will be used
  #          to clip the SNODAS data prior to saving to geotiff
  require(raster)
  require(parallel)
  require(doParallel)
  require(foreach)
  require(lubridate)
  require(purrr)
  require(rgdal)
  #'
  #SNODAS info:
  #File format: Flat binary, 16-bit signed integer (big-endian)
  #From NSDIC Special report (http://nsidc.org/pubs/documents/special/nsidc_special_report_11.pdf)
  # byteorder M  big-endian
  # layout bil   raster is formatted using Band Interleave by Line (BIL): X[col,band,row]
  # nbands 1     Only one band for SWE
  # nbits 16     2-byte
  # ncols 6935  lines
  # nrows 3351  samples
  ulxmap <- -124.729583333331703
  ulymap <- 52.871249516804028
  llxmap <- -66.9420833564
  llymap <- 24.9504161946
  xdim <- 0.00833333333
  ydim <- 0.00833333333
  
  # SNODAS spatial reference and scaling
  nCol <- 6935
  nRow <- 3351 #columns and rows number:masked version of contiguous US
  snodasCRS <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  ####For testing: #####
  # dates = seq(as.Date("2017-10-01"), as.Date("2018-09-30"), "d")
  saveDir <- "raster"
  # parameters = c("swe","snow depth", "accumulated total precipitation",
  #                "accumulated solid precipitation",
  #                "accumulated liquid precipitation")
  # maxDlAttempts=5
  # refExtentObject <- readOGR(dsn = "vector\\DistrictBoundary_100mi_Buffer.shp")

  #provides a pairing between the parameter and the file
  #  that gets downloaded from the ftp.  Only need
  #  to add in a date in the form of YYYYMMDD 
  #  to the file template using sprintf
  snodasMeta <- data.frame(
    par=c("swe",
          "snow depth",
          "snow melt runoff at the base of the snow pack",
          "sublimination from the snow pack", 
          "sublimination from blowing snow",
          "solid precipitation",
          "liquid precipitation",
          "snow pack average temperature",
          "accumulated total precipitation",
          "accumulated solid precipitation",
          "accumulated liquid precipitation"
    ),
    code = c(1034,1036,1044,1050,1039,1025, 1025,1038,NA,NA,NA),
    #https://nsidc.org/data/g02158: To convert integers in files to model output values, divide integers in files by scale factor
    scaleFactor =c(1000,1000,100000,100000,100000,10,10,1,NA,NA,NA), 
    dirName=c("swe","depth","snow_melt","sublim_snow_pack",
              "sublim_blowing_snow","precip_solid","precip_liquid","snow_temp",
              "acc_precip_total","acc_precip_solid","acc_precip_liquid"),
    vcode= c("",  "",  "",  "",  "", "IL01","IL00",  "",NA,NA,NA),
    fileTemplate=c(
      "us_ssmv11034tS__T0001TTNATS%s05HP001.dat.gz",
      "us_ssmv11036tS__T0001TTNATS%s05HP001.dat.gz",
      "us_ssmv11044bS__T0024TTNATS%s05DP000.dat.gz",
      "us_ssmv11050lL00T0024TTNATS%s05DP000.dat.gz",
      "us_ssmv11039lL00T0024TTNATS%s05DP000.dat.gz",
      "us_ssmv01025SlL01T0024TTNATS%s05DP001.dat.gz",
      "us_ssmv01025SlL00T0024TTNATS%s05DP001.dat.gz",
      "us_ssmv11038wS__A0024TTNATS%s05DP001.dat.gz", NA,NA,NA
    ),
    units=c("m","m","m","m","m","kg_m2","kg_m2","K", NA,NA,NA),
    stringsAsFactors = F)
  
  monthLabels = c("01_Jan", "02_Feb", "03_Mar", "04_Apr", "05_May", "06_Jun", 
                   "07_Jul", "08_Aug", "09_Sep", "10_Oct", "11_Nov", "12_Dec")
  
  ### Functions ###############
  
  waterYear <-  function(dateFun){
    require(lubridate)
    wy <- as.numeric(format(x=dateFun, format="%Y"))
    wy[month(dateFun) >= 10] <- wy[month(dateFun) >= 10] + 1
    return(wy)
  }
  
  #Extracts from right side of character (similar to Excel function)
  right <- function(text, nchars)   substring(text, nchar(text)-nchars+1, nchar(text))
  
  #Extracts from left side of character(similar to Excel function)
  left <- function(text, nchars)  substring(text, 1, nchars)
  
  #Extracts %Y-%m-%d as string from all string elements
  datesFromString <- function(strings)str_extract_all(strings,"\\d{4}-\\d{2}-\\d{2}")

  #Gets the parameter name from the parameter directory
  parFromParDir <- function(parDir) snodasMeta$par[match(parDir,snodasMeta$dirName)]
  parDirFromPar <- function(par) snodasMeta$dirName[match(par,snodasMeta$par)]
  
  is.consecutive <- function(vec,allowGaps=T){
    #returns logical indicating which elements in the input vector are consecutive
    #  with adjacent vector elements. If allowGaps=F, then will only return true for the
    #  first series of consecutive numbers
    # e.g., allowGaps=T; vec=c(1,2,3,5,6,7) -> c(T,T,T,T,T,T)
    #       allowGaps=F; vec=c(1,2,3,5,6,7) -> c(T,T,T,F,F,F)
    out <- rep(T,length(vec))
    vecDiff <- diff(vec)==1
    foundGap <- F
    for( k in 1:length(vec) )
      if(allowGaps){
        if( !any(vecDiff[ pmin(pmax(c(k,k-1),1),length(vecDiff)) ]) ) out[k] <- F
      }else{
        if( !any(vecDiff[ pmin(pmax(c(k,k-1),1),length(vecDiff)) ]) | foundGap ) {
          out[k] <- F
        }
        if( any(!vecDiff[ pmin(pmax(c(k,k-1),1),length(vecDiff))]) ) foundGap=T
      }
    return(out)
  }
  
  datFileToRaster <- function(datFileName, par){
    #Performs operations required to read in SNODAS dat file to 
    #Read in binary data
    dataScaleFactor <- snodasMeta$scaleFactor[snodasMeta$par==par]
    
    rasCon <- gzcon(file(datFileName,"rb"))
    rasData <- readBin(rasCon, integer(), n=nRow*nCol, size=2,
                       signed=TRUE, endian='big')/dataScaleFactor
    close(rasCon)
    
    #assign data to raster
    r <- raster(x = matrix(rasData, nrow=nRow, ncol=nCol, byrow = T))
    
    #Assign extents, crs
    extent(r) <- extent(c(ulxmap, llxmap, llymap,  ulymap) )
    crs(r) = snodasCRS
    return(r)
  }
  
  readUntarredSNODAS <- function(unpackedDir,parameters){
    #Inputs are 1) the fully qualified directory path to an untarred 
    #  SNODAS dataset for a single day, and 2) the SNODAS
    #  parameters to read.  See snodasMeta$par for exact names
    r <- list()
    dateTag <- right(unpackedDir,8)
    for(par in parameters){
      
      rasFileTemplate <- snodasMeta$fileTemplate[snodasMeta$par==par]
      rasFile <- sprintf(sprintf("%s/%s",unpackedDir,rasFileTemplate), dateTag)
      
      r[[par]] <- datFileToRaster(rasFile, par)
      
    } #end par loop
    return(r)
  }
  
  convertRasterToInchesFromMeters <- function(ras){
    #https://nsidc.org/data/g02158
    #SWE, snow depth, sublimination in meters
    ras[ras == -0.09999] <- NA  #Setting NA; scale factor of 100,000
    ras[ras == -9.999] <- NA  #Setting NA; scale factor of 1,000
    ras <- ras*1000/25.4     #to inches
    return(ras)
  }
  
  #assume 1000kg/m^3 of water
  #[kg/m^2] = [kg/m^2]*[1m^3/1000kg][1000mm/m] -> [m]
  convertRasterToInchesFromKGM2 <- function(ras){
    #https://nsidc.org/data/g02158
    #precipitation in kg/m^2
    ras[ras == -999.9] <- NA  #Setting NA; used scale factor of 10
    ras <- ras/25.4     #to inches from kg/m^2
    return(ras)
  }
  
  convertRasterToDegF <- function(ras){
    #https://nsidc.org/data/g02158
    #snow pack temp in Kelvin
    ras[ras == -9999] <- NA      #Setting NA; used scale factor of 1
    ras <- (ras-273.15)*9/5 + 32 #Conversion to deg F from Kelvin
    return(ras)
  }
  
  cleanSNODASRasterList <- function(rList){
    #Performs unit conversions and removal of undefined grid cells in SNODAS raster
    #  datasets.  Input is a list of raster objects, as returned from
    #  the 'dlSNODAStoRaster' function
    for( rName in names(rList)){
      if(rName %in% c("swe", "snow depth","sublimination from the snow pack", "sublimination from blowing snow" ))
        rList[[rName]] <- convertRasterToInchesFromMeters(rList[[rName]])
      
      if(rName %in% c("solid precipitation", "liquid precipitation"))
        rList[[rName]] <- convertRasterToInchesFromKGM2(rList[[rName]])
      
      if(rName %in% c("snow pack average temperature"))
        rList[[rName]] <- convertRasterToDegF(rList[[rName]])
    }
    return(rList)
  }
  
  dlSNODAStoRasterList <- function(date,parameters, cleanRaster=T){
    
    date <- as.Date(date, origin = "1970-01-01")
    dateTag <- format(date,"%Y%m%d") #Creating a tag that will be used in the ftp address
    
    ###___Retrieving Files from FTP ###########################
    #Creating urls/directory references for files
    tarFileName <- sprintf("SNODAS_%s.tar", dateTag)
    
    tarURL <- sprintf("ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/%g/%s/%s",
                      year(date), monthLabels[month(date)], tarFileName)
    
    
    ###___Downloading SNODAS tar file ####################
    
    dlFilePath <- paste0(tempDir,"/",tarFileName) #download location
    dateExists <- F; attempts <- 1 #initializing loop parameters
    while( !dateExists ){
      tryUrl <- try(expr=download.file(url=tarURL, destfile = dlFilePath, quiet = F, mode="wb"))
      if( class(tryUrl) != "try-error" ) dateExists <- T
      if(attempts > maxDlAttempts) break
    }
    
    #If data doesn't exist for that date, then removing the generated file and continuing
    #  to next iteration
    if(!dateExists){file.remove(dlFilePath); return(NA) }
    
    #If statement for if date exists in SNODAS FTP
    if( dateExists ){
      
      #Establishing directory for where to untar the results
      unpackedDir <- paste0(tempDir,"/",gsub(".tar","",tarFileName))
      tarError <- tryCatch(expr= untar(tarfile=dlFilePath,exdir = unpackedDir),
                           error=function(e) return(T))
      
      file.remove(dlFilePath) #Getting rid of data download
      #Because there are occassionally errors with the download, this code will
      #  then check that the tar file can actually be read without breaking the 
      #  foreach loop
      if( class(tarError) == "try-error" ){unlink(unpackedDir,recursive = T); return(NA)}
      
      r <- readUntarredSNODAS(unpackedDir,parameters)
      
      unlink(unpackedDir,recursive = T)
      if(cleanRaster) r <- cleanSNODASRasterList(r)

    } #End date exists if
    return(r)
  }
  
  saveRasters <- function(r,date){
    #Iterate through each list element
    #save to a folder corresponding to same name as list element (e.g., 'swe')
    #  in 'saveDir'.  Save file name is prefixed with parameter name
    for(par in names(r)){
      dirName <- snodasMeta$dirName[snodasMeta$par==par] #directory name
      parSaveDir <- sprintf("%s\\%s",saveDir,dirName)    #full directory path
      if(!dir.exists(parSaveDir)) dir.create(parSaveDir,recursive = T) #creating directory if needed
      tifFile <- sprintf("%s_%s.tif",dirName,format(date,"%Y-%m-%d"))
      saveFileName <- sprintf("%s/%s",parSaveDir,tifFile) #full save file path
      writeRaster(x=r[[par]], filename=saveFileName,overwrite=TRUE) #saving to raster
    }
  }
  
  "%!in%" <- function(x,y) !(x %in% y)
  
  sumRasters <- function(rList,rIndex=NULL){
    #inputs are a list of rasters, and an optional index
    #  indicating which of the raster elements should be summed
    if(is.null(rIndex)) rIndex <- 1:length(rList) #sum all if not specified
    Reduce("+",rList[rIndex]) #applying summation to each specified element
  }
  
  accumulateRasterList <- function(rSubList, pars){
    #Takes a raster list and accumulates from smalles to largest date
    parDirs <- parDirFromPar(pars) #named in raster list as 
    wyDates <- names(rSubList)
    wyDates <- sort(wyDates) #sort low-high
    accR <- list()
    accR[[wyDates[1]]] <- rSubList[[wyDates[1]]] #Initializing accumulation raster
    for( k in 2:length(wyDates)){ #Iterating through each day and required accumulation parameter
      accR[[wyDates[k]]] <- list()
      for( parDir in parDirs ){
        accR[[wyDates[k]]][[parDir]] <- 
          sumRasters( c(accR[[wyDates[k-1]]][[parDir]], rSubList[[wyDates[k]]][[parDir]]) )
        }
      }
    return(accR)
  }
  
  computeTotalAccPrecip <- function(accRList){
    #Sums each raster by date
    wyDates <- unique(sapply(accRList,names,simplify=F))
    wyDates <- sort(wyDates) #sort low-high
    out <- list()
    for( k in 1:length(wyDates))
      out[[k]] <- sumRasters( c(accRList[[wyDates[k]]], accRList[[wyDates[k]]]) )
    return(out)
  }
  
  getAvailableDates <- function(rasDir){
    #From the directory containing the processed tif files, returns
    #  a vector of all the available dates
    require(stringr)
    as.character(str_extract_all(dir(rasDir,"tif"),"\\d{4}-\\d{2}-\\d{2}",T))
  }
  
  getWYDateExtents <- function(dates){
    #From a vector of dates, returns a dataframe tabulating the minimum and
    #  maximum date within each water year
    #e.g.,
    # > str(out)
    # 'data.frame':	3 obs. of  3 variables:
    # $ wy       : num  2017 2018 2019
    # $ startDate: Date, format: "2017-09-29" "2017-10-01" "2018-10-01"
    # $ endDates : Date, format: "2017-09-30" "2018-09-30" "2019-08-16"
    out <- NULL
    for( wy in unique(waterYear(dates)) ){
      wyDates <- dates[waterYear(dates) %in% wy] #extracting dates in wy
      wyDates <- wyDates[is.consecutive(wyDates,F)] #reducing to consecutive dates
      out <- rbind(out,
                   data.frame(wy=wy,
                              startDate=min(dates[waterYear(dates) %in% wy]),
                              endDates = max(dates[waterYear(dates) %in% wy])))
    }
    return(out)
  }
  
  getAvailablePrecipData <- function(saveDir){
    #Checks the save directory for both the available precipitation data.
    #Returns a nested list that provides the start and end dates of each
    #  precip parameter data by water year
    
    #Getting directories associated with solid and liquid precip
    precipDirs <- snodasMeta$dirName[grepl("precip",snodasMeta$par) & !grepl("acc",snodasMeta$par)]
    precipPars <- snodasMeta$par[grepl("precip",snodasMeta$par) & !grepl("acc",snodasMeta$par)]
    
    #within each directory, extract dates
    allDates <- lapply( sprintf("%s\\%s",saveDir,precipDirs), getAvailableDates)
    allDates <- lapply(allDates,as.Date)
    
    dateExtents <- lapply(allDates, getWYDateExtents)
    names(dateExtents) <- precipPars
    
    return(dateExtents)
  }
  
  
  getRequiredPrecipPars <- function(accParameters){
    #Determining the required precipitation parameters from the input accumulation parameters
    requiredPrecipPars <- NULL
    if( "accumulated total precipitation" %in% accParameters ) requiredPrecipPars <- precipPars
    if(  "accumulated solid precipitation" %in% accParameters )
      requiredPrecipPars <- unique(requiredPrecipPars, grep("solid",precipPars,value = T))
    if(  "accumulated liquid precipitation" %in% accParameters )
      requiredPrecipPars <- unique(requiredPrecipPars, grep("liquid",precipPars,value = T))
    return(requiredPrecipPars)
  }
  
  getWYsWithData <- function(availablePrecipData){
    #Extracts the water years that start on 10-01
    wys <- sapply(availablePrecipData,function(x) x$wy[month(x$startDate)==10 & day(x$startDate)==1])
    return(unique(as.vector(wys)))
  }
  
  loadRasters <- function(wyDates, pars, saveDir){
    #Returns a list of rasters named after 'pars'
    #looks inside directories associated with pars and returns
    #  those corresponding to dates specified in 'wyDates'
    
    parDirs <- snodasMeta$dirName[match(pars,snodasMeta$par) ]
    
    #Nested list of raster file names: date > parameter
    rasFileNames <- 
      lapply(wyDates, function(x) sapply(sprintf("%s/%s",saveDir,parDirs),dir,pattern = as.character(x),full.names=T) )

    rList <- lapply(rasFileNames, function(x) lapply(x, raster))
    names(rList) <- wyDates #Naming lists after dates
    #naming nested lists after parameters
    for( k in 1:length(rList)) names(rList[[k]]) <- basename(names(rList[[k]]))
    return(rList)
  }

  ### QC checks #############

  # - all parameter strings in snodasMeta$par
  if( !all(parameters %in% snodasMeta$par) )
    stop( sprintf("\nAllowable 'parameters' arguments are:\n\t%s\n\nUnallowed parameter passed to function:\n\t%s\n\n",
                  paste0(snodasMeta$par,collapse="\n\t"),
                  paste0(parameters[parameters %!in% c(snodasMeta$par,"precip_accumulated")],collapse="\n\t")) )
  
  # - dates are convertible to date objects
  dDates <- as.Date(dates) #date strings as date objects
  if( any(is.na(dDates)) )
    stop(sprintf( paste0("\nAll 'dates' need to be convertible to date objects",
                         " and passed as unambiguous format (%%Y-%%m-%%d).\nBad dates:\n\t%s"),
                  paste0(dates[is.na(dDates)],collapse="\n\t")) )
  
  # - dates are after 2003-09-30 (when SNODAS started) and before tomorrow
  if( any(dDates < "2003-09-30") )
    warning( sprintf(paste0("\nDates need to be after beginning of available SNODAS",
                            " data (2003-09-30).\nSkipping %g dates:\n\t%s"),
                     sum(dDates < "2003-09-30"),
                     paste0(dDates[dDates< "2003-09-30"],collapse="\n\t")) )
  dDates <- dDates[dDates >= "2003-09-30" & dDates <= Sys.Date()]
  
  # - if creating archive, save directory needs to be defined
  if( createArchive & is.null(saveDir))
    stop("\nIf 'createArchive' is true, need defined a directory path for 'saveDir'.")

  
  ### Additional checks ###########
  # - data have already been processed - prompt if they have and remove date from 'dDates'
  ### TODO ####
  
  
  # - if processing the accumulated precip, check if have enough data
  #   stop if not enough data and indicate which dates need processing and how to 
  #   process (i.e., 'set parameters to liquid and solid precip for these dates...')
  
  #Separating the accumulation parameters (computed after download) from those parameters
  #  avaialable directly from SNODAS tar file
  precipPars <- snodasMeta$par[grepl("precip",snodasMeta$par) & !grepl("acc",snodasMeta$par)]
  accParameters <- parameters[is.na(snodasMeta$code[match(parameters,snodasMeta$par)])]
  accumulatePrecip <- length(accParameters) > 0
  #removing from accumulation parameters list for download configuration
  parameters <- parameters[ parameters %!in% accParameters] 

  # - Add readme file defining units for each raster and extents used for clipping
  

  
  ### Main ############################

  #Creating temporary directory to work in (if needed)
  tempDir <- paste0(ifelse(is.null(saveDir),tempdir(),dirname(saveDir)),"/temp")
  cat(sprintf("\nUsing '%s' as temporary directory",tempDir))
  if(!dir.exists(tempDir)) dir.create(tempDir)
  
  nCores <- min(length(dDates),detectCores()-1)
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl, cores = nCores)
  
  rList <-
    foreach::foreach(date = dDates, .packages = c("raster","rgdal", "lubridate")) %dopar% {
      try({
        r <- dlSNODAStoRasterList(date, parameters,T) #download, create raster list

        #Clip raster using reference object if provided
        if( !is.null(refExtentObject) ) r <- lapply(r,crop,y=refExtentObject)
        
        if(createArchive) saveRasters(r,date) #Saving to disk
        gc() #running garbage collection
        return(r)
      })
    } #End parallel foreach loop
  stopCluster(cl)
  names(rList) <- as.character(dDates) #naming after dates
  

  ### _Accumulating Precipitation ##################
  ### TODO #########
  #compute accumulated precip by WY
  
  # determining which water years have necessary data to perform accumulation calc
  if( accumulatePrecip ){ 
    cat("\nChecking if enough available data to compute accumulated precip")
    #Precipitation Parameters required to compute the desired accumulation parameters
    requiredPrecipPars <- getRequiredPrecipPars(accParameters)
    #dataframe providing available date ranges for previously downloaded data
    print(availablePrecipData <- getAvailablePrecipData(saveDir)) 
    availablePrecipData <- availablePrecipData[requiredPrecipPars] #reducing only to required parameters
    
    allWYs <- unique(unlist(purrr::map(availablePrecipData,"wy"))) #all available WYs with any data
    wysToComputeAccumulation <- getWYsWithData(availablePrecipData) #available WYs with data starting Oct-01
    
    if( any(allWYs %!in% wysToComputeAccumulation) )
      warning(sprintf(paste0("\nUnable to compute accumalted precipitation for the following water years:\n\t%s",
                             "\n\nWill need to run this function to archive precip data for all",
                             " necessary dates in water year(s)."),
                      paste0(allWYs[allWYs %!in% wysToComputeAccumulation],collapse="\n\t")))
    cat(sprintf("\nComputing accumulated precip for the following water year(s):\n\t%s\n",
                paste0(wysToComputeAccumulation,collapse="\n\t")))

    #for wy in wysWithAvailableData
    for( wy in wysToComputeAccumulation){
      cat(sprintf("\nAccumulating water year:\t%g",wy))
      wyDates <- seq(as.Date(sprintf("%g-10-01",wy-1)), as.Date(sprintf("%g-09-30",wy)),"d")
      
      #Loading all available dates
      precipRList <- loadRasters(wyDates, requiredPrecipPars,saveDir)
      
      
      computeTotalPrecip <- function(precipRList){
        #Takes a raster list of precip (solid and liquid) and computes total
        parDirs <- parDirFromPar(pars) #named in raster list as 
        wyDates <- names(rSubList)
        wyDates <- sort(wyDates) #sort low-high
        for( k in 2:length(wyDates)){ #Iterating through each day and summing liquid and solid precip
          accR[[wyDates[k]]][["total precipitation"]] <- sumRasters( accR[[wyDates[k]]] )
        }
      }
      
      ### CONTINUE HERE ###########
      #QC - check this works
      #
      
      #Computing total accumulated precipitation
      if( "accumulated total precipitation" %in% accParameters )
        # accRList[["accumulated total precipitation"]] <- computeTotalAccPrecip(accRList)
        accRList <- computeTotalPrecip(precipRList)
      
      accRList <- accumulateRasterList(accRList, accParameters) #Accumulating specified parameters
      


      #Writing out by date
      for(wyDate in wyDates){
        saveRasters(accRList,wyDate)
      }

      
      
      
    }
    
    warning("\n\nHaven't yet developed ability to compute accumulated precip")
  }
  return(rList)
}


#### TESTING #######


startDate <- as.Date("2015-09-29")
endDate <- as.Date("2018-09-29")
# # endDate <- Sys.Date()
dates <- seq(startDate,endDate,by="d")
# saveDir <- "D:\\RunoffVolumeForecasting\\raster"
# parameters = c("swe","liquid precipitation","solid precipitation")
# 
# refExtentObject <-
#   rgdal::readOGR(dsn = "D:\\RunoffVolumeForecasting\\vector\\DistrictBoundary_100mi_Buffer.shp")
# 
# rList <- downloadSNODASparallel(dates = dates,
#                                 saveDir = saveDir,
#                                 createArchive = F,
#                                 parameters = parameters,
#                                 refExtentObject = refExtentObject)

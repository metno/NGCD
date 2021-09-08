`nc4out` <- 
function(
         grid.list, #list of grids: (x,y,t) or (x,y,ens,t)
         grid.type="utm33",
         x,
         y,
         file.name="dummy.nc",
         var.name="DUMMY",
         var.longname = var.name, 
         var.standardname = var.name, 
         var.unit = NULL,
         var.mv = -999.99,
         var.version = "no version",
         times = NULL,
         times.unit,
         times.ref = "190001010000",
         prod.date = substr(Sys.time(),1,10),
         proj4.string = "",
         # global attributes
         source.string="",
         reference = "unknown",
         summary = "",
         title = "",
         license = NULL,
         comment="none",
         # add forecast reference time as scalar variable
         for_ref_time.out=F,
         for_ref_time.val=NULL,
         # additional variable attributes
         atts.var.add=NULL,
         # add lon-lat as 2D variables
         lonlat.out=T,
         # round digit
         round.dig=NULL
         ) {
#==============================================================================
# input:
# grid.list: A list og grids object, i.e. a matrix/array (possibly for several times, i.e. 3D, and/or for several ensemble members,i.e. 4D)
# grid.type: The type of the grid: either "utm33", "lonlat", "radlcc", "ngcd"
# x,y: The coordinate vectors (x: easting or longitude, y: northing or latitude)
# file.name: The file name to write. A character string
# var.name: The variable name to write. A character string
# var.longname: The long name of the variable. Check out long-names
# var.unit: The variable's units, a character string
# var.mv: The missval to use (NAs in R are replaced by this number )
# var.version: A version indicator
# times: The time(s) of the grid.A character vector of convention "yyyymmddhhmm"
# If not specified, the function searches for a "time" attribute in grid
# to get the information
# times.unit: The time unit to use for the time coordinate in the file. 
# Possible options are: "Y" for annual, "M" for months, "D" for days,
# times.ref: "H" for hours and "T" for minutes. 
# The reference time to use for the time coordinate. 
# Time units in the file are given as the elapsed time since \code{t.ref}.
# prod.date: A production date/time to be written into the file
# proj4.string: proj4.string
# source.string: source
# reference: A reference to a paper or documentation of the dataset.
# lonlat.out: should we add the lat/lon fields? (case of grid.type!="latlon")
# round.dig: rounds the grid values to the specified number of decimal places.
#==============================================================================
  # license
  if (is.null(license)) 
    license<-"https://www.met.no/en/free-meteorological-data/Licensing-and-crediting"
#-------------------------------------------------------------------------------
# set compression level for writing netcdf
  compr.lev<-9
#-----------------------------------------
# load netcdf library
  if (!requireNamespace("ncdf4",quietly=TRUE)) {
    stop("Package \'ncdf\' required for this function.") 
  }
  if (!requireNamespace("rgdal",quietly=TRUE)) {
    stop("Package \'rgdal\' required for this function.") 
  }
# -----------------------------------------
# checks and default handling
# -----------------------------------------
  dimv<-length(grid.list)
  # check all the vectors are consistent 
  if (length(var.name)!=dimv) 
    stop(paste("Length of \"var.name\" inconsistent with length of \"grid.list\"."))
  if (length(var.longname)!=dimv) 
    stop(paste("Length of \"var.longname\" inconsistent with length of \"grid.list\"."))
  if (length(var.standardname)!=dimv) 
    stop(paste("Length of \"var.standardname\" inconsistent with length of \"grid.list\"."))
  if (length(var.unit)!=dimv) 
    stop(paste("Length of \"var.unit\" inconsistent with length of \"grid.list\"."))
  # check time specifications
  if (missing(times.unit)) {
     stop("Argument \"times.unit\" must be specified.")
  }
  if (is.null(times)) {
     times <- attr(grid.list[[1]],"time")
     if (is.null(times)) {
        stop("No argument \"times\" and no attribute \"time\" found in \"grid\".")
     }
  }
  if (!is.character(times)) {
     stop("Inappropriate class of argument/attribute \"times\".")
  }
  # all the grid should have the same length of "times"
  if (is.matrix(grid.list[[1]])) { 
    dimt <- 1 
  } else {
    dimt <- dim(grid.list[[1]])[length(dim(grid.list[[1]]))] 
  }
  for (v in 1:dimv) {
    if (is.matrix(grid.list[[v]])) {
      dimt.v<-1
    } else {
      dimt.v<-dim(grid.list[[v]])[length(dim(grid.list[[v]]))]
    }
    if (dimt.v!=dimt) {
      stop(paste("Length of \"times\" for ",v,"-th variable inconsistent with length of \"times\" for the 1-st variable. (",dimt.v,"instead of",dimt,")"))
    }
  }
  if (length(times) != dimt) {
     stop("Length of \"times\" inconsistent with dimensions of \"grid\".")
  }
  # check x and y dimensions
  for (v in 1:dimv) {
    if (length(x) != dim(grid.list[[v]])[1]) {
         stop(paste(v,"-th variable: Length of \"x\" inconsistent with dimensions of \"grid\"."))
    }
    if (length(y) != dim(grid.list[[v]])[2]) {
       stop(paste(v,"-th variable: Length of \"y\" inconsistent with dimensions of \"grid\"."))
    }
  }
  if (!is.character(var.unit)) {
      warning("No units specified for grid. <dimensionless> chosen instead.")
      var.unit <- "no dim"
  }
  # ensemble member numbers
  isEns<-vector(length=dimv,mode="logical")
  dime<-NA
  for (v in 1:dimv) {
    if (is.matrix(grid.list[[v]])) { 
      isEns[v]<-F
    } else {
      if (length(dim(grid.list[[v]]))==3) {
        isEns[v]<-F 
      } else if (length(dim(grid.list[[v]]))==4) {
        isEns[v]<-T 
        if (is.na(dime)) {
          dime<-dim(grid.list[[v]])[3]
        } else if (dime!=dim(grid.list[[v]])[3]) {
          stop(paste(v,"-th variable: inconsistent number of ensemble members."))
        }
      } else {
       stop(paste(v,"-th variable: incorrect number of grid dimensions."))
      }
    }
  }
  # check if grid.type is available
  if (!(grid.type %in% c("lonlat","utm33","radlcc","yrlcc","ngcd"))) {
     stop("Writing netcdf for grid.type= ",attr(grid,"grid.type")," not implemented.")
  }
#------------------------------------------------------------------------------
# ==> prepare the variables for writing <==
#------------------------------------------------------------------------------
  # convert grid into array if it is a matrix
  for (v in 1:dimv) {
    if (is.matrix(grid.list[[v]])) {
      if (isEns[v]) { 
        grid.list[[v]] <- array(grid.list[[v]],dim=c(dim(grid.list[[v]]),dime,dimt))
      } else {
        grid.list[[v]] <- array(grid.list[[v]],dim=c(dim(grid.list[[v]]),dimt))
      }
    }
  }
#..............................................................................
#==>> determine times for netcdf file
`date.element` <-
function(ts,format="%Y-%m-%d %H:%M:%S",nam) {
# ===========================================================================
# retrieve an element of a date representation ts as numeric
# ts can be either of class character (then format defines the format)
# or of class POSIXt, POSIXct, POSIXlt
# this is a generic function for year(), month(), day(), dofy(), dofm(), etc.
   if (!(class(ts)[1] %in% c("character","POSIXt","POSIXct","POSIXlt"))) {
      stop("** ERROR ** input is not of class character or POSIX")
   }

   doit4one <- function(ts) {
       if ("character" %in% class(ts)) {
           Xlt <- as.POSIXlt(str2Rdate(ts,format=format))
       } else {
           Xlt <- as.POSIXlt(ts)
       }
       unlist(Xlt)[nam]
   }

   qq <- sapply(ts,FUN=doit4one)
   names(qq) <- c()
   qq
}

`year` <-
function(ts,format="%Y-%m-%d") {
  date.element(ts,format=format,nam="year")+1900
}

`month` <-
function(ts,format="%Y-%m-%d") {date.element(ts,format=format,nam="mon")+1}

`str2Rdate` <-
function(ts,format="%Y-%m-%d %H:%M:%S") {
# ===========================================================================
# converts a string into an R (POSIXt,POSIXct) date object
# date objects can be used with arithmetic operations +/-
# ts is a character or a vector of characters with date information
# in the format specified in format
# Output is a date object

     # the lengthy bunch of testing is necessary because strptime needs
     # explicit specifications for month and day, otherwise it returns NA.
     # this extension allows inputs of the format "%Y-%m" in which case the
     # the first day of the month is taken as a reference.

     #Êcheck if year/month/day is specified
     ysp <- length(c(grep("%Y",format,fixed=TRUE),
                     grep("%y",format,fixed=TRUE)))
     msp <- length(c(grep("%m",format,fixed=TRUE),
                     grep("%b",format,fixed=TRUE),
                     grep("%B",format,fixed=TRUE)))
     jsp <- length(c(grep("%j",format,fixed=TRUE)))
     dsp <- length(c(grep("%d",format,fixed=TRUE)))
     if (ysp > 1) { stop("ERROR: Multiple specification of year in 
                         date format.") }
     if (ysp == 0) { stop("ERROR: No year specification in 
                         date format.") }
     if (msp > 1) { stop("ERROR: Multiple specification of month in 
                         date format.") }
     if (dsp > 1) { stop("ERROR: Multiple specification of day in 
                         date format.") }

     # append month or day if not specified
     tss <- ts
     formati <- format
     if (jsp == 0) {
     if (msp == 0) {
        tss <- paste(tss,"01",sep="")
        formati <- paste(formati,"%m",sep="")
     }
     if (dsp == 0) {
        tss <- paste(tss,"01",sep="")
        formati <- paste(formati,"%d",sep="")
     }
     }

     # this is necessary because strptime() returns NA otherwise
     as.POSIXct(strptime(tss,format=formati),tz="GMT")
}

`Rdate2str` <-
function(date,format="%Y-%m-%d %H:%M:%S") {
# ===========================================================================
# converts a R (POSIXt,POSIXct) date object into a string
# date is one or a vector of R date-time objects
# format specifies the desired character format analogous to str2Rdate
# Output is one or a vector of characters
     format(date,format=format,tz="GMT")
}


`str2nc4t` <-
function(str,t.unit,format="%Y%m%d%H%M",
         t.ref="190001010000",format.ref="%Y%m%d%H%M") {

      # seconds elapsed
      ss <- str2Rdate(str,format=format)
      ss.ref <- str2Rdate(t.ref,format=format.ref)
      secs.elapsed <- unclass(ss)-unclass(ss.ref)
      names(secs.elapsed) <- c()

      # date indicator in netcdf attribute
      fmt.nc <- "%Y-%m-%d %H:%M:00"
      ref.date <- Rdate2str(ss.ref,format=fmt.nc)

      # years, months elapsed
      yy <- year(str,format=format)
      yy.ref <- year(t.ref,format=format.ref)
      mm <- month(str,format=format)
      mm.ref <- month(t.ref,format=format.ref)
      yeas.elapsed <- yy-yy.ref
      names(yeas.elapsed) <- c()
      mons.elapsed <- (12*yy+(mm-1))-(12*yy.ref+(mm.ref-1))
      names(mons.elapsed) <- c()

      tt <- switch(t.unit,
          "S" = secs.elapsed,
          "T" = secs.elapsed/60,
          "H" = secs.elapsed/3600,
          "D" = secs.elapsed/86400,
          "M" = mons.elapsed,
          "Y" = yeas.elapsed)
      tt.attr <- switch(t.unit,
          "S" = paste("seconds since",ref.date,sep=" "),
          "T" = paste("minutes since",ref.date,sep=" "),
          "H" = paste("hours since",ref.date,sep=" "),
          "D" = paste("days since",ref.date,sep=" "),
          "M" = paste("months since",ref.date,sep=" "),
          "Y" = paste("years since",ref.date,sep=" "))
      attr(tt,"since.lab") <- tt.attr
      attr(tt,"tzone") <- NULL

      tt
}

  nctims <- str2nc4t(str=times,t.unit=times.unit,t.ref=times.ref)
#..............................................................................
#==>> define space coordinates and grid mapping variable
  if (grid.type == "lonlat") {  
    #  - case for normal lon-lat grid  -------------------------
    grid.mapping <- "longitude_latitude"
    grid.mapping.name <- "latitude_longitude"
    xx <- ncdim_def("lon","degrees",x)
    yy <- ncdim_def("lat","degrees",y)
    atts.gmv <- list(list("grid_mapping_name", grid.mapping.name, "text"),
                     list("longitude_of_prime_meridian", 0.0, "double"),
                     list("semi_major_axis", 6378137.0, "double"),
                     list("inverse_flattening", 298.257223563, "double")) 
    atts.var <- list(list("grid_mapping", grid.mapping, "text"))
    atts.xx <-  list(list("standard_name","longitude","text"),
                     list("long_name", "grid longitude", "text"))
    atts.yy <-  list(list("standard_name","latitude","text"),
                     list("long_name", "grid latitude", "text"))
  }
  if (grid.type == "utm33") {  
    #  - case for Norwegian UTM zone 33 grid  ------------------------
    grid.mapping <- "UTM_Zone_33"
    grid.mapping.name <- "transverse_mercator"
    xx <- ncdim_def("X","meters",x)
    yy <- ncdim_def("Y","meters",y)
    atts.gmv <- list(list("grid_mapping_name", grid.mapping.name, "text"),
                     list("utm_zone_number", 33, "integer"),
                     list("inverse_flattening", 298.257222101, "double"),
                     list("semi_major_axis", 6378137, "double"),
                     list("proj4", proj4.string, "text"),
                     list("_CoordinateTransformType", "Projection", "text"),
                     list("_CoordinateAxisType", "GeoX GeoY", "text") )
    atts.var <- list(list("grid_mapping", grid.mapping, "text"))
    atts.xx <-  list(list("standard_name","projection_x_coordinate","text"),
                     list("long_name", "x coordinate of projection", "text"))
    atts.yy <-  list(list("standard_name","projection_y_coordinate","text"),
                     list("long_name", "y coordinate of projection", "text"))
  }
  if (grid.type == "radlcc") {  
    #  - case for radlcc  ------------------------
    grid.mapping <- "projection_lambert"
    grid.mapping.name <- "lambert_conformal_conic"
    xx <- ncdim_def("X","meters",x)
    yy <- ncdim_def("Y","meters",y)
    atts.gmv <- list(list("grid_mapping_name", grid.mapping.name, "text"),
                     list("standard_parallel", c(56.6, 69.7), "double"),
                     list("longitude_of_central_meridian", 15., "double"),
                     list("semi_major_axis", 6378137, "double"),
                     list("semi_minor_axis", 6356752.3142, "double"),
                     list("proj4", proj4.string, "text"),
                     list("x0","-790000", "text"),
                     list("y0","7540000", "text"),
                     list("dx","1000", "text"),
                     list("dy","1000", "text"),
                     list("_CoordinateTransformType", "Projection", "text"),
                     list("_CoordinateAxisType", "GeoX GeoY", "text") )
    atts.var <- list(list("grid_mapping", grid.mapping, "text"))
    atts.xx <-  list(list("standard_name","projection_x_coordinate","text"),
                     list("long_name", "x coordinate of projection", "text"))
    atts.yy <-  list(list("standard_name","projection_y_coordinate","text"),
                     list("long_name", "y coordinate of projection", "text"))
  }
  if (grid.type == "yrlcc") {  
    #  - case for yr lcc  ------------------------
    grid.mapping <- "projection_lambert"
    grid.mapping.name <- "lambert_conformal_conic"
    xx <- ncdim_def("x","m",x)
    yy <- ncdim_def("y","m",y)
    atts.gmv <- list(list("grid_mapping_name", grid.mapping.name, "text"),
                     list("standard_parallel", c(63.,63.), "double"),
                     list("longitude_of_central_meridian", 15., "double"),
                     list("latitude_of_projection_origin", 63., "double"),
                     list("earth_radius", 6371000, "double"),
                     list("proj4", proj4.string, "text") )
    atts.var <- list(list("grid_mapping", grid.mapping, "text"))
    atts.xx <-  list(list("standard_name","projection_x_coordinate","text"),
                     list("long_name", "x-coordinate in Cartesian system", "text"))
    atts.yy <-  list(list("standard_name","projection_y_coordinate","text"),
                     list("long_name", "y-coordinate in Cartesian system", "text"))
  }
  if (grid.type == "ngcd") {  
    #  - case for ngcd 
    # http://spatialreference.org/ref/epsg/etrs89-etrs-laea/
    grid.mapping <- "projection_laea"
    grid.mapping.name <- "lambert_azimuthal_equal_area"
    xx <- ncdim_def("X","meters",x)
    yy <- ncdim_def("Y","meters",y)
    atts.gmv <- list(list("grid_mapping_name", grid.mapping.name, "text"),
                     list("datum","European_Terrestrial_Reference_System_1989","text"),
                     list("spheroid","GRS 1980","text"),
                     list("semi_major_axis", 6378137, "double"),
                     list("reference","http://spatialreference.org/ref/epsg/etrs89-etrs-laea/","text"),
                     list("proj4", proj4.string, "text"),
                     list("latitude_of_center",52,"double"),
                     list("longitude_of_center",10,"double"),
                     list("false_easting",4321000,"double"),
                     list("false_northing",3210000,"double"),
                     list("EPSG", "3035", "text") )
    atts.var <- list(list("grid_mapping", grid.mapping, "text"))
    atts.xx <-  list(list("standard_name","projection_x_coordinate","text"),
                     list("long_name", "x coordinate of projection", "text"))
    atts.yy <-  list(list("standard_name","projection_y_coordinate","text"),
                     list("long_name", "y coordinate of projection", "text"))
  }
  # calculate longitute and latitude fields if the grid is not in lonlat
  # and define variables and attributes for these additional fields
  if (lonlat.out & (grid.type != "lonlat")) {
    xy<-expand.grid(x,y)
    names(xy)<-c("x","y")
    coordinates(xy)<-c("x","y")
    proj4string(xy)<-proj4.string
    # transformation is actually done by using the gdal library (www.gdal.org)
    gg<-spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
    lon<-attributes(gg)$coords[,1]
    lon<-round(lon,6)
    lat<-attributes(gg)$coords[,2]
    lat<-round(lat,6)
#     lon <- gg$x 
#     lat <- gg$y
    lon.var <- ncvar_def(name="lon",units="degrees_east",missval=NA,
                         dim=list(xx,yy),prec="single",
                         compression=compr.lev)
    atts.lon <- list(list("long_name","longitude coordinate","text"),
                     list("standard_name","longitude","text")) 
    lat.var <- ncvar_def(name="lat",units="degrees_north",missval=NA,
                         dim=list(xx,yy),prec="single",
                         compression=compr.lev)
    atts.lat <- list(list("long_name","latitude coordinate","text"),
                     list("standard_name","latitude","text"))
    atts.var <- c(atts.var,list(list("coordinates", "lon lat", "text")))
  }
#..............................................................................
#==>> define time coordinate
  tim <- ncdim_def("time",attr(nctims,"since.lab"),nctims,unlim=TRUE)
  atts.tim <- list(list("axis", "T", "text"), 
                   list("calendar", "standard", "text"), 
                   list("standard_name", "time", "text"), 
                   list("long_name", "time", "text"))
#..............................................................................
#==>> define ensemble number as a coordinate
  if (any(isEns)) {
    ens <- ncdim_def("ensemble_member","",1:dime)
    atts.ens <- list( list("_CoordinateAxisType","Ensemble","text"),
                      list("long_name", "ensemble run number", "text"),
                      list("standard_name", "realization", "text"),
                      list("_ChunkSizes", "10", "integer"))
  }
#..............................................................................
#==>> define the grid mapping variable
  # gmvdim <- ncdim_def(name="dummy",units="",vals=c(1))
  gmvdim<-list()
  gmv <- ncvar_def(grid.mapping,"",dim=gmvdim,missval=-1,prec="double",
                   compression=compr.lev)
#..............................................................................
#==>> define the forecast reference time variable
  if (for_ref_time.out) { 
    forecast.reference.time<-"forecast_reference_time"
    atts.frtv <- list(list("standard_name","forecast_reference_time","text"))
    frtvdim<-list()
    frtv <- ncvar_def(forecast.reference.time,attr(nctims,"since.lab"),
                      dim=frtvdim,missval=-1,
                      prec="double",compression=compr.lev)
  }
#..............................................................................
#==>> define the variables and their attributes
  var.list<-list()
#  atts.var.list<-list()
  for (v in 1:dimv) {
    if (isEns[v]) {
      pre<-ncvar_def(name=var.name[v],units=var.unit[v],
                     dim=list(xx,yy,ens,tim),missval=var.mv,prec="single",
                     compression=compr.lev)
    } else {
      pre<-ncvar_def(name=var.name[v],units=var.unit[v],
                     dim=list(xx,yy,tim),missval=var.mv,prec="single",
                     compression=compr.lev)
    }
    var.list[[v]]<-pre
  }
#..............................................................................
#==>> define global attributes
  atts.glob <- list(list("institution",
                         "Norwegian Meteorological Institute, MET Norway", 
                         "text"),
                    list("creator_url","met.no","text"),
                    list("summary",summary,"text"),
                    list("title",title,"text"),
                    list("license",license,"text"),
                    list("source",source.string,"text"),
                    list("References", reference, "text"),
                    list("comment", comment, "text"),
                    list("Conventions", "CF-1.6", "text"))
#------------------------------------------------------------------------------
# ==> write now <==
#------------------------------------------------------------------------------
# function to add attributes
  add.att.to.nc <-
  function(list,var,...) {
     ncatt_put(nc=nc.con,varid=var,
               attname=list[[1]],attval=list[[2]],prec=list[[3]],...)
     return()
  }
#..............................................................................
#==>> create file, write variables and add additional attributes, close the file
  var.list[[dimv+1]]<-gmv
  dimv1<-dimv+1
  if (exists("lon.var") & exists("lat.var")) {
    var.list[[dimv+2]]<-lon.var
    var.list[[dimv+3]]<-lat.var
    dimv1<-dimv+3
  } 
  if (for_ref_time.out) {
    var.list[[dimv1+1]]<-frtv
    dimv1<-dimv1+1
  }
#..............................................................................
  nc.con <- nc_create(filename=file.name,vars=var.list,force_v4=T)
  for (v in 1:dimv) {
    atts.var1<-c(atts.var,
                 list(list("long_name", var.longname[v], "text"),
                      list("standard_name", var.standardname[v], "text"),
#                      list("_FillValue", var.mv, "single"),
                      list("version", var.version, "text"),
                      list("prod_date", as.character(prod.date), "text")))
    if (length(atts.var.add)>=v) {
      if (!is.null(atts.var.add[[v]])) {
        atts.var1<-c(atts.var1,list(atts.var.add[[v]]))
      }
    }
    # variable attributes
    if (length(atts.var1)>0) 
      hhh <- lapply(atts.var1,FUN=add.att.to.nc,var=var.list[[v]])
  }  
  # time attributes
  if (length(atts.tim)>0) 
    hhh <- lapply(atts.tim,FUN=add.att.to.nc,var="time")
  # ensemble attributes
  if (any(isEns)) {
    if (length(atts.ens)>0) 
      hhh <- lapply(atts.ens,FUN=add.att.to.nc,var="ensemble_member")
  }
  # attributes for grid mapping variable
  if (length(atts.gmv)>0) 
    hhh <- lapply(atts.gmv,FUN=add.att.to.nc,var=gmv)
  # global attributes
  if (length(atts.glob)>0) 
    hhh <- lapply(atts.glob,FUN=add.att.to.nc,var=0)
  # x-coordinate attributes
  if (length(atts.xx)>0) 
    hhh <- lapply(atts.xx,FUN=add.att.to.nc,var=xx$name)
  # y-coordinate attributes
  if (length(atts.yy)>0) 
    hhh <- lapply(atts.yy,FUN=add.att.to.nc,var=yy$name)
  # forecast reference time attributes
  if (for_ref_time.out)  
    hhh <- lapply(atts.frtv,FUN=add.att.to.nc,var=frtv)
  #
  for (v in 1:dimv) {
    if (!is.null(round.dig)) 
      grid.list[[v]]<-round(grid.list[[v]],round.dig)
    ncvar_put(nc=nc.con,varid=var.list[[v]],vals=grid.list[[v]])
  }
  #
  ncvar_put(nc=nc.con,varid=gmv,vals=1)
  #
  if (exists("lon.var") & exists("lat.var")) {
    ncvar_put(nc=nc.con,varid=lon.var,vals=lon)
    ncvar_put(nc=nc.con,varid=lat.var,vals=lat)
    # attributes for lon coordinate
    if (length(atts.lon)>0) 
      hhh <- lapply(atts.lon,FUN=add.att.to.nc,var=lon.var)  
    # attributes for lat coordinate
    if (length(atts.lat)>0) 
      hhh <- lapply(atts.lat,FUN=add.att.to.nc,var=lat.var)  
  }
  if (for_ref_time.out) {
    ncfrtims <- str2nc4t(str=for_ref_time.val,
                         t.unit=times.unit,t.ref=times.ref)
    ncvar_put(nc=nc.con,varid=frtv,vals=ncfrtims)
  }
  #
  nc_close(nc.con)
  # return NULL
  NULL
}

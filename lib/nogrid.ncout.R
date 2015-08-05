`nogrid.ncout` <- 
function(
         # A grid object, i.e. a matrix/array (possibly for several times, i.e. 3D)
         grid,
         
         # The type of the grid: either "utm33" or "lonlat"
         grid.type="utm33",

         # The coordinate vectors (x: easting or longitude, y: northing or latitude)
         x,
         y,

         # The file name to write. A character string
         file.name="dummy.nc",
         
         # The variable name to write. A character string
         var.name="DUMMY",
         
         # The long name of the variable. Check out long-names
         var.longname = var.name, 
         
         # The variable's units, a character string
         var.unit = NULL,
         
         # The missval to use (NAs in R are replaced by this number )
         var.mv = -999.99,
         
         # A version indicator
         var.version = "no version",
         
         # The time(s) of the grid. A character vector of convention "yyyymmddhhmm"
         # If not specified, the function searches for a "time" attribute in grid
         # to get the information
         times = NULL,
         
         # The time unit to use for the time coordinate in the file. 
         # Possible options are: "Y" for annual, "M" for months, "D" for days,
         # "H" for hours and "T" for minutes. 
         times.unit,
         # The reference time to use for the time coordinate. 
         # Time units in the file are given as the elapsed time since \code{t.ref}.
         times.ref = "190001010000",

         # A production date/time to be written into the file
         prod.date = substr(Sys.time(),1,10),

         # proj4.string
         proj4.string = "",

         # source
         source.string="",

         # A reference to a paper or documentation of the dataset.
         reference = "") {


   # load netcdf library
   require(ncdf)


   # -----------------------------------------
   # functions for date/time handling (stripped from package "datefuns" at MCH)
   # -----------------------------------------
   `str2nct` <- 
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
             "T" = secs.elapsed/60,
             "H" = secs.elapsed/3600,
             "D" = secs.elapsed/86400,
             "M" = mons.elapsed,
             "Y" = yeas.elapsed)
         tt.attr <- switch(t.unit,
             "T" = paste("minutes since",ref.date,sep=" "),
             "H" = paste("hours since",ref.date,sep=" "),
             "D" = paste("days since",ref.date,sep=" "),
             "M" = paste("months since",ref.date,sep=" "),
             "Y" = paste("years since",ref.date,sep=" "))
         attr(tt,"since.lab") <- tt.attr
         attr(tt,"tzone") <- NULL

         tt
   }

   `str2Rdate` <-
   function(ts,format="%Y-%m-%d %H:%M:%S") {
   # converts a string into an R (POSIXt,POSIXct) date object
   # date objects can be used with arithmetic operations +/-
   # ts is a character or a vector of characters with date information
   # in the format specified in format
   # Output is a date object

        # the lengthy bunch of testing is necessary because strptime needs
        # explicit specifications for month and day, otherwise it returns NA.
        # this extension allows inputs of the format "%Y-%m" in which case the
        # the first day of the month is taken as a reference.

        # check if year/month/day is specified
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
   # converts a R (POSIXt,POSIXct) date object into a string
   # date is one or a vector of R date-time objects
   # format specifies the desired character format analogous to str2Rdate
   # Output is one or a vector of characters
        format(date,format=format,tz="GMT")
   }

   `year` <-
   function(ts,format="%Y-%m-%d") {date.element(ts,format=format,nam="year")+1900}
   `month` <-
   function(ts,format="%Y-%m-%d") {date.element(ts,format=format,nam="mon")+1}
   `day` <-
   function(ts,format="%Y-%m-%d") {date.element(ts,format=format,nam="mday")}

   `date.element` <-
   function(ts,format="%Y-%m-%d %H:%M:%S",nam) {
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


   # -----------------------------------------
   # checks and default handling
   # -----------------------------------------
   if (missing(times.unit)) {
      stop("Argument \"times.unit\" must be specified.")
   }
   if (is.null(times)) {
      times <- attr(grid,"time")
      if (is.null(times)) {
         stop("No argument \"times\" and no attribute \"time\" found in \"grid\".")
      }
   }
   if (!is.character(times)) {
      stop("Inappropriate class of argument/attribute \"times\".")
   }
   if (is.matrix(grid)) { dim3 <- 1 } else { dim3 <- dim(grid)[3] }
   if (length(times) != dim3) {
      stop("Length of \"times\" inconsistent with dimensions of \"grid\".")
   }
   if (length(x) != dim(grid)[1]) {
      stop("Length of \"x\" inconsistent with dimensions of \"grid\".")
   }
   if (length(y) != dim(grid)[2]) {
      stop("Length of \"y\" inconsistent with dimensions of \"grid\".")
   }
   if (!is.character(var.unit)) {
       warning("No units specified for grid. <dimensionless> chosen instead.")
       var.unit <- "no dim"
   }
   # check if grid.type is available
   if (!(grid.type %in% c("lonlat","utm33","laea"))) {
      stop("Writing netcdf for grid.type= ",attr(grid,"grid.type")," not implemented.")
   }

   
   
   
   # -----------------------------------------
   # prepare the variables for writing 
   # -----------------------------------------
   
   # convert grid into array if it is a matrix
   if (dim3 == 1) {
      grid <- array(grid,dim=c(dim(grid),dim3))
   }

   # determine times for netcdf file
   nctims <- str2nct(str=times,t.unit=times.unit,t.ref=times.ref)


   # define space coordinates and grid mapping variable
   if (grid.type == "lonlat") {  
     #  - case for normal lon-lat grid  -------------------------
     grid.mapping <- "longitude_latitude"
     grid.mapping.name <- "latitude_longitude"
     xx <- dim.def.ncdf("lon","degrees",x)
     yy <- dim.def.ncdf("lat","degrees",y)
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
     xx <- dim.def.ncdf("X","meters",x)
     yy <- dim.def.ncdf("Y","meters",y)
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

   if (grid.type == "laea") {  
     #  - case for Norwegian UTM zone 33 grid  ------------------------
     grid.mapping <- "Lambert_Azimuthal_Equal_Area"
     grid.mapping.name <- "lambert_azimuthal_equal_area"
     xx <- dim.def.ncdf("X","meters",x)
     yy <- dim.def.ncdf("Y","meters",y)
     atts.gmv <- list(list("grid_mapping_name", grid.mapping.name, "text"),
                      list("projection", proj4.string, "text"),
                      list("projection_format", "PROJ.4", "text") )
     atts.var <- list(list("grid_mapping", grid.mapping, "text"))
     atts.xx <-  list(list("standard_name","projection_x_coordinate","text"),
                      list("long_name", "x coordinate of projection", "text"))
     atts.yy <-  list(list("standard_name","projection_y_coordinate","text"),
                      list("long_name", "y coordinate of projection", "text"))
   }
   
   # calculate longitute and latitude fields if the grid is not in lonlat
   # and define variables and attributes for these additional fields
   # (*** the following would need to be completed ***)
   #if (grid.type != "lonlat") {
   #  lon <- 
   #  lat <- 
   #  lon.var <- var.def.ncdf(name="lon",units="degrees_east",missval=NA,
   #                          dim=list(xx,yy),prec="single")
   #  atts.lon <- list(list("long_name","longitude coordinate","text"),
   #                   list("standard_name","longitude","text")) 
   #  lat.var <- var.def.ncdf(name="lat",units="degrees_north",missval=NA,
   #                          dim=list(xx,yy),prec="single")
   #  atts.lat <- list(list("long_name","latitude coordinate","text"),
   #                   list("standard_name","latitude","text"))
   #  atts.var <- c(atts.var,list(list("coordinates", "lon lat", "text")))
   #}


   # define time coordinate
   tim <- dim.def.ncdf("time",attr(nctims,"since.lab"),nctims,unlim=TRUE)
   atts.tim <- list(list("axis", "T", "text"), 
                    list("calendar", "standard", "text"), 
                    list("long_name", "time", "text"))

   # define the grid mapping variable
   gmvdim <- dim.def.ncdf(name="dummy",units="",vals=c(1))
   gmv <- var.def.ncdf(grid.mapping,"",dim=gmvdim,missval=-1,prec="double")

   # define the variable and its attributes
   pre <- var.def.ncdf(name=var.name,units=var.unit,
                       dim=list(xx,yy,tim),missval=var.mv,prec="single")
   atts.var <- c(atts.var,
                 list(list("long_name", var.longname, "text"),
                      list("_FillValue", var.mv, "single"),
                      list("version", var.version, "text"),
                      list("prod_date", as.character(prod.date), "text")))

   # define global attributes
   atts.glob <- list(list("Conventions", "CF-1.4", "text"),
                     list("institution",
                          "Norwegian Meteorological Institute, met.no", 
                          "text"),
                     list("source",source.string,"text"),
                     list("References", reference, "text"))

   
   # -----------------------------------------
   # write now 
   # -----------------------------------------

   # function to add attributes
   add.att.to.nc <- function(list,var,...) { 
                            att.put.ncdf(nc=nc.con,varid=var,
                                         attname=list[[1]],attval=list[[2]],prec=list[[3]],...)
                            return() }
   
   # create file, write variables and add additional attributes, close the file
   if (exists("lon.var") & exists("lat.var")) {
      nc.con <- create.ncdf(filename=file.name,vars=list(pre,gmv,lon.var,lat.var))
   } else {
      nc.con <- create.ncdf(filename=file.name,vars=list(pre,gmv))
   }
   put.var.ncdf(nc=nc.con,varid=pre,vals=grid)
   put.var.ncdf(nc=nc.con,varid=gmv,vals=1)
   if (length(atts.var)>0) {hhh <- lapply(atts.var,FUN=add.att.to.nc,var=pre)}  # variable attributes
   if (length(atts.tim)>0) {hhh <- lapply(atts.tim,FUN=add.att.to.nc,var="time")}  # time attributes
   if (length(atts.gmv)>0) {hhh <- lapply(atts.gmv,FUN=add.att.to.nc,var=gmv)}  # attributes for grid mapping variable
   if (length(atts.glob)>0) {hhh <- lapply(atts.glob,FUN=add.att.to.nc,var=0)}   # global attributes
   if (length(atts.xx)>0) {hhh <- lapply(atts.xx,FUN=add.att.to.nc,var=xx$name)}   # x-coordinate attributes
   if (length(atts.yy)>0) {hhh <- lapply(atts.yy,FUN=add.att.to.nc,var=yy$name)}   # y-coordinate attributes
   if (exists("lon.var") & exists("lat.var")) {
      put.var.ncdf(nc=nc.con,varid=lon.var,vals=lon)
      put.var.ncdf(nc=nc.con,varid=lat.var,vals=lat)
      if (length(atts.lon)>0) {hhh <- lapply(atts.lon,FUN=add.att.to.nc,var=lon.var)}  # attributes for lon coordinate
      if (length(atts.lat)>0) {hhh <- lapply(atts.lat,FUN=add.att.to.nc,var=lat.var)}  # attributes for lat coordinate
   }
   close.ncdf(nc.con)
      
   # return NULL
   NULL

}

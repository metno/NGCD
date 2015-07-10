
load(file="example.rda")


grid.type <- "utm33"
file.name <- "example.nc"
var.name <- "thefield"
var.longname <- "the_field_that_shows_the_distribution"
var.unit <- "mm"
var.mv <- -999.99
var.version <- "no version"
times <- c("201312010300")
times.unit <- "H"
times.ref <- "200101010000"
prod.date <- substr(Sys.time(),1,10)
reference <- "Lussana C. and O.E. Tveito, 2014: The new hourly grids for Norway. J. ..."

source("nogrid.ncout.R")
nogrid.ncout(grid=grid,x=x,y=y,file.name="example.nc",var.name="thefield",
             var.longname="the_field",var.unit="mm",var.mv=-999.9,
             var.version="no_version",times=c("201312010300"),times.unit="H",
             times.ref="200101010000",
             reference="Lussana C. and O.E. Tveito, 2014: The new ..")


# try reducing file size by rounding numbers
grid1 <- round(grid,digits=1)
nogrid.ncout(grid=grid1,x=x,y=y,file.name="example_rounded1.nc",var.name="thefield",
             var.longname="the_field",var.unit="mm",var.mv=-999.9,
             var.version="no_version",times=c("201312010300"),times.unit="H",
             times.ref="200101010000",
             reference="Lussana C. and O.E. Tveito, 2014: The new ..")
grid2 <- round(grid,digits=2)
nogrid.ncout(grid=grid2,x=x,y=y,file.name="example_rounded2.nc",var.name="thefield",
             var.longname="the_field",var.unit="mm",var.mv=-999.9,
             var.version="no_version",times=c("201312010300"),times.unit="H",
             times.ref="200101010000",
             reference="Lussana C. and O.E. Tveito, 2014: The new ..")

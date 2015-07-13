# << vedo_TEMP1d.R >>
# Create maps for daily mean temperature.
#==============================================================================
rm(list=ls())
# Libraries
library(raster)
library(rgdal)
library(ncdf)
#------------------------------------------------------------------------------
# Graphical parameter
xlim.sw<-3975220
xlim.ne<-5672220
ylim.sw<-3576520
ylim.ne<-5462520
# South Norway only
#xlim.sw<--75000
#xlim.ne<-500000
#ylim.sw<-6450000
#ylim.ne<-7150000
#
scale.TEMP1d<-c(-1000,-35, -30,-25,-20,-15,-10,
                 -7.5, -5,-2.5,  0,  1,  3,  5,
                  7.5, 10, 12, 14, 16, 18,20,22,
                   24, 26, 28, 30, 32, 34,36,38,1000)
#==============================================================================
# read command line
arguments <- commandArgs()
arguments
yyyy.mm.dd<-arguments[3]
config_file<-arguments[4]
config_par<-arguments[5]
if (length(arguments)!=5) {
  print("Error in command line arguments:")
  print("R --vanilla yyyy.mm.dd config_file config_par")
  quit(status=1)
}
# define/check paths
if (!file.exists(config_file)) 
  ext<-error_exit("Fatal Error: configuration file not found")
source(config_file)
for (p in 1:length(config_list)) {
  if (config_list[[p]]$pname==config_par) break
}
if (p==length(config_list) & (config_list[[p]]$pname!=config_par) )  
  ext<-error_exit("Fatal Error: configuration parameter not in configuration file")
main.path<-config_list[[p]]$opt$main.path
main.path.output<-config_list[[p]]$opt$main.path.output
testmode<-config_list[[p]]$opt$testmode
if ( !(file.exists(main.path)) | !(file.exists(main.path.output)) ) 
  ext<-error_exit("Fatal Error: path not found")
# geographical information
main.path.geoinfo<-paste(main.path,"/geoinfo",sep="")
filenamedem<-paste(main.path.geoinfo,"/fennodem_ETRS89_ETRS-LAEA.nc",sep="")
fileborders<-paste(main.path.geoinfo,"/TM_WORLD_BORDERS_ETRS89_ETRS-LAEA/TM_WORLD_BORDERS_ETRS89_ETRS-LAEA-0.2.shp",sep="")
if (!file.exists(paste(main.path.geoinfo,"/fennodem_ETRS89_ETRS-LAEA.nc",sep=""))) 
  ext<-error_exit(paste("File not found:",main.path.geoinfo,"/fennodem_ETRS89_ETRS-LAEA.nc"))
# common libs and etcetera
path2lib.com<-paste(main.path,"/lib",sep="")
path2etc.com<-paste(main.path,"/etc",sep="")
if (!file.exists(paste(path2lib.com,"/Bspat_plot.R",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/Bspat_plot.R"))
# call to external function
source(paste(path2lib.com,"/Bspat_plot.R",sep=""))
# set Time-related variables
yyyy<-substr(yyyy.mm.dd,1,4)
mm<-substr(yyyy.mm.dd,6,7)
dd<-substr(yyyy.mm.dd,9,10)
yyyymm<-paste(yyyy,mm,sep="")
yyyymmdd<-paste(yyyy,mm,dd,sep="")
# Color tables
banded<-read.table(file=paste(path2etc.com,"/color_table/NCV_banded.rgb",sep=""),skip=2,stringsAsFactors=F)
banded.r<-as.numeric(banded$V1)
banded.g<-as.numeric(banded$V2)
banded.b<-as.numeric(banded$V3)
banded.col<-rgb(banded.r,banded.g,banded.b,maxColorValue = 256)
rainbow<-read.table(file=paste(path2etc.com,"/color_table/NCV_rainbow2.rgb",sep=""),skip=2,stringsAsFactors=F)
rainbow.r<-as.numeric(rainbow$V1)
rainbow.g<-as.numeric(rainbow$V2)
rainbow.b<-as.numeric(rainbow$V3)
rainbow.col<-rgb(rainbow.r,rainbow.g,rainbow.b,maxColorValue = 256)
stepseq<-read.table(file=paste(path2etc.com,"/color_table/MPL_StepSeq.rgb",sep=""),skip=2,stringsAsFactors=F)
stepseq.r<-as.numeric(stepseq$V1)
stepseq.g<-as.numeric(stepseq$V2)
stepseq.b<-as.numeric(stepseq$V3)
stepseq.col<-rev(rgb(stepseq.r,stepseq.g,stepseq.b,maxColorValue = 1))
stepseq25<-read.table(file=paste(path2etc.com,"/color_table/StepSeq25.rgb",sep=""),skip=7,stringsAsFactors=F)
stepseq25.r<-as.numeric(stepseq25$V1)
stepseq25.g<-as.numeric(stepseq25$V2)
stepseq25.b<-as.numeric(stepseq25$V3)
stepseq25.col<-rev(rgb(stepseq25.r,stepseq25.g,stepseq25.b,maxColorValue = 256))
aux<-stepseq25.col[11:15]
stepseq25.col[11:15]<-stepseq25.col[1:5]
stepseq25.col[1:5]<-aux
t2m_29lev<-read.table(file=paste(path2etc.com,"/color_table/t2m_29lev.rgb",sep=""),skip=6,stringsAsFactors=F)
t2m_29lev.r<-as.numeric(t2m_29lev$V1)
t2m_29lev.g<-as.numeric(t2m_29lev$V2)
t2m_29lev.b<-as.numeric(t2m_29lev$V3)
t2m_29lev.col<-rgb(t2m_29lev.r,t2m_29lev.g,t2m_29lev.b,maxColorValue = 256)
# input directories
path2input.1d.main<-paste(main.path.output,"/NGCD/TEMP1d",sep="")
path2input.1d.main.stn<-paste(path2input.1d.main,"/station_dataset",sep="")
path2input.1d.main.grd<-paste(path2input.1d.main,"/gridded_dataset",sep="")
path2input.1d.add<-paste(main.path.output,"/NGCD_addInfo/TEMP1d",sep="")
path2input.1d.add.grd<-paste(path2input.1d.add,"/gridded_dataset",sep="")
path2input.1d.add.eve<-paste(path2input.1d.add,"/event_dataset",sep="")
in.1d.file.stn<- paste(path2input.1d.main.stn,"/",yyyymm,
                       "/NGCD_TEMP1d_station_",
                       yyyymmdd,".txt",sep="")
in.1d.file.grd.ana<- paste(path2input.1d.main.grd,"/",yyyymm,
                           "/NGCD_TEMP1d_grid_",
                           yyyymmdd,".nc",sep="")
in.1d.file.grd.idi<- paste(path2input.1d.add.grd,"/",yyyymm,
                           "/NGCD_TEMP1d_grid_idi_",
                           yyyymmdd,".nc",sep="")
# output directories
dir.create(file.path(main.path.output,"NGCD_scratch"), showWarnings = FALSE)
path2output.main<-paste(main.path.output,"/NGCD_scratch/TEMP1d",sep="")
path2output.main.grd<-paste(path2output.main,"/TEMP",sep="")
path2output.add<-paste(main.path.output,"/NGCD_scratch/TEMP1d",sep="")
path2output.add.grd<-paste(path2output.add,"/IDI",sep="")
if (!(file.exists(path2output.main)))     dir.create(path2output.main,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) dir.create(path2output.main.grd,showWarnings=F) 
if (!(file.exists(path2output.add)))      dir.create(path2output.add,showWarnings=F) 
if (!(file.exists(path2output.add.grd)))  dir.create(path2output.add.grd,showWarnings=F) 
# output files 
dir.create(paste(path2output.main.grd,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm,sep=""),showWarnings=F)
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm,
                         "/NGCD_TEMP1d_grid_",
                         yyyymmdd,".png",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm,
                         "/NGCD_TEMP1d_grid_idi_",
                         yyyymmdd,".png",sep="")
#
print("Output files:")
print("analysis on the grid (netcdf)")
print(out.file.grd.ana)
print("event-normalized idi on the grid (netcdf)")
print(out.file.grd.idi)
#print("station outputs (text)")
#print(out.file.stn)
#print("event outputs (text)")
#print(out.file.eve)
#------------------------------------------------------------------------------
# read geographical info
orog<-raster(filenamedem)
borders<-readOGR(fileborders,"TM_WORLD_BORDERS_ETRS89_ETRS-LAEA-0.2")
# open/read/close netcdf file
nc <- open.ncdf(in.1d.file.grd.ana)
nc.idi <- open.ncdf(in.1d.file.grd.idi)
data <- get.var.ncdf(nc)
data.idi <- get.var.ncdf(nc.idi)
aux<-att.get.ncdf(nc,"Lambert_Azimuthal_Equal_Area","projection")
projstr<-aux$value
dx<-nc$dim$X$vals[2]-nc$dim$X$vals[1]
ex.xmin<-min(nc$dim$X$vals)-dx/2
ex.xmax<-max(nc$dim$X$vals)+dx/2
dy<-nc$dim$Y$vals[2]-nc$dim$Y$vals[1]
ex.ymin<-min(nc$dim$Y$vals)-dy/2
ex.ymax<-max(nc$dim$Y$vals)+dy/2
nx<-nc$dim$X$len
ny<-nc$dim$Y$len
close.ncdf(nc)
close.ncdf(nc.idi)
# Define raster variable "r"
r <-raster(ncol=nx, nrow=ny,
           xmn=ex.xmin, xmx=ex.xmax,
           ymn=ex.ymin, ymx=ex.ymax,
           crs=projstr)
r[]<-NA
# put data on raster variable (t=transpose)
r[]<-t(data)
data<-extract(r,1:ncell(r))
aux<-which(data==0)
r[aux]<-rep(NA,length(aux))
# read station data
y.data<-read.table(file=in.1d.file.stn,sep=";",header=T,stringsAsFactors=F)
y.data$dqcflag<-as.integer(y.data$dqcflag)
#year;month;day;nday;stid;x;y;z;eve.lab;yo;yb;ya;yav;yidi;yidiv;dqcflag;
#------------------------------------------------------------------------------
# Plot analysis 
par.TEMP1d<-list(col.scale=t2m_29lev.col,scale=scale.TEMP1d,
                 main=paste(yyyy.mm.dd,"TEMP1d","daily mean temperature [06-06 UTC]",sep=" - "),
                 xlab="",ylab="",
                 xl=c(xlim.sw,xlim.ne),yl=c(ylim.sw,ylim.ne))
#
plot<-TEMPplot(namefileout=out.file.grd.ana,
               y.data=y.data,
               r.data=r,
               orog=orog,
               bound=borders,
               par=par.TEMP1d)
# Plot IDI
scale.TEMP1h.IDI<-seq(0,110,length=257)
scale.TEMP1h.IDI[length(scale.TEMP1h.IDI)]<-1000
#
r[]<-NA
# put data on raster variable (t=transpose)
r[]<-t(data.idi)
data.idi<-extract(r,1:ncell(r))
aux<-which(data.idi==0)
r[aux]<-rep(NA,length(aux))
par.TEMP1h.IDI<-list(col.scale=banded.col,scale=scale.TEMP1h.IDI,
                 main=paste(yyyy.mm.dd,"TEMP1d/IDI","daily mean temperature [UTC]",sep=" - "),
                 xlab="",ylab="",
                 xl=c(xlim.sw,xlim.ne),yl=c(ylim.sw,ylim.ne))
#
plot<-TEMPplot.IDI(namefileout=out.file.grd.idi,
                   y.data=y.data,
                   r.data=r,
                   orog=orog,
                   bound=borders,
                   par=par.TEMP1h.IDI)
#------------------------------------------------------------------------------
# Exit
quit(status=0)

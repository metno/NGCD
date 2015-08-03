# --~- Bspat_PREC1d.R -~--
# Bayesian Spatial Interpolation of daily cumulated precipitation
#..............................................................................
# == Command line ==
#  + $>R --vanilla yyyy.mm.dd yyyy.mm.dd blacklist errobs cfg_file cfg_par
#  + time stamps mark begin/end of the accumulation period 
#  + blacklist: station blacklist
#  + errobs: suspect/erroneous observations (external DQC)
#  + cfg_file: configuration file
#  + cfg_par: configuration parameter name within the config file
# == Time specification ==
# Timezone is UTC: hour [0,23]. timestamp = end of accumulation period.
#  (i.e. 2014/09/01 12 -> precipitation sum 2014/09/01 11:01 2014/09/01 12:00)
# == Grid specifications ==
# 2 Grids:
# 1. FG -> Fine Grid, which is the high-resolution (1Km) output grid masked 
#          on Norwegian mainland
# 2. CG -> Coarse Grid, which is the low-resolution (7Km) grid larger than the 
#          Norwegian mainland
# == OUTPUT:  netCDF files + text files ==
# @@@@ in the code marks output sessions
# directory tree:
# +analysis/prediction at gridpoints (FG)
# | seNorge2->PREC1d->gridded_dataset->yyyymm
# | seNorge_v2_0_PREC1d_grid_"yyyymmdd.b"_"yyyymmdd.e".txt
# +analysis/prediction at station locations
# | seNorge2->PREC1d->station_dataset->yyyymm
# | seNorge_v2_0_PREC1d_station_"yyyymmdd.b"_"yyyymmdd.e".txt
# +IDI Integral Data Influence
# | seNorge2_addInfo->PREC1d->gridded_dataset->yyyymm
# | seNorge_v2_0_PREC1d_grid_normIDI_"yyyymmdd.b"_"yyyymmdd.e".txt
# +additional informations on events
# | seNorge2_addInfo->PREC1d->event_dataset->yyyymm
# | seNorge_v2_0_PREC1d_event_"yyyymmdd.b"_"yyyymmdd.e".txt
#
# 1. analysis/prediction at station locations
# see the output sessione marked with  
# HEADER HEADER HEADER HEADER HEADER HEADER HEADER HEADER
# ASCII file with data in columns separated by ";"
# description:
# year/month/day hour => timestamp
# nhour => number of hours accumulated
# stid => station identifier
# x,y => geographical coordinates (UTM 33, Km)
# z => elevation [m]
# eve => Contiguous Rain Area
# yo => observation [mm]
# yb => Background value used in the last iteration [mm]
# ya => analysis [mm] (prediction)
# yav => Cross-Validated (CV) analysis [mm]
# yidi => Integral Data Influence (IDI) (accumulated and normalized)
# yidiv => CV IDI (accumulated and normalized)
# dqcflag => Data Quality Control (DQC) Flag
#   NA  observation is NA
#   -1  missing DQC info
#   0   good observation
#   100 bad  observation: KDVH flag > 2 | observation not good in external DQC | 
#                         observed value not plausible | station in blacklist/s  
#                         # available observations != # time steps (in case of accumulation) 
#   200 bad  observation: dry-station surrounded only by wet-stations (close enough)
#   300 bad  observation: wet-stations surrounded only by dry-stations (close enough)
#   400 bad  observation: dry observation is (1) not included in a dry area
#                         (2) is in Norway 
#   500 bad  observation: wet observation is (1) not included in an event (2) in Norway
#
# 2. additional informations on events
# see the output sessione marked with  
# HEADER HEADER HEADER HEADER HEADER HEADER HEADER HEADER
# ASCII file with data in columns separated by ";"
# description:
# year/month/day hour => timestamp
# nhour => number of hours accumulated
# eve.lab => eve label
# nobs => number of observations in the eve
# area => eve areas (Km**2 == number of gridpoints if gridcell dimensions are 1x1 Km)
# volume => eve areas * eve rain field [mm * Km**2]
# mean.idi.x => average of the IDI field
# mean.idi.y => average of IDI at station locations
# mean.idiv.y => average of CV IDI at station locations
# mean.rain => average of rain field
# max.rain.x => maximum rain intensity in the rain field
# max.rain.yo => maximum rain intensity in the observations
# max.rain.ya => maximum rain intensity in the analysis 
# max.rain.yav => maximum rain intensity in the CV analysis 
# x.loc => x coordinate of the ellipsoid hull 
# “ellipsoid hull” or “spanning ellipsoid”, i.e. the ellipsoid of minimal volume
#  (‘area’ in 2D) such that all given points lie just inside or on the boundary
#  of the ellipsoid) center
# y.loc => y coordinate of the ellipsoid hull center
# s.maj.ax => semi-major axis of the ellipsoid hull
# s.min.ax => semi-minor axis of the ellipsoid hull
# dir.s.maj.ax => semi-major axis direction
#   dir=0 N-S orientation; dir=45 NE-SW; dir=90 E-W; dir=135 NW-SE; dir=180 N-S
# cv.rel.all  => mean(Cross-Validated relative error) set=all observation in that OAP
# cv.bias.all => mean(CV bias) set=all observation in that OAP
# cv.rmse.all => mean(CV root mean square error) set=all observation in that OAP
# cv.made.all => 1.4826 * Median Absolut Deviation set=all observation in that OAP
# cv.rel.q50  => mean(Cross-Validated relative error) set=observation in OAP >=q50.daily
# cv.bias.q50 => mean(CV bias) set=observation in OAP >=q50.daily
# cv.rmse.q50 => mean(CV root mean square error) set=observation in OAP >=q50.daily
# mean.idiv.y.q50 => mean(CV IDI for observation in OAP >=q50.daily)
# n.q50       => number of observation in OAP >=q50.daily
# cv.made.q50 => 1.4826 * Median Absolut Deviation set=observation in OAP >=q50.daily
# cv.rel.q75  => mean(Cross-Validated relative error) set=observation in OAP >=q75.daily
# cv.bias.q75 => mean(CV bias) set=observation in OAP >=q75.daily
# cv.rmse.q75 => mean(CV root mean square error) set=observation in OAP >=q75.daily
# cv.made.q75 => 1.4826 * Median Absolut Deviation set=observation in OAP >=q75.daily
# mean.idiv.y.q50 => mean(CV IDI for observation in OAP >=q75.daily)
# n.q50       => number of observation in OAP >=q75.daily
# idi.norm.fac => IDI normalization factor for that event 
# 
# History:
# 26.02.2015 - Cristian Lussana. original code from Bspat_TAMRR_v1_0.R
# -----------------------------------------------------------------------------
rm(list=ls())
# Libraries
library(igraph)
library(raster)
#rasterOptions(chunksize=2e+07,maxmemory=3e+08)
library(rgdal)
library(ncdf)
require(tripack)
require(cluster)
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#..............................................................................
# Functions
# + st.log:
st.log<-function(x,lc=NULL) {
  if (is.null(lc)) lc<-1.5
  if (is.vector(x)) {
    res<-vector(mode="numeric",length=length(x))
    aux<-which(x>lc)
    if (length(aux)>0) {
      res[aux]<-log10(x[aux])
      res[-aux]<-log10(lc)+(x[-aux]-lc)/(lc*log(10))
    } else {
      res<-log10(lc)+(x-lc)/(lc*log(10))
    }
  } else {
    if (x>lc) {
      res<-log10(x)
    } else {
      res<-log10(lc)+(x-lc)/(lc*log(10))
    }
  }
  return(res)
}
# + manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}
#+
read.cg2fg<-function(fname) {
  to.read<-file(fname, "rb")
  n<-readBin(to.read,integer(),size=4,n=1)
  print("n")
  print(n)
#  fglab<-vector(mode="integer",length=n)
#  nc<-vector(mode="integer",length=n)
  c1<-vector(mode="integer",length=n)
  c2<-vector(mode="integer",length=n)
  c3<-vector(mode="integer",length=n)
  c4<-vector(mode="integer",length=n)
  w1<-vector(mode="numeric",length=n)
  w2<-vector(mode="numeric",length=n)
  w3<-vector(mode="numeric",length=n)
  w4<-vector(mode="numeric",length=n)
#  d1<-vector(mode="numeric",length=n)
#  d2<-vector(mode="numeric",length=n)
#  d3<-vector(mode="numeric",length=n)
#  d4<-vector(mode="numeric",length=n)
  aux<-matrix(readBin(to.read,numeric(),size=4,n*8),ncol=8,nrow=n,byrow=T)
  print(aux[1:5,])
  print(ncol(aux))
#  fglab[]<-aux[,1]
#  nc[]<-aux[,2]
  c1[]<-aux[,1]
  c2[]<-aux[,2]
  c3[]<-aux[,3]
  c4[]<-aux[,4]
  w1[]<-aux[,5]
  w2[]<-aux[,6]
  w3[]<-aux[,7]
  w4[]<-aux[,8]
#  d1[]<-aux[,11]
#  d2[]<-aux[,12]
#  d3[]<-aux[,13]
#  d4[]<-aux[,14]
  close(to.read)
  return(list(c1=c1,c2=c2,c3=c3,c4=c4,
              w1=w1,w2=w2,w3=w3,w4=w4))
}

CG2FG_bilinear<-function(cg2fg,x.FG,x.CG,na.value) {
  n<-length(x.FG)
  out<-.C("cg2fg_bilinear",nfg=as.integer(n),
                           fg_mask=as.double(x.FG),
                           fg_field=as.double(x.FG),
                           cg_field=as.double(x.CG),
                           cg_c1=as.integer(cg2fg$c1),
                           cg_c2=as.integer(cg2fg$c2),
                           cg_c3=as.integer(cg2fg$c3),
                           cg_c4=as.integer(cg2fg$c4),
                           cg_w1=as.double(cg2fg$w1),
                           cg_w2=as.double(cg2fg$w2),
                           cg_w3=as.double(cg2fg$w3),
                           cg_w4=as.double(cg2fg$w4),
                           na=as.double(na.value) )
  return(out$fg_field)
}

CG2FG_ngb<-function(cg2fg,x.FG,x.CG,na.value,na.rm) {
  n<-length(x.FG)
  out<-.C("cg2fg_ngb",nfg=as.integer(n),
                      fg_mask=as.double(x.FG),
                      fg_field=as.double(x.FG),
                      cg_field=as.double(x.CG),
                      cg_c1=as.integer(cg2fg$c1),
                      cg_c2=as.integer(cg2fg$c2),
                      cg_c3=as.integer(cg2fg$c3),
                      cg_c4=as.integer(cg2fg$c4),
                      na=as.double(na.value),
                      na_rm=as.integer(na.rm) )
  return(out$fg_field)
}

CG2FG_cgweights<-function(cg2fg,x.FG,x.CG=x.eve.CG.1,na.value=na1) {
  nfg<-length(x.FG)
  ncg<-length(x.CG)
  out<-.C("cg2fg_cgweights",nfg=as.integer(nfg),
                            ncg=as.integer(ncg),
                            fg_weight=as.double(x.FG),
                            cg_weight=as.double(x.CG),
                            cg_mask=as.double(x.CG),
                            cg_c1=as.integer(cg2fg$c1),
                            cg_c2=as.integer(cg2fg$c2),
                            cg_c3=as.integer(cg2fg$c3),
                            cg_c4=as.integer(cg2fg$c4),
                            na=as.double(na.value) )
  return(out$cg_weight)
}

#G_fast<-function(xg,yg,zg,xo,yo,zo,Dh,Dz,g) {
#  ng<-length(xg)
#  no<-length(xo)
#  out<-.C("g_fast",ng=as.integer(ng),no=as.integer(no),
#                   xg=as.double(xg),yg=as.double(yg),zg=as.double(zg),
#                   xo=as.double(xo),yo=as.double(yo),zo=as.double(zo),
#                   Dh=as.double(Dh),Dz=as.double(Dz),g=as.double(g) )
#  return(out$g)
#}

OI_fast<-function(yo.sel,yb.sel,xb.sel,
                  xg,yg,zg,
                  xo,yo,zo,
                  xsmooth,xinf,
                  Dh,Dz) {
  no<-length(yo.sel)
  ng<-length(xb.sel)
  xa.sel<-vector(mode="numeric",length=ng)
  xidi.sel<-vector(mode="numeric",length=ng)
  vec<-vector(mode="numeric",length=no)
  vec1<-vector(mode="numeric",length=no)
  xa.sel[]<-0
  xidi.sel[]<-0
  d<-yo.sel-yb.sel
#  print(InvD[1,1:10])
#  print(d[1:10])
  out<-.C("oi_first",no=as.integer(no), 
                     innov=as.double(d),
                     SRinv=as.numeric(InvD),
                     vec=as.double(vec), vec1=as.double(vec1) )
  vec[]<-as.numeric(out$vec[])
  vec1[]<-as.numeric(out$vec1[])
  rm(out)
#  print(any(is.na(xg)))
#  print(any(is.na(yg)))
#  print(any(is.na(zg)))
#  print(any(is.na(xo)))
#  print(any(is.na(yo.sel)
#  print(any(is.na(zo)))
#  print(any(is.na(xb.sel)))
#  print(any(is.na(xsmooth)))
#  print(any(is.na(yo.sel)
#  print(any(is.na(yb.sel)))
#  print(any(is.na(d)))
#  print(any(is.na(vec)))
#  print(any(is.na(vec1)))
#  print(vec)
#  print(vec1)
  out<-.C("oi_fast",ng=as.integer(ng),no=as.integer(no),
                    xg=as.double(xg),yg=as.double(yg),zg=as.double(zg),
                    xo=as.double(xo),yo=as.double(yo),zo=as.double(zo),
                    Dh=as.double(Dh),Dz=as.double(Dz),
                    xb=as.double(xb.sel),
                    vec=as.double(vec),vec1=as.double(vec1),
                    smooth=as.double(xsmooth),xinf=as.double(xinf),
                    xa=as.double(xa.sel),xidi=as.double(xidi.sel) )
#  print(out$xa)
  return(list(xa=as.numeric(out$xa[]),xidi=as.numeric(out$xidi[])))
}
#-------------------------------------------------------------------
# [] Setup parameters
# eps2 meaning: believe in observations
eps2<-0.1
# rain/norain threshold
rr.inf<-0.1
# Contiguous NO-Rain and Rain Area (first step) 
Lsubsample.max<-20
Lsubsample.DHmax<-150
# identify adjacent nodes - Establish (local-) triangulation (Delauney)
n.sector<-16
sector.angle<-360/n.sector
# DQC thresholds 
yo.dqc.plausible.min<-0 # mm/h
yo.dqc.plausible.max<-500 # mm/h
# triangulation based chacks (flags=200 or 300)
DQC.min.dist.allowed<-20 #Km
# evaluation thresholds
q50.hourly<-0.4 #mm/h (weak/normal)
q75.hourly<-1.0 #mm/h (normal/moderate-strong)
q50.daily<-3 #mm/daily (weak/normal)
q75.daily<-8 #mm/daily (normal/moderate-strong)
# definition of Dh sequence - horizontal de-correlation length scale
min.Dh.seq.allowed<-10 # Km
#Dh.seq.ref<-100 # Km
Dh.seq.ref<-50 # Km
#Dh.seq.reference<-c(5000,4000,3000,2500,2000,1750,1500,1400,1300,1200,
#                    1100,1000, 900, 800, 700, 650, 600, 550, 500, 450,
#                     400, 350, 300, 260, 230, 200, 180, 160, 140, 120,
#                     100,  90,  80,  70,  60,  50,  40,  30,  20,  10)
Dh.seq.reference<-c(5000,4000,3000,2000,1000, 900, 800, 700, 600, 500,
                     400, 300, 200, 150, 100,  80,  60,  40,  20,  10)
#Dh.seq.reference<-c(500, 100,  40,  20,  10)
Dz.seq.gt.ref<-c(10000,1000) # m
Dz.seq.le.ref<-c(1000,500,250) # m
# analysis cycle on FG
ndim.FG.iteration<-10000
#ndim.FG.iteration<-5000
# stations outside Norway
max.Km.stnINdomain<-300
# MAIN ========================================================================
# [] Read command line arguments
arguments <- commandArgs()
print("Arguments")
print(arguments)
date.b.string<-arguments[3]
date.e.string<-arguments[4]
file_blacklist<-arguments[5]
file_errobs<-arguments[6]
config_file<-arguments[7]
config_par<-arguments[8]
if (length(arguments)!=8) 
  ext<-error_exit(paste("Error in command line arguments: \n",
  " R --vanilla yyyy.mm.dd yyyy.mm.dd blacklist errobs cfgFILE cfgPAR \n",
  sep=""))
# [] define/check paths
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
filenamedem<-paste(main.path.geoinfo,"/fennodem1000_ETRS89_mask.nc",sep="")
filenamedem.CG<-paste(main.path.geoinfo,"/fennodem10000_ETRS89_mask.nc",sep="")
if (!file.exists(filenamedem)) 
  ext<-error_exit(paste("File not found:",filenamedem))
if (!file.exists(filenamedem.CG)) 
  ext<-error_exit(paste("File not found:",filenamedem.CG))
# common libs and etcetera
path2lib.com<-paste(main.path,"/lib",sep="")
path2etc.com<-paste(main.path,"/etc",sep="")
if (!file.exists(paste(path2lib.com,"/nogrid.ncout.R",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/nogrid.ncout.R"))
if (!file.exists(paste(path2lib.com,"/ncout.spec.list.r",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/ncout.spec.list.r"))
if (!file.exists(paste(path2lib.com,"/getStationData.R",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/getStationData.R"))
source(paste(path2lib.com,"/nogrid.ncout.R",sep=""))
source(paste(path2lib.com,"/ncout.spec.list.r",sep=""))
source(paste(path2lib.com,"/getStationData.R",sep=""))
#
dyn.load(paste(main.path,"/Bspat_PREC1d/src/cg2fg_ngb.so",sep=""))
dyn.load(paste(main.path,"/Bspat_PREC1d/src/cg2fg_bilinear.so",sep=""))
dyn.load(paste(main.path,"/Bspat_PREC1d/src/cg2fg_cgweights.so",sep=""))
dyn.load(paste(main.path,"/Bspat_PREC1d/src/oi_first.so",sep=""))
dyn.load(paste(main.path,"/Bspat_PREC1d/src/oi_fast.so",sep=""))
# CG to FG, faster! 
file.cg2fg<-paste(main.path.geoinfo,"/fennodem10000to1000_easy.dat",sep="")
cg2fg<-read.cg2fg(fname=file.cg2fg)
# test mode
print(testmode)
if (testmode) {
  print("TESTMODE TESTMODE TESTMODE")
  if (file.exists(paste(main.path,"/Bspat_PREC1d/testbed",sep=""))) {
    testbed<-paste(main.path,"/Bspat_PREC1d/testbed",sep="")
    station.info<-paste(testbed,"/station_data.csv",sep="")
    observed.data<-paste(testbed,"/observed_data.csv",sep="")
  } else {
    ext<-error_exit(paste("testbed not found"))
  }
}
# netcdf fixed parameters
grid.type <- "laea"
prod.date <- substr(Sys.time(),1,10)
xa.source.nc<-"daily precipitation from station data"
xa.var.version <- "1.0"
xa.pname<-"PREC1d"
for (p in 1:length(ncout.spec.list)) {
  if (ncout.spec.list[[p]]$pname==xa.pname) break
}
if (p==length(ncout.spec.list) & (ncout.spec.list[[p]]$pname!=xa.pname) ) {
  xa.flag.write<-F
} else {
  xa.flag.write<-T
  xa.var.name<-ncout.spec.list[[p]]$opts$var.name
  xa.var.longname<-ncout.spec.list[[p]]$opts$var.longname
  xa.var.unit<-ncout.spec.list[[p]]$opts$var.unit
  xa.var.mv<-as.numeric(ncout.spec.list[[p]]$opts$var.mv)
  xa.times.unit <-ncout.spec.list[[p]]$opts$t.unit
  xa.times.ref <-ncout.spec.list[[p]]$opts$t.ref
  xa.reference <- ncout.spec.list[[p]]$opts$reference
}
xidi.source.nc<-"IDIms for daily precipitation from station data"
xidi.var.version <- "1.0"
xidi.pname<-"IDIms"
for (p in 1:length(ncout.spec.list)) {
  if (ncout.spec.list[[p]]$pname==xidi.pname) break
}
if (p==length(ncout.spec.list) & (ncout.spec.list[[p]]$pname!=xidi.pname) ) {
  xidi.flag.write<-F
} else {
  xidi.flag.write<-T
  xidi.var.name<-ncout.spec.list[[p]]$opts$var.name
  xidi.var.longname<-ncout.spec.list[[p]]$opts$var.longname
  xidi.var.unit<-ncout.spec.list[[p]]$opts$var.unit
  xidi.var.mv<-as.numeric(ncout.spec.list[[p]]$opts$var.mv)
  xidi.times.unit <-ncout.spec.list[[p]]$opts$t.unit
  xidi.times.ref <-ncout.spec.list[[p]]$opts$t.ref
  xidi.reference <- ncout.spec.list[[p]]$opts$reference
}
#
if (xa.pname=="PREC1d") {
  q50<-q50.daily
  q75<-q75.daily
} else {
  q50<-q50.hourly
  q75<-q75.hourly
}
# set Time-related variables
start.string<-paste(date.b.string,sep="")
end.string<-paste(date.e.string,sep="")
start.string.day<-paste(substr(start.string,1,10),".01",sep="")
end.string.day  <-paste(substr(end.string,1,10),".24",sep="")
start <- strptime(start.string,"%Y.%m.%d","UTC")
end <- strptime(end.string,"%Y.%m.%d","UTC")
timeseq<-as.POSIXlt(seq(as.POSIXlt(start),as.POSIXlt(end),by="1 hour"),"UTC")
nhour<-length(timeseq)
start.day <- strptime(start.string.day,"%Y.%m.%d.%H","UTC")
end.day <- strptime(end.string.day,"%Y.%m.%d.%H","UTC")
dayseq<-as.POSIXlt(seq(as.POSIXlt(start.day),as.POSIXlt(end.day),by="1 day"),"UTC")
nday<-length(dayseq)
yyyy.b<-timeseq$year[1]+1900
mm.b<-timeseq$mon[1]+1
dd.b<-timeseq$mday[1]
yyyy.e<-timeseq$year[nhour]+1900
mm.e<-timeseq$mon[nhour]+1
dd.e<-timeseq$mday[nhour]
date.b<-paste(dd.b,".",mm.b,".",yyyy.b,sep="")
date.e<-paste(dd.e,".",mm.e,".",yyyy.e,sep="")
yyyymmdd.b<-paste(yyyy.b,formatC(mm.b,width=2,flag="0"),formatC(dd.b,width=2,flag="0"),sep="")
yyyymm.b<-paste(yyyy.b,formatC(mm.b,width=2,flag="0"),sep="")
datestring.b<-paste(yyyy.b,".",mm.b,".",dd.b," ",sep="")
yyyymmdd.e<-paste(yyyy.e,formatC(mm.e,width=2,flag="0"),formatC(dd.e,width=2,flag="0"),sep="")
yyyymm.e<-paste(yyyy.e,formatC(mm.e,width=2,flag="0"),sep="")
datestring.e<-paste(yyyy.e,".",mm.e,".",dd.e," ",sep="")
print(timeseq)
print(nhour)
print(dayseq)
print(nday)
# output directories
dir.create(file.path(main.path.output,"NGCD"), showWarnings = FALSE)
dir.create(file.path(main.path.output,"NGCD_addInfo"), showWarnings = FALSE)
path2output.main<-paste(main.path.output,"/NGCD/PREC1d",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/NGCD_addInfo/PREC1d",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
path2output.add.eve<-paste(path2output.add,"/event_dataset",sep="")
if (!(file.exists(path2output.main)))     dir.create(path2output.main,showWarnings=F) 
if (!(file.exists(path2output.main.stn))) dir.create(path2output.main.stn,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) dir.create(path2output.main.grd,showWarnings=F) 
if (!(file.exists(path2output.add)))      dir.create(path2output.add,showWarnings=F) 
if (!(file.exists(path2output.add.grd)))  dir.create(path2output.add.grd,showWarnings=F) 
if (!(file.exists(path2output.add.eve)))  dir.create(path2output.add.eve,showWarnings=F) 
# Setup output files 
dir.create(paste(path2output.main.stn,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.main.grd,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.add.eve,"/",yyyymm.b,sep=""),showWarnings=F)
out.file.stn<- paste(path2output.main.stn,"/",yyyymm.b,
                     "/NGCD_PREC1d_station_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.eve<- paste(path2output.add.eve,"/",yyyymm.b,
                     "/NGCD_PREC1d_event_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm.b,
                         "/NGCD_PREC1d_grid_",
                         yyyymmdd.b,"_",yyyymmdd.e,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm.b,
                         "/NGCD_PREC1d_grid_normIDI_",
                         yyyymmdd.b,"_",yyyymmdd.e,".nc",sep="")
#
print("Output files:")
print("analysis on the grid (netcdf)")
print(out.file.grd.ana)
print("event-normalized idi on the grid (netcdf)")
print(out.file.grd.idi)
print("station outputs (text)")
print(out.file.stn)
print("event outputs (text)")
print(out.file.eve)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# HEADER HEADER HEADER HEADER HEADER HEADER HEADER HEADER
# define header for the station data output file
cat(paste("year","month","day","nday","stid",
          "x","y","z","eve.lab","yo",
          "yb","ya","yav","yidi","yidiv","dqcflag","\n",sep=";"),
          file=out.file.stn,append=F)
cat(paste("year","month","day","nday","eve.lab","nobs",
          "area","volume",
          "mean.idi.x","mean.idi.y","mean.idiv.y",
          "mean.rain","max.rain.x",
          "max.rain.yo","max.rain.ya","max.rain.yav",
          "x.loc","y.loc","s.maj.ax","s.min.ax","dir.s.maj.ax",
          "cv.rel.all","cv.bias.all","cv.rmse.all","cv.made.all",
          "cv.rel.q50","cv.bias.q50","cv.rmse.q50","cv.made.q50",
          "mean.idiv.y.q50","n.q50",
          "cv.rel.q75","cv.bias.q75","cv.rmse.q75","cv.made.q75",
          "mean.idiv.y.q75","n.q75",
          "idi.norm.fac",
          "\n",sep=";"),
          file=out.file.eve,append=F)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [] Grid
# CRS Coordinate Reference System
r.orog.FG<-raster(filenamedem)
nx.FG<-ncol(r.orog.FG)
ny.FG<-nrow(r.orog.FG)
dx.FG<-xres(r.orog.FG)
dy.FG<-yres(r.orog.FG)
area.1cell.FG<-dx.FG*dy.FG/(1000*1000) # Area Km**2
# 4 borders point (SW corner (xmn,ymn); NE corner (xmx,ymx))
xmn.FG<-xmin(r.orog.FG)
xmx.FG<-xmax(r.orog.FG)
ymn.FG<-ymin(r.orog.FG)
ymx.FG<-ymax(r.orog.FG)
#
r.orog.CG<-raster(filenamedem.CG)
nx.CG<-ncol(r.orog.CG)
ny.CG<-nrow(r.orog.CG)
dx.CG<-xres(r.orog.CG)
dy.CG<-yres(r.orog.CG)
xmn.CG<-xmin(r.orog.CG)
xmx.CG<-xmax(r.orog.CG)
ymn.CG<-ymin(r.orog.CG)
ymx.CG<-ymax(r.orog.CG)
# extract all the cell values: zvalues[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues.FG<-getValues(r.orog.FG)
storage.mode(zvalues.FG)<-"numeric"
xy.FG<-xyFromCell(r.orog.FG,1:ncell(r.orog.FG))
x.FG<-sort(unique(xy.FG[,1]))
y.FG<-sort(unique(xy.FG[,2]),decreasing=T)
mask.FG<-which(!is.na(zvalues.FG))
zgrid<-zvalues.FG[mask.FG]
xgrid<-xy.FG[mask.FG,1]
ygrid<-xy.FG[mask.FG,2]
Lgrid.FG<-length(mask.FG)
# CG
zvalues.CG<-getValues(r.orog.CG)
storage.mode(zvalues.CG)<-"numeric"
xy.CG<-xyFromCell(r.orog.CG,1:ncell(r.orog.CG))
mask.CG<-which(!is.na(zvalues.CG))
zgrid.CG<-zvalues.CG[mask.CG]
xgrid.CG<-xy.CG[mask.CG,1]
ygrid.CG<-xy.CG[mask.CG,2]
Lgrid.CG<-length(mask.CG)
# clean memory
rm(zvalues.CG,zvalues.FG)
# debug info
print("grid parameters")
print(paste("nx ny dx dy",as.integer(nx.FG),as.integer(ny.FG),round(dx.FG,2),round(dy.FG,2)))
print(paste("xmn xmx ymn ymx",round(xmn.FG,2),round(xmx.FG,2),round(ymn.FG,2),round(ymx.FG,2)))
print(paste("Lgrid.FG",as.integer(Lgrid.FG)))
print("grid.CG parameters")
print(paste("nx.CG ny.CG dx.CG dy.CG",as.integer(nx.CG),as.integer(ny.CG),round(dx.CG,2),round(dy.CG,2)))
print(paste("xmn.CG xmx.CG ymn.CG ymx.CG",round(xmn.CG,2),round(xmx.CG,2),round(ymn.CG,2),round(ymx.CG,2)))
print(paste("Lgrid.CG",as.integer(Lgrid.CG)))
#
#
na1<-(-999.)
#x.mask.FG<-vector(mode="numeric",length=Lgrid.FG)
#x.ngb.FG<-vector(mode="numeric",length=Lgrid.FG)
#x.bil.FG<-vector(mode="numeric",length=Lgrid.FG)
#r.ngb1.FG<-r.orog.FG
#r.bil1.FG<-r.orog.FG
#x.mask.FG[]<-1
#zgrid.CG[is.na(zgrid.CG)]<-na1
#print("ngb")
#x.ngb.FG<-CG2FG_ngb(cg2fg,x.FG=x.mask.FG,x.CG=zgrid.CG,na.value=na1,na.rm=1)
#print("ngb ok")
#r.ngb1.FG[mask.FG]<-x.ngb.FG
#rnc<-writeRaster(r.ngb1.FG,filename=paste("rn.nc",sep=""),format="CDF",overwrite=TRUE)
#print("bilinear")
#x.bil.FG<-CG2FG_bilinear(cg2fg,x.FG=x.mask.FG,x.CG=zgrid.CG,na.value=na1)
#print("bilinear ok")
#r.bil1.FG[mask.FG]<-x.bil.FG
#rnc<-writeRaster(r.bil1.FG,filename=paste("rb.nc",sep=""),format="CDF",overwrite=TRUE)
#print("Tests")
#r.ngb0.FG<-resample(r.orog.CG,r.orog.FG,method="ngb")
#r.bil0.FG<-resample(r.orog.CG,r.orog.FG,method="bilinear")
#r.ngbdiff<-r.ngb1.FG-r.ngb0.FG
#r.bildiff<-r.bil1.FG-r.bil0.FG
#rnc<-writeRaster(r.ngbdiff,filename=paste("rndiff.nc",sep=""),format="CDF",overwrite=TRUE)
#rnc<-writeRaster(r.bildiff,filename=paste("rbdiff.nc",sep=""),format="CDF",overwrite=TRUE)
#print("end")
#q()
#------------------------------------------------------------------------------
# [] Read Station Information 
# conditions:
# 1. stations in KDVH having: lat, lon and elevation. Note that UTM33 is 
#    obtained from lat,lon. Furthermore, the location must be in Norway or on
#    the border (less than max.Km.stnINdomain)
# 2. stations in CG
stations.tmp<-getStationMetadata_NGCD(from.year=yyyy.b,to.year=yyyy.b,path=main.path.output,
                                      max.Km=max.Km.stnINdomain,var="PREC1d")
n.tmp<-length(stations.tmp$stnr)
stn.out.CG<-vector(length=n.tmp)
stn.out.CG[1:n.tmp]<-F
# avoid duplicates
dup.indx<-which(duplicated(cbind(stations.tmp$x,stations.tmp$y)))
if (length(dup.indx)>0) stn.out.CG[dup.indx]<-T
# check if in CG
aux<-extract(r.orog.CG,cbind(stations.tmp$x,stations.tmp$y),na.rm=T)
stn.out.CG[which(is.na(aux))]<-T
# definitive station list
L.y.tot<-length(which(!stn.out.CG))
stations<-data.frame(matrix(nrow=L.y.tot,ncol=5))
names(stations)<-c("stnr","z","x","y","NO")
stations$stnr<-stations.tmp$stnr[which(!stn.out.CG)]
stations$cn<-stations.tmp$cn[which(!stn.out.CG)]
stations$z<-stations.tmp$z[which(!stn.out.CG)]
stations$y<-stations.tmp$y[which(!stn.out.CG)]
stations$x<-stations.tmp$x[which(!stn.out.CG)]
# [] define Vectors and Matrices
VecX<-vector(mode="numeric",length=L.y.tot)
VecY<-vector(mode="numeric",length=L.y.tot)
VecZ<-vector(mode="numeric",length=L.y.tot)
VecS<-vector(mode="numeric",length=L.y.tot)
Disth<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
Distz<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
n.yo<-vector(mode="numeric",length=L.y.tot)
yo<-vector(mode="numeric",length=L.y.tot)
y.eve<-vector(mode="numeric",length=L.y.tot)
#yo.n<-vector(mode="numeric",length=L.y.tot)
xb<-vector(mode="numeric",length=Lgrid.FG)
ya<-vector(mode="numeric",length=L.y.tot)
yav<-vector(mode="numeric",length=L.y.tot)
ydqc<-vector(mode="numeric",length=L.y.tot)
ydqc.flag<-vector(mode="numeric",length=L.y.tot)
y.mask<-vector(mode="numeric",length=L.y.tot)
Lsubsample.vec<-vector(mode="numeric",length=L.y.tot)
D<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
S<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
# fill VecX/Y/Z/S
VecX<-as.numeric(as.vector(stations$x))
VecY<-as.numeric(as.vector(stations$y))
VecZ<-as.numeric(as.vector(stations$z))
VecS<-as.numeric(as.vector(stations$stnr))
# identify stations outside FG
z.FG<-extract(r.orog.FG,cbind(VecX,VecY),na.rm=T)
stn.out.FG<-vector(length=L.y.tot)
stn.out.FG[1:L.y.tot]<-F
stn.out.FG[which(is.na(z.FG))]<-T


#Test
#i<-0
#while ((i*ndim.FG.iteration)<=Lgrid.FG) {
#  start<-min( i*ndim.FG.iteration+1, Lgrid.FG)
#  end<-min( (i+1)*ndim.FG.iteration, Lgrid.FG)
#  ndimaux.FG<-end-start+1
#  print(paste("start end ndimaux.FG",start,end,ndimaux.FG,Lgrid.FG))
#  pos<-start:end
#  start.time <- Sys.time()
#  print("G elaboration")
#  aux.FG<-matrix(ncol=length(VecX),nrow=ndimaux.FG,data=0.)
#  auxz.FG<-matrix(ncol=length(VecX),nrow=ndimaux.FG,data=0.)
#  print("outer h")
#  aux.FG<-(outer(ygrid[pos],VecY,FUN="-")**2. +
#           outer(xgrid[pos],VecX,FUN="-")**2.)**0.5/1000.
#  print("outer z")
#  auxz.FG<-abs(outer(zgrid[pos],VecZ,FUN="-"))
#  print("realG ")
#  G.FG<-matrix(ncol=length(VecX),nrow=ndimaux.FG,data=0.)
#  G<-matrix(ncol=length(VecX),nrow=ndimaux.FG,data=0.)
#  G.FG<-exp(-0.5*(aux.FG/20000)**2.-0.5*(auxz.FG/500)**2.)
#  rm(aux.FG,auxz.FG)
#  print("G elaboration end")
#  end.time <- Sys.time()
#  time.taken <- end.time - start.time
#  print(time.taken)
#  print("G fast")
#  start.time <- Sys.time()
#  a<-length(VecX)
#  out<-.C("g_fast",ng=as.integer(ndimaux.FG),no=as.integer(length(VecX)),
#                   xg=as.double(xgrid[pos]),yg=as.double(ygrid[pos]),zg=as.double(zgrid[pos]),
#                   xo=as.double(VecX),yo=as.double(VecY),zo=as.double(VecZ),
#                   Dh=as.double(20000.),Dz=as.double(500.),g=as.double(G) )
##  G<-as.matrix(out$g,nrow=ndimaux.FG,ncol=a,byrow=T)
#  G[]<-out$g
#  end.time <- Sys.time()
#  time.taken <- end.time - start.time
#  print(time.taken)
#  print("G fast end")
#  print(G[1,1:10])
#  print(G.FG[1,1:10])
#  A<-G-G.FG
#  print(max(A))
#  i<-i+1
#}
#q()








#------------------------------------------------------------------------------
# [] compute Disth and Distz (symmetric) matrices: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
#  Distz(i,j)=elevation difference between i-th station and j-th station [m]
Disth<-(outer(VecY,VecY,FUN="-")**2.+outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
Distz<-abs(outer(VecZ,VecZ,FUN="-"))
print("list of station ids")
print(cbind(1:L.y.tot,VecS))
#------------------------------------------------------------------------------
# [] Get observations from KDVH
n.yo[]<-NA
yo[]<-NA
# set DQC flag to NA for every station/observation
ydqc.flag[]<-NA
# cycle over days to get data from KDVH
#if (!testmode) {
#  if (nday==1) {
#    data<-getStationData(var="RR", from.yyyy=yyyy.b, from.mm=mm.b, from.dd=dd.b,
#                         to.yyyy=yyyy.e, to.mm=mm.e, to.dd=dd.e,
#                         qa=NULL, statlist=stations, outside.Norway=T,
#                         err.file=file_errobs, blist=file_blacklist,
#                         val.min.allowed=yo.dqc.plausible.min,
#                         val.max.allowed=yo.dqc.plausible.max,verbose=T)
#  } else {
#    data<-getStationData(var="RR", from.yyyy=yyyy.b, from.mm=mm.b, from.dd=dd.b,
#                         to.yyyy=yyyy.e, to.mm=mm.e, to.dd=dd.e,
#                         qa=NULL, statlist=stations, outside.Norway=T,
#                         err.file=file_errobs, blist=file_blacklist, 
#                         val.min.allowed=yo.dqc.plausible.min, 
#                         val.max.allowed=yo.dqc.plausible.max,fun="sum",verbose=T)
#  }
#} else {
#  data<-read.csv(file=observed.data)
#}
data<-getStationData_NGCD(var="PREC1d", yyyy=yyyy.b, mm=formatC(mm.b,width=2,flag="0"), dd=formatC(dd.b,width=2,flag="0"), statlist=stations,path=main.path.output)
# debug: start
print("Observed data (alligned with station info):")
print(data)
# debug: end
yo<-data$value
y.notNA<-which(!is.na(yo))
L.y.notNA<-length(y.notNA)
for (i in 1:L.y.tot) {
  if (is.na(data$value[i])) next
  ydqc.flag[i]<-0
#  ydqc.flag[i]<--1
#  if (!is.na(data$KDVHflag[i])) if (data$KDVHflag[i]>2) ydqc.flag[i]<-100
#  if (!is.na(data$KDVHflag[i])) if (data$KDVHflag[i]==0 | data$KDVHflag[i]==1 | data$KDVHflag[i]==2) ydqc.flag[i]<-0
#  if (data$ntime[i]!=data$nvalue[i] |
#      !data$plausible[i] |
#      data$err.ext[i] |
#      data$blist[i]) ydqc.flag[i]<-100
}
# Stations on the output file 
stn.output<-which((!stn.out.FG) | (!is.na(yo)) )
n.stn.output<-length(stn.output)
#data.names<-c("stnr","year","month","day","hour","ntime",
# "value","nvalue","DQC","KDVHflag","plausible","err.ext","blist.perm","blist.curr")
print(paste("number of station having at least one available observation =",L.y.notNA,"(tot stations in KDVH=",L.y.tot,")"))
print(paste("number of station having good observation           (so far)=",length(which(ydqc.flag<=0 & !is.na(yo)))))
#print(paste("  # station having at least one erroneous observation (plausibility check) =",length(which(!data$plausible & !is.na(yo)))))
#print(paste("  # station having at least one erroneous observation           (KDVH DQC) =",length(which(data$KDVHflag>2 & !is.na(yo)))))
#print(paste("  # station having at least one erroneous observation       (external DQC) =",length(which(data$err.ext & !is.na(yo)))))
#print(paste("  # station blacklisted for the time period considered =",length(which(data$blist & !is.na(yo)))))
print(paste("  # station in output file (on Norwegian mainland or within CG and not NA) =",length(stn.output)))
#------------------------------------------------------------------------------
# Elaborations
# loop for DQC
yo.ok.pos<-which(ydqc.flag<=0 & !is.na(yo))
L.yo.ok<-length(yo.ok.pos)
while (L.yo.ok>0) {
# define Vectors and Matrices
  if (exists("xidi.CG.wet")) rm(xidi.CG.wet)
  if (exists("xidi.CG.dry")) rm(xidi.CG.dry)
  xidi.CG.wet<-vector(mode="numeric",length=Lgrid.CG)
  xidi.CG.dry<-vector(mode="numeric",length=Lgrid.CG)
# vector with the positions (pointers to VecS) of good observations 
  yo.ok.pos<-which(ydqc.flag<=0 & !is.na(yo))
  L.yo.ok<-length(yo.ok.pos)
  yo.ok.wet<-ydqc.flag<=0 & !is.na(yo) & yo>=rr.inf
  yo.ok.dry<-ydqc.flag<=0 & !is.na(yo) & yo<rr.inf
  yo.ok.pos.wet<-which(yo.ok.wet)
  yo.ok.pos.dry<-which(yo.ok.dry)
  L.yo.ok.wet<-length(yo.ok.pos.wet)
  L.yo.ok.dry<-length(yo.ok.pos.dry)
  print(paste("observation not NAs, presumably good =",
              L.yo.ok,"(wet=",L.yo.ok.wet,"dry=",L.yo.ok.dry,")"))
# NO-RAIN OVER THE WHOLE DOMAIN
  if (L.yo.ok.wet==0) {
    print("no rain over the whole domain")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    cat(paste(yyyy.b,mm.b,dd.b,nday,
              formatC(VecS[stn.output],format="f",digits=0),
              formatC(VecX[stn.output],format="f",digits=0),
              formatC(VecY[stn.output],format="f",digits=0),
              formatC(VecZ[stn.output],format="f",digits=0),
              rep(NA,n.stn.output), #eve.lab
              formatC(yo[stn.output],format="f",digits=1),
              rep(NA,n.stn.output), #yb
              rep(0,n.stn.output), #ya
              rep(NA,n.stn.output), #yav
              rep(NA,n.stn.output), #yidi.eve
              rep(NA,n.stn.output), #yidiv.eve
              formatC(ydqc.flag[stn.output],format="f",digits=0),
              "\n",sep=";"),file=out.file.stn,append=T)
    # Figures
    r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                  ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    xa.FG<-vector(mode="numeric",length=Lgrid.FG)
    xa.FG[]<-0
    r.aux.FG[mask.FG]<-round(xa.FG,1)
    nogrid.ncout(file.name=out.file.grd.ana,
                 grid=t(as.matrix(r.aux.FG)),
                 x=x.FG,y=y.FG,grid.type=grid.type,
                 times=c(paste(yyyymmdd.b,"0000",sep="")),
                 prod.date=prod.date,
                 proj4.string=proj4.utm33,
                 var.name=xa.var.name,
                 var.longname=xa.var.longname,
                 var.unit=xa.var.unit,
                 var.mv=xa.var.mv,
                 var.version=xa.var.version,
                 times.unit=xa.times.unit,
                 times.ref=xa.times.ref,
                 reference=xa.reference,
                 source.string=xa.source.nc)
    nogrid.ncout(file.name=out.file.grd.idi,
                 grid=t(as.matrix(r.aux.FG)),
                 x=x.FG,y=y.FG,grid.type=grid.type,
                 times=c(paste(yyyymmdd.b,"0000",sep="")),
                 prod.date=prod.date,
                 proj4.string=proj4.utm33,
                 var.name=xidi.var.name,
                 var.longname=xidi.var.longname,
                 var.unit=xidi.var.unit,
                 var.mv=xidi.var.mv,
                 var.version=xidi.var.version,
                 times.unit=xidi.times.unit,
                 times.ref=xidi.times.ref,
                 reference=xidi.reference,
                 source.string=xidi.source.nc)
    print("Success exit")
    quit(status=0)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  } # end case of no-rain over the whole domain
#------------------------------------------------------------------------------
# [] Contiguous NO-Rain (nor) Area identification
# First Guess
  Lsubsample.vec[]<-NA
  # vectors:
  # +nnor.1st.vec>number of no-rain areas
  # +Lnor.1st.vec[i]>(i=1,..,nnor.1st.vec) number of stns in i-th no-rain area 
  # +nor.1st.vec[i,n]-> n-th stn (i.e. pointers to VecS) in i-th no-rain area
  #  matrix(i,j) i=1,nnor.1st.vec j=1,Lnor.1st.vec[i]
  #  note: an no-rain areas could contain both wet and dry stations
  nor.1st.vec<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Lnor.1st.vec<-vector(mode="numeric",length=L.yo.ok)
  Lnor.1st.vec[]<-NA
  nnor.1st.vec<-0
  nor.1st.aux<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Lnor.1st.aux<-vector(mode="numeric",length=L.yo.ok)
  Lnor.1st.aux[]<-NA
  nnor.1st.aux<-0
  for (b in yo.ok.pos.dry) {  # START: Cycle over dry stations 
    # a. identify the closest stations to the b-th dry station
    #  outputs: -close2b-> pointer to VecS of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -Lsubsample.vec[b]-> actual number of stations used
    close2b.aux<-order(Disth[b,],decreasing=F)
    close2b.au1<-close2b.aux[which(close2b.aux%in%yo.ok.pos)][1:Lsubsample.max]
    Lsubsample.vec[b]<-Lsubsample.max
    if (Disth[b,close2b.au1[Lsubsample.max]]>Lsubsample.DHmax) {
      if (any(Disth[b,close2b.au1]<=Lsubsample.DHmax)) {
        Lsubsample.vec[b]<-max(which(Disth[b,close2b.au1]<=Lsubsample.DHmax))
      } else {
        Lsubsample.vec[b]<-0
      }
    }
    if (Lsubsample.vec[b]<3) Lsubsample.vec[b]<-Lsubsample.max
    close2b<-close2b.au1[1:Lsubsample.vec[b]]
    rm(close2b.aux,close2b.au1)
    # b. setup vectors 
##    VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
##    VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
##    yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
##    VecX.b<-VecX[close2b]
##    VecY.b<-VecY[close2b]
##    wet.b.aux<-close2b %in% yo.ok.pos.wet
##    dry.b.aux<-close2b %in% yo.ok.pos.dry
##    wet.b<-close2b[wet.b.aux]
##    dry.b<-close2b[dry.b.aux]
##    nwet<-length(wet.b)
##    ndry<-length(dry.b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
##    tri.rr<-tri.mesh(VecX.b,VecY.b)
    tri.rr<-tri.mesh(VecX[close2b],VecY[close2b])
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th dry station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th dry station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,VecX[b],VecY[b])
    # note: due to the fact that (VecX[b],VecY[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-VecX[b]+sin(sector.angle*(cont.clock-1)*pi/180.)*1000
      y.aux<-VecY[b]+cos(sector.angle*(cont.clock-1)*pi/180.)*1000
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,nodes.aux[which( (!(nodes.aux %in% nodes)) & (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify nor 
    #    merge in a (new) nor all the (old) temporary nor (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary nor merged
    aux<-vector()
    if (nnor.1st.aux>0) {
      for (nor.1st in 1:nnor.1st.aux) {
        if (Lnor.1st.aux[nor.1st]==0) next
        if ( any(nodes %in% nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]]) ) {
          aux<-c(aux[which(!(aux %in% nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]]))],nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]])
          nor.1st.aux[nor.1st,]<-NA
          Lnor.1st.aux[nor.1st]<-0
        }
      }
    }
    nnor.1st.aux<-nnor.1st.aux+1
    Lnor.1st.aux[nnor.1st.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    nor.1st.aux[nnor.1st.aux,1:Lnor.1st.aux[nnor.1st.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the dry stations
  rm(aux,close2b,tri.rr,nodes)
# Reorganise no-rain area labels
  y.nor.1st<-vector(length=L.y.tot)
  y.nor.1st[1:L.y.tot]<-NA
  nnor.1st.vec<-0
  for (nor.1st in 1:nnor.1st.aux) {
    if (Lnor.1st.aux[nor.1st]==0) next
    nnor.1st.vec<-nnor.1st.vec+1
    Lnor.1st.vec[nnor.1st.vec]<-Lnor.1st.aux[nor.1st]
    nor.1st.vec[nnor.1st.vec,1:Lnor.1st.vec[nnor.1st.vec]]<-nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]]
    y.nor.1st[nor.1st.vec[nnor.1st.vec,1:Lnor.1st.vec[nnor.1st.vec]]]<-nnor.1st.vec
  }
  rm(nor.1st.aux,Lnor.1st.aux,nnor.1st.aux)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# identify dry-stations surrounded only by wet-stations (close enough)
  flag.dry.dqc<-F
  for (j in 1:nnor.1st.vec) {
    nor.1st.j<-nor.1st.vec[j,1:Lnor.1st.vec[j]]
    indx.dry.stns<-which(nor.1st.j %in% yo.ok.pos.dry)
    n.dry<-length(indx.dry.stns)
    if (n.dry>1) next
    nor.1st.j.dry<-nor.1st.j[indx.dry.stns]
    aux.dist<-Disth[nor.1st.j.dry,nor.1st.j]
    min.dist.from.dry.stn<-min(aux.dist[aux.dist>0])
    if (min.dist.from.dry.stn<DQC.min.dist.allowed) {
      flag.dry.dqc<-T
      ydqc.flag[nor.1st.j.dry]<-200
    }
  }
  rm(nor.1st.j,indx.dry.stns)
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.dry.stn,nor.1st.j.dry)
  if (flag.dry.dqc) next
#------------------------------------------------------------------------------
# [] Precipitation events identification (eve) (contiguous rain areas)
  Lsubsample.vec[]<-NA
  # vectors:
  # +neve.1st.vec>number of events 
  # +Leve.1st.vec[i]>(i=1,..,neve.1st.vec) number of stns in i-th event
  # +eve.1st.vec[i,n]-> n-th stn (i.e. pointers to VecS) in i-th event
  #  matrix(i,j) i=1,neve.1st.vec j=1,Leve.1st.vec[i]
  #  note: an event could contain both wet and dry stations
  eve.1st.vec<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Leve.1st.vec<-vector(mode="numeric",length=L.yo.ok)
  Leve.1st.vec[]<-NA
  neve.1st.vec<-0
  eve.1st.aux<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Leve.1st.aux<-vector(mode="numeric",length=L.yo.ok)
  Leve.1st.aux[]<-NA
  neve.1st.aux<-0
  for (b in yo.ok.pos.wet) {  # START: Cycle over wet stations 
    # a. identify the closest stations to the b-th wet station
    #  outputs: -close2b-> position (i.e. pointer to VecS) of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -Lsubsample.vec[b]-> actual number of stations used
    close2b.aux<-order(Disth[b,],decreasing=F)
    close2b.au1<-close2b.aux[which(close2b.aux%in%yo.ok.pos)][1:Lsubsample.max]
    Lsubsample.vec[b]<-Lsubsample.max
    if (Disth[b,close2b.au1[Lsubsample.max]]>Lsubsample.DHmax) {
      if (any(Disth[b,close2b.au1]<=Lsubsample.DHmax)) {
        Lsubsample.vec[b]<-max(which(Disth[b,close2b.au1]<=Lsubsample.DHmax))
      } else {
        Lsubsample.vec[b]<-0
      }
    }
    if (Lsubsample.vec[b]<3) Lsubsample.vec[b]<-Lsubsample.max
    close2b<-close2b.au1[1:Lsubsample.vec[b]]
    rm(close2b.aux,close2b.au1)
    # b. setup vectors 
##    VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
##    VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
##    yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
##    VecX.b<-VecX[close2b]
##    VecY.b<-VecY[close2b]
##    wet.b.aux<-close2b %in% yo.ok.pos.wet
##    dry.b.aux<-close2b %in% yo.ok.pos.dry
##    wet.b<-close2b[wet.b.aux]
##    dry.b<-close2b[dry.b.aux]
##    nwet<-length(wet.b)
##    ndry<-length(dry.b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
##    tri.rr<-tri.mesh(VecX.b,VecY.b)
    tri.rr<-tri.mesh(VecX[close2b],VecY[close2b])
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th wet station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th wet station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,VecX[b],VecY[b])
    # note: due to the fact that (VecX[b],VecY[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-VecX[b]+sin(sector.angle*(cont.clock-1)*pi/180.)*1000
      y.aux<-VecY[b]+cos(sector.angle*(cont.clock-1)*pi/180.)*1000
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,nodes.aux[which( (!(nodes.aux %in% nodes)) & (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify eve 
    #    merge in a (new) eve all the (old) temporary eve (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary eve merged
    aux<-vector()
    if (neve.1st.aux>0) {
      for (eve.1st in 1:neve.1st.aux) {
        if (Leve.1st.aux[eve.1st]==0) next
        if ( any(nodes %in% eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]]) ) {
          aux<-c(aux[which(!(aux %in% eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]]))],eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]])
          eve.1st.aux[eve.1st,]<-NA
          Leve.1st.aux[eve.1st]<-0
        }
      }
    }
    neve.1st.aux<-neve.1st.aux+1
    Leve.1st.aux[neve.1st.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    eve.1st.aux[neve.1st.aux,1:Leve.1st.aux[neve.1st.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the wet stations
  rm(aux,close2b,tri.rr,nodes)
# reorganise event labels
  y.eve.1st<-vector(length=L.y.tot)
  y.eve.1st[1:L.y.tot]<-NA
  neve.1st.vec<-0
  for (eve.1st in 1:neve.1st.aux) {
    if (Leve.1st.aux[eve.1st]==0) next
    neve.1st.vec<-neve.1st.vec+1
    Leve.1st.vec[neve.1st.vec]<-Leve.1st.aux[eve.1st]
    eve.1st.vec[neve.1st.vec,1:Leve.1st.vec[neve.1st.vec]]<-eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]]
    y.eve.1st[eve.1st.vec[neve.1st.vec,1:Leve.1st.vec[neve.1st.vec]]]<-neve.1st.vec
  }
  rm(eve.1st.aux,Leve.1st.aux,neve.1st.aux)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# wet-stations surrounded only by dry-stations (close enough)
  flag.wet.dqc<-F
  for (j in 1:neve.1st.vec) {
    eve.1st.j<-eve.1st.vec[j,1:Leve.1st.vec[j]]
    indx.wet.stns<-which(eve.1st.j %in% yo.ok.pos.wet)
    n.wet<-length(indx.wet.stns)
    if (n.wet>1) next
    eve.1st.j.wet<-eve.1st.j[indx.wet.stns]
    aux.dist<-Disth[eve.1st.j.wet,eve.1st.j]
    min.dist.from.wet.stn<-min(aux.dist[aux.dist>0])
    if (min.dist.from.wet.stn<DQC.min.dist.allowed) {
      flag.wet.dqc<-T
      ydqc.flag[eve.1st.j.wet]<-300
    }
  }
  rm(eve.1st.j,indx.wet.stns)
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.wet.stn,eve.1st.j.wet)
  if (flag.wet.dqc) next
#--------------------------------------------------------------------------------
# GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID
# [] Contiguous NO-Rain Area on the grid.CG
# A gridpoint is assigned to the same area/event of the closest station
  print("++ define no-rain areas and precipitation events on the grid.CG: triangulation")
  lab.nor.CG.1st<-vector(mode="numeric",length=Lgrid.CG)
  lab.nor.CG.1st[]<-NA
  # do the same for Contiguous Rain Area [to be used later]
  x.eve.CG.1st<-vector(mode="numeric",length=Lgrid.CG)
  x.eve.CG.1st[]<-NA
  for (i in 1:Lgrid.CG) {
    Disth.i<-( (ygrid.CG[i]-VecY[yo.ok.pos])**2+(xgrid.CG[i]-VecX[yo.ok.pos])**2 )**0.5 / 1000.
    close2i<-yo.ok.pos[which.min(Disth.i)]
    if (!is.na(y.nor.1st[close2i])) {
      lab.nor.CG.1st[i]<-y.nor.1st[close2i]
    } else {
      lab.nor.CG.1st[i]<-0
    }
    if (!is.na(y.eve.1st[close2i])) {
      x.eve.CG.1st[i]<-y.eve.1st[close2i]
    } else {
      x.eve.CG.1st[i]<-0
    }
  }
  rm(Disth.i,close2i)
##  r.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
##                        ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
##  r.CG[]<-NA
##  r.CG[mask.CG]<-lab.nor.CG.1st
  print("++ define no-rain areas on the grid.CG: IDI")
  xa.Dh.nor.1st<-vector(mode="numeric",length=Lgrid.CG)
  xa.Dz.nor.1st<-vector(mode="numeric",length=Lgrid.CG)
  x.dry.CG<-vector(mode="numeric",length=Lgrid.CG)
  n.nor<-0
# x.dry.CG -> =1 gridpoint not classified as dry; =0 grid-point classified as wet
  x.dry.CG[]<-0
  for (nor.1st in 1:nnor.1st.vec) { # START: cycle over nor.1st
#    print(paste(nor.1st,"nor.1st +++++++++++++++++++++++++++++++++"))
    # define vectors / variables / indexes 
    Lnor.1st<-Lnor.1st.vec[nor.1st]
    VecX.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecY.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecZ.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecS.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    yo.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecX.nor.1st<-VecX[nor.1st.vec[nor.1st,1:Lnor.1st]]
    VecY.nor.1st<-VecY[nor.1st.vec[nor.1st,1:Lnor.1st]]
    VecZ.nor.1st<-VecZ[nor.1st.vec[nor.1st,1:Lnor.1st]]
    VecS.nor.1st<-VecS[nor.1st.vec[nor.1st,1:Lnor.1st]]
    yo.nor.1st<-yo[nor.1st.vec[nor.1st,1:Lnor.1st]]
    yoindx.nor.1st.wet<-which(yo.nor.1st>=rr.inf & !is.na(yo.nor.1st))
    yoindx.nor.1st.dry<-which(yo.nor.1st<rr.inf & !is.na(yo.nor.1st))
    yoindx.wet<-nor.1st.vec[nor.1st,yoindx.nor.1st.wet]
    yoindx.dry<-nor.1st.vec[nor.1st,yoindx.nor.1st.dry]
    n.nor.1st.wet<-length(yoindx.nor.1st.wet)
    n.nor.1st.dry<-length(yoindx.nor.1st.dry)
    # the next line should not happen...
    if (n.nor.1st.dry==0) next
    # case of an isolated dry station...
    if (n.nor.1st.wet==0) {
      x.dry.CG[which(lab.nor.CG.1st==nor.1st)]<-1
      next
    }
    xindx.nor.CG<-which(lab.nor.CG.1st==nor.1st)
    Lgrid.nor.1st<-length(xindx.nor.CG)
    # define Dz for the current nor (min allowed value is 500 m)
    xa.Dz.nor.1st[xindx.nor.CG]<-max(500,max(zgrid.CG[xindx.nor.CG]))
    for (i in xindx.nor.CG) {    #START: cycle over xindx.nor.CG
      Disth.nor.1st<-( (ygrid.CG[i]-VecY.nor.1st)**2+(xgrid.CG[i]-VecX.nor.1st)**2 )**0.5 / 1000.
      Distz.nor.1st<-abs(zgrid.CG[i]-VecZ.nor.1st)
      close2i.ord<-order(Disth.nor.1st,decreasing=F)
      xa.Dh.nor.1st[i] <- Disth.nor.1st[close2i.ord[1]]
      notsofar.wet.i<-which(Disth.nor.1st[yoindx.nor.1st.wet]<(5*xa.Dh.nor.1st[i]))
      notsofar.dry.i<-which(Disth.nor.1st[yoindx.nor.1st.dry]<(5*xa.Dh.nor.1st[i]))
      n.notsofar.wet.i<-length(notsofar.wet.i)
      n.notsofar.dry.i<-length(notsofar.dry.i)
      if (n.notsofar.wet.i==0) {
        x.dry.CG[i]<-1
      } else if (n.notsofar.dry.i>0) {
# note: if (length(notsofar.dry.i)==0) next is implicit
        # IDI wet stations
        G.wet<-exp(-0.5*(Disth.nor.1st[yoindx.nor.1st.wet[notsofar.wet.i]]/xa.Dh.nor.1st[i])**2.
                   -0.5*(Distz.nor.1st[yoindx.nor.1st.wet[notsofar.wet.i]]/xa.Dz.nor.1st[i])**2.)
        if (n.notsofar.wet.i==1) {
          K.wet<-G.wet/(1+eps2)
        } else {
          D.wet<-as.matrix(exp(-0.5*(Disth[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dh.nor.1st[i])**2.
                               -0.5*(Distz[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dz.nor.1st[i])**2. ))
          D.wet[row(D.wet)==col(D.wet)]<-D.wet[row(D.wet)==col(D.wet)]+eps2
          InvD.wet<-solve(D.wet)
          K.wet<-tcrossprod(G.wet,InvD.wet)
          rm(D.wet,InvD.wet)
        }
        rm(G.wet)
        xidi.CG.wet[i]<-sum(K.wet)
        rm(K.wet)
        # IDI dry stations
        G.dry<-exp(-0.5*(Disth.nor.1st[yoindx.nor.1st.dry[notsofar.dry.i]]/xa.Dh.nor.1st[i])**2.-
                    0.5*(Distz.nor.1st[yoindx.nor.1st.dry[notsofar.dry.i]]/xa.Dz.nor.1st[i])**2.)
        if (n.notsofar.dry.i==1) {
          K.dry<-G.dry/(1+eps2)
        } else {
          D.dry<-as.matrix(exp(-0.5*(Disth[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dh.nor.1st[i])**2.+
                               -0.5*(Distz[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dz.nor.1st[i])**2. ))
          D.dry[row(D.dry)==col(D.dry)]<-D.dry[row(D.dry)==col(D.dry)]+eps2
          InvD.dry<-solve(D.dry)
          K.dry<-tcrossprod(G.dry,InvD.dry)
          rm(D.dry,InvD.dry)
        }
        rm(G.dry)
        xidi.CG.dry[i]<-sum(K.dry)
        rm(K.dry)
        if (xidi.CG.dry[i]>xidi.CG.wet[i]) x.dry.CG[i]<-1
      } # END: compare IDIs
    } # END: cycle over xindx.nor.CG
  } # END: cycle over nor.1st
  rm(xa.Dh.nor.1st,xa.Dz.nor.1st)
  rm(Disth.nor.1st,Distz.nor.1st)
  rm(VecX.nor.1st,VecY.nor.1st,VecZ.nor.1st,VecS.nor.1st,yo.nor.1st)
  rm(yoindx.nor.1st.wet,yoindx.nor.1st.dry,yoindx.wet,yoindx.dry,n.nor.1st.wet,n.nor.1st.dry)
  rm(close2i.ord)
  rm(notsofar.wet.i,notsofar.dry.i,n.notsofar.wet.i,n.notsofar.dry.i)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# if a dry observation is (1) not included in a dry area and (2) is in Norway then
# flag it as presumably bad observation and repeat spatial interpolation
  r.CG.nor.dry <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                        ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
  r.CG.nor.dry[]<-NA
  r.CG.nor.dry[mask.CG]<-round(x.dry.CG,0)
  stn.nor.dry<-extract(r.CG.nor.dry,cbind(VecX,VecY),na.rm=T)
  yo.ok.dry.fail<-!(yo.ok.pos.dry %in% which(stn.nor.dry==1 | is.na(stn.nor.dry)))
  if (any(yo.ok.dry.fail)) {
    ydqc.flag[yo.ok.pos.dry[which(yo.ok.dry.fail)]]<-400
    next
  }
  rm(r.CG.nor.dry,stn.nor.dry,yo.ok.dry.fail,x.dry.CG)
#------------------------------------------------------------------------------
# [] Precipitation events on the grid.CG
  print("++ define precipitation events on the grid.CG: IDI")
  xa.Dh.eve.1st<-vector(mode="numeric",length=Lgrid.CG)
  xa.Dz.eve.1st<-vector(mode="numeric",length=Lgrid.CG)
  x.wet.CG<-vector(mode="numeric",length=Lgrid.CG)
  n.eve<-0
# x.wet.CG -> =0 gridpoint not classified as wet; =1 grid-point classified as wet
# x.eve.CG -> =0 gridpoint "dry"; =-1 gridpoint "wet", waiting to labelled;
#             n>0 gridpoint "wet" belonging to eve n (only for isolated wet station)
  x.wet.CG[]<-0
  for (eve.1st in 1:neve.1st.vec) { # START: cycle over eve.1st
    # define vectors / variables / indexes 
    Leve.1st<-Leve.1st.vec[eve.1st]
    VecX.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecY.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecZ.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecS.eve.1st<-vector(mode="numeric",length=Leve.1st)
    yo.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecX.eve.1st<-VecX[eve.1st.vec[eve.1st,1:Leve.1st]]
    VecY.eve.1st<-VecY[eve.1st.vec[eve.1st,1:Leve.1st]]
    VecZ.eve.1st<-VecZ[eve.1st.vec[eve.1st,1:Leve.1st]]
    VecS.eve.1st<-VecS[eve.1st.vec[eve.1st,1:Leve.1st]]
    yo.eve.1st<-yo[eve.1st.vec[eve.1st,1:Leve.1st]]
    yoindx.eve.1st.wet<-which(yo.eve.1st>=rr.inf & !is.na(yo.eve.1st))
    yoindx.eve.1st.dry<-which(yo.eve.1st<rr.inf & !is.na(yo.eve.1st))
    yoindx.wet<-eve.1st.vec[eve.1st,yoindx.eve.1st.wet]
    yoindx.dry<-eve.1st.vec[eve.1st,yoindx.eve.1st.dry]
    n.eve.1st.wet<-length(yoindx.eve.1st.wet)
    n.eve.1st.dry<-length(yoindx.eve.1st.dry)
    # the next line should not happen...
    if (n.eve.1st.wet==0) next
    # case of an isolated wet station...
    if (n.eve.1st.dry==0) {
      x.wet.CG[which(x.eve.CG.1st==eve.1st)]<-1
      next
    }
    xindx.eve.CG<-which(x.eve.CG.1st==eve.1st)
    Lgrid.eve.1st<-length(xindx.eve.CG)
    # more than one single wet-station
    # define Dz for the current eve (min allowed value is 500 m)
    xa.Dz.eve.1st[xindx.eve.CG]<-max(500,max(zgrid.CG[xindx.eve.CG]))
    for (i in xindx.eve.CG) {    #START: cycle over xindx.eve.CG
      Disth.eve.1st<-( (ygrid.CG[i]-VecY.eve.1st)**2+(xgrid.CG[i]-VecX.eve.1st)**2 )**0.5 / 1000.
      Distz.eve.1st<-abs(zgrid.CG[i]-VecZ.eve.1st)
      close2i.ord<-order(Disth.eve.1st,decreasing=F)
      xa.Dh.eve.1st[i] <- Disth.eve.1st[close2i.ord[1]]
      notsofar.wet.i<-which(Disth.eve.1st[yoindx.eve.1st.wet]<(5*xa.Dh.eve.1st[i]))
      notsofar.dry.i<-which(Disth.eve.1st[yoindx.eve.1st.dry]<(5*xa.Dh.eve.1st[i]))
      n.notsofar.wet.i<-length(notsofar.wet.i)
      n.notsofar.dry.i<-length(notsofar.dry.i)
      if (n.notsofar.dry.i==0) {
        x.wet.CG[i]<-1
      } else if (n.notsofar.wet.i>0) {
# note: if (length(notsofar.wet.i)==0) next is implicit
        # IDI wet stations
        G.wet<-exp(-0.5*(Disth.eve.1st[yoindx.eve.1st.wet[notsofar.wet.i]]/xa.Dh.eve.1st[i])**2.
                   -0.5*(Distz.eve.1st[yoindx.eve.1st.wet[notsofar.wet.i]]/xa.Dz.eve.1st[i])**2.)
        if (n.notsofar.wet.i==1) {
          K.wet<-G.wet/(1+eps2)
        } else {
          D.wet<-as.matrix(exp(-0.5*(Disth[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dh.eve.1st[i])**2.
                               -0.5*(Distz[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dz.eve.1st[i])**2. ))
          D.wet[row(D.wet)==col(D.wet)]<-D.wet[row(D.wet)==col(D.wet)]+eps2
          InvD.wet<-solve(D.wet)
          K.wet<-tcrossprod(G.wet,InvD.wet)
          rm(D.wet,InvD.wet)
        }
        rm(G.wet)
        xidi.CG.wet[i]<-sum(K.wet)
        rm(K.wet)
        # IDI dry stations
        G.dry<-exp(-0.5*(Disth.eve.1st[yoindx.eve.1st.dry[notsofar.dry.i]]/xa.Dh.eve.1st[i])**2.-
                    0.5*(Distz.eve.1st[yoindx.eve.1st.dry[notsofar.dry.i]]/xa.Dz.eve.1st[i])**2.)
        if (n.notsofar.dry.i==1) {
          K.dry<-G.dry/(1+eps2)
        } else {
          D.dry<-as.matrix(exp(-0.5*(Disth[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dh.eve.1st[i])**2.+
                               -0.5*(Distz[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dz.eve.1st[i])**2. ))
          D.dry[row(D.dry)==col(D.dry)]<-D.dry[row(D.dry)==col(D.dry)]+eps2
          InvD.dry<-solve(D.dry)
          K.dry<-tcrossprod(G.dry,InvD.dry)
          rm(D.dry,InvD.dry)
        }
        rm(G.dry)
        xidi.CG.dry[i]<-sum(K.dry)
        rm(K.dry)
        if (xidi.CG.wet[i]>=xidi.CG.dry[i]) x.wet.CG[i]<-1
      }
    } # END: cycle over xindx.eve.CG
  } # END: cycle over eve.1st
  rm(xa.Dh.eve.1st,xa.Dz.eve.1st)
  rm(xidi.CG.wet,xidi.CG.dry)
  rm(Disth.eve.1st,Distz.eve.1st)
  rm(VecX.eve.1st,VecY.eve.1st,VecZ.eve.1st,VecS.eve.1st,yo.eve.1st)
  rm(yoindx.eve.1st.wet,yoindx.wet,n.eve.1st.wet)
  rm(close2i.ord)
  rm(notsofar.wet.i,notsofar.dry.i,n.notsofar.wet.i,n.notsofar.dry.i)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# if a wet observation is (1) not included in an event and (2) is in Norway then
# flag it as presumably bad observation and repeat spatial interpolation
  r.CG.eve.wet <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                        ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
  r.CG.eve.wet[]<-NA
  r.CG.eve.wet[mask.CG]<-round(x.wet.CG,0)
  stn.eve.wet<-extract(r.CG.eve.wet,cbind(VecX,VecY),na.rm=T)
  yo.ok.wet.fail<-!(yo.ok.pos.wet %in% which(stn.eve.wet==1 | is.na(stn.eve.wet)))
  if (any(yo.ok.wet.fail)) {
    ydqc.flag[yo.ok.pos.wet[which(yo.ok.wet.fail)]]<-500
    next
  }
  rm(stn.eve.wet,yo.ok.wet.fail,x.wet.CG)
#------------------------------------------------------------------------------
# [] Precipitation events labeling
# y.eve[i] i-th station association to precip event
# eve.lables[i] (i=1,..,n.eve) i-th precip event label value
# x.eve.CG[i] i-th CG gridpoint association to precip events
# r.eve.CG[i] correspondent raster data structure
# y.bwgt.CG[i] border-weight (between 0-1) for smoothing along the border
# x.eve.FG[i] i-th FG gridpoint association to precip events
# r.eve.FG[i] correspondent raster data structure
# *.bwgt.* border-weight (between 0-1) for smoothing along the border
#    CG-> minimum dist from borders/(2*mean(c(dx.CG,dy.CG))) 
#    FG-> resample CG events with bilinear interpolation (events on CG=1, outside events=0)
  x.eve.CG<-vector(mode="numeric",length=Lgrid.CG)
  x.eve.CG[]<-NA
  print("++ Precipitation events labeling: CG")
  # Detect clumps (patches) of connected cells. Each clump gets a unique ID.
  r.eve.CG<-clump(r.CG.eve.wet,directions=8,gaps=F)
  f.lab<-freq(r.eve.CG)
  aux<-which(!is.na(f.lab[,1]))
  f.lab.val<-f.lab[aux,1]
  f.lab.n<-f.lab[aux,2]
  x.eve.CG<-getValues(r.eve.CG)[mask.CG]
  storage.mode(x.eve.CG)<-"integer"
  rm(f.lab,aux,r.CG.eve.wet)
  # Detect (outer) edges - used for smoothing along the event borders 
  r.edge.CG<-edge(r.eve.CG,type="outer")
  x.edge.CG<-getValues(r.edge.CG)[mask.CG]
  x.edge.CG[x.edge.CG!=1]<-NA
  r.edge.CG[]<-NA
  r.edge.CG[mask.CG]<-x.edge.CG
  edge.CG<-which(!is.na(getValues(r.edge.CG)))
  xedge.CG<-xy.CG[edge.CG,1]
  yedge.CG<-xy.CG[edge.CG,2]
  # compute border-weights
#  r.dist<-distanceFromPoints(r.eve.CG,cbind(xedge.CG,yedge.CG))
#  aux<-getValues(r.dist)[mask.CG]/(2*mean(c(dx.CG,dy.CG)))
#  aux[is.na(x.eve.CG)|x.eve.CG<1]<-NA
#  aux[aux>1]<-1
#  x.bwgt.CG<-exp(-0.5*((aux-1)/0.25)**2)
#  r.bwgt.CG<-r.orog.CG
#  r.bwgt.CG[]<-NA
#  r.bwgt.CG[mask.CG]<-x.bwgt.CG
#  y.bwgt.CG<-extract(r.bwgt.CG,cbind(VecX,VecY),method="bilinear")
#  storage.mode(y.bwgt.CG)<-"numeric"
#  y.bwgt.CG[is.na(y.bwgt.CG)]<-1
##  rm(r.edge.CG,x.edge.CG,edge.CG,xedge.CG,yedge.CG,r.dist,aux,r.bwgt.CG)
#  rm(r.edge.CG,x.edge.CG,edge.CG,xedge.CG,yedge.CG,r.dist,aux)
#  rnc<-writeRaster(r.bwgt.CG,filename=paste("rbwgtcg0.nc",sep=""),format="CDF",overwrite=TRUE)
# assign stations at precipitation events
  print("++ Precipitation events labeling: stations")
  y.eve<-extract(r.eve.CG,cbind(VecX,VecY),na.rm=T)
  eve.labels<-as.vector(na.omit(unique(y.eve[yo.ok.pos.wet])))
  n.eve<-length(eve.labels)
  eve.nostn.pos<-which(!(f.lab.val %in% eve.labels))
  # events including 0 stations are re-labeled
  if (length(eve.nostn.pos)>0) {
    r <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
               ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
    x.aux<-vector(mode="numeric",length=Lgrid.CG)
    for (i in eve.nostn.pos) {
      aux<-which(x.eve.CG==f.lab.val[i])
      x.aux[]<-NA
      x.aux[aux]<-x.eve.CG[aux]
      r[]<-NA
      r[mask.CG]<-round(x.aux,0)
      dist.min<-NA
      for (e in eve.labels) {
        dist<-cellStats(distanceFromPoints(r,cbind(VecX[which(y.eve==e)],VecY[which(y.eve==e)])),mean)
        if (is.na(dist.min)) {
          dist.min<-dist
          new.lab<-e
        } else {
          if (dist<dist.min) {
            dist.min<-dist
            new.lab<-e
          }
        }
      }
      new.lab.pos<-which(f.lab.val==new.lab)
      x.eve.CG[aux]<-new.lab
      f.lab.val[i]<-NA
      f.lab.n[new.lab.pos]<-f.lab.n[new.lab.pos]+f.lab.n[i]
      f.lab.n[i]<-NA
      r.eve.CG[mask.CG[aux]]<-x.eve.CG[aux]
    }
    y.eve<-extract(r.eve.CG,cbind(VecX,VecY),na.rm=T)
    eve.labels<-as.vector(na.omit(unique(y.eve[yo.ok.pos.wet])))
    n.eve<-length(eve.labels)
  }
  if (exists("r")) rm(r,x.aux,aux,dist)
  #----------------------------------------------------------------------------
  print("++ Precipitation events labeling: FG")
# x.pnop.CG <- prec/no-prec distinction over CG (1=prec;0=no)
# x.pnop.FG <- prec/no-prec distinction over FG (1=prec;0=no)
#              --> smooth transition 0-1 by using bilinear int. <--
# x.eve.CG <- event labels on CG (NA = no precip)
# x.eve.CG.1 <- event labels on CG (na1 = no precip)
# x.eve.FG <- event labels on FG (NA = no precip)
# x.eve.FG.1 <- event labels on FG (na1 = no precip)
# --> given an event, a FG point is always surrounded by at least one CG point of the same event
# x.bwgt.FG <- background weights for more realistic prec/no-prec transition (0-1)
# y.bwgt.FG <- background weights for more realistic prec/no-prec transition (0-1)
# x.lngb.FG <- labels from CG to FG (na1 = no precip)
#  x.lngb.FG<-vector(mode="numeric",length=Lgrid.FG)
#  out<-.C("cg2fg_ngb",nsel=as.integer(length(cg2fg$fglab)),nlab=as.integer(length(cg2fg$fglab)),
#                      fglab=as.integer(cg2fg$fglab),fglabsel=as.integer(cg2fg$fglab),
#                      fgf=as.double(x.lngb.FG),cgf=as.numeric(x.eve.CG.1),cgngb=as.integer(cg2fg$cngb))
#  x.lngb.FG<-out$fgf
#  rm(out)
  # CG to FG resample event area with bilinear
  x.eve.CG<-getValues(r.eve.CG)[mask.CG]
  x.eve.CG.1<-x.eve.CG
  x.eve.CG.1[is.na(x.eve.CG.1)]<-na1
  x.pnop.CG<-vector(mode="numeric",length=Lgrid.CG)
  x.pnop.FG<-vector(mode="numeric",length=Lgrid.FG)
  x.mask.FG<-vector(mode="numeric",length=Lgrid.FG)
  x.pnop.CG[x.eve.CG.1!=na1]<-1
  x.pnop.CG[x.eve.CG.1==na1]<-0
  x.mask.FG[]<-1
  x.pnop.FG<-CG2FG_bilinear(cg2fg,x.FG=x.mask.FG,x.CG=x.pnop.CG,na.value=na1)
  x.pnop.FG[which(x.pnop.FG<0.01)]<-rep(0,length(which(x.pnop.FG<0.01)))
  # FG clump 
  x.eve.FG<-vector(mode="numeric",length=Lgrid.FG)
  x.eve.FG.1<-vector(mode="numeric",length=Lgrid.FG)
  x.mask.FG[which(x.pnop.FG==0)]<-na1
  x.eve.FG.1<-CG2FG_ngb(cg2fg,x.FG=x.mask.FG,x.CG=x.eve.CG.1,na.value=na1,na.rm=1)
  x.eve.FG[]<-x.eve.FG.1
  x.eve.FG[which(x.eve.FG.1==na1)]<-NA
  # compute the weights to be used in the transition area between rain and no-rain 
  x.bwgt.FG<-vector(mode="numeric",length=Lgrid.FG)
  x.bwgt.FG[]<-exp(-0.5*((x.pnop.FG-1)/0.25)**2)
  r.bwgt.FG<-r.orog.FG
  r.bwgt.FG[]<-NA
  r.bwgt.FG[mask.FG]<-x.bwgt.FG
  y.bwgt.FG<-extract(r.bwgt.FG,cbind(VecX,VecY),method="bilinear")
  storage.mode(y.bwgt.FG)<-"numeric"
  y.bwgt.FG[is.na(y.bwgt.FG)]<-1
  # compute CG weights (given FG weights)
  x.pnop.areal.CG<-vector(mode="numeric",length=Lgrid.CG)
  x.pnop.areal.CG<-CG2FG_cgweights(cg2fg,x.FG=x.pnop.FG,x.CG=x.eve.CG.1,na.value=na1)
  x.bwgt.CG<-vector(mode="numeric",length=Lgrid.CG)
  x.bwgt.CG[]<-exp(-0.5*((x.pnop.areal.CG-1)/0.25)**2)
  r.bwgt.CG<-r.orog.CG
  r.bwgt.CG[]<-NA
  r.bwgt.CG[mask.CG]<-x.bwgt.CG
  y.bwgt.CG<-extract(r.bwgt.CG,cbind(VecX,VecY),method="bilinear")
  storage.mode(y.bwgt.CG)<-"numeric"
  y.bwgt.CG[is.na(y.bwgt.CG)]<-1
# clean
#  rm(r.lngb.FG,x.lngb.FG,r.aux.CG,r.CGtoFG,r.clu.FG)
#  rm(f.lab,f.lab.val,f.lab.n,x.CGtoFG)
# debug
#  r<-r.orog.FG
#  r[mask.FG]<-x.pnop.FG
#  rnc<-writeRaster(r,filename=paste("xpnopfg.nc",sep=""),format="CDF",overwrite=TRUE)
#  r[mask.FG]<-x.eve.FG
#  rnc<-writeRaster(r,filename=paste("xevefg.nc",sep=""),format="CDF",overwrite=TRUE)
#  rnc<-writeRaster(r.eve.CG,filename=paste("revecg.nc",sep=""),format="CDF",overwrite=TRUE)
#  rnc<-writeRaster(r.bwgt.FG,filename=paste("rbwgtfg.nc",sep=""),format="CDF",overwrite=TRUE)
#  rnc<-writeRaster(r.bwgt.CG,filename=paste("rbwgtcg.nc",sep=""),format="CDF",overwrite=TRUE)
#------------------------------------------------------------------------------
# finally, we're interested only in events occurring in Norway
  aux<-which(eve.labels %in% x.eve.FG)
  if (length(aux)<n.eve) {
    tmp<-eve.labels
    eve.labels<-NA
    eve.labels<-tmp[aux]
    n.eve<-length(eve.labels)
    rm(tmp)
  }
  rm(aux)
# NO-RAIN OVER THE WHOLE DOMAIN
  if (n.eve==0) {
    print("no rain over the whole domain")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    cat(paste(yyyy.b,mm.b,dd.b,nday,
              formatC(VecS[stn.output],format="f",digits=0),
              formatC(VecX[stn.output],format="f",digits=0),
              formatC(VecY[stn.output],format="f",digits=0),
              formatC(VecZ[stn.output],format="f",digits=0),
              rep(NA,n.stn.output), #eve.lab
              formatC(yo[stn.output],format="f",digits=1),
              rep(NA,n.stn.output), #yb
              rep(0,n.stn.output), #ya
              rep(NA,n.stn.output), #yav
              rep(NA,n.stn.output), #yidi.eve
              rep(NA,n.stn.output), #yidiv.eve
              formatC(ydqc.flag[stn.output],format="f",digits=0),
              "\n",sep=";"),file=out.file.stn,append=T)
    # Figures
    r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                  ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    xa.FG<-vector(mode="numeric",length=Lgrid.FG)
    xa.FG[]<-0
    r.aux.FG[mask.FG]<-round(xa.FG,1)
    nogrid.ncout(file.name=out.file.grd.ana,
             grid=t(as.matrix(r.aux.FG)),
             x=x.FG,y=y.FG,grid.type=grid.type,
             times=c(paste(yyyymmdd.b,"0000",sep="")),
             prod.date=prod.date,
             proj4.string=proj4.utm33,
             var.name=xa.var.name,
             var.longname=xa.var.longname,
             var.unit=xa.var.unit,
             var.mv=xa.var.mv,
             var.version=xa.var.version,
             times.unit=xa.times.unit,
             times.ref=xa.times.ref,
             reference=xa.reference,
             source.string=xa.source.nc)
    print("Success exit")
    quit(status=0)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  } # end case of no-rain over the whole domain
#------------------------------------------------------------------------------
# ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS
# [] Analysis
  print("++ Analysis on Events")
# + data structures definition
  xa.FG<-vector(mode="numeric",length=Lgrid.FG)
  xb.FG<-vector(mode="numeric",length=Lgrid.FG)
  xa.CG<-vector(mode="numeric",length=Lgrid.CG)
  xb.CG<-vector(mode="numeric",length=Lgrid.CG)
  xidi.CG<-vector(mode="numeric",length=Lgrid.CG)
  xidi.FG<-vector(mode="numeric",length=Lgrid.FG)
  xtmp.FG<-vector(mode="numeric",length=Lgrid.FG)
# Y vectors (station locations)
  ya<-vector(mode="numeric",length=L.y.tot)
  yb<-vector(mode="numeric",length=L.y.tot)
  yb1.tmp<-vector(mode="numeric",length=L.y.tot)
  yb2.tmp<-vector(mode="numeric",length=L.y.tot)
  ybv.mat<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
  yidi.eve<-vector(mode="numeric",length=L.y.tot)
  yidiv.eve<-vector(mode="numeric",length=L.y.tot)
  yo.n<-vector(mode="numeric",length=L.y.tot)
  yav<-vector(mode="numeric",length=L.y.tot)
  ya.tmp<-vector(mode="numeric",length=L.y.tot)
  yav.tmp<-vector(mode="numeric",length=L.y.tot)
# events descriptive vectors
  n.y.eve<-vector(mode="numeric",length=n.eve)
  n.ya.eve<-vector(mode="numeric",length=n.eve)
  area.eve<-vector(mode="numeric",length=n.eve)
  volume.eve<-vector(mode="numeric",length=n.eve)
  meanrain.eve<-vector(mode="numeric",length=n.eve)
  maxrain.x.eve<-vector(mode="numeric",length=n.eve)
  maxrain.yo.eve<-vector(mode="numeric",length=n.eve)
  maxrain.ya.eve<-vector(mode="numeric",length=n.eve)
  maxrain.yav.eve<-vector(mode="numeric",length=n.eve)
  meanidi.x.eve<-vector(mode="numeric",length=n.eve)
  meanidi.y.eve<-vector(mode="numeric",length=n.eve)
  meanidiv.y.eve<-vector(mode="numeric",length=n.eve)
  meanidiv.y.eve.q50<-vector(mode="numeric",length=n.eve)
  meanidiv.y.eve.q75<-vector(mode="numeric",length=n.eve)
  # ell -> ellipse
  ell.locx.eve<-vector(mode="numeric",length=n.eve)
  ell.locy.eve<-vector(mode="numeric",length=n.eve)
  ell.smajor.eve<-vector(mode="numeric",length=n.eve)
  ell.sminor.eve<-vector(mode="numeric",length=n.eve)
  ell.smadir.eve<-vector(mode="numeric",length=n.eve)
  # cv -> cross-validated
  cv.rel.eve.all<-vector(mode="numeric",length=n.eve)
  cv.bias.eve.all<-vector(mode="numeric",length=n.eve)
  cv.rmse.eve.all<-vector(mode="numeric",length=n.eve)
  cv.made.eve.all<-vector(mode="numeric",length=n.eve)
  cv.rel.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.bias.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.rmse.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.made.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.rel.eve.q75<-vector(mode="numeric",length=n.eve)
  cv.bias.eve.q75<-vector(mode="numeric",length=n.eve)
  cv.rmse.eve.q75<-vector(mode="numeric",length=n.eve)
  cv.made.eve.q75<-vector(mode="numeric",length=n.eve)
  idi.norm.fac<-vector(mode="numeric",length=n.eve)
  n.q50<-vector(mode="numeric",length=n.eve)
  n.q75<-vector(mode="numeric",length=n.eve)
# + initialization
  xa.FG[]<-NA
  xb.FG[]<-NA
  xidi.FG[]<-0
  xidi.CG[]<-0
  xtmp.FG[]<-0
  xa.CG[]<-NA
  xb.CG[]<-NA
  yo.n[]<-NA
  yo.n[yo.ok.pos]<-0
  yb[]<-NA
  yb[yo.ok.pos]<-0
  yb1.tmp[]<-NA
  yb1.tmp[yo.ok.pos]<-0
  yb2.tmp[]<-NA
  yb2.tmp[yo.ok.pos]<-0
  yidi.eve[]<-0
  yidiv.eve[]<-0
  ya[]<-NA
  ya[yo.ok.pos]<-0
  yav[]<-NA
  yav[yo.ok.pos]<-0
  ell.locx.eve[]<-NA
  ell.locy.eve[]<-NA
  ell.smajor.eve[]<-NA
  ell.sminor.eve[]<-NA
  ell.smadir.eve[]<-NA
#------------------------------------------------------------------------------
# ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS 
# + Ellipsoid hulls
# ellipsoid hulls computation could fail (i.e. small events) and return NAs. 
  # debug: start
#  png(file="../../seNorge2_scratch/Bspat_PREC1d/elli.png",width=1200,height=1200)
#  mx<-as.integer(max(eve.labels,na.rm=T))
#  cols<-c("gray",rainbow((mx-1)))
#  plot(r.eve.CG,breaks=seq(0.5,(mx+0.5),by=1),col=cols)
  # debug: end
  for (n in 1:n.eve) {
#    print(paste("n eve.labels[n]",n,eve.labels[n]))
    xindx.eve.CG<-which(x.eve.CG==eve.labels[n])
    Lgrid.eve<-length(xindx.eve.CG)
#    print(Lgrid.eve)
    if (Lgrid.eve<=(min.Dh.seq.allowed)**2.) { 
      ell.smajor.eve[n] <- 2*min.Dh.seq.allowed
      ell.sminor.eve[n] <- 2*min.Dh.seq.allowed
      ell.smadir.eve[n] <- NA
      next
    }
    xy<-cbind(xgrid.CG[xindx.eve.CG],ygrid.CG[xindx.eve.CG])
    ell<-ellipsoidhull(xy)
    ell.locx.eve[n]<-ell$loc[1]
    ell.locy.eve[n]<-ell$loc[2]
    if ( is.na(ell$cov[1,1]) | is.na(ell$cov[2,2]) | is.na(ell$cov[1,2]) | 
         ell$cov[1,1]==0 | ell$cov[2,2]==0 | ell$cov[1,2]==0) {
      ell.smajor.eve[n] <- 2*min.Dh.seq.allowed
      ell.sminor.eve[n] <- 2*min.Dh.seq.allowed
      ell.smadir.eve[n] <- NA
      next
    }
#   debug plot: lines(predict(ell), col="blue")
#    lines(predict(ell), col="blue")
#    points(ell.locx.eve[n],ell.locy.eve[n],pch=19,col="black",cex=1.5)
#    eigen -> values sorted in decreasing order; 
#             vectors matrix whose columns contain the eigenvectors
    eigenval<-eigen(ell$cov)$values
    eigenvec<-eigen(ell$cov)$vectors
    e <- sqrt(eigenval)
    ell.smajor.eve[n] <- sqrt(ell$d2) * e[1] /1000  # semi-major axis [Km]
    ell.sminor.eve[n] <- sqrt(ell$d2) * e[2] /1000  # semi-minor axis [Km]
    # orientation of the ellipse: angle between the major axis and the y-axis
    # dir=0 N-S orientation; dir=45 NE-SW; dir=90 E-W; dir=135 NW-SE; dir=180 N-S
    ell.smadir.eve[n]<-atan2(eigenvec[1,1],eigenvec[2,1])
    if (ell.smadir.eve[n]<0) ell.smadir.eve[n]<-ell.smadir.eve[n]+pi
    ell.smadir.eve[n]<-ell.smadir.eve[n]/pi*180.
#    print(paste(eve.labels[n]))
#    print(paste(ell.locx.eve[n],ell.locy.eve[n]))
#    print(paste(ell.smajor.eve[n],ell.sminor.eve[n]))
#    print(paste(ell.smadir.eve[n]))
  } # end cycle: compute ellipsoid hulls
# debug: close plot session
#  dev.off()
  if (exists("xindx.eve.CG")) rm(xindx.eve.CG)
  if (exists("eigenval")) rm(xy,ell,eigenval,eigenvec,e)
#------------------------------------------------------------------------------
# ANALYSIS OVER EVENTS  ANALYSIS OVER EVENTS  ANALYSIS OVER EVENTS  ANALYSIS OVER EVENTS  
# + Analysis cycle over events 
  print("++ Analysis cycle over events")
  for (n in 1:n.eve) {
    # initialization
    r.xb.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                     ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
    r.xb.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                     ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    #
    r.xb.CG[]<-NA
    r.xb.FG[]<-NA
    xb.CG[]<-NA
    xidi.CG[]<-0
    # identify event extension on the CG grid
    xindx.eve.CG<-which(x.eve.CG==eve.labels[n])
    xindx.eve.FG<-which(x.eve.FG==eve.labels[n])
    Lgrid.eve.CG<-length(xindx.eve.CG)
    Lgrid.eve.FG<-length(xindx.eve.FG)
    # stations located in the event area without valid observation
    ya.indx<-which(y.eve==eve.labels[n] & !(yo.ok.wet))
    # stations located in the event area having valid observation
    yindx<-which(y.eve==eve.labels[n] & yo.ok.wet)
    n.y.eve[n]<-length(yindx)
    n.ya.eve[n]<-length(ya.indx)
    # debug
    eve.rx<-range(xgrid.CG[xindx.eve.CG])
    eve.ry<-range(ygrid.CG[xindx.eve.CG])
    print(paste("+ eve lab >",eve.labels[n],
                "(",n,"/",n.eve,")",
                " #obs=",n.y.eve[n]," #stns=",(n.y.eve[n]+n.ya.eve[n]),
                sep=""))
    print(paste("event range xmin,xmax (diff[Km]) / ymin,ymax (diff[Km])= ",
                round(eve.rx[1],0)," ",round(eve.rx[2],0),
                " (",round((eve.rx[2]-eve.rx[1])/1000,0),")/ ",
                round(eve.ry[1],0)," ",round(eve.ry[2],0),
                " (",round((eve.ry[2]-eve.ry[1])/1000,0),")",sep=""))
    print(paste("area [Km**2]=",Lgrid.eve.FG))
    print(paste("smaj smin [Km] dir[degree]=",
                round(ell.smajor.eve[n],0),
                round(ell.sminor.eve[n],0),
                round(ell.smadir.eve[n],2)))
#   ---------------------------------------------------------------------------
#   + event observed by a single (isolated) station
    if (n.y.eve[n]==1) {
      # background=observation -> analysis field =observation smoothed at the borders
      # note: xidi=-1; yidi/v=-1; yb=NA; yav=0
      xa.FG[xindx.eve.FG]<-yo[yindx]*x.bwgt.FG[xindx.eve.FG]
      xidi.FG[xindx.eve.FG]<-(-1)
      ya[yindx]<-yo[yindx]*y.bwgt.FG[yindx]
      yb[yindx]<-NA
      yav[yindx]<-0
      yidi.eve[yindx]<-(-1)
      yidiv.eve[yindx]<-(-1)
      if (n.ya.eve[n]>0) {
        ya[ya.indx]<-yo[yindx]*y.bwgt.FG[ya.indx] 
        yb[ya.indx]<-NA
        yav[ya.indx]<-0
        yidi.eve[ya.indx]<-(-1)
        yidiv.eve[ya.indx]<-(-1)
      }
      # event properties 
      area.eve[n]<-Lgrid.eve.FG*area.1cell.FG # Area Km**2
      volume.eve[n]<-sum(xa.FG[xindx.eve.FG])
      meanidi.x.eve[n]<-NA
      meanidi.y.eve[n]<-NA
      meanidiv.y.eve[n]<-NA
      meanidiv.y.eve.q50[n]<-NA
      meanidiv.y.eve.q75[n]<-NA
      meanrain.eve[n]<-mean(xa.FG[xindx.eve.FG])
      maxrain.x.eve[n]<-max(xa.FG[xindx.eve.FG])
      maxrain.yo.eve[n]<-max(yo[yindx])
      maxrain.ya.eve[n]<-max(ya[yindx])
      maxrain.yav.eve[n]<-NA
      cv.rel.eve.all[n]<-NA
      cv.bias.eve.all[n]<-NA
      cv.rmse.eve.all[n]<-NA
      cv.made.eve.all[n]<-NA
      cv.rel.eve.q50[n]<-NA
      cv.bias.eve.q50[n]<-NA
      cv.rmse.eve.q50[n]<-NA
      cv.made.eve.q50[n]<-NA
      cv.rel.eve.q75[n]<-NA
      cv.bias.eve.q75[n]<-NA
      cv.rmse.eve.q75[n]<-NA
      cv.made.eve.q75[n]<-NA
      idi.norm.fac[n]<-NA
      n.q50[n]<-NA
      n.q75[n]<-NA
      next
    }
#   ---------------------------------------------------------------------------
#   + event including more than one observation
#   + Larger-scale Background = mode(distribution of observations)
    histog<-hist(yo[yindx],breaks=seq(0,round((max(yo[yindx])+1),0)),plot=F)
    mode<-which.max(histog$counts)-0.5
    if (mode<rr.inf) mode<-rr.inf
    xb.CG[]<-NA
    xb.CG[xindx.eve.CG]<-mode
    yb[yindx]<-mode
    ybv.mat[yindx,yindx]<-mode
    if (n.ya.eve[n]>0) yb[ya.indx]<-mode
#   + define Dh.seq: sequence of horizontal decorrelation length scales
#     general case: from ellipsoid hull major semiaxis (coarser scale) 
#                   to a minimum predefined value (finer scale)
    Disth.aux<-Disth[yindx,yindx]
    min.Dh.seq<-as.integer(min.Dh.seq.allowed)
    max.Dh.seq<-as.integer(max(7*min.Dh.seq.allowed,ceiling(ell.smajor.eve[n])))
    aux<-which(Dh.seq.reference>=min.Dh.seq & Dh.seq.reference<=max.Dh.seq)
    Dh.seq<-c(Dh.seq.reference[aux],min.Dh.seq.allowed)
    n.Dh.seq<-length(Dh.seq)
    # debug
    print("Dh.seq [Km]:")
    print(Dh.seq)
#   + pre-elaboration to obtain err-cov matrices on CG
#    G.CG<-matrix(ncol=n.y.eve[n],nrow=Lgrid.eve.CG,data=0.)
#    aux.CG<-(outer(ygrid.CG[xindx.eve.CG],VecY[yindx],FUN="-")**2. +
#             outer(xgrid.CG[xindx.eve.CG],VecX[yindx],FUN="-")**2.)**0.5/1000.
#    auxz.CG<-abs(outer(zgrid.CG[xindx.eve.CG],VecZ[yindx],FUN="-"))
    InvD<-matrix(data=0,ncol=n.y.eve[n],nrow=n.y.eve[n])
#   + cycle on horizontal decorrelation length scales
    flag.CGtoFG<-T
    n.iter<-0
    for (Dh.test in Dh.seq) {
      n.iter<-n.iter+1
      # background
      if (Dh.test>Dh.seq.ref) {
        y.bwgt<-y.bwgt.CG
      } else {
        y.bwgt<-y.bwgt.FG
        # from CG to FG
        if (flag.CGtoFG) {
          flag.CGtoFG<-F
          #
          r.xb.CG[mask.CG]<-xb.CG
          xb.CG[is.na(xb.CG)]<-0
          xidi.CG[is.na(xidi.CG)]<-0
          #
          xtmp.FG[]<-na1
          xtmp.FG[xindx.eve.FG]<-1
          xtmp.FG<-CG2FG_bilinear(cg2fg,x.FG=xtmp.FG,x.CG=xb.CG,na.value=na1)
          xb.FG[]<-NA
          xb.FG[xindx.eve.FG]<-xtmp.FG[xindx.eve.FG]
          r.xb.FG[]<-NA
          r.xb.FG[mask.FG]<-xb.FG
          #
          xtmp.FG[]<-na1
          xtmp.FG[xindx.eve.FG]<-1
          xtmp.FG<-CG2FG_bilinear(cg2fg,x.FG=xtmp.FG,x.CG=xidi.CG,na.value=na1)
          xidi.FG[xindx.eve.FG]<-xtmp.FG[xindx.eve.FG]
          #
          yb1.tmp[yindx]<-extract(r.xb.CG,cbind(VecX[yindx],VecY[yindx]))
          yb2.tmp[yindx]<-extract(r.xb.FG,cbind(VecX[yindx],VecY[yindx]))
          aux<-which(!is.na(yb2.tmp[yindx]) & !is.na(yb1.tmp[yindx]))
          if (length(aux)>0) yb[yindx][aux]<-yb[yindx][aux]+(yb2.tmp[yindx][aux]-yb1.tmp[yindx][aux])
          if (n.ya.eve[n]>0) {
            yb1.tmp[ya.indx]<-extract(r.xb.CG,cbind(VecX[ya.indx],VecY[ya.indx]))
            yb2.tmp[ya.indx]<-extract(r.xb.FG,cbind(VecX[ya.indx],VecY[ya.indx]))
            aux<-which(!is.na(yb2.tmp[ya.indx]) & !is.na(yb1.tmp[ya.indx]))
            if (length(aux)>0) yb[ya.indx][aux]<-yb[ya.indx][aux]+(yb2.tmp[ya.indx][aux]-yb1.tmp[ya.indx][aux])
          }
        }
      }
      #
      D.test.part<-exp(-0.5*(Disth[yindx,yindx]/Dh.test)**2.)
      J<-1000000
      if (Dh.test>=(5*Dh.seq.ref)) Dz.seq<-c(10000)
      if (Dh.test>Dh.seq.ref & Dh.test<(5*Dh.seq.ref)) Dz.seq<-Dz.seq.gt.ref
      if (Dh.test<=Dh.seq.ref) Dz.seq<-Dz.seq.le.ref
      # last iteration= no filter on observation
      if (n.iter==n.Dh.seq) {
          yo.n.tmp<-vector(mode="numeric",length=L.y.tot)
          yo.n.tmp[yindx]<-yo[yindx]
      }
#     + cycle on vertical decorrelation length scales
      print("cycle on vertical decorrelation length scales")
      for (Dz.test in Dz.seq) {
        # compute superobservations according to decorrelation lenght-scales
        if (n.iter!=n.Dh.seq) {
          yo.n.tmp<-vector(mode="numeric",length=L.y.tot)
          for (b in yindx) {
            aux.yo<-which(Disth[b,yindx]<=(Dh.test) & Distz[b,yindx]<=(Dz.test))
            yo.n.tmp[b]<-mean(yo[yindx[aux.yo]])
          }
        }
        # analysis at eve station locations
        D.test<-D.test.part*exp(-0.5*(Distz[yindx,yindx]/Dz.test)**2.)
        S.test<-D.test
        D.test[row(D.test)==col(D.test)]<-D.test[row(D.test)==col(D.test)]+eps2
        InvD.test<-solve(D.test)
        W.test<-tcrossprod(S.test,InvD.test)
        aux<-t(yo.n.tmp[yindx]-yb[yindx])
        ya.tmp[yindx]<-yb[yindx] +tcrossprod(W.test,aux)
        rm(aux)
        yav.tmp[yindx]<-yo.n.tmp[yindx] +
                           1./(1.-diag(W.test)) * (ya.tmp[yindx]-yo.n.tmp[yindx])
        # cost function
        # st.log to avoid small observations to influence too much the result 
        J.tmp<-(1/n.y.eve[n]*sum((st.log(yav.tmp[yindx])-st.log(yo.n.tmp[yindx]))**2))**0.5
        # J.tmp<-(1/n.y.eve[n]*sum((yav.tmp[yindx]-yo.n.tmp[yindx])**2))**0.5
        if (Dz.test==Dz.seq[1] | J.tmp<J) {
          # save the best values (minimum value of the cost function)
          J<-J.tmp
          InvD<-InvD.test
          D<-D.test
          S<-S.test
          W<-W.test
          Dz.choice<-Dz.test
# + ya! analysis at station locations
          ya[yindx]<-ya.tmp[yindx] * y.bwgt[yindx]
          aux<-which(ya[yindx]<rr.inf & y.bwgt[yindx]>=0.9)
          if (length(aux)>0) ya[yindx][aux]<-rr.inf
          aux<-which(ya[yindx]<0 | is.na(ya[yindx]))
          if (length(aux)>0) ya[yindx][aux]<-0
          yo.n[yindx]<-yo.n.tmp[yindx]
# + yidi and yidiv are accumulated
          yidi.eve[yindx]<-yidi.eve[yindx]+rowSums(W)
          yidiv.eve[yindx]<-yidiv.eve[yindx]+
                            rep(1,n.y.eve[n])+ 1./(1.-diag(W)) * (rowSums(W)-rep(1,n.y.eve[n]))
        }
      } # end cycle on vertical decorrelation lenght scales
      print(paste("Dh=",Dh.test,"Km ==> Dz=",Dz.choice,"m"))
      rm(D.test.part,D.test,S.test,InvD.test,W.test)
      t.d<-t(yo.n[yindx]-yb[yindx])
      uno<-t(rep(1,length(yo.n[yindx])))
# optimization: if innovation is very close to zero go to next iteration
#      if ( !(any(abs(t.d)>0.01)) & n.iter!=n.Dh.seq) next
#     + CrossValidated-analysis 
      print("CrossValidated-analysis...NO")
#      i<-0
#      for (b in yindx) {
#        i<-i+1
#        InvD.1<-InvD[-i,-i]-1/InvD[i,i]*outer(InvD[-i,i],InvD[i,-i])
#        W.1<-tcrossprod(S[,-i],t(InvD.1))
#        aux<-t(yo.n[yindx[-i]]-ybv.mat[yindx[-i],b])
#        ya.tmp[yindx]<-ybv.mat[yindx,b]+tcrossprod(W.1,aux)
#        ya.tmp[yindx]<-ya.tmp[yindx] * y.bwgt[yindx]
#        aux<-which(ya.tmp[yindx]<rr.inf & y.bwgt[yindx]>=0.9)
#        if (length(aux)>0) ya.tmp[yindx][aux]<-rr.inf
#        aux<-which(ya.tmp[yindx]<0 | is.na(ya.tmp[yindx]))
#        if (length(aux)>0) ya.tmp[yindx][aux]<-0
#        yav[yindx[i]]<-ya.tmp[yindx[i]]
#        ybv.mat[yindx,b]<-ya.tmp[yindx]
#      }
      yav[yindx]<-ya[yindx]
#      rm(W.1,InvD.1,D,W,S,aux)
#     + Analysis on other station points within the event areas
      print("Analysis on other station points within the event areas...NO")
#      if (n.ya.eve[n]>0) {
#        G.ya<-matrix(ncol=n.y.eve[n],nrow=n.ya.eve[n],data=0.)
#        G.ya<-exp(-0.5*(Disth[ya.indx,yindx]/Dh.test)**2.
#                  -0.5*(Distz[ya.indx,yindx]/Dz.choice)**2.)
#        K<-tcrossprod(G.ya,InvD)
#        yidi.eve[ya.indx]<-yidi.eve[ya.indx]+rowSums(K)
#        yidiv.eve[ya.indx]<-yidi.eve[ya.indx]
#        ya[ya.indx]<-tcrossprod(K,t.d)
#        ya[ya.indx]<-ya[ya.indx] + yb[ya.indx]
#        aux<-which(ya[ya.indx]<rr.inf & y.bwgt[ya.indx]>=0.9)
#        if (length(aux)>0) ya[ya.indx][aux]<-rr.inf
#        aux<-which(ya[ya.indx]<0 | is.na(ya[ya.indx]))
#        if (length(aux)>0) ya[ya.indx][aux]<-0
#        yav[ya.indx]<-ya[ya.indx]
#        rm(K,G.ya)
#      }
      yav[ya.indx]<-ya[ya.indx]
#     optimization: 
#     + Dh>ref.value ->analysis on CG; Dh<ref.value ->analysis on FG
#     + Analysis on CG
      if (Dh.test>Dh.seq.ref) {
        print("Analysis on CG")
        oi<-OI_fast(yo.sel=yo.n[yindx],yb.sel=yb[yindx],
                    xb.sel=xb.CG[xindx.eve.CG],
                    xg=xgrid[xindx.eve.CG],yg=ygrid[xindx.eve.CG],zg=zgrid[xindx.eve.CG],
                    xo=VecX[yindx],yo=VecY[yindx],zo=VecZ[yindx],
                    xsmooth=x.bwgt.CG[xindx.eve.CG],xinf=rr.inf,
                    Dh=Dh.test,Dz=Dz.choice) 
        xa.CG[xindx.eve.CG]<-oi$xa
        xidi.CG[xindx.eve.CG]<-xidi.CG[xindx.eve.CG]+oi$xidi
#        print("xa.CG[xindx.eve.CG]")
#        print(xa.CG[xindx.eve.CG])
        rm(oi)
#        q()
        # next background is the current analysis, only a little bit smoothed
        # smoothing is for realistic border effects
        xb.CG[]<-NA
        xb.CG[xindx.eve.CG]<-xa.CG[xindx.eve.CG]
        print("end Analysis on CG")
      } else { 
#     + Analysis on FG
        print("Analysis on FG")
        oi<-OI_fast(yo.sel=yo.n[yindx],yb.sel=yb[yindx],
                    xb.sel=xb.FG[xindx.eve.FG],
                    xg=xgrid[xindx.eve.FG],yg=ygrid[xindx.eve.FG],zg=zgrid[xindx.eve.FG],
                    xo=VecX[yindx],yo=VecY[yindx],zo=VecZ[yindx],
                    xsmooth=x.bwgt.FG[xindx.eve.FG],xinf=rr.inf,
                    Dh=Dh.test,Dz=Dz.choice) 
        xa.FG[xindx.eve.FG]<-oi$xa
        xidi.FG[xindx.eve.FG]<-xidi.FG[xindx.eve.FG]+oi$xidi
        rm(oi)
        print("end Analysis on FG")
        # current analysis is next iteration background
        xb.FG[]<-NA
        xb.FG[xindx.eve.FG]<-xa.FG[xindx.eve.FG]
      } # end: analysis on FG
      # current analysis is next iteration background
      yb[yindx]<-ya[yindx]
      yb[ya.indx]<-ya[ya.indx]
    } # end cycle on horizontal decorrelation lenght scales
    if (exists("K.CG")) rm(K.CG,G.CG)
    #
    idi.norm.fac[n]<-max(xidi.FG[xindx.eve.FG],yidi.eve[yindx],yidiv.eve[yindx])
    xidi.FG[xindx.eve.FG]<-xidi.FG[xindx.eve.FG]/idi.norm.fac[n]
    yidi.eve[yindx]<-yidi.eve[yindx]/idi.norm.fac[n]
    yidiv.eve[yindx]<-yidiv.eve[yindx]/idi.norm.fac[n]
    if (n.ya.eve[n]>0) {
      yidi.eve[ya.indx]<-yidi.eve[ya.indx]/idi.norm.fac[n]
      yidiv.eve[ya.indx]<-yidiv.eve[ya.indx]/idi.norm.fac[n]
    }
    area.eve[n]<-Lgrid.eve.FG*area.1cell.FG # Area Km**2
    volume.eve[n]<-sum(xa.FG[xindx.eve.FG])
    meanidi.x.eve[n]<-mean(xidi.FG[xindx.eve.FG])
    meanidi.y.eve[n]<-mean(yidi.eve[yindx])
    meanidiv.y.eve[n]<-mean(yidiv.eve[yindx])
    meanrain.eve[n]<-mean(xa.FG[xindx.eve.FG])
    maxrain.x.eve[n]<-max(xa.FG[xindx.eve.FG])
    yo.aux<-yo[yindx]
    ya.aux<-ya[yindx]
    yav.aux<-yav[yindx]
    yidiv.aux<-yidiv.eve[yindx]
    maxrain.yo.eve[n]<-max(yo.aux)
    maxrain.ya.eve[n]<-max(ya.aux)
    maxrain.yav.eve[n]<-max(yav.aux)
    cv.rel.eve.all[n]<-mean(yav.aux/yo.aux)
    aux<-yav.aux-yo.aux
    cv.bias.eve.all[n]<-mean(aux)
    cv.rmse.eve.all[n]<-sqrt(mean(aux**2.))
    cv.made.eve.all[n]<-1.4826*median(abs(aux))
    q50.pos<-which(yo.aux>=q50)
    q75.pos<-which(yo.aux>=q75)
    n.q50[n]<-length(q50.pos)
    n.q75[n]<-length(q75.pos)
    if (n.q50[n]>0) {
      cv.rel.eve.q50[n]<-mean(yav.aux[q50.pos]/yo.aux[q50.pos])
      cv.bias.eve.q50[n]<-mean(aux[q50.pos])
      cv.rmse.eve.q50[n]<-sqrt(mean(aux[q50.pos]**2.))
      cv.made.eve.q50[n]<-1.4826*median(abs(aux[q50.pos]))
      meanidiv.y.eve.q50[n]<-mean(yidiv.aux[q50.pos])
    } else {
      cv.rel.eve.q50[n]<-NA
      cv.bias.eve.q50[n]<-NA
      cv.rmse.eve.q50[n]<-NA
      cv.made.eve.q50[n]<-NA
      meanidiv.y.eve.q50[n]<-NA
    }
    if (n.q75[n]>0) {
      cv.rel.eve.q75[n]<-mean(yav.aux[q75.pos]/yo.aux[q75.pos])
      cv.bias.eve.q75[n]<-mean(aux[q75.pos])
      cv.rmse.eve.q75[n]<-sqrt(mean(aux[q75.pos]**2.))
      cv.made.eve.q75[n]<-1.4826*median(abs(aux[q75.pos]))
      meanidiv.y.eve.q75[n]<-mean(yidiv.aux[q75.pos])
    } else {
      cv.rel.eve.q75[n]<-NA
      cv.bias.eve.q75[n]<-NA
      cv.rmse.eve.q75[n]<-NA
      cv.made.eve.q75[n]<-NA
      meanidiv.y.eve.q75[n]<-NA
    }
  } # END CYCLE over events
  if (exists("r.xb.CG")) rm(r.xb.CG)
  if (exists("r.xb.FG")) rm(r.xb.FG)
  break
} # end of DQC loop
#
xa.FG[is.na(xa.FG)]<-0.
#
r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                  ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
r.aux.FG[]<-NA
r.aux.FG[mask.FG]<-xa.FG
#r.aux.FG<-trim(r.aux.FG)
ya.tmp<-extract(r.aux.FG,cbind(VecX,VecY))
y.eve[is.na(ya)]<-NA
ya[is.na(ya)]<-round(ya.tmp[is.na(ya)],1)
yb[is.na(y.eve)]<-NA
yav[is.na(y.eve)]<-NA
yidi.eve[is.na(y.eve)]<-NA
yidiv.eve[is.na(y.eve)]<-NA
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
print("++ Output")
# Station Points - Write output on file 
cat(paste(yyyy.b,mm.b,dd.b,nday,
          formatC(VecS[stn.output],format="f",digits=0),
          formatC(VecX[stn.output],format="f",digits=0),
          formatC(VecY[stn.output],format="f",digits=0),
          formatC(VecZ[stn.output],format="f",digits=0),
          formatC(y.eve[stn.output],format="f",digits=0),
          formatC(yo[stn.output],format="f",digits=1),
          formatC(yb[stn.output],format="f",digits=2),
          formatC(ya[stn.output],format="f",digits=2),
          formatC(yav[stn.output],format="f",digits=2),
          formatC(yidi.eve[stn.output],format="f",digits=4),
          formatC(yidiv.eve[stn.output],format="f",digits=4),
          formatC(ydqc.flag[stn.output],format="f",digits=0),
          "\n",sep=";"),file=out.file.stn,append=T)
cat(paste(yyyy.b,mm.b,dd.b,nday,
          eve.labels,
          n.y.eve,
          formatC(area.eve,format="f",digits=0),
          formatC(volume.eve,format="f",digits=2),
          formatC(meanidi.x.eve,format="f",digits=4),
          formatC(meanidi.y.eve,format="f",digits=4),
          formatC(meanidiv.y.eve,format="f",digits=4),
          formatC(meanrain.eve,format="f",digits=2),
          formatC(maxrain.x.eve,format="f",digits=2),
          formatC(maxrain.yo.eve,format="f",digits=2),
          formatC(maxrain.ya.eve,format="f",digits=2),
          formatC(maxrain.yav.eve,format="f",digits=2),
          formatC(ell.locx.eve,format="f",digits=1),
          formatC(ell.locy.eve,format="f",digits=1),
          formatC(ell.smajor.eve,format="f",digits=1),
          formatC(ell.sminor.eve,format="f",digits=1),
          formatC(ell.smadir.eve,format="f",digits=1),
          formatC(cv.rel.eve.all,format="f",digits=4),
          formatC(cv.bias.eve.all,format="f",digits=4),
          formatC(cv.rmse.eve.all,format="f",digits=4),
          formatC(cv.made.eve.all,format="f",digits=4),
          formatC(cv.rel.eve.q50,format="f",digits=4),
          formatC(cv.bias.eve.q50,format="f",digits=4),
          formatC(cv.rmse.eve.q50,format="f",digits=4),
          formatC(cv.made.eve.q50,format="f",digits=4),
          formatC(meanidiv.y.eve.q50,format="f",digits=4),
          formatC(n.q50,format="f",digits=0),
          formatC(cv.rel.eve.q75,format="f",digits=4),
          formatC(cv.bias.eve.q75,format="f",digits=4),
          formatC(cv.rmse.eve.q75,format="f",digits=4),
          formatC(cv.made.eve.q75,format="f",digits=4),
          formatC(meanidiv.y.eve.q75,format="f",digits=4),
          formatC(n.q75,format="f",digits=0),
          formatC(idi.norm.fac,format="f",digits=5),
          "\n",sep=";"),file=out.file.eve,append=T)
# Figure: idi
r.aux.FG[]<-NA
r.aux.FG[mask.FG]<-round(xidi.FG,1)
r.aux.FG[mask.FG[xidi.FG>0]]<-round(100*xidi.FG[xidi.FG>0],1)
if (xidi.flag.write) {
  nogrid.ncout(file.name=out.file.grd.idi,
               grid=t(as.matrix(r.aux.FG)),
               x=x.FG,y=y.FG,grid.type=grid.type,
               times=c(paste(yyyymmdd.b,"0000",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.utm33,
               var.name=xidi.var.name,
               var.longname=xidi.var.longname,
               var.unit=xidi.var.unit,
               var.mv=xidi.var.mv,
               var.version=xidi.var.version,
               times.unit=xidi.times.unit,
               times.ref=xidi.times.ref,
               reference=xidi.reference,
               source.string=xidi.source.nc)
}
# Figure: Analysis on high-resolution grid
r.aux.FG[mask.FG]<-round(xa.FG,1)
if (xa.flag.write) {
  nogrid.ncout(file.name=out.file.grd.ana,
               grid=t(as.matrix(r.aux.FG)),
               x=x.FG,y=y.FG,grid.type=grid.type,
               times=c(paste(yyyymmdd.b,"0000",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.utm33,
               var.name=xa.var.name,
               var.longname=xa.var.longname,
               var.unit=xa.var.unit,
               var.mv=xa.var.mv,
               var.version=xa.var.version,
               times.unit=xa.times.unit,
               times.ref=xa.times.ref,
               reference=xa.reference,
               source.string=xa.source.nc)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Success: exit!
print("Success exit")
quit(status=0)

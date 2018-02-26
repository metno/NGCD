#!/usr/bin/env Rscript
# --~- Bspat_T.R -~--
# Bayesian Spatial Interpolation of daily aggregated temperature
# See the software repository here: https://github.com/metno/NGCD
#..............................................................................
#Copyright and license
# Copyright (C) 2018 MET Norway. NGCD software is licensed under GPL version 3 
# or (at your option) any later version.
# https://www.gnu.org/licenses/gpl-3.0.en.html
# 
# History:
# 22.02.2018 - Cristian Lussana. Original code.
# -----------------------------------------------------------------------------
#
rm(list=ls())
#
#-------------------------------------------------------------------
# Libraries
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("raster"))
suppressPackageStartupMessages(library("rgdal"))
suppressPackageStartupMessages(library("ncdf4"))
options(warn = 2, scipen = 999)
# 
# -----------------------------------------------------------------------------
# Constants
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# netcdf fixed parameters
varunit<-"Celsius degree"
varversion<-"1.0"
reference<-"Lussana, C., Tveito, O. E. and Uboldi, F. (2018), Three-dimensional spatial interpolation of 2 m temperature over Norway. Q.J.R. Meteorol. Soc.. doi:10.1002/qj.3208"
diground<-1
sourcestring<-"MET Norway"
comment<-"Our open data are licensed under Norwegian Licence for Open Government Data (NLOD) or a Creative Commons Attribution 4.0 International License at your preference. Credit should be given to The Norwegian Meteorological institute, shortened “MET Norway”, as the source of data."
#
#-------------------------------------------------------------------
# FUNCTIONS 
# + manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}
#
#==============================================================================
# MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN -
#==============================================================================
t0<-Sys.time()
# [] Read command line arguments and/or set parameters to default
# create parser object
p <- arg_parser("ngcdt")
# specify our desired options 
# by default ArgumentParser will add an help option 
p <- add_argument(p, "date",
                  help="date (format YYYY.MM.DD)",
                  type="character")
p <- add_argument(p, "var",
                  help="variable (TG mean daily T, TX max, TN min)",
                  type="character")
# options 
p <- add_argument(p, "--verbose",help="debug mode",flag=T,short="-v")
#
p <- add_argument(p, "--blacklist.file",
                  help="blacklist file",
                  type="character",
                  default="none")
p <- add_argument(p, "--errobs.file",
                  help="file with erroneous observations",
                  type="character",
                  default="none")
# Spatial interpolation parameters
p <- add_argument(p, "--eps2",
                  help="ratio of observation to background error covariance",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--sig2o",
                  help="observation error variance",
                  type="numeric",
                  default=1.3)
p <- add_argument(p, "--Dh",
                  help="background error covariance matrix, horizontal decorrelation length scale (km)",
                  type="numeric",
                  default=60)
p <- add_argument(p, "--Dz",
                  help="background error covariance matrix, vertical decorrelation length scale (m)",
                  type="numeric",
                  default=600)
p <- add_argument(p, "--lafmin",
                  help="land area fraction influence in the OI correlation function",
                  type="numeric",
                  default=0.5)

p <- add_argument(p, "--T2",
                  help="SCT threshold",
                  type="numeric",
                  default=15)
# Background parameters
p <- add_argument(p, "--z.range.min",
                  help="(m)",
                  type="numeric",
                  default=50)
p <- add_argument(p, "--nstat.min.inv",
                  help="",
                  type="numeric",
                  default=20)
p <- add_argument(p, "--Dh.b",
                  help="background elaboration. background error covariance matrix, horizontal decorrelation length scale (km)",
                  type="numeric",
                  default=70)
p <- add_argument(p, "--Dz.b",
                  help="background elaboration. background error covariance matrix, vertical decorrelation length scale (m)",
                  type="numeric",
                  default=100000)
p <- add_argument(p, "--eps2.b",
                  help="background elaboration. ratio of observation to background error covariance",
                  type="numeric",
                  default=0.5)
# Background parameters: selection of a centroid station 
p <- add_argument(p, "--nstat.min",
                  help="minimum number of station in a subdomain",
                  type="numeric",
                  default=5)
p <- add_argument(p, "--Lsubsample",
                  help="maximum number of stations in a subdomain",
                  type="numeric",
                  default=100)
p <- add_argument(p, "--Lsubsample.DHmax",
                  help="maximum distance between a centroid station and the other stations of a subdomain (km)",
                  type="numeric",
                  default=400)
p <- add_argument(p, "--centroid.DHmin",
                  help="minimum distance between two centroid stations (km)",
                  type="numeric",
                  default=50)
# paths
p <- add_argument(p, "--main.path",
                  help="path of the main program",
                  type="character",
                  default="none")
p <- add_argument(p, "--main.path.output",
                  help="path where to save the output",
                  type="character",
                  default="none")
p <- add_argument(p, "--case.path",
                  help="path to CASE",
                  type="character",
                  default="none")
p <- add_argument(p, "--geoinfo.path",
                  help="path to digital elevation model and grid mask",
                  type="character",
                  default="none")
#
argv <- parse_args(p)
#
#------------------------------------------------------------------------------
# set variable dependent netcdf fields 
if (argv$var=="TG") {
  varname<-c("tg")
  varlongname<-c("daily average air temperature (from 06 UTC prev day to 06 UTC day)")
  varstandardname<-c("air temperature")
  summary<-c("Nordic Gridded Climate Dataset (NGCD) daily mean temperature dataset (variable TG). The aggregation time interval ranges from 0600 UTC previous day to 0600 UTC day reported as timestamp. NGCD is an observational gridded dataset and the input data comes from the Norwegian Meteorological Institute Climate Database and the European Climate Assessment & Dataset (ecad.eu). For more information on the datasets see https://github.com/metno/CASE. NGCD is based on the software https://github.com/metno/NGCD.") 
  title<-"NGCD_TG"
} else if (argv$var=="TX") {
  varname<-c("tx")
  varunit<-c("Celsius degree")
  varlongname<-c("daily max air temperature (from 06 UTC prev day to 06 UTC day)")
  varstandardname<-c("air temperature")
  summary<-c("Nordic Gridded Climate Dataset (NGCD) daily maximum temperature dataset (variable TX). The aggregation time interval ranges from 0600 UTC previous day to 0600 UTC day reported as timestamp. NGCD is an observational gridded dataset and the input data comes from the Norwegian Meteorological Institute Climate Database and the European Climate Assessment & Dataset (ecad.eu). For more information on the datasets see https://github.com/metno/CASE. NGCD is based on the software https://github.com/metno/NGCD.") 
  title<-"NGCD_TX"
} else if (argv$var=="TN") {
  varname<-c("tn")
  varunit<-c("Celsius degree")
  varlongname<-c("daily min air temperature (from 06 UTC prev day to 06 UTC day)")
  varstandardname<-c("air temperature")
  summary<-c("Nordic Gridded Climate Dataset (NGCD) daily minimum temperature dataset (variable TN). The aggregation time interval ranges from 0600 UTC previous day to 0600 UTC day reported as timestamp. NGCD is an observational gridded dataset and the input data comes from the Norwegian Meteorological Institute Climate Database and the European Climate Assessment & Dataset (ecad.eu). For more information on the datasets see https://github.com/metno/CASE. NGCD is based on the software https://github.com/metno/NGCD.") 
  title<-"NGCD_TN"
} else {
  ext<-error_exit("Fatal Error: variable must be TG or TX or TN")
}
#
#------------------------------------------------------------------------------
# [] define/check paths and load external functions
if ( !(file.exists(argv$main.path))        | 
     !(file.exists(argv$main.path.output)) |
     !(file.exists(argv$case.path))        | 
     !(file.exists(argv$geoinfo.path)) ) 
  ext<-error_exit("path not found")
# geographical information
filenamedem<-file.path(argv$geoinfo.path,"fennodem_NGCD.nc")
if (!file.exists(filenamedem)) 
  ext<-error_exit(paste("File not found:",filenamedem))
filenamemask<-file.path(argv$geoinfo.path,"fennomask_NGCD.nc")
if (!file.exists(filenamemask)) 
  ext<-error_exit(paste("File not found:",filenamemask))
filenamelaf<-file.path(argv$geoinfo.path,"ngcd_1km_max_laf.nc")
if (!file.exists(filenamelaf)) 
  ext<-error_exit(paste("File not found:",filenamelaf))
# load functions
# load functions
path2lib.com<-paste(argv$main.path,"/Bspat/utilities",sep="")
source(paste(path2lib.com,"/Bspat_PseudoBackground.R",sep=""))
source(paste(path2lib.com,"/nc4out.R",sep=""))
source(paste(path2lib.com,"/OI_T_fast.R",sep=""))
source(paste(path2lib.com,"/OI_T_xb_fast.R",sep=""))
# load external C functions
dyn.load(paste(argv$main.path,"/Bspat/src/T/oi_first.so",sep=""))
dyn.load(paste(argv$main.path,"/Bspat/src/T/oi_fast.so",sep=""))
dyn.load(paste(argv$main.path,"/Bspat/src/T/oi_xb_fast.so",sep=""))
#
# set Time-related variables
Rdate <- strptime(argv$date,"%Y.%m.%d")
yyyy<-Rdate$year[1]+1900
mm<-Rdate$mon[1]+1
dd<-Rdate$mday[1]
nday<-1
yyyymmdd<-paste(yyyy,
                formatC(mm,width=2,flag="0"),
                formatC(dd,width=2,flag="0"),sep="")
yyyymm<-paste(yyyy,
              formatC(mm,width=2,flag="0"),sep="")
yyyymmdd0600<-paste(yyyymmdd,"0600",sep="")
#
#------------------------------------------------------------------------------
# Grid - it is defined by the DEM file
# CRS Coordinate Reference System
# laf and dem must be available over a larger region than the mask
rtmp1<-raster(filenamedem)
rtmp2<-raster(filenamemask)
rtmp3<-raster(filenamelaf)
r.laf<-mask(rtmp3,mask=rtmp2)/100
r.orog<-mask(rtmp1,mask=rtmp2)
rm(rtmp1,rtmp2,rtmp3)
nx<-ncol(r.orog)
ny<-nrow(r.orog)
dx<-xres(r.orog)
dy<-yres(r.orog)
# 4 borders point
xmn<-xmin(r.orog)
xmx<-xmax(r.orog)
ymn<-ymin(r.orog)
ymx<-ymax(r.orog)
# extract all the cell values: cells[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
xy<-xyFromCell(r.orog,1:ncell(r.orog))
x.G<-sort(unique(xy[,1]))
y.G<-sort(unique(xy[,2]),decreasing=T)
aux<-as.vector(extract(r.orog,1:ncell(r.orog)))
mask<-which(!is.na(aux))
zgrid<-aux[mask]
aux<-as.vector(extract(r.laf,1:ncell(r.laf)))
lgrid<-aux[mask]
xgrid<-xy[mask,1]
ygrid<-xy[mask,2]
Lgrid<-length(xgrid)
rm(xy,aux)
# debug info
if (argv$verbose) {
  print("+----------------------------+")
  print("+ grid parameters")
  print(paste("nx ny dx dy",
    as.integer(nx),as.integer(ny),round(dx,2),round(dy,2)))
  print(paste("xmn xmx ymn ymx",
    round(xmn,2),round(xmx,2),round(ymn,2),round(ymx,2)))
  print(paste("# grid points=",as.integer(Lgrid)))
}
#
#------------------------------------------------------------------------------
# [] Read data from CASE 
ffin<-file.path(argv$case.path,
                paste(argv$var,"_date_dqc",sep=""),
                yyyymm,
                paste("case_",argv$var,"_date_",yyyymmdd,".txt",sep="")) 
if (!file.exists(ffin))
  ext<-error_exit(paste("File not found:",ffin))
data<-read.table(file=ffin,header=T,sep=";",
                         stringsAsFactors=F,strip.white=T)
n.tmp<-length(data$staid)
# select observations to use
aux<-extract(r.orog,
     cbind(data$etrs_laea_x,data$etrs_laea_y),na.rm=T)
stn.output<-which(!is.na(aux) & 
                  !is.na(data$value) &
                  data$dqc %in% c(0) &
                  data$qcode %in% c(0,1,2))
n.stn.output<-length(stn.output)
rm(aux)
# definitive station list
L.y.tot<-length(stn.output)
VecX<-data$etrs_laea_x[stn.output]
VecY<-data$etrs_laea_y[stn.output]
VecLat<-data$lat[stn.output]
VecLon<-data$lon[stn.output]
VecZ<-data$elev[stn.output]
VecS<-data$staid[stn.output]
yo<-data$value[stn.output]
ydqc.flag<-rep(0,length=L.y.tot)
VecLaf<-extract(r.laf,cbind(VecX,VecY))
rm(data)
if (argv$verbose) {
  print("+----------------------------+")
  print(paste("# observations =",L.y.tot))
}
#
#------------------------------------------------------------------------------
# [] compute Disth (symmetric) matrix: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
Disth<-(outer(VecY,VecY,FUN="-")**2.+
        outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
#D.b<-exp(-0.5*(Disth/argv$Dh.b)**2.)
#S.b<-D.b
#diag(D.b)<-diag(D.b)+argv$eps2.b
#------------------------------------------------------------------------------
# Elaborations
yo.h.pos<-which(!is.na(yo))
# BACKGROUND AT STATION LOCATIONS
# For each station, compute a non-linear vertical profile using 
# the Lsubsample (closest) surrounding stations.
# The non-linear profile allows for a lat-lon dependance and
# for an inversion in the vertical profile.
#
# for each iteration, the yb is computed using only Lsubsample stations
# but over all stations (also the distant ones). A weight is then assigned
# for all the stations, giving more weight to the stations closest to the 
# actually used stations. This procedure should ensure a continuous background
yb.set<-matrix(data=NA,ncol=L.y.tot,nrow=L.y.tot)
ybweights.set<-matrix(data=0,ncol=L.y.tot,nrow=L.y.tot)
yb.param<-matrix(data=NA,ncol=12,nrow=L.y.tot)
yb.param1.mat<-matrix(data=NA,ncol=12,nrow=L.y.tot)
yb.param2.mat<-matrix(data=NA,ncol=12,nrow=L.y.tot)
VecS.set<-matrix(data=NA,ncol=argv$Lsubsample,nrow=L.y.tot)
VecS.set.pos<-matrix(data=NA,ncol=argv$Lsubsample,nrow=L.y.tot)
yb.h.pos<-vector()
Lsubsample.vec<-vector()
# cycle over the available stations
b.inc<-0
for (b in yo.h.pos) {
  # A. selection of the centroid station
  #  A1. it must not be too close to a previous centroid station
  if (b.inc>0) {
    finoa<-as.integer(argv$Lsubsample*1/10)
    if (VecS[b] %in% VecS.set[yb.h.pos[1:b.inc],2:finoa]) next
    distaux<-min(sqrt((VecX[b]-VecX[yb.h.pos[1:b.inc]])**2+
                      (VecY[b]-VecY[yb.h.pos[1:b.inc]])**2)
                 /1000.)
    if (distaux<argv$centroid.DHmin) next 
  }
  # select the closest argv$Lsubsample stations
  close2b.aux<-order(Disth[b,],decreasing=F)
  close2b.au1<-close2b.aux[which(close2b.aux%in%yo.h.pos)][1:argv$Lsubsample]
  #  A2. if the closest station is too far away then b is not a good candidate 
  if (Disth[b,close2b.au1[2]]>(2*argv$Dh.b)) next
  aux<-argv$Lsubsample
  if (Disth[b,close2b.au1[argv$Lsubsample]]>argv$Lsubsample.DHmax) 
    aux<-max(which(Disth[b,close2b.au1]<=argv$Lsubsample.DHmax))
  #  A3. if the number of stations within a pres-set radius is less than
  #      a pre-def threshold then b is not a good candidate 
  if (aux<argv$nstat.min) next
  # ==> b has been chosen as a centroid station <==
  print("+-------------------------------------------------------------------")
  # B. compute the background for the subdomain having b as centroid
  Lsubsample.vec[b]<-aux
  close2b<-close2b.au1[1:Lsubsample.vec[b]]
  b.inc<-b.inc+1
  yb.h.pos[b.inc]<-b
  yb.set0<-matrix(data=NA,nrow=L.y.tot)
  yb.set1<-matrix(data=NA,nrow=L.y.tot)
  yb.set2<-matrix(data=NA,nrow=L.y.tot)
  VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  VecX.b<-VecX[close2b]
  VecY.b<-VecY[close2b]
  VecZ.b<-VecZ[close2b]
  VecS.set[b,]<-NA
  VecS.set.pos[b,]<-NA
  VecS.set[b,1:Lsubsample.vec[b]]<-VecS[close2b]
  VecS.set.pos[b,1:Lsubsample.vec[b]]<-close2b
  yo.b<-yo[close2b]
  # ++ B1. Background 0  
  yb.param0<-XYZnoinv_step0(VecX.b,VecY.b,VecZ.b,yo.b)
  # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
  # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
  yb.set0[]<-yb.param0[3]+ yb.param0[6]*(VecX-yb.param0[10])+
                           yb.param0[8]*(VecY-yb.param0[11])+
                           yb.param0[4]* VecZ
  J0<-sum((yb.set0[close2b]-yo.b)**2.)
  z.range<-max(VecZ.b)-min(VecZ.b)
  z.range.90<-round((quantile(VecZ.b,probs=0.9)-quantile(VecZ.b,probs=0.1)),0)
  flag.only0<-(z.range.90<=argv$z.range.min | 
               Lsubsample.vec[b]<argv$nstat.min.inv)
  if (!flag.only0) {
   # ++ B2. Background 1  
    back1.flag<-T
    par.aux<-c(mean(VecZ.b),50,yb.param0[3],
               yb.param0[4],yb.param0[4],
               yb.param0[6],yb.param0[6],
               yb.param0[8],yb.param0[8])
    yb.param1<-XYZinv_step1(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
    yb.param1.mat[b,1:11]<-yb.param1
    if (!is.na(yb.param1[11])) {
      # 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
      # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
      zinv<-yb.param1[1]
      zabov<-zinv+yb.param1[2]
      zbelo<-zinv-yb.param1[2]
      bfabov<-yb.param1[3]+ yb.param1[6]*(VecX-yb.param1[10])+
                            yb.param1[8]*(VecY-yb.param1[11])+
                            yb.param1[4]*(VecZ-zinv)
      bfbelo<-yb.param1[3]+ yb.param1[7]*(VecX-yb.param1[10])+
                            yb.param1[9]*(VecY-yb.param1[11])+ 
                            yb.param1[5]*(VecZ-zinv)
      yb.set1[VecZ>zabov]<-bfabov[VecZ>zabov]
      yb.set1[VecZ<=zbelo]<-bfbelo[VecZ<=zbelo]
      aux<-(VecZ>zbelo)&(VecZ<=zabov)
      yb.set1[aux]<-(bfabov[aux]*(VecZ[aux]-zbelo) + 
                     bfbelo[aux]*(zabov-VecZ[aux]) ) / (zabov-zbelo)
      J1<-sum((yb.set1[close2b]-yo.b)**2.)
    } else {
      back1.flag<-F
      J1<-J0+1
    }
    # ++ B3. Background 2  
    # NOTE: yb.set* is defined for all the LOBS stations
    # 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
    # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
    back2.flag<-T
    if (back1.flag) {
      par.aux<-c(zinv,
                 (max(VecZ.b)-min(VecZ.b))/10,
                 yb.param1[3]+yb.param1[4]*(-zinv),
                 yb.param1[4],
                 0,
                 yb.param1[6],yb.param1[7],
                 yb.param1[8],yb.param1[9])
    } else {
      par.aux<-c(mean(VecZ.b),
                 (max(VecZ.b)-min(VecZ.b))/10,
                 yb.param0[3],
                 yb.param0[4],
                 0,
                 yb.param0[6],yb.param0[6],
                 yb.param0[8],yb.param0[8])
    }
    yb.param2<-XYZinv_step2(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
    yb.param2.mat[b,1:11]<-yb.param2
    if (!is.na(yb.param2[11])) {
      yb.param2[2]<-abs(yb.param2[2])
      ab<-which( VecZ>=(yb.param2[1]+yb.param2[2]) )
      bl<-which( VecZ<=yb.param2[1] )
      bw<-which( (VecZ>yb.param2[1]) & (VecZ<(yb.param2[1]+yb.param2[2])) )
      yb.set2[ab]<-yb.param2[3] + yb.param2[4]*VecZ[ab] + 
                     yb.param2[6]*(VecX[ab]-yb.param2[10]) + 
                     yb.param2[8]*(VecY[ab]-yb.param2[11])
      if (length(bl)>0)
        yb.set2[bl]<-yb.param2[3] + yb.param2[4]*VecZ[bl] - 
                       yb.param2[5] +
                       yb.param2[7]*(VecX[bl]-yb.param2[10]) + 
                       yb.param2[9]*(VecY[bl]-yb.param2[11])
      if (length(bw)>0)
        yb.set2[bw]<-yb.param2[3] + 
                     yb.param2[4]*VecZ[bw] - 
                     yb.param2[5]/2.*
                      (1+cos(pi*(VecZ[bw]-yb.param2[1])/yb.param2[2])) + 
                     ( (yb.param2[6]*(VecX[bw]-yb.param2[10])+
                        yb.param2[8]*(VecY[bw]-yb.param2[11]))*
                                     (VecZ[bw]-yb.param2[1]) + 
                        (yb.param2[7]*(VecX[bw]-yb.param2[10])+
                         yb.param2[9]*(VecY[bw]-yb.param2[11]))*
                                      (yb.param2[1]+yb.param2[2]-VecZ[bw]) ) / 
                      yb.param2[2]
      J2<-sum((yb.set2[close2b]-yo.b)**2.)
    } else {
      back2.flag<-F
      J2<-J0+1
    }
  } else {
    J1<-NA
    J2<-NA
    yb.param1<-rep(NA,12)
    yb.param2<-rep(NA,12)
  }
  # B4. select the Best background
  if (flag.only0) {
    best<-0
  } else {
    aux<-order(c(J0,J1,J2))
    best<-aux[1]-1
  }
  if (best==0) {
    yb.set[,b]<-yb.set0[]
    yb.param[b,]<-c(yb.param0,0)
  }
  if (best==1) {
    yb.set[,b]<-yb.set1[]
    yb.param[b,]<-c(yb.param1,1)
  }
  if (best==2) {
    yb.set[,b]<-yb.set2[]
    yb.param[b,]<-c(yb.param2,2)
  }
# DEBUG: plot maps showing the stations used for each sub-domain
#  png(file=paste("backg_",
#                 formatC(b.inc,width=4,flag="0"),
#                 ".png",sep=""),
#      width=800,height=800)
#  plot(VecX,VecY)
#  points(VecX.b,VecY.b,col="red",pch=19)
#  points(VecX[b],VecY[b],col="purple",pch=19,cex=2)
#  dev.off()
#  png(file=paste("backg_vert_",
#                 formatC(b.inc,width=4,flag="0"),
#                 ".png",sep=""),
#      width=800,height=800)
#  mnn<-min(c(yo.b,yb.set[close2b,b]))
#  mxx<-max(c(yo.b,yb.set[close2b,b]))
#  plot(yb.set[close2b,b],VecZ.b,col="blue",pch=19,
#       ylim=c(0,2100),xlim=c(mnn,mxx))
#  points(yo.b,VecZ.b,col="red",pch=19)
#  dev.off()
# DEBUG end
# DEBUG start
  zero.str<-"0."
  uno.str<-"1."
  due.str<-"2."
  if (best==0) zero.str<-"+0."
  if (best==1) uno.str<-"+1."
  if (best==2) due.str<-"+2."
  print(paste("@@",b.inc,
              ". id pos Zmn/x (Zrange) [Z.q90-Z.q10] ",
              "DisthMAX / J0 J1 J2 / #stn:",
              VecS[b],
              b,
              round(min(VecZ.b),0),
              round(max(VecZ.b),0),
              paste("(",z.range,")",sep=""),
              paste("[",z.range.90,"]",sep=""),
              round(Disth[b,close2b[Lsubsample.vec[b]]],0),"/",
              round(J0,0),round(J1,0),round(J2,0),"/",
              Lsubsample.vec[b]))
  print(paste(zero.str,
              "alpha beta gamma:",
               round(yb.param0[6],6),
               round(yb.param0[8],6),
               round(yb.param0[4],4)))
  print(paste(uno.str,"zinv dz tinv aA aB bA bB ga gb:",
        round(yb.param1[1],1), round(yb.param1[2],0), round(yb.param1[3],1),
        round(yb.param1[6],6), round(yb.param1[7],6), round(yb.param1[8],6),
        round(yb.param1[9],6), round(yb.param1[4],4), round(yb.param1[5],4)))
  print(paste(due.str,"h0 h1-h0 t0 a aA aB bA bB gamma:",
        round(yb.param2[1],1), round(yb.param2[2],0), round(yb.param2[3],1),
        round(yb.param2[5],6), round(yb.param2[6],6), round(yb.param2[7],6),
        round(yb.param2[8],6), round(yb.param2[9],6), round(yb.param2[4],4)))
# DEBUG stop
# B5. compute weights used to blend the local backgrounds into a global one
  G<-exp(-0.5*(Disth[,close2b]/argv$Dh.b)**2.)
  D<-G[close2b,]
  diag(D)<-diag(D)+argv$eps2.b
  InvD<-chol2inv(chol(D))
  K<-tcrossprod(G,InvD)
  rm(G)
  diag(D)<-diag(D)-argv$eps2.b
  W<-tcrossprod(D,InvD)
  rm(D,InvD)
  # this is idi (sort of)
  ybweights.set[,b]<-rowSums(K)
  if (any(is.na(ybweights.set[,b]))) {
    print("Error: NA value presents in ybweights.set vector. Unexpected.")
    print(paste(round(yb.set[,b],2),
                      round(yo,2),
                      round(ybweights.set[,b],3),"\n"))
    q(status=1)
  }
  rm(K,W)
} # end cycle LOBS
if (argv$verbose) {
  print(paste("# centroid stations or subdomains=",b.inc))
  # print sub-domain diagnostics: 
  #  the average number of stations in a sub-domain
  #  
}
# 
# At this point I've this 5 outputs
# 1. yb.set<-matrix(data=NA,ncol=btimes,nrow=LOBSt)
# 2. ybweights.set<-matrix(data=NA,ncol=btimes,nrow=LOBSt)
# 3. yb.param<-matrix(data=NA,ncol=12,nrow=LOBSt)
# 4. VecS.set<-matrix(data=NA,ncol=Lsubsample,nrow=LOBSt)
# 5. VecS.set.pos[b,]<-close2b
#====================================================================
# Station Analysis cycle (with SCT!)
D<-exp(-0.5*(Disth/argv$Dh)**2.-
        0.5*(abs(outer(VecZ,VecZ,FUN="-"))/argv$Dz)**2.) *
   (1-(1-argv$lafmin)*abs(outer(VecLaf,VecLaf,FUN="-")))
S<-D
diag(D)<-diag(D)+argv$eps2
#
yoBad.id<-NA
yoBad.pos<-NA
isct<-0
ydqc.flag[]<-rep(0,L.y.tot)
ydqc.flag[-yo.h.pos]<-(-1)
yo.OKh.pos<-which(!is.na(yo) & ydqc.flag!=1)
LOBStOK<-length(yo.OKh.pos)
#
yb<-vector(mode="numeric",length=length(yo))
ya<-vector(mode="numeric",length=length(yo))
yav<-vector(mode="numeric",length=length(yo))
yidi<-vector(mode="numeric",length=length(yo))
yidiv<-vector(mode="numeric",length=length(yo))
ydqc<-vector(mode="numeric",length=length(yo))
yb[]<-NA
ya[]<-NA
yav[]<-NA
yidi[]<-NA
yidiv[]<-NA
ydqc[]<-NA
safe.cont<-0
safe.limit<-100000
while (safe.cont<safe.limit) {
  safe.cont<-safe.cont+1
  print(paste(">> Total number of observations [not NA & not ERR(so far)] =",
              LOBStOK))
  # In case of a rejected station, re-compute those local backgrounds that
  #  include the rejected stations
  if (!is.na(yoBad.id)) {
    isct<-isct+1
    for (b in yb.h.pos) {
      if (any(VecS.set[b,1:Lsubsample.vec[b]]==yoBad.id)) {
        yb.set0<-matrix(data=NA,nrow=L.y.tot)
        yb.set1<-matrix(data=NA,nrow=L.y.tot)
        yb.set2<-matrix(data=NA,nrow=L.y.tot)
        Lsubsample.vec[b]<-Lsubsample.vec[b]-1
        VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        close2b.aux<-order(Disth[b,],decreasing=F)
        close2b<-close2b.aux[which(close2b.aux%in%yo.OKh.pos)][1:Lsubsample.vec[b]]
        VecX.b<-VecX[close2b]
        VecY.b<-VecY[close2b]
        VecZ.b<-VecZ[close2b]
        VecS.set[b,]<-NA
        VecS.set.pos[b,]<-NA
        VecS.set[b,1:Lsubsample.vec[b]]<-VecS[close2b]
        VecS.set.pos[b,1:Lsubsample.vec[b]]<-close2b
        yo.b<-yo[close2b]
        # Background 0  
        yb.param0<-XYZnoinv_step0(VecX.b,VecY.b,VecZ.b,yo.b)
        # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
        # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
        yb.set0[]<-yb.param0[3]+ yb.param0[6]*(VecX-yb.param0[10])+
                                 yb.param0[8]*(VecY-yb.param0[11])+
                                 yb.param0[4]*VecZ
        J0<-sum((yb.set0[close2b]-yo.b)**2.)
        z.range<-max(VecZ.b)-min(VecZ.b)
        z.range.90<-round((quantile(VecZ.b,probs=0.9)-
                           quantile(VecZ.b,probs=0.1)),0)
        flag.only0<-(z.range.90<=argv$z.range.min | 
                     Lsubsample.vec[b]<argv$nstat.min.inv)
        if (!flag.only0) {
          # Background 1  
          back1.flag<-T
          if (!is.na(yb.param1.mat[b,11])) {
            par.aux<-yb.param1.mat[b,1:9]
          } else {
            par.aux<-c(mean(VecZ.b),50,yb.param0[3],
                       yb.param0[4],yb.param0[4],
                       yb.param0[6],yb.param0[6],
                       yb.param0[8],yb.param0[8])
          }
          yb.param1<-XYZinv_step1(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
          yb.param1.mat[b,1:11]<-yb.param1
          if (!is.na(yb.param1[11])) {
            # 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
            # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
            zinv<-yb.param1[1]
            zabov<-zinv+yb.param1[2]
            zbelo<-zinv-yb.param1[2]
            bfabov<-yb.param1[3]+ yb.param1[6]*(VecX-yb.param1[10])+
                                  yb.param1[8]*(VecY-yb.param1[11])+
                                  yb.param1[4]*(VecZ-zinv)
            bfbelo<-yb.param1[3]+ yb.param1[7]*(VecX-yb.param1[10])+
                                  yb.param1[9]*(VecY-yb.param1[11])+ 
                                  yb.param1[5]*(VecZ-zinv)
            yb.set1[VecZ>zabov]<-bfabov[VecZ>zabov]
            yb.set1[VecZ<=zbelo]<-bfbelo[VecZ<=zbelo]
            aux<-(VecZ>zbelo)&(VecZ<=zabov)
            yb.set1[aux]<-(bfabov[aux]*(VecZ[aux]-zbelo) + bfbelo[aux]*(zabov-VecZ[aux]) ) / (zabov-zbelo)
            J1<-sum((yb.set1[close2b]-yo.b)**2.)
          } else {
            back1.flag<-F
            J1<-J0+1
          }
          # Background 2  
          # NOTE: yb.set* is defined for all the LOBSt stations
          # 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
          # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
          back2.flag<-T
          if (!is.na(yb.param2.mat[b,11])) {
            par.aux<-yb.param2.mat[b,1:9]
          } else {
            if (back1.flag) {
              par.aux<-c(zinv,
                     (max(VecZ.b)-min(VecZ.b))/10,
                     yb.param1[3]+yb.param1[4]*(-zinv),
                     yb.param1[4],
                     0,
                     yb.param1[6],yb.param1[7],
                     yb.param1[8],yb.param1[9])
            } else {
              par.aux<-c(mean(VecZ.b),
                         (max(VecZ.b)-min(VecZ.b))/10,
                         yb.param0[3],
                         yb.param0[4],
                         0,
                         yb.param0[6],yb.param0[6],
                         yb.param0[8],yb.param0[8])
            }
          }
          yb.param2<-XYZinv_step2(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
          yb.param2.mat[b,1:11]<-yb.param2
          if (!is.na(yb.param2[11])) {
            yb.param2[2]<-abs(yb.param2[2])
            ab<-which( VecZ>=(yb.param2[1]+yb.param2[2]) )
            bl<-which( VecZ<=yb.param2[1] )
            bw<-which( (VecZ>yb.param2[1]) & (VecZ<(yb.param2[1]+yb.param2[2])) )
            yb.set2[ab]<-yb.param2[3] + yb.param2[4]*VecZ[ab] + 
                           yb.param2[6]*(VecX[ab]-yb.param2[10]) + 
                           yb.param2[8]*(VecY[ab]-yb.param2[11])
            if (length(bl)>0)
              yb.set2[bl]<-yb.param2[3] + yb.param2[4]*VecZ[bl] - 
                             yb.param2[5] +
                             yb.param2[7]*(VecX[bl]-yb.param2[10]) + 
                             yb.param2[9]*(VecY[bl]-yb.param2[11])
            if (length(bw)>0)
              yb.set2[bw]<-yb.param2[3] + yb.param2[4]*VecZ[bw] - 
                             yb.param2[5]/2.*
                             (1+cos(pi*(VecZ[bw]-yb.param2[1])/yb.param2[2])) + 
                            ( (yb.param2[6]*(VecX[bw]-yb.param2[10])+
                               yb.param2[8]*(VecY[bw]-yb.param2[11]))*
                               (VecZ[bw]-yb.param2[1]) + 
                              (yb.param2[7]*(VecX[bw]-yb.param2[10])+
                               yb.param2[9]*(VecY[bw]-yb.param2[11]))*
                               (yb.param2[1]+yb.param2[2]-VecZ[bw]) ) / 
                              yb.param2[2]
            J2<-sum((yb.set2[close2b]-yo.b)**2.)
          } else {
            back2.flag<-F
            J2<-J0+1
          }
        }
        # Best background
        if (flag.only0) {
          best<-0
        } else {
          aux<-order(c(J0,J1,J2))
          best<-aux[1]-1
        }
        if (best==0) {
          yb.set[,b]<-yb.set0[]
          yb.param[b,]<-c(yb.param0,0)
        }
        if (best==1) {
          yb.set[,b]<-yb.set1[]
          yb.param[b,]<-c(yb.param1,1)
        }
        if (best==2) {
          yb.set[,b]<-yb.set2[]
          yb.param[b,]<-c(yb.param2,2)
        }
        # DEBUG start
        zero.str<-"0."
        uno.str<-"1."
        due.str<-"2."
        if (best==0) zero.str<-"*0."
        if (best==1) uno.str<-"*1."
        if (best==2) due.str<-"*2."
        print(paste("@@",isct,
                    ". id pos Zmn/x (Zrange) [Z.q90-Z.q10] ",
                    "DisthMAX / J0 J1 J2 / #stn:",VecS[b],b,
                    round(min(VecZ.b),0),round(max(VecZ.b),0),
                    paste("(",z.range,")",sep=""),
                    paste("[",z.range.90,"]",sep=""),
                    round(Disth[b,close2b[Lsubsample.vec[b]]],0),"/",
                    round(J0,0),round(J1,0),round(J2,0),"/",Lsubsample.vec[b]))
        print(paste(zero.str,"alpha beta gamma:",
                    round(yb.param0[6],6),round(yb.param0[8],6),
                    round(yb.param0[4],4)))
        print(paste(uno.str,"zinv dz tinv aA aB bA bB ga gb:",
           round(yb.param1[1],1), round(yb.param1[2],0), round(yb.param1[3],1),
           round(yb.param1[6],6), round(yb.param1[7],6), round(yb.param1[8],6),
           round(yb.param1[9],6), round(yb.param1[4],4), round(yb.param1[5],4)))
        print(paste(due.str,"h0 h1-h0 t0 a aA aB bA bB gamma:",
           round(yb.param2[1],1), round(yb.param2[2],0), round(yb.param2[3],1),
           round(yb.param2[5],6), round(yb.param2[6],6), round(yb.param2[7],6),
           round(yb.param2[8],6), round(yb.param2[9],6), round(yb.param2[4],4)))
        # DEBUG stop
        # update weights
        # G1 is LOBStOK x Lsubsample matrix to interpolate the Lsubsample values
        # on the whole LOBStOK station dataset
        G<-exp(-0.5*(Disth[,close2b]/argv$Dh.b)**2.)
        Da<-G[close2b,]
        diag(Da)<-diag(Da)+argv$eps2.b
        InvDa<-chol2inv(chol(Da))
        K<-tcrossprod(G,InvDa)
        rm(G)
        diag(Da)<-diag(Da)-argv$eps2.b
        W<-tcrossprod(Da,InvDa)
        rm(Da,InvDa)
        # this is idi (sort of)
        ybweights.set[,b]<-rowSums(K)
        if (any(is.na(ybweights.set[,b]))) {
          print("Error: NA value presents in ybweights.set vector. Unexpected.")
          print(paste(round(yb.set[,b],2),round(yo,2),
                      round(ybweights.set[,b],3),"\n"))
          q(status=1)
        }
        rm(K,W)
      }
    } # end cycle btimes (no optimization: over the station number)
  } # end cycle to recompute the background in case of a rejected station
# normalization of the yb weights such that their sum equals to one
  ybweights.norm<-ybweights.set[,yb.h.pos] / 
                  rowSums(ybweights.set[,yb.h.pos])
# background on station points
  for (m in 1:L.y.tot) {
    yb[m]<-tcrossprod(ybweights.norm[m,],t(yb.set[m,yb.h.pos]))
  }
# deallocate memory
#      rm(D.b,S.b,yb.set,ybweights.set,ybweights.norm)
#      rm(VecY.b,VecX.b,VecZ.b,yo.b,close2b,y1.b)
# Station (CV)Analysis/(CV)IDI
  InvD<-chol2inv(chol(D[yo.OKh.pos,yo.OKh.pos]))
  W<-tcrossprod(S[yo.OKh.pos,yo.OKh.pos],InvD)
  K<-tcrossprod(S[,yo.OKh.pos],InvD)
  rm(InvD)
  ya<-yb + tcrossprod(K,t(yo[yo.OKh.pos]-yb[yo.OKh.pos]))
  yav<-ya
  yav[yo.OKh.pos]<-yo[yo.OKh.pos] + 
                   1./(1.-diag(W)) * (ya[yo.OKh.pos]-yo[yo.OKh.pos])
  yidi<-rowSums(K)
  yidiv<-yidi
  yidiv[yo.OKh.pos]<-rep(1,LOBStOK) + 
                     1./(1.-diag(W)) * (yidi[yo.OKh.pos]-rep(1,LOBStOK))
# DQC CHECK - Spatial Continuity Check
  ydqc[]<--9999.
  ydqc[yo.OKh.pos]<-(yo[yo.OKh.pos]-yav[yo.OKh.pos]) *
                    (yo[yo.OKh.pos]-ya[yo.OKh.pos]) / argv$sig2o
  rm(W,K)
# DQC test: T2 is the SCT threshold, if any(ydqc>T2) is true then reject the 
# station with highest ydqc value. This mean that the station is not used to
# compute the analysis
  if (max(ydqc)<=argv$T2) break # EXIT the SCT loop
  aux<-which.max(ydqc)
  yoBad.id<-VecS[aux]
  yoBad.pos<-aux
  ydqc.flag[yoBad.pos]<-1
  yo.OKh.pos<-which(!is.na(yo) & ydqc.flag!=1)
  LOBStOK<-length(yo.OKh.pos)
# DQC flag "1" means erroneous observation
##      OBS$DQC[indx[aux]]<-1
  print(paste("SCT found GE-> id yo yb ya yav ydqc",VecS[aux],
              round(yo[aux],1),
              round(yb[aux],2),
              round(ya[aux],2),
              round(yav[aux],2),
              round(ydqc[aux],2),"\n"))
# write output file
} # End of Station Analysis cycle (with SCT!)
if (argv$verbose)
  print(paste(">>>>Total number of observations [not NA & good] =",LOBStOK))
# if something strange take place in the Analysis/SCT cycle then exit
if (safe.cont==safe.limit) error_exit("error during the SCT/analysis cycle")
if (LOBStOK<=0) error_exit("no valid observations found after the SCT")
#------------------------------------------------------------------------------
# Grid - Background
xb<-vector(mode="numeric",length=Lgrid)
xb.tmp<-vector(mode="numeric",length=Lgrid)
xidi.tmp<-vector(mode="numeric",length=Lgrid)
xidi.norm<-vector(mode="numeric",length=Lgrid)
xb[]<-0
xb.tmp[]<-0
xidi.tmp[]<-0
xidi.norm[]<-0
b.aux<-0
if (argv$verbose) print("++ Grid - Background elaborations\n")
for (b in yb.h.pos) {
  b.aux<-b.aux+1
  D.b<-exp(-0.5*(Disth[VecS.set.pos[b,1:Lsubsample.vec[b]],
                       VecS.set.pos[b,1:Lsubsample.vec[b]]]/argv$Dh.b)**2.)
  diag(D.b)<-diag(D.b)+argv$eps2.b
  InvD.b<-chol2inv(chol(D.b))
  rm(D.b)
  oi<-OI_T_xb_fast(b.param=yb.param[b,1:12],
                   xgrid.sel=xgrid,
                   ygrid.sel=ygrid,
                   zgrid.sel=zgrid,
                   VecX.sel=VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]],
                   VecY.sel=VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]],
                   VecZ.sel=VecZ[VecS.set.pos[b,1:Lsubsample.vec[b]]],
                   Dh.cur=argv$Dh.b,
                   Dz.cur=argv$Dz.b)
  rm(InvD.b)
  xb.tmp[]<-oi$xb
  xidi.tmp[]<-oi$xidi
  rm(oi)
  xb<-xb+xidi.tmp*xb.tmp
  xidi.norm<-xidi.norm+xidi.tmp
}
rm(Disth)
aux<-which(xidi.norm!=0)
if (length(aux!=0)) xb[aux]<-xb[aux]/xidi.norm[aux]
aux<-which(xidi.norm==0)
if (length(aux!=0)) xb[aux]<-NA
rm(xb.tmp,xidi.tmp,xidi.norm,aux)
# Gridded Analysis/IDI
InvD<-chol2inv(chol(D[yo.OKh.pos,yo.OKh.pos]))
rm(D)
oi<-OI_T_fast(yo.sel=yo[yo.OKh.pos],
              yb.sel=yb[yo.OKh.pos],
              xb.sel=xb,
              xgrid.sel=xgrid,
              ygrid.sel=ygrid,
              zgrid.sel=zgrid,
              lgrid.sel=lgrid,
              VecX.sel=VecX[yo.OKh.pos],
              VecY.sel=VecY[yo.OKh.pos],
              VecZ.sel=VecZ[yo.OKh.pos],
              VecLaf.sel=VecLaf[yo.OKh.pos],
              Dh.cur=argv$Dh,
              Dz.cur=argv$Dz,
              lafmin=argv$lafmin) 
xa<-vector(mode="numeric",length=Lgrid)
xa[]<-oi$xa
rm(oi,InvD)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (argv$verbose) print("++ Output")
# set output
dir.create(file.path(argv$main.path.output,"NGCD"), showWarnings = FALSE)
path2output.main<-file.path(argv$main.path.output,"NGCD",argv$var)
path2output.main.stn<-file.path(path2output.main,"station_dataset",yyyymm)
path2output.main.grd<-file.path(path2output.main,"gridded_dataset",yyyymm)
if (!(file.exists(path2output.main)))     
  dir.create(path2output.main,recursive=T,showWarnings=F) 
if (!(file.exists(path2output.main.stn))) 
  dir.create(path2output.main.stn,recursive=T,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) 
  dir.create(path2output.main.grd,recursive=T,showWarnings=F) 
# output files 
out.file.stn<-file.path(path2output.main.stn,
                 paste("NGCD_",argv$var,"_station_",yyyymmdd,".txt",sep=""))
out.file.verbg<-file.path(path2output.main.stn,
                 paste("NGCD_",argv$var,"_verbg_",yyyymmdd,".txt",sep=""))
out.file.veran<-file.path(path2output.main.stn,
                 paste("NGCD_",argv$var,"_veran_",yyyymmdd,".txt",sep=""))
out.file.vercv<-file.path(path2output.main.stn,
                 paste("NGCD_",argv$var,"_vercv_",yyyymmdd,".txt",sep=""))
out.file.grd<-file.path(path2output.main.grd,
                 paste("NGCD_",argv$var,"_grid_",yyyymmdd,".nc",sep=""))
# Station Points
# background for verif
cat("# variable: Temperature\n",file=out.file.verbg,append=F)
cat("# units: $^oC$\n",file=out.file.verbg,append=T)
cat("date     leadtime location  lat     lon      altitude obs      fcst\n",
    file=out.file.verbg,append=T)
ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(yb))
cat(paste( yyyymmdd," ",
           rep(0,length(ix))," ",
           VecS[ix]," ",
           VecLat[ix]," ",
           VecLon[ix]," ",
           VecZ[ix]," ",
           round(yo[ix],1)," ",
           round(yb[ix],2),"\n",
           sep=""),
    file=out.file.verbg,append=T)
print(paste("data saved on file",out.file.verbg))
# analysis for verif
cat("# variable: Temperature\n",file=out.file.veran,append=F)
cat("# units: $^oC$\n",file=out.file.veran,append=T)
cat("date     leadtime location  lat     lon      altitude obs      fcst\n",
    file=out.file.veran,append=T)
ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya))
cat(paste( yyyymmdd," ",
           rep(0,length(ix))," ",
           VecS[ix]," ",
           VecLat[ix]," ",
           VecLon[ix]," ",
           VecZ[ix]," ",
           round(yo[ix],1)," ",
           round(ya[ix],2),"\n",
           sep=""),
    file=out.file.veran,append=T)
print(paste("data saved on file",out.file.veran))
# CVanalysis for verif 
cat("# variable: Temperature\n",file=out.file.vercv,append=F)
cat("# units: $^oC$\n",file=out.file.vercv,append=T)
cat("date     leadtime location  lat     lon      altitude obs      fcst\n",
    file=out.file.vercv,append=T)
ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(yav))
cat(paste( yyyymmdd," ",
           rep(0,length(ix))," ",
           VecS[ix]," ",
           VecLat[ix]," ",
           VecLon[ix]," ",
           VecZ[ix]," ",
           round(yo[ix],1)," ",
           round(yav[ix],2),"\n",
           sep=""),
    file=out.file.vercv,append=T)
print(paste("data saved on file",out.file.vercv))
# output at station locations
cat("yyyy;mm;dd;stid;x;y;z;yo;yb;ya;yav;yidi;yidiv;dqc;\n",
    file=out.file.stn,append=F)
cat(paste(yyyy,mm,dd,
          formatC(VecS,format="f",digits=0),
          formatC(VecX,format="f",digits=0),
          formatC(VecY,format="f",digits=0),
          formatC(VecZ,format="f",digits=0),
          formatC(yo,format="f",digits=1),
          formatC(yb,format="f",digits=2),
          formatC(ya,format="f",digits=2),
          formatC(yav,format="f",digits=2),
          formatC(yidi,format="f",digits=4),
          formatC(yidiv,format="f",digits=4),
          formatC(ydqc.flag,format="f",digits=0),
          "\n",sep=";"),
    file=out.file.stn,append=T)
print(paste("data saved on file",out.file.stn))
# grid
ra<-r.orog
ra[]<-NA
ra[mask]<-round(xa,1)
r.list<-list()
r.list[[1]]<-matrix(data=getValues(ra),
                    ncol=length(y.G),
                    nrow=length(x.G))
out<-nc4out(grid.list=r.list,
            times=yyyymmdd0600,
            file.name=out.file.grd,
            grid.type="ngcd",
            x=x.G,
            y=y.G,
            var.name=varname,
            var.longname=varlongname,
            var.standardname=varstandardname,
            var.version=varversion,
            times.unit="D",
            reference=reference,
            proj4.string=proj4.ETRS_LAEA,
            var.unit=varunit,
            lonlat.out=T,
            round.dig=diground,
            summary=summary,
            source.string=sourcestring,
            title=title,
            comment=comment,
            atts.var.add=NULL)
print(paste("data saved on file",out.file.grd))
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# normal exit
t1<-Sys.time()
print(paste("Normal exit, time=",round(t1-t0,1),attr(t1-t0,"unit")))
quit(status=0)

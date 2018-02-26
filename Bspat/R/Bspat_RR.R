#!/usr/bin/env Rscript
# --~- Bspat_RR.R -~--
# Bayesian Spatial Interpolation of daily cumulated precipitation
# See the software repository here: https://github.com/metno/NGCD
#..............................................................................
#Copyright and license
# Copyright (C) 2018 MET Norway. NGCD software is licensed under GPL version 3 
# or (at your option) any later version.
# https://www.gnu.org/licenses/gpl-3.0.en.html
# 
# History:
# 20.02.2018 - Cristian Lussana. Original code.
# -----------------------------------------------------------------------------
#
rm(list=ls())
#
# -----------------------------------------------------------------------------
# Libraries
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("raster"))
suppressPackageStartupMessages(library("rgdal"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("tripack"))
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("ncdf4"))
options(warn = 2, scipen = 999)
# 
# -----------------------------------------------------------------------------
# Constants
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# DQC - identify adjacent nodes - Establish (local-) triangulation (Delauney)
n.sector<-16
sector.angle<-360/n.sector
# netcdf fixed parameters
varname<-c("rr")
varunit<-c("mm")
varlongname<-c("daily total precipitation (from 06 UTC prev day to 06 UTC day)")
varstandardname<-c("daily total precipitation")
varversion<-c("1.0")
reference<-c("Lussana, C., Saloranta, T., Skaugen, T., Magnusson, J., Tveito, O. E., and Andersen, J.: seNorge2 daily precipitation, an observational gridded dataset over Norway from 1957 to the present day, Earth Syst. Sci. Data, 10, 235-249, https://doi.org/10.5194/essd-10-235-2018, 2018.")
diground<-1
summary<-c("Nordic Gridded Climate Dataset (NGCD) daily precipitation dataset (variable RR). The aggregation time interval ranges from 0600 UTC previous day to 0600 UTC day reported as timestamp. NGCD is an observational gridded dataset and the input data comes from the Norwegian Meteorological Institute Climate Database and the European Climate Assessment & Dataset (ecad.eu). For more information on the datasets see https://github.com/metno/CASE. NGCD is based on the software https://github.com/metno/NGCD.") 
sourcestring<-"MET Norway"
title<-"NGCD_RR"
comment<-"Our open data are licensed under Norwegian Licence for Open Government Data (NLOD) or a Creative Commons Attribution 4.0 International License at your preference. Credit should be given to The Norwegian Meteorological institute, shortened “MET Norway”, as the source of data."
#
#..............................................................................
# Functions
# + manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}

# + write output files in case of no-precipitation over the whole domain
writeIfNoPrec<-function(){
  cat("# variable: Precipitation\n",file=out.file.ver,append=F)
  cat("# units: $mm$\n",file=out.file.ver,append=T)
  cat("date     leadtime location  lat     lon      altitude obs      fcst\n",
      file=out.file.ver,append=T)
  ix<-which(ydqc.flag<=0 & !is.na(yo))
  cat(paste( yyyymmdd," ",
             rep(0,length(ix))," ",
             VecS[ix]," ",
             VecLat[ix]," ",
             VecLon[ix]," ",
             VecZ[ix]," ",
             round(yo[ix],1)," ",
             rep(0,length(ix)),"\n",
             sep=""),
      file=out.file.ver,append=T)
  print(paste("data saved on file",out.file.ver))
  cat("yyyy;mm;dd;stid;x;y;z;yo;ya;dqc;\n",
      file=out.file.stn,append=F)
  cat(paste(yyyy,mm,dd,
            formatC(VecS,format="f",digits=0),
            formatC(VecX,format="f",digits=0),
            formatC(VecY,format="f",digits=0),
            formatC(VecZ,format="f",digits=0),
            formatC(yo,format="f",digits=1),
            rep(0,length(yo)), #ya
            formatC(ydqc.flag,format="f",digits=0),
            "\n",sep=";"),
      file=out.file.stn,append=T)
  print(paste("data saved on file",out.file.stn))
  # Figures
  r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG,
                    xmn=xmn.FG, xmx=xmx.FG,
                    ymn=ymn.FG, ymx=ymx.FG,
                    crs=proj4.ETRS_LAEA)
  xa.FG<-vector(mode="numeric",length=Lgrid.FG)
  xa.FG[]<-0
  r.aux.FG[mask.FG]<-round(xa.FG,1)
  r.list<-list()
  r.list[[1]]<-matrix(data=getValues(r.aux.FG),
                      ncol=length(y.FG),
                      nrow=length(x.FG))
  out<-nc4out(grid.list=r.list,
              times=yyyymmdd0600,
              file.name=out.file.grd,
              grid.type="ngcd",
              x=x.FG,
              y=y.FG,
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
}

#+ remove clumps determined by less than n observations
remove_ltNobsClumps<-function(r,n,reverse=F,obs) {
  obsflag<-obs$x
  obsflag[]<-0
  d<-getValues(r)
  if (reverse) {
    dd<-d
    d[dd==0]<-1
    d[dd>0]<-0
    r[]<-d
  }
  rc<-clump(r)
  dc<-getValues(rc)
  f<-freq(rc)
  ytmp<-extract(rc,cbind(obs$x,obs$y),na.rm=T)
  sel<-which(!is.na(ytmp))
  for (i in 1:length(f[,1])) {
    if (is.na(f[i,1])) next
    if (any(ytmp[sel]==f[i,1])) {
      ix<-which(ytmp[sel]==f[i,1])
      if (length(ix)<n) {
        dc[dc==f[i,1]]<-NA
        obsflag[sel[ix]]<-1
      }
    }
  }
  dc[is.na(dc)]<-(-1)
  if (any(dc==0)) d[which(dc==0)]<-0
  if (reverse) {
    dd<-d
    d[dd==0]<-1
    d[dd>0]<-0
  }
  r[]<-d
  return(list(r=r,obsflag=obsflag))
}
#
#==============================================================================
# MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN -
#==============================================================================
t0<-Sys.time()
# [] Read command line arguments and/or set parameters to default
# create parser object
p <- arg_parser("ngcdrr")
# specify our desired options 
# by default ArgumentParser will add an help option 
p <- add_argument(p, "date",
                  help="date (format YYYY.MM.DD)",
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
#
p <- add_argument(p, "--eps2",
                  help="ratio of observation to background error covariance",
                  type="numeric",
                  default=0.1)
#
p <- add_argument(p, "--rrinf",
                  help="precipitation yes/no threshold",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--fact",
                  help="aggregation factor to speed up the elaboration",
                  type="numeric",
                  default=2)
# DQC
p <- add_argument(p, "--Lsubsample.max",
                  help="check for dry(wet) station surrounded only by wet(dry) stations: max number of stations defining a neighbourhood to use",
                  type="numeric",
                  default=20)
p <- add_argument(p, "--Lsubsample.DHmax",
                  help="check for dry(wet) station surrounded only by wet(dry) stations: max distance (km) defining a neighbourhood to use",
                  type="numeric",
                  default=150)
p <- add_argument(p, "--DQC.min.dist.allowed",
                  help="check for dry(wet) station surrounded only by wet(dry) stations: max distance (km) form a neighbour for an observation to be flagged",
                  type="numeric",
                  default=150)
# multi-scale OI
p <- add_argument(p, "--min.Dh.seq.allowed",
                  help="smaller length-scale used for interpolation (km)",
                  type="numeric",
                  default=5)
p <- add_argument(p, "--Dh.seq.ref",
                  help="auxiliary length-scale used for interpolation (km)",
                  type="numeric",
                  default=50000)
p <- add_argument(p, "--Dh.seq.reference",
                  help="sequence of decreasing length-scales for interpolation (km)",
                  type="numeric",
                  default=c(5000,4000,3000,2000,1000,
                            900,800,700,600,500,
                            400,300,200,150,100,
                            80,60,50,40,30,20,10,5),
                  nargs=Inf)
p <- add_argument(p, "--Dz.seq.gt.ref",
                  help="sequence of decreasing vertical length-scales for interpolation (km). Case of horizontal length scale greater than Dh.seq.ref",
                  type="numeric",
                  default=c(100000),
                  nargs=Inf)
p <- add_argument(p, "--Dz.seq.le.ref",
                  help="sequence of decreasing vertical length-scales for interpolation (km). Case of horizontal length scale smaller or equal to Dh.seq.ref",
                  type="numeric",
                  default=c(10000),
                  nargs=Inf)
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
# load functions
path2lib.com<-paste(argv$main.path,"/Bspat/utilities",sep="")
source(paste(path2lib.com,"/RR_dqc_isolatedDry.R",sep=""))
source(paste(path2lib.com,"/RR_dqc_isolatedWet.R",sep=""))
source(paste(path2lib.com,"/nc4out.R",sep=""))
source(paste(path2lib.com,"/OI_RR_fast.R",sep=""))
# load external C functions
dyn.load(paste(argv$main.path,"/Bspat/src/RR/oi_first.so",sep=""))
dyn.load(paste(argv$main.path,"/Bspat/src/RR/oi_fast.so",sep=""))
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
# set output
dir.create(file.path(argv$main.path.output,"NGCD"), showWarnings = FALSE)
path2output.main<-file.path(argv$main.path.output,"NGCD","RR")
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
                 paste("NGCD_RR_station_",yyyymmdd,".txt",sep=""))
out.file.ver<-file.path(path2output.main.stn,
                 paste("NGCD_RR_verif_",yyyymmdd,".txt",sep=""))
out.file.grd<-file.path(path2output.main.grd,
                 paste("NGCD_RR_grid_",yyyymmdd,".nc",sep=""))
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [] Grid
# CRS Coordinate Reference System
rtmp1<-raster(filenamedem)
rtmp2<-raster(filenamemask)
r.orog.FG<-mask(rtmp1,mask=rtmp2)
rm(rtmp1,rtmp2)
nx.FG<-ncol(r.orog.FG)
ny.FG<-nrow(r.orog.FG)
dx.FG<-xres(r.orog.FG)
dy.FG<-yres(r.orog.FG)
# 4 borders point (SW corner (xmn,ymn); NE corner (xmx,ymx))
xmn.FG<-xmin(r.orog.FG)
xmx.FG<-xmax(r.orog.FG)
ymn.FG<-ymin(r.orog.FG)
ymx.FG<-ymax(r.orog.FG)
#
r.orog.CG<-aggregate(r.orog.FG,fact=argv$fact,expand=T)
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
#zgrid<-zvalues.FG[mask.FG]
#xgrid<-xy.FG[mask.FG,1]
#ygrid<-xy.FG[mask.FG,2]
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
rm(zvalues.CG,zvalues.FG,xy.FG,xy.CG)
# debug info
if (argv$verbose) {
  print("+----------------------------+")
  print("+ fine grid parameters")
  print(paste("nx ny dx dy",
    as.integer(nx.FG),as.integer(ny.FG),round(dx.FG,2),round(dy.FG,2)))
  print(paste("xmn xmx ymn ymx",
    round(xmn.FG,2),round(xmx.FG,2),round(ymn.FG,2),round(ymx.FG,2)))
  print(paste("# grid points=",as.integer(Lgrid.FG)))
  print("+ coarse grid parameters")
  print(paste("nx.CG ny.CG dx.CG dy.CG",
    as.integer(nx.CG),as.integer(ny.CG),round(dx.CG,2),round(dy.CG,2)))
  print(paste("xmn.CG xmx.CG ymn.CG ymx.CG",
    round(xmn.CG,2),round(xmx.CG,2),round(ymn.CG,2),round(ymx.CG,2)))
  print(paste("# grid points=",as.integer(Lgrid.CG)))
}
#
#------------------------------------------------------------------------------
# [] Read data from CASE 
ffin<-file.path(argv$case.path,
                "RR_date_dqc",
                yyyymm,
                paste("case_RR_date_",yyyymmdd,".txt",sep="")) 
if (!file.exists(ffin))
  ext<-error_exit(paste("File not found:",ffin))
data<-read.table(file=ffin,header=T,sep=";",
                         stringsAsFactors=F,strip.white=T)
n.tmp<-length(data$staid)
# select observations to use
aux1<-extract(r.orog.CG,
      cbind(data$etrs_laea_x,data$etrs_laea_y),na.rm=T)
aux2<-extract(r.orog.FG,
      cbind(data$etrs_laea_x,data$etrs_laea_y),na.rm=T)
stn.output<-which(!is.na(aux1) & 
                  !is.na(aux2) &
                  !is.na(data$value) &
                  data$dqc %in% c(0) &
                  data$qcode %in% c(0,1,2))
n.stn.output<-length(stn.output)
rm(aux1,aux2)
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
rm(data)
#
#------------------------------------------------------------------------------
# [] compute Disth (symmetric) matrix: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
Disth<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
Disth<-(outer(VecY,VecY,FUN="-")**2.+
        outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
if (argv$verbose) {
  print("+----------------------------+")
  print(paste("# observations =",L.y.tot))
}
#
#------------------------------------------------------------------------------
# Elaborations
# loop for DQC
yo.ok.pos<-which(ydqc.flag<=0 & !is.na(yo))
L.yo.ok<-length(yo.ok.pos)
while (L.yo.ok>0) {
# vector with the positions (pointers to VecS) of good observations 
  yo.ok.pos<-which(ydqc.flag<=0 & !is.na(yo))
  L.yo.ok<-length(yo.ok.pos)
  yo.ok.wet<-ydqc.flag<=0 & !is.na(yo) & yo>=argv$rrinf
  yo.ok.dry<-ydqc.flag<=0 & !is.na(yo) & yo<argv$rrinf
  yo.ok.pos.wet<-which(yo.ok.wet)
  yo.ok.pos.dry<-which(yo.ok.dry)
  L.yo.ok.wet<-length(yo.ok.pos.wet)
  L.yo.ok.dry<-length(yo.ok.pos.dry)
  if (argv$verbose) {
    print("+ -- (RE)START elaboration -------------------------------------------")
    print(paste("observation not NAs, presumably good =",
                L.yo.ok,"(wet=",L.yo.ok.wet,"dry=",L.yo.ok.dry,")"))
  }
# NO-RAIN OVER THE WHOLE DOMAIN
  if (L.yo.ok.wet==0) {
    if (argv$verbose) print("no rain over the whole domain")
    writeIfNoPrec()
    if (argv$verbose) print("Success exit")
    quit(status=0)
  } # end case of no-rain over the whole domain
#
#------------------------------------------------------------------------------
# DQC 200: remove isolated dry observations
  ytmp<-RR_dqc_isolatedDry(obs=data.frame(stid=VecS[yo.ok.pos],
                                          x=VecX[yo.ok.pos],
                                          y=VecY[yo.ok.pos],
                                          yo=yo[yo.ok.pos]),
                           rrinf=argv$rrinf,
                           pmax=argv$Lsubsample.max,
                           dmax=argv$Lsubsample.DHmax,
                           n.sector=n.sector,
                           dmin.dry=argv$DQC.min.dist.allowed) 
  if (any(ytmp!=0)) {
    ydqc.flag[yo.ok.pos[which(ytmp!=0)]]<-200
    next
  }
#
#------------------------------------------------------------------------------
# DQC 300: remove isolated WET observations
  ytmp<-RR_dqc_isolatedWet(obs=data.frame(stid=VecS[yo.ok.pos],
                                          x=VecX[yo.ok.pos],
                                          y=VecY[yo.ok.pos],
                                          yo=yo[yo.ok.pos]),
                           rrinf=argv$rrinf,
                           pmax=argv$Lsubsample.max,
                           dmax=argv$Lsubsample.DHmax,
                           n.sector=n.sector,
                           dmin.wet=argv$DQC.min.dist.allowed) 
  if (any(ytmp!=0)) {
    ydqc.flag[yo.ok.pos[which(ytmp!=0)]]<-300
    next
  }
#
#------------------------------------------------------------------------------
# assign grid points to precipitation yes/no
  Dh<-25
  t00<-Sys.time()
  D0<-Disth[yo.ok.dry,yo.ok.dry]
  D<-exp(-0.5*(D0/Dh)**2.)
  diag(D)<-diag(D)+argv$eps2
  InvD<-chol2inv(chol(D))
  rm(D,D0)
  xidi.dry<-OI_RR_fast(yo.sel=rep(1,L.yo.ok.dry),
                       yb.sel=rep(0,L.yo.ok.dry),
                       xb.sel=rep(0,Lgrid.CG),
                       xgrid.sel=xgrid.CG,
                       ygrid.sel=ygrid.CG,
                       zgrid.sel=zgrid.CG,
                       VecX.sel=VecX[yo.ok.dry],
                       VecY.sel=VecY[yo.ok.dry],
                       VecZ.sel=VecZ[yo.ok.dry],
                       Dh.cur=Dh,
                       Dz.cur=100000) 
  t11<-Sys.time()
  if (argv$verbose) 
    print(paste("xIDI dry, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  t00<-Sys.time()
  D0<-Disth[yo.ok.wet,yo.ok.wet]
  D<-exp(-0.5*(D0/Dh)**2.)
  diag(D)<-diag(D)+argv$eps2
  InvD<-chol2inv(chol(D))
  rm(D,D0)
  xidi.wet<-OI_RR_fast(yo.sel=rep(1,L.yo.ok.wet),
                       yb.sel=rep(0,L.yo.ok.wet),
                       xb.sel=rep(0,Lgrid.CG),
                       xgrid.sel=xgrid.CG,
                       ygrid.sel=ygrid.CG,
                       zgrid.sel=zgrid.CG,
                       VecX.sel=VecX[yo.ok.wet],
                       VecY.sel=VecY[yo.ok.wet],
                       VecZ.sel=VecZ[yo.ok.wet],
                       Dh.cur=Dh,
                       Dz.cur=100000) 
  t11<-Sys.time()
  if (argv$verbose) 
    print(paste("xIDI wet, time=",round(t11-t00,1),attr(t11-t00,"unit")))
# xp=1 for grid points where precipitation occurs
  xp<-xidi.wet
  xp[]<-0
#  if (any(xidi.wet>=(0.5*xidi.dry))) {
#    xp[which(xidi.wet>=(0.5*xidi.dry))]<-1
  if (any(xidi.wet>=(xidi.dry))) {
    xp[which(xidi.wet>=(xidi.dry))]<-1
  } else {
    if (argv$verbose) print("no rain over the whole domain")
    writeIfNoPrec()
    if (argv$verbose) print("Success exit")
    quit(status=0)
  }
#
#------------------------------------------------------------------------------
# DQC 400: remove dry(wet) observations close to a lot of wet(dry) ones
# observations might be correct but they will create unrealistic patterns, 
# either a hole or a isolated max, in the analysis field because the local 
# observation density is not sufficient to properly detect such a small scale
# process
  rp<-r.orog.CG
  rp[]<-NA
  rp[mask.CG]<-xp
  ytmp<-extract(rp,cbind(VecX[yo.ok.pos],VecY[yo.ok.pos]),na.rm=T)
  cond<-(ytmp==0 & yo[yo.ok.pos]>=argv$rrinf) |
        (ytmp==1 & yo[yo.ok.pos]< argv$rrinf)
  if (any(cond)) {
    ydqc.flag[yo.ok.pos[which(cond)]]<-400
    next
  }
  out<-remove_ltNobsClumps(r=rp,
                           n=2,
                           reverse=F,
                           obs=data.frame(x=VecX[yo.ok.pos],
                                          y=VecY[yo.ok.pos]))  
  if (any(out$obsflag==1)) {
    ydqc.flag[yo.ok.pos[which(out$obsflag==1)]]<-500
    next
  }
  rm(out)
  out<-remove_ltNobsClumps(r=rp,
                           n=2,
                           reverse=T,
                           obs=data.frame(x=VecX[yo.ok.pos],
                                          y=VecY[yo.ok.pos]))  
  if (any(out$obsflag==1)) {
    ydqc.flag[yo.ok.pos[which(out$obsflag==1)]]<-500
    next
  }
  rm(out)
#
# END of DQC 
  break
}
rm(r.orog.CG,mask.CG,zgrid.CG,xgrid.CG,ygrid.CG)
#
#==============================================================================
# ANALYSIS
ix<-which(ydqc.flag==0)
if (length(yo.ok.wet)==0) {
  if (argv$verbose) print("no rain over the whole domain")
  writeIfNoPrec()
  if (argv$verbose) print("Success exit")
  quit(status=0)
}

D0<-Disth[ix,ix]
rm(Disth)
#rm(x.FG,y.FG,mask.FG,Lgrid.FG)
nl<-length(argv$Dh.seq.reference)
Dhnl<-argv$Dh.seq.reference[nl]
for (l in 1:nl) {
  Dh<-argv$Dh.seq.reference[l]
  if (argv$verbose) 
    print(paste("scale # ",l," of ",length(argv$Dh.seq.reference),
                " (",round(Dh,0),"km)",sep=""))
  D<-exp(-0.5*(D0/Dh)**2.)
  diag(D)<-diag(D)+argv$eps2
  InvD<-chol2inv(chol(D))
  if (l==nl) {
    r<-r.orog.FG
    rm(r.orog.FG)
  } else {
    r<-aggregate(r.orog.FG,
                 fact=round(Dh/Dhnl,0),
                 expand=T,na.rm=T)
  }
  zvalues.l<-getValues(r)
  storage.mode(zvalues.l)<-"numeric"
  xy.l<-xyFromCell(r,1:ncell(r))
  x.l<-sort(unique(xy.l[,1]))
  y.l<-sort(unique(xy.l[,2]),decreasing=T)
  mask.l<-which(!is.na(zvalues.l))
  zgrid.l<-zvalues.l[mask.l]
  xgrid.l<-xy.l[mask.l,1]
  ygrid.l<-xy.l[mask.l,2]
  rm(xy.l,zvalues.l)
  if (l==1) {
    yb<-rep(mean(yo[ix]),length=length(ix))
    xb<-rep(mean(yo[ix]),length=length(xgrid.l))
  } else {
    rb<-resample(ra,r,method="bilinear")
    xb<-getValues(rb)[mask.l]
    if (any(is.na(xb))) {
      ib<-which(is.na(xb))
      aux<-extract(rb,cbind(xgrid.l[ib],ygrid.l[ib]),na.rm=T,buffer=(round(Dh/5,0)*1000))
      for (ll in 1:length(aux)) xb[ib[ll]]<-mean(aux[[ll]],na.rm=T)
      rb[mask.l]<-xb
      rm(aux,ib)
    }
    yb<-extract(rb,cbind(VecX[ix],VecY[ix]),method="bilinear")
    rm(rb)
    if (any(is.na(xb))) print("xb is NA")
    if (any(is.na(yb))) print("yb is NA")
  }
  xa.l<-OI_RR_fast(yo.sel=yo[ix],
                   yb.sel=yb,
                   xb.sel=xb,
                   xgrid.sel=xgrid.l,
                   ygrid.sel=ygrid.l,
                   zgrid.sel=zgrid.l,
                   VecX.sel=VecX[ix],
                   VecY.sel=VecY[ix],
                   VecZ.sel=VecZ[ix],
                   Dh.cur=Dh,
                   Dz.cur=10000) 
  ra<-r
  ra[]<-NA
  ra[mask.l]<-xa.l
}
rm(r,xb,yb,xa.l,mask.l)
ya<-extract(ra,cbind(VecX,VecY),method="bilinear")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (argv$verbose) print("++ Output")
# Station Points
cat("# variable: Precipitation\n",file=out.file.ver,append=F)
cat("# units: $mm$\n",file=out.file.ver,append=T)
cat("date     leadtime location  lat     lon      altitude obs      fcst\n",
    file=out.file.ver,append=T)
ix<-which(ydqc.flag<=0 & !is.na(yo))
cat(paste( yyyymmdd," ",
           rep(0,length(ix))," ",
           VecS[ix]," ",
           VecLat[ix]," ",
           VecLon[ix]," ",
           VecZ[ix]," ",
           round(yo[ix],1)," ",
           round(ya[ix],1),"\n",
           sep=""),
    file=out.file.ver,append=T)
print(paste("data saved on file",out.file.ver))
cat("yyyy;mm;dd;stid;x;y;z;yo;ya;dqc;\n",
    file=out.file.stn,append=F)
cat(paste(yyyy,mm,dd,
          formatC(VecS,format="f",digits=0),
          formatC(VecX,format="f",digits=0),
          formatC(VecY,format="f",digits=0),
          formatC(VecZ,format="f",digits=0),
          formatC(yo,format="f",digits=1),
          formatC(ya,format="f",digits=1),
          formatC(ydqc.flag,format="f",digits=0),
          "\n",sep=";"),
    file=out.file.stn,append=T)
print(paste("data saved on file",out.file.stn))
# grid
r.list<-list()
r.list[[1]]<-matrix(data=getValues(ra),
                    ncol=length(y.l),
                    nrow=length(x.l))
out<-nc4out(grid.list=r.list,
            times=yyyymmdd0600,
            file.name=out.file.grd,
            grid.type="ngcd",
            x=x.l,
            y=y.l,
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

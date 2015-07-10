# --~- Bspat_TEMP1d.R  -~--
# Bayesian spatial interpolation of daily mean air temperature (06-06 UTC)
# Spatial Consistency Test (SCT) is included.
# input: TAMRR from KDVH
#
# Outputs: look for "@@@@@@@@@" in the code and you'll get the ouput formats
#
# History:
# 02.12.2014 - Cristian Lussana. Original code.
# 25.02.2015 - Cristian Lussana. version 1.0.
#  change log: 
#  - allow for the use of observations outside Norway
#  - revisied queries to KDVH
#  - geographical information on seNorge2_dem_UTM33.nc
#  - definition of TEMP1d
#  - definition of a new directory tree
# ==~==========================================================================
rm(list=ls())
# Libraries
library(raster)
library(rgdal)
#library(colorspace)
library(ncdf)
#-------------------------------------------------------------------
# FUNCTIONS 
# manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}
#-------------------------------------------------------------------
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#
max.Km.stnINdomain<-300
#-----------------------------------------------------------------------------
# [] Setup OI parameters
sig2o<-1.3
eps2<-0.5
Dh<-60
Dz<-600
T2<-15
# Background - parameters
Dh.b<-70.
eps2.b<-0.5
Lsubsample<-50
Lsubsample.max<-50
Lsubsample.DHmax<-200
Lsubsample.vec<-vector()
#
print("ANALYSIS and DQC parameters")
print(paste("EPS2 Dh[Km] Dz[m] > ", round(eps2,3),
                                    round(Dh,3),
                                    round(Dz,3),
                                    sep=" "))
print(paste("VARobs[C^2] T^2 = ", round(sig2o,3), round(T2,3),sep=" "))
# MAIN ========================================================================
# Read command line arguments
arguments <- commandArgs()
arguments
date.string<-arguments[3]
config_file<-arguments[4]
config_par<-arguments[5]
if (length(arguments)!=5) 
  ext<-error_exit(paste("Error in command line arguments: \n",
  " R --vanilla yyyy.mm.dd configFILE configPAR \n",
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
filenamedem<-paste(main.path.geoinfo,"/fennodem_ETRS89_ETRS-LAEA.nc",sep="")
if (!file.exists(paste(main.path.geoinfo,"/fennodem_ETRS89_ETRS-LAEA.nc",sep=""))) 
  ext<-error_exit(paste("File not found:",main.path.geoinfo,"/fennodem_ETRS89_ETRS-LAEA.nc"))
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
source(paste(path2lib.com,"/Bspat_PseudoBackground.R",sep=""))
# test mode
print(testmode)
if (testmode) {
  print("TESTMODE TESTMODE TESTMODE")
  if (file.exists(paste(main.path,"/Bspat_TEMP1d/testbed",sep=""))) {
    testbed<-paste(main.path,"/Bspat_TEMP1d/testbed",sep="")
    station.info<-paste(testbed,"/station_data.csv",sep="")
    observed.data<-paste(testbed,"/observed_data.csv",sep="")
  } else {
    ext<-error_exit(paste("testbed not found"))
  }
}
# netcdf fixed parameters
grid.type <- "laea"
prod.date <- substr(Sys.time(),1,10)
xa.source.nc<-"daily mean air temperature from station data (06-06 UTC)"
xa.var.version <- "1.0"
xa.pname<-"TEMP1d"
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
xidi.source.nc<-"IDI for daily mean air temperature from station data (06-06 UTC)"
xidi.var.version <- "1.0"
xidi.pname<-"IDI"
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
# set Time-related variables
yyyy<-substr(date.string,1,4)
mm<-substr(date.string,6,7)
dd<-substr(date.string,9,10)
date.dot<-paste(dd,".",mm,".",yyyy,sep="")
yyyymm<-paste(yyyy,mm,sep="")
yyyymmdd<-paste(yyyymm,dd,sep="")
# output directories
dir.create(file.path(main.path.output,"NGCD"), showWarnings = FALSE)
dir.create(file.path(main.path.output,"NGCD_addInfo"), showWarnings = FALSE)
path2output.main<-paste(main.path.output,"/NGCD/TEMP1d",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/NGCD_addInfo/TEMP1d",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
if (!(file.exists(path2output.main)))     dir.create(path2output.main,showWarnings=F) 
if (!(file.exists(path2output.main.stn))) dir.create(path2output.main.stn,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) dir.create(path2output.main.grd,showWarnings=F) 
if (!(file.exists(path2output.add)))      dir.create(path2output.add,showWarnings=F) 
if (!(file.exists(path2output.add.grd)))  dir.create(path2output.add.grd,showWarnings=F) 
# Setup output files 
dir.create(paste(path2output.main.stn,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.main.grd,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm,sep=""),showWarnings=F)
#
out.file.stn<- paste(path2output.main.stn,"/",yyyymm,
                "/NGCD_TEMP1d_station_",yyyymmdd,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm,
                   "/NGCD_TEMP1d_grid_",yyyymmdd,".nc",sep="")
out.file.grd.bck<- paste(path2output.add.grd,"/",yyyymm,
                   "/NGCD_TEMP1d_grid_background_",yyyymmdd,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm,
                   "/NGCD_TEMP1d_grid_idi_",yyyymmdd,".nc",sep="")
#
print("Output files:")
print("analysis on the grid (netcdf)")
print(out.file.grd.ana)
print("background on the grid (netcdf)")
print(out.file.grd.bck)
print("idi on the grid (netcdf)")
print(out.file.grd.idi)
print("station outputs (text)")
print(out.file.stn)
#------------------------------------------------------------------------------
# Grid - it is defined by the DEM file
# CRS Coordinate Reference System
stackGeoGrid<-raster(filenamedem)
nx<-ncol(stackGeoGrid)
ny<-nrow(stackGeoGrid)
dx<-xres(stackGeoGrid)
dy<-yres(stackGeoGrid)
# 4 borders point
xmn<-xmin(stackGeoGrid)
xmx<-xmax(stackGeoGrid)
ymn<-ymin(stackGeoGrid)
ymx<-ymax(stackGeoGrid)
# South-West Point Coordinates
#Xnodesw<-xmn+dx/2
#Ynodesw<-ymn+dy/2
Xnodesw<-xmn
Ynodesw<-ymn
Xnode<-Xnodesw+(0:nx-1)*dx
Ynode<-Ynodesw+(0:ny-1)*dy
# Extract orography on unmasked gridpoints only
orog<-stackGeoGrid
# extract all the cell values: cells[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
cells<-extract(orog,1:ncell(orog))
xy<-xyFromCell(orog,1:ncell(orog))
x.G<-sort(unique(xy[,1]))
y.G<-sort(unique(xy[,2]),decreasing=T)
rc<-rowColFromCell(orog,1:ncell(orog))
aux<-as.vector(cells)
mask<-which(!is.na(aux))
zgrid<-aux[mask]
xgrid<-xy[mask,1]
ygrid<-xy[mask,2]
rowgrid<-rc[mask,1]
colgrid<-rc[mask,2]
Lgrid<-length(xgrid)
print(Lgrid)
rm(xy,rc,rowgrid,colgrid,cells,aux,stackGeoGrid)
#------------------------------------------------------------------------------
# [] Read Station Information 
# conditions:
# 1. stations in KDVH having: lat, lon and elevation. Note that UTM33 is 
#    obtained from lat,lon. Furthermore, the location must be in Norway or on
#    the border (less than max.Km.stnINdomain)
# 2. stations in CG
stations<-getStationMetadata_NGCD(from.year=yyyy,to.year=yyyy,
                                  max.Km=max.Km.stnINdomain,var="TEMP1d")
LOBS<-length(stations$stnr)
print(LOBS)
# define Vectors and Matrices
VecX<-vector(mode="numeric",length=LOBS)
VecY<-vector(mode="numeric",length=LOBS)
VecZ<-vector(mode="numeric",length=LOBS)
VecS<-vector(mode="numeric",length=LOBS)
Disth<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
Distz<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
yo<-vector(mode="numeric",length=LOBS)
yb<-vector(mode="numeric",length=LOBS)
ya<-vector(mode="numeric",length=LOBS)
yav<-vector(mode="numeric",length=LOBS)
yidi<-vector(mode="numeric",length=LOBS)
yidiv<-vector(mode="numeric",length=LOBS)
ydqc<-vector(mode="numeric",length=LOBS)
ydqc.flag<-vector(mode="numeric",length=LOBS)
D<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
S<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
VecX<-as.numeric(as.vector(stations$x))
VecY<-as.numeric(as.vector(stations$y))
VecZ<-as.numeric(as.vector(stations$z))
VecS<-as.numeric(as.vector(stations$stnr))
# compute S and D=S+R matrices (R is assumed to be diagonal = sigma_obs**2*I)
# Disth and Distz are the (symmetric) matrices where 
#      Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
#      Distz(i,j)=elevation difference between i-th station and j-th station [m]
Disth<-(outer(VecY,VecY,FUN="-")**2.+outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
Distz<-abs(outer(VecZ,VecZ,FUN="-"))
D<-exp(-0.5*(Disth/Dh)**2.-0.5*(Distz/Dz)**2.)
S<-D
D[row(D)==col(D)]<-D[row(D)==col(D)]+eps2
#
D.b<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
S.b<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
# NOTE: for the background Dh.b replaces Dh and Dz.b replaces Dz
# CL: is it eps2.b needed??
#D.b<-exp(-0.5*(Disth/Dh.b)**2.-0.5*(Distz/Dz.b)**2.)
D.b<-exp(-0.5*(Disth/Dh.b)**2.)
S.b<-D.b
D.b[row(D.b)==col(D.b)]<-D.b[row(D.b)==col(D.b)]+eps2.b
#
# DEBUG
#print(stations)
#write.table(file="stations.txt",stations)
#stations<-read.table(file="stations.txt")
#print(stations[1,])
#stations<-read.table("stations.txt",header=TRUE)
#print(stations[1,])
#  stnr z      x       y fennomean fenno_min4 fenno_long fenno_lat
# number of days
#print(ndays)
#------------------------------------------------------------------------------
# Elaborations
# define header for the station data output file
cat(paste("year","month","day","stid","x","y","z","yo",
          "yb","ya","yav","yidi","yidiv","dqcflag","\n",sep=";"),
    file=out.file.stn,append=F)
# r is the raster structure used for map production
r <-raster(ncol=nx, nrow=ny, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx,
            crs=proj4.ETRS_LAEA)
r[]<-NA
# get data 
data<-getStationData_NGCD(var="TEMP1d", yyyy=yyyy, mm=mm, dd=dd, statlist=stations)
yo<-as.numeric(data$value)
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
yb.set<-matrix(data=NA,ncol=LOBS,nrow=LOBS)
ybweights.set<-matrix(data=0,ncol=LOBS,nrow=LOBS)
yb.param<-matrix(data=NA,ncol=12,nrow=LOBS)
yb.param1.mat<-matrix(data=NA,ncol=12,nrow=LOBS)
yb.param2.mat<-matrix(data=NA,ncol=12,nrow=LOBS)
VecS.set<-matrix(data=NA,ncol=Lsubsample.max,nrow=LOBS)
VecS.set.pos<-matrix(data=NA,ncol=Lsubsample.max,nrow=LOBS)
yb.h.pos<-vector()
# cycle over the available stations
b.inc<-0
for (b in yo.h.pos) {
  if (b.inc>0) {
    finoa<-as.integer(Lsubsample.max*1/4)
    if (VecY[b]<7750000 & VecZ[b]<800) finoa<-as.integer(Lsubsample.max*2/4)
    if (VecY[b]<7100000 & VecZ[b]<800) finoa<-as.integer(Lsubsample.max*3/4)
    if (VecS[b]%in%VecS.set[yb.h.pos[1:b.inc],2:finoa]) next
  }
# select the closest (respect to b-th station) Lsubsample stations
  close2b.aux<-order(Disth[b,],decreasing=F)
  close2b.au1<-close2b.aux[which(close2b.aux%in%yo.h.pos)][1:Lsubsample]
  if (Disth[b,close2b.au1[2]]>(2*Dh.b)) next
  aux<-Lsubsample.max
  if (Disth[b,close2b.au1[Lsubsample.max]]>Lsubsample.DHmax) 
    aux<-max(which(Disth[b,close2b.au1]<=Lsubsample.DHmax))
  if (aux<10) next
  Lsubsample.vec[b]<-aux
  close2b<-close2b.au1[1:Lsubsample.vec[b]]
  b.inc<-b.inc+1
  yb.h.pos[b.inc]<-b
  yb.set0<-matrix(data=NA,nrow=LOBS)
  yb.set1<-matrix(data=NA,nrow=LOBS)
  yb.set2<-matrix(data=NA,nrow=LOBS)
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
# NOTE: this is D.b! not D...
  ide.b<-matrix(data=0,ncol=Lsubsample.vec[b],nrow=Lsubsample.vec[b])
  ide.b[row(ide.b)==col(ide.b)]=1
  InvD.b<-solve(D.b[close2b,close2b],ide.b)
# Background 0  
  yb.param0<-XYZnoinv_step0(VecX.b,VecY.b,VecZ.b,yo.b)
  # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
  # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
  yb.set0[]<-yb.param0[3]+ yb.param0[6]*(VecX-yb.param0[10])+
                           yb.param0[8]*(VecY-yb.param0[11])+
                           yb.param0[4]* VecZ
  J0<-sum((yb.set0[close2b]-yo.b)**2.)
# Background 1  
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
      yb.set2[bw]<-yb.param2[3] + yb.param2[4]*VecZ[bw] - 
                     yb.param2[5]/2.*(1+cos(pi*(VecZ[bw]-yb.param2[1])/yb.param2[2])) + 
                    ( (yb.param2[6]*(VecX[bw]-yb.param2[10])+yb.param2[8]*(VecY[bw]-yb.param2[11]))*(VecZ[bw]-yb.param2[1]) + 
                      (yb.param2[7]*(VecX[bw]-yb.param2[10])+yb.param2[9]*(VecY[bw]-yb.param2[11]))*(yb.param2[1]+yb.param2[2]-VecZ[bw]) ) / yb.param2[2]
    J2<-sum((yb.set2[close2b]-yo.b)**2.)
  } else {
    back2.flag<-F
    J2<-J0+1
  }
# DEBUG start      
#      png(file=paste("pro_",VecS[b],".png",sep=""),width=1200,height=1200)
#      plot(yo.b,VecZ.b,pch=19,col="black",cex=2.5,
#           xlim=c(min(c(yo.b,yb.set0[close2b],yb.set1[close2b],yb.set2[close2b]),na.rm=T),
#                  max(c(yo.b,yb.set0[close2b],yb.set1[close2b],yb.set2[close2b]),na.rm=T)))
#      points(yb.set0[close2b],VecZ[close2b],col="gray",cex=1.8)
#      points(yb.set1[close2b],VecZ[close2b],pch=19,col="blue",cex=1.8)
#      points(yb.set2[close2b],VecZ[close2b],pch=19,col="red",cex=1.8)
##      mtext(,side=3,cex=1.6)
#      dev.off()
# DEBUG stop      
# Best background
  aux<-order(c(J0,J1,J2))
  best<-aux[1]-1
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
  if (best==0) zero.str<-"+0."
  if (best==1) uno.str<-"+1."
  if (best==2) due.str<-"+2."
  print(paste("@@",b.inc,". id pos Zmn/x DisthMAX / J0 J1 J2 / #stn:",VecS[b],b,
                           round(min(VecZ.b),0),round(max(VecZ.b),0),
                           round(Disth[b,close2b[Lsubsample.vec[b]]],0),"/",
                           round(J0,0),round(J1,0),round(J2,0),"/",Lsubsample.vec[b]))
  print(paste(zero.str,"alpha beta gamma:",round(yb.param0[6],6),round(yb.param0[8],6),round(yb.param0[4],4)))
  print(paste(uno.str,"zinv dz tinv aA aB bA bB ga gb:",
        round(yb.param1[1],1), round(yb.param1[2],0), round(yb.param1[3],1),
        round(yb.param1[6],6), round(yb.param1[7],6), round(yb.param1[8],6),
        round(yb.param1[9],6), round(yb.param1[4],4), round(yb.param1[5],4)))
  print(paste(due.str,"h0 h1-h0 t0 a aA aB bA bB gamma:",
        round(yb.param2[1],1), round(yb.param2[2],0), round(yb.param2[3],1),
        round(yb.param2[5],6), round(yb.param2[6],6), round(yb.param2[7],6),
        round(yb.param2[8],6), round(yb.param2[9],6), round(yb.param2[4],4)))
# DEBUG stop
# G1 is LOBStOK x Lsubsample matrix to interpolate the Lsubsample values
# on the whole LOBStOK station dataset
  G1.b<-S.b[,close2b]
  K.b<-tcrossprod(G1.b,InvD.b)
  rm(G1.b)
  W.b<-tcrossprod(S.b[close2b,close2b],InvD.b)
  # this is idi (sort of)
  ybweights.set[,b]<-rowSums(K.b)
  if (any(is.na(ybweights.set[,b]))) {
    print("Error: NA value presents in ybweights.set vector. Unexpected.")
    print(paste(round(yb.set[,b],2),round(yo,2),round(ybweights.set[,b],3),"\n"))
    q(status=1)
  }
  rm(K.b,W.b)
} # end cycle LOBSt
LBAKh<-b.inc
print(paste("# stations used in background elaborations=",LBAKh))
# At this point I've this 5 outputs
# 1. yb.set<-matrix(data=NA,ncol=btimes,nrow=LOBSt)
# 2. ybweights.set<-matrix(data=NA,ncol=btimes,nrow=LOBSt)
# 3. yb.param<-matrix(data=NA,ncol=12,nrow=LOBSt)
# 4. VecS.set<-matrix(data=NA,ncol=Lsubsample,nrow=LOBSt)
# 5. VecS.set.pos[b,]<-close2b
#====================================================================
# Station Analysis cycle (with SCT!)
yoBad.id<-NA
yoBad.pos<-NA
isct<-0
ydqc.flag[]<-rep(0,LOBS)
ydqc.flag[-yo.h.pos]<-(-1)
yo.OKh.pos<-which(!is.na(yo) & ydqc.flag!=1)
LOBStOK<-length(yo.OKh.pos)
while (TRUE) {
  print(paste(">> Total number of observations [not NA & not ERR(so far)] =",
              LOBStOK))
# Station points - Background
  if (!is.na(yoBad.id)) {
    isct<-isct+1
    for (b in yb.h.pos) {
      if (any(VecS.set[b,1:Lsubsample.vec[b]]==yoBad.id)) {
        yb.set0<-matrix(data=NA,nrow=LOBS)
        yb.set1<-matrix(data=NA,nrow=LOBS)
        yb.set2<-matrix(data=NA,nrow=LOBS)
        Lsubsample.vec[b]<-Lsubsample.vec[b]-1
        VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
# select the closest (respect to indxb stations) Lsubsample stations
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
# NOTE: this is D.b! not D...
        ide.b<-matrix(data=0,ncol=Lsubsample.vec[b],nrow=Lsubsample.vec[b])
        ide.b[row(ide.b)==col(ide.b)]<-1
        InvD.b<-solve(D.b[close2b,close2b],ide.b)
# Background 0  
        yb.param0<-XYZnoinv_step0(VecX.b,VecY.b,VecZ.b,yo.b)
        # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
        # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
        yb.set0[]<-yb.param0[3]+ yb.param0[6]*(VecX-yb.param0[10])+
                                 yb.param0[8]*(VecY-yb.param0[11])+
                                 yb.param0[4]*VecZ
        J0<-sum((yb.set0[close2b]-yo.b)**2.)
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
                           yb.param2[5]/2.*(1+cos(pi*(VecZ[bw]-yb.param2[1])/yb.param2[2])) + 
                          ( (yb.param2[6]*(VecX[bw]-yb.param2[10])+yb.param2[8]*(VecY[bw]-yb.param2[11]))*(VecZ[bw]-yb.param2[1]) + 
                            (yb.param2[7]*(VecX[bw]-yb.param2[10])+yb.param2[9]*(VecY[bw]-yb.param2[11]))*(yb.param2[1]+yb.param2[2]-VecZ[bw]) ) / yb.param2[2]
          J2<-sum((yb.set2[close2b]-yo.b)**2.)
        } else {
          back2.flag<-F
          J2<-J0+1
        }
# DEBUG start
#            png(file=paste("pro_",VecS[b],".png",sep=""),width=1200,height=1200)
#            plot(yo.b,VecZ.b,pch=19,col="black",cex=2.5,xlim=c(min(c(yo.b,yb.set0[close2b],yb.set1[close2b],yb.set2[close2b])),
#                                                               max(yo.b,yb.set0[close2b],yb.set1[close2b],yb.set2[close2b])))
#            points(yb.set0[close2b],VecZ[close2b],col="gray",cex=1.8)
#            points(yb.set1[close2b],VecZ[close2b],pch=19,col="blue",cex=1.8)
#            points(yb.set2[close2b],VecZ[close2b],pch=19,col="red",cex=1.8)
#            dev.off()
# DEBUG stop
# Best background
        aux<-order(c(J0,J1,J2))
        best<-aux[1]-1
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
        print(paste("@@",isct,". id pos Zmn/x DisthMAX / J0 J1 J2 / #stn:",VecS[b],b,
                                 round(min(VecZ.b),0),round(max(VecZ.b),0),
                                 round(Disth[b,close2b[Lsubsample.vec[b]]],0),"/",
                                 round(J0,0),round(J1,0),round(J2,0),"/",Lsubsample.vec[b]))
        print(paste(zero.str,"alpha beta gamma:",round(yb.param0[6],6),round(yb.param0[8],6),round(yb.param0[4],4)))
        print(paste(uno.str,"zinv dz tinv aA aB bA bB ga gb:",
              round(yb.param1[1],1), round(yb.param1[2],0), round(yb.param1[3],1),
              round(yb.param1[6],6), round(yb.param1[7],6), round(yb.param1[8],6),
              round(yb.param1[9],6), round(yb.param1[4],4), round(yb.param1[5],4)))
        print(paste(due.str,"h0 h1-h0 t0 a aA aB bA bB gamma:",
              round(yb.param2[1],1), round(yb.param2[2],0), round(yb.param2[3],1),
              round(yb.param2[5],6), round(yb.param2[6],6), round(yb.param2[7],6),
              round(yb.param2[8],6), round(yb.param2[9],6), round(yb.param2[4],4)))
# DEBUG stop
# G1 is LOBStOK x Lsubsample matrix to interpolate the Lsubsample values
# on the whole LOBStOK station dataset
        G1.b<-S.b[,close2b]
        K.b<-tcrossprod(G1.b,InvD.b)
        rm(G1.b)
        W.b<-tcrossprod(S.b[close2b,close2b],InvD.b)
        # this is idi (sort of)
        ybweights.set[,b]<-rowSums(K.b)
        if (any(is.na(ybweights.set[,b]))) {
          print("Error: NA value presents in ybweights.set vector. Unexpected.")
          print(paste(round(yb.set[,b],2),round(yo,2),round(ybweights.set[,b],3),"\n"))
          q(status=1)
        }
        rm(K.b,W.b)
      }
    } # end cycle btimes (no optimization: over the station number)
  }
# normalization of the yb weights such that their sum equals to one
  ybweights.norm<-ybweights.set[,yb.h.pos] / 
                  rowSums(ybweights.set[,yb.h.pos])
# background on station points
  for (m in 1:LOBS) {
    yb[m]<-tcrossprod(ybweights.norm[m,],t(yb.set[m,yb.h.pos]))
  }
# deallocate memory
#      rm(D.b,S.b,yb.set,ybweights.set,ybweights.norm)
#      rm(VecY.b,VecX.b,VecZ.b,yo.b,close2b,y1.b)
# Station (CV)Analysis/(CV)IDI
  ide<-matrix(data=0,ncol=LOBStOK,nrow=LOBStOK)
  ide[row(ide)==col(ide)]=1
  InvD<-solve(D[yo.OKh.pos,yo.OKh.pos],ide)
  W<-tcrossprod(S[yo.OKh.pos,yo.OKh.pos],InvD)
  G1<-S[,yo.OKh.pos]
  K<-tcrossprod(G1,InvD)
  rm(G1)
#  ya<-yb + K %*% (yo[yo.OKh.pos]-yb[yo.OKh.pos])
  ya<-yb + tcrossprod(K,t(yo[yo.OKh.pos]-yb[yo.OKh.pos]))
  yav<-ya
  yav[yo.OKh.pos]<-yo[yo.OKh.pos] + 1./(1.-diag(W)) * (ya[yo.OKh.pos]-yo[yo.OKh.pos])
  yidi<-rowSums(K)
  yidiv<-yidi
  yidiv[yo.OKh.pos]<-rep(1,LOBStOK) + 1./(1.-diag(W)) * (yidi[yo.OKh.pos]-rep(1,LOBStOK))
# DQC CHECK - Spatial Continuity Check
  ydqc[]<--9999.
  ydqc[yo.OKh.pos]<-(yo[yo.OKh.pos]-yav[yo.OKh.pos]) *
                               (yo[yo.OKh.pos]-ya[yo.OKh.pos]) / sig2o
# DQC test: T2 is the SCT threshold, if any(ydqc>T2) is true then reject the 
# station with highest ydqc value. This mean that the station is not used to
# compute the analysis
  if (max(ydqc)<=T2) {
    break
  }
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
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Station Points - Write output on file 
cat(paste(yyyy,mm,dd,
          round(VecS,0),round(VecX,0),round(VecY,0),round(VecZ,0),
          round(yo,1),round(yb,2),round(ya,2),round(yav,2),
          round(yidi,3),round(yidiv,3),round(ydqc.flag,2),
          "\n",sep=";"),
          file=out.file.stn,append=T)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# if something strange take place in the Analysis/SCT cycle the pass to the next iteration
if (LOBStOK<=0) next
# DQC flag "0" means good observations
print(paste(">>>>Total number of observations [not NA & good] =",LOBStOK))
#------------------------------------------------------------------------------
# Grid - Background
xb<-vector(mode="numeric",length=Lgrid)
xb.set<-matrix(data=0,ncol=LBAKh,nrow=Lgrid)
xbweights.set<-matrix(data=NA,ncol=LBAKh,nrow=Lgrid)
b.aux<-0
print("++ Grid - Background elaborations\n")
print("#/tot stnid #grid.points/#grid.points.tot")
for (b in yb.h.pos) {
  b.aux<-b.aux+1
  xindx<-which( (xgrid-min(VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]]))>(-3*Dh.b*1000) & (xgrid-max(VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]]))<(3*Dh.b*1000) &
                (ygrid-min(VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]]))>(-3*Dh.b*1000) & (ygrid-max(VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]]))<(3*Dh.b*1000) )
  Lgrid.b<-length(xindx)
  if (Lgrid.b==0) {
    xb.set[,b.aux]<-0
    xbweights.set[,b.aux]<-0.
    next
  }
  ide.b<-matrix(data=0,ncol=Lsubsample.vec[b],nrow=Lsubsample.vec[b])
  ide.b[row(ide.b)==col(ide.b)]<-1
  print("Inv")
  InvD.b<-solve(D.b[VecS.set.pos[b,1:Lsubsample.vec[b]],VecS.set.pos[b,1:Lsubsample.vec[b]]],ide.b)
#
  xbweights.set[,b.aux]<-0
  ndim<-10000
  i<-0
  while ((i*ndim)<Lgrid.b) {
    start<-i*ndim+1
    end<-(i+1)*ndim
    if (end>Lgrid.b) {
      end<-Lgrid.b
    }
    ndimaux<-end-start+1
    print(paste(round(i),round(start,0),round(end,0),round(Lgrid,0)))
    aux<-matrix(ncol=Lsubsample.vec[b],nrow=ndimaux,data=0.)
    G.b<-matrix(ncol=Lsubsample.vec[b],nrow=ndimaux,data=0.)
    aux<-(outer(ygrid[xindx[start:end]],VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2. +
          outer(xgrid[xindx[start:end]],VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.)**0.5/1000.
    G.b<-exp(-0.5*(aux/Dh)**2.)
    rm(aux)
    K.b<-tcrossprod(G.b,InvD.b)
    rm(G.b)
    xbweights.set[xindx[start:end],b.aux]<-rowSums(K.b)
    rm(K.b)
    i<-i+1
  }
## G matrix
#  print("G")
#  Disth.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
#  G.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
#  Disth.b<-(outer(ygrid[xindx],VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.+
#            outer(xgrid[xindx],VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.)**0.5/1000.
#  G.b<-exp(-0.5*(Disth.b/Dh.b)**2.)
#  rm(Disth.b)
##  compute analysis/idi over grid 
#  print("tcross")
#  K.b<-tcrossprod(G.b,InvD.b)
#  print("tcross end")
#  xbweights.set[,b.aux]<-0
#  xbweights.set[xindx,b.aux]<-rowSums(K.b)
  xindx1<-which(xbweights.set[,b.aux]>0.05)
  Lgrid.b1<-length(xindx1)
  print(paste(b.aux,"/",LBAKh," ",VecS[b]," ",Lgrid.b,"-->",Lgrid.b1,"/",Lgrid,"\n",sep=""))
  rm(G.b,K.b)
  if (yb.param[b,12]==0) {
    # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
    # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
    xb.set[xindx1,b.aux]<-yb.param[b,3]+ yb.param[b,6]*(xgrid[xindx1]-yb.param[b,10])+
                                        yb.param[b,8]*(ygrid[xindx1]-yb.param[b,11])+
                                        yb.param[b,4]* zgrid[xindx1]
  }
  if (yb.param[b,12]==1) {
    # 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
    # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
    zinv<-yb.param[b,1]
    zabov<-zinv+yb.param[b,2]
    zbelo<-zinv-yb.param[b,2]
    bfabov<-yb.param[b,3]+ yb.param[b,6]*(xgrid[xindx1]-yb.param[b,10])+
                           yb.param[b,8]*(ygrid[xindx1]-yb.param[b,11])+
                           yb.param[b,4]*(zgrid[xindx1]-zinv)
    bfbelo<-yb.param[b,3]+ yb.param[b,7]*(xgrid[xindx1]-yb.param[b,10])+
                           yb.param[b,9]*(ygrid[xindx1]-yb.param[b,11])+ 
                           yb.param[b,5]*(zgrid[xindx1]-zinv)
    aux.ab<-which(zgrid[xindx1]>zabov)
    aux.bl<-which(zgrid[xindx1]<=zbelo)
    aux.bw<-which((zgrid[xindx1]>zbelo)&(zgrid[xindx1]<=zabov))
    xb.set[xindx1[aux.ab],b.aux]<-bfabov[aux.ab]
    xb.set[xindx1[aux.bl],b.aux]<-bfbelo[aux.bl]
    xb.set[xindx1[aux.bw],b.aux]<-(bfabov[aux.bw]*(zgrid[xindx1[aux.bw]]-zbelo) + bfbelo[aux.bw]*(zabov-zgrid[xindx1[aux.bw]]) ) / (zabov-zbelo)
    rm(aux.ab,aux.bl,aux.bw)
  }
  if (yb.param[b,12]==2) {
# 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
    h0<-yb.param[b,1]
    h1<-yb.param[b,1]+yb.param[b,2]
    aux.ab<-which(zgrid[xindx1]>=h1)
    aux.bl<-which(zgrid[xindx1]<=h0)
    aux.bw<-which((zgrid[xindx1]>h0) & (zgrid[xindx1]<h1))
    xb.set[xindx1[aux.ab],b.aux]<-yb.param[b,3] + 
                             yb.param[b,4]* zgrid[xindx1[aux.ab]] + 
                             yb.param[b,6]*(xgrid[xindx1[aux.ab]]-yb.param[b,10]) + 
                             yb.param[b,8]*(ygrid[xindx1[aux.ab]]-yb.param[b,11])
    xb.set[xindx1[aux.bl],b.aux]<-yb.param[b,3] + 
                             yb.param[b,4]*zgrid[xindx1[aux.bl]] - 
                             yb.param[b,5] +
                             yb.param[b,7]*(xgrid[xindx1[aux.bl]]-yb.param[b,10]) + 
                             yb.param[b,9]*(ygrid[xindx1[aux.bl]]-yb.param[b,11])
    xb.set[xindx1[aux.bw],b.aux]<-yb.param[b,3] +
                             yb.param[b,4]*zgrid[xindx1[aux.bw]] - 
                             yb.param[b,5]/2.*(1+cos(pi*(zgrid[xindx1[aux.bw]]-h0)/yb.param[b,2])) + 
     ( (yb.param[b,6]*(xgrid[xindx1[aux.bw]]-yb.param[b,10])+yb.param[b,8]*(ygrid[xindx1[aux.bw]]-yb.param[b,11]))*(zgrid[xindx1[aux.bw]]-h0) + 
       (yb.param[b,7]*(xgrid[xindx1[aux.bw]]-yb.param[b,10])+yb.param[b,9]*(ygrid[xindx1[aux.bw]]-yb.param[b,11]))*(h1-zgrid[xindx1[aux.bw]]) ) / yb.param[b,2]
    rm(aux.ab,aux.bl,aux.bw)
  }
}
xbweights.norm<-xbweights.set/rowSums(xbweights.set)
# DEBUG start
#    for (b in 1:length(yb.h.pos)) {
#      b.pos<-yb.h.pos[b]
#      xindx<-which( (xgrid-min(VecX[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]]))>(-3*Dh.b*1000) & 
#                    (xgrid-max(VecX[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]]))<(3*Dh.b*1000) &
#                    (ygrid-min(VecY[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]]))>(-3*Dh.b*1000) & 
#                    (ygrid-max(VecY[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]]))<(3*Dh.b*1000) )
#      r[mask]<-NA
#      r[mask[xindx]]<-xbweights.norm[xindx,b]
#      png(file=paste("xb_norm_",VecS[b.pos],".png",sep=""),width=1400,height=1400)
#      image(r,breaks=seq(0,1,by=0.01),col=rainbow(100))
#      points(VecX[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]],VecY[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]],pch=19,col="black")
#      points(VecX[b.pos],VecY[b.pos],pch=print(paste(round(i),round(start,0),round(end,0),round(Lgrid,0)))19,col="black",cex=2)
#      legend(x="bottomright",legend=seq(0,1,by=0.01),fill=rainbow(100),cex=0.9)
#      dev.off()
#      r[mask[xindx]]<-xbweights.set[xindx,b]
#      png(file=paste("xb_set_",VecS[b.pos],".png",sep=""),width=1400,height=1400)
#      image(r,breaks=seq(0,1,by=0.01),col=rainbow(100))
#      points(VecX[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]],VecY[VecS.set.pos[b.pos,1:Lsubsample.vec[b.pos]]],pch=19,col="black")
#      points(VecX[b.pos],VecY[b.pos],pch=19,col="black",cex=2)
#      legend(x="bottomright",legend=seq(0,1,by=0.01),fill=rainbow(100),cex=0.9)
#      dev.off()
#    }
# DEBUG end
#    print(paste(round(rowSums(ybweights.set),9),"\n"))
for (i in 1:Lgrid) {
  xb[i]<-tcrossprod(xbweights.norm[i,],t(xb.set[i,]))
}
rm(xb.set,xbweights.norm,xbweights.set,xindx,xindx1)
# Gridded Analysis/IDI    
xa<-vector(mode="numeric",length=Lgrid)
xidi<-vector(mode="numeric",length=Lgrid)
ndim<-10000
i<-0
print("++ Grid - Analysis and IDI elaborations\n")
print("# gridpoint.start gridpoint.end gridpoint.total\n")
#to.write <- file("aux.bin", "wb")
#to.read<-file("aux.bin", "rb")
while ((i*ndim)<Lgrid) {
  start<-i*ndim+1
  end<-(i+1)*ndim
  if (end>Lgrid) {
    end<-Lgrid
  }
  ndimaux<-end-start+1
  print(paste(round(i),round(start,0),round(end,0),round(Lgrid,0)))
  aux<-matrix(ncol=LOBStOK,nrow=ndimaux,data=0.)
  auxz<-matrix(ncol=LOBStOK,nrow=ndimaux,data=0.)
  G<-matrix(ncol=LOBStOK,nrow=ndimaux,data=0.)
  aux<-(outer(ygrid[start:end],VecY[yo.OKh.pos],FUN="-")**2. +
        outer(xgrid[start:end],VecX[yo.OKh.pos],FUN="-")**2.)**0.5/1000.
#  for (j in 1:Linfo) writeBin(as.numeric(aux[,j]),to.write,size=4)
#  aux<-matrix(readBin(to.read,numeric(),size=4,ndimaux*Linfo),ndimaux,Linfo)
# vertical distance
  auxz<-abs(outer(zgrid[start:end],VecZ[yo.OKh.pos],FUN="-"))
  G<-exp(-0.5*(aux/Dh)**2.-0.5*(auxz/Dz)**2.)
  rm(aux,auxz)
  K<-tcrossprod(G,InvD)
  rm(G)
  xa[start:end]<-xb[start:end]+tcrossprod(K,t(yo[yo.OKh.pos]-yb[yo.OKh.pos]))
  xidi[start:end]<-rowSums(K)
  rm(K)
  i<-i+1
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# output - ANALYSIS
r[]<-NA
r[mask]<-round(xa,1)
if (xa.flag.write) {
  nogrid.ncout(file.name=out.file.grd.ana,
               grid=t(as.matrix(r)),
               x=x.G,y=y.G,grid.type=grid.type,
               times=c(paste(yyyymmdd,"0000",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.ETRS_LAEA,
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
# output - background
r[]<-NA
r[mask]<-round(xb,1)
if (xa.flag.write) {
  nogrid.ncout(file.name=out.file.grd.bck,
               grid=t(as.matrix(r)),
               x=x.G,y=y.G,grid.type=grid.type,
               times=c(paste(yyyymmdd,"0000",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.ETRS_LAEA,
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
# output - IDI
r[]<-NA
r[mask]<-round(100*xidi,1)
if (xidi.flag.write) {
  nogrid.ncout(file.name=out.file.grd.idi,
               grid=t(as.matrix(r)),
               x=x.G,y=y.G,grid.type=grid.type,
               times=c(paste(yyyymmdd,"0000",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.ETRS_LAEA,
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
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Exit - Success
q(status=0)


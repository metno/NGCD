require(gstat)
require(foreign)

require(sp)

#source("plotgrd.R")

days<-c(31,28,31,30,31,30,31,31,30,31,30,31)
nr.stations<-0
nr.stations.date<-0
kl <- 1
fyear<-2009
lyear<-2009
for (year in lyear:fyear){
 Jday <- 0
 iar<-year
 lmnd<-12
 fmnd<-ifelse(year==lyear,12,1)
 
 for (imnd in fmnd:fmnd){
# Check for leap year
  ifelse(((iar/4)-as.integer(iar/4))==0,days[2]<-29,days[2]<-28)
  for (idag in 1:days[imnd]){
   ifelse(imnd<10,mnd<-paste("0",imnd,sep=""),mnd<-imnd)
   ifelse(idag<10,day<-paste("0",idag,sep=""),day<-idag)
   date<-paste(day,".",mnd,".",year,sep="")
   Jday<-Jday+1


# Station coordinates
   normal.koord<-read.dbf("RR_stations_f.dbf")
   normal.koord$X <-normal.koord$X
   normal.koord$Y <-normal.koord$Y
   normal.koord$Stnr<-normal.koord$SOUID
   stations<-normal.koord

# Select observations from daily data files in ../NGCD_dataset/data/TG

   date <- paste(day,".",mnd,".",year,sep="")
   print(date)

# Define CASE datafile
   case.path<-c("/lustre/storeB/project/metkl/senorge2/case/case_20180112/RR_date_dqc/")
   case.dir<-paste(case.path,year,mnd,"/",sep="")
   case.file<-paste(case.dir,"case_RR_date_",year,mnd,day,".txt",sep="")
#   infile<-paste("/lustre/NGCD_dataset/data/RR/",year,"/",mnd,"/RRNGCD_",year,mnd,day,".txt",sep="")
   o.nor<-read.table(case.file,header=TRUE,sep=";")
   o.nor$id<-ifelse(is.na(o.nor$metnostaid)==FALSE,o.nor$metnostaid,o.nor$souid)
   n <- dim(o.nor)[1]

   o.i <- data.frame(matrix(nrow=n,ncol=9))
   names(o.i)<-c("Stnr","TAM","X","Y","Z","Fennomean","Fennomin","Lat","Long")
   o.i$Stnr<-o.nor$id
   o.i$RR<-o.nor$value
   j<-0
   for (i in 1:n){
#   print(paste("X1",i))
# Find coordinates
    if (length(stations$X[stations$Stnr==o.i$Stnr[i]]) == 1){
     if (is.na(o.i$RR)[i]==FALSE){
      j<-j+1
      o.i$RR[j]<-o.i$RR[i] 
      o.i$X[j]<-stations$X[stations$Stnr==o.i$Stnr[i]]
      o.i$Y[j]<-stations$Y[stations$Stnr==o.i$Stnr[i]]
      o.i$Z[j]<-stations$HGHT[stations$Stnr==o.i$Stnr[i]]
      o.i$Stnr2[j]<-o.i$Stnr[i]
     }
    }
   }

   o.i.bkp <-o.i
   o.i$Stnr<-o.i$Stnr2
   o.i<-o.i[1:j,]
#     Check for duplicate points

   if (anyDuplicated(o.i$X)){o.i$Stnr[duplicated(o.i$X)]}->ax
  

   while (anyDuplicated(o.i$X)>0){
#     print(anyDuplicated(o.i$X))
     kk<-anyDuplicated(o.i$X)
     m<-dim(o.i)[1] - 1
       o.i.2<-rbind(o.i[1:(kk-1),],o.i[(kk+1):dim(o.i)[1],]) 
     o.i<-o.i.2 
   }
   nst<-dim(o.i)[1]
   print(paste("Number of stations:",nst))
   nr.stations[kl+1]<-nst
   nr.stations.date[kl+1]<-date
   o.i$it<-2
   

   o.i.out<-cbind(o.i$X,o.i$Y,o.i$RR, o.i$it, o.i$Stnr, o.i$Z)

# XVAL-loop
  o.i.res<-data.frame(matrix(ncol=6,nrow=nst))

  for (i in 1:nst){
#    print(paste("Validating station",i,"of",nst)) 
    if (i==1){o.i.xval<-o.i.out[(i+1):nst,]}
    if (i>1){
     if (i<nst){o.i.xval<-rbind(o.i.out[1:(i-1),],o.i.out[(i+1):nst,])}
     if (i==nst){o.i.xval<-o.i.out[1:(nst-1),]}
    }
    o.i.x<-rbind(o.i.out[i,])  

# drop     

   write.table(o.i.xval,file="rr_x.par", row.names=FALSE,col.names=FALSE)
   write.table(o.i.x,file="x_station.par", row.names=FALSE,col.names=FALSE)
 
#   print("Triangulation")
   system("./test_temp_stationsX")
#   print("Gridding")
   system("./point_in_triangleX")

   
    o.i.res[i,]<-read.table("rr1km_out_X.txt")
   }   
   
   names(o.i.res)<-c("Stnr","RRobs","RRest","Z","Zest","RRestZ")
   plot(o.i.res$RRobs,o.i.res$RRestZ,main=date,pch=19,cex=.7)
   lines(c(0,200),c(0,200),lty=3)
   write.table(o.i.res,file=paste("RR.xv/rr_",date,".txt",sep=""),row.names=FALSE)
   kl <- kl+1
  }
 }
}



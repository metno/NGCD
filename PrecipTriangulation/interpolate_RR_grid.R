require(gstat)
require(foreign)

require(sp)

source("plotgrd.R")

days<-c(31,28,31,30,31,30,31,31,30,31,30,31)
nr.stations<-0
nr.stations.date<-0
kl <- 1
fyear<-1981
lyear<-2003
for (year in lyear:fyear){
 Jday <- 0
 iar<-year
 lmnd<-1
 fmnd<-ifelse(year==lyear,1,1)
 
 for (imnd in fmnd:lmnd){
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

# Select observations from daily data files (CASE)

   date <- paste(day,".",mnd,".",year,sep="")
   print(date)

# Define datafile
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
   write.table(o.i.out,file="rr.par", row.names=FALSE,col.names=FALSE)


# Triangulation
   print("Triangulation")
   system("./test_temp_stations2")
   print("Gridding")
   system("./point_in_triangle2")

# Plot result
   
   plot.rrgrd.asc("rr1km_out_2.asc",date)
   points(o.i$X,o.i$Y,pch=19,cex=.2)

# Output
   if(!file.exists(paste("RR.grd/",iar,sep=""))){ 
       dir.create(paste("RR.grd/",iar,sep=""))
       print(paste("Directory RR.grd/",iar," created",sep=""))
   }
   system(paste("cp rr1km_out_2.asc RR.grd/",iar,"/rr_",date,".asc",sep=""))

   if(!file.exists(paste("RR.dbf/",iar,sep=""))){ 
       dir.create(paste("RR.dbf/",iar,sep=""))
       print(paste("Directory RR.dbf/",iar," created",sep=""))
   }
   o.i.out<-data.frame(o.i.out)
   names(o.i.out)<-c("X","Y","RR","it","Stnr","Z")
   write.dbf(o.i.out,file=paste("RR.dbf/",iar,"/rr_",date,".dbf",sep=""))
   kl <- kl+1
  }
 }
}



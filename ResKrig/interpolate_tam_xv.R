require(gstat)
require(foreign)
require(sp)
days<-c(31,28,31,30,31,17,31,31,30,31,30,31)
nr.stations<-0
nr.stations.date<-0
kl <- 1
fyear<-1971
lyear<-2016

iplot<-FALSE

# Read external predicor grids
dem<-read.asciigrid("dem1km.asc")
fennomean<-read.asciigrid("dem_mean20.asc")
fennomin<-read.asciigrid("dem_min20.asc")
lat<-read.asciigrid("ngcd_lat.asc")
long<-read.asciigrid("ngcd_long.asc")
fennomask<-read.asciigrid("fennomask.asc")



for (year in fyear:lyear){
 Jday <- 0
 iar<-year
 lmnd<-12
 fmnd<-ifelse(year==fyear,1,1)
 
 for (imnd in fmnd:lmnd){
# Check for leap year
  ifelse(((iar/4)-as.integer(iar/4))==0,days[2]<-29,days[2]<-28)
  for (idag in 1:days[imnd]){
   ifelse(imnd<10,mnd<-paste("0",imnd,sep=""),mnd<-imnd)
   ifelse(idag<10,day<-paste("0",idag,sep=""),day<-idag)
   date<-paste(day,".",mnd,".",year,sep="")
   Jday<-Jday+1
 
# Station coordinates and metadata
   normal.koord<-read.dbf("ngcd_laea_c3s_meta_tg.dbf")
   normal.koord$X <-normal.koord$X/1000.
   normal.koord$Y <-normal.koord$Y/1000.
   normal.koord$Stnr<-normal.koord$Stnr
   stations<-normal.koord

# Select observations from daily data files in /ProdData/TG_<date>.txt

   date <- paste(day,".",mnd,".",year,sep="")
   print(date)

# Define TAX datafile
# Define datafile
    case.path<-c("/lustre/storeB/project/metkl/senorge2/case/case_20180112/TG_date_dqc/")
   case.dir<-paste(case.path,year,mnd,"/",sep="")
   case.file<-paste(case.dir,"case_TG_date_",year,mnd,day,".txt",sep="")
   o.nor<-read.table(case.file,header=TRUE,sep=";")
   o.nor$id<-ifelse(is.na(o.nor$metnostaid)==FALSE,o.nor$metnostaid,o.nor$souid)

   n <- dim(o.nor)[1]

   o.i <- data.frame(matrix(nrow=n,ncol=9))
   names(o.i)<-c("Stnr","TAM","X","Y","Z","Fennomean","Fennomin","Lat","Long")
   
   o.i$Stnr<-o.nor$id
   o.i$TAM<-o.nor$value
   j<-0
   for (i in 1:n){
#   print(paste("X1",i))
# Find coordinates
    if (length(stations$X[stations$Stnr==o.i$Stnr[i]]) == 1){
     if (o.i$Stnr[i]!=76900){      
      if (is.na(o.i$TAM)[i]==FALSE){
       j<-j+1
       o.i$TAM[j]<-o.i$TAM[i] 
       o.i$X[j]<-stations$X[stations$Stnr==o.i$Stnr[i]]
       o.i$Y[j]<-stations$Y[stations$Stnr==o.i$Stnr[i]]
       o.i$Z[j]<-stations$Z[stations$Stnr==o.i$Stnr[i]]
       o.i$Fennomean[j]<-stations$dem_mean20[stations$Stnr==o.i$Stnr[i]]
       o.i$Fennomin[j]<-stations$dem_min20[stations$Stnr==o.i$Stnr[i]]
       o.i$Lat[j]<-stations$Lat[stations$Stnr==o.i$Stnr[i]]
       o.i$Long[j]<-stations$Long[stations$Stnr==o.i$Stnr[i]]
       o.i$Stnr2[j]<-o.i$Stnr[i]
      }
     } 
    }
   }

   o.i.bkp <-o.i
   o.i$Stnr<-o.i$Stnr2
   o.i<-o.i[1:j,]
#     Check for duplicate points

   if (anyDuplicated(o.i$X)){o.i$Stnr[duplicated(o.i$X)]}->ax
  

   while (anyDuplicated(o.i$X)>0){
     print(anyDuplicated(o.i$X))
     kk<-anyDuplicated(o.i$X)
     m<-dim(o.i)[1] - 1
       o.i.2<-rbind(o.i[1:(kk-1),],o.i[(kk+1):dim(o.i)[1],]) 
     o.i<-o.i.2 
   }
   nst<-dim(o.i)[1]
   print(paste("Number of stations:",nst))
   nr.stations[kl+1]<-nst
   nr.stations.date[kl+1]<-date

# Detrending

# Define parameters and constants.
   k <- c(19.4879,14.7045,19.8752,26.0346,35.1963,39.4180,35.9240,36.8990,34.4318, 29.1849,25.9764,18.0301)
   v1 <-c(-.0012,-.0019, -.0046,-.0061,-.0063,-.0063,-.0061, -.0057,-.0055,-.0046, -.0032, -.0016)
   v2 <-c(-.0051,-.0043, -.0021, -.0006, 0.0007,0.0011, 0.0009,0.0000,-.0011,-.0016, -.0031, -.0041)
   v3 <-c(-.0083,-.0062,-.0033,-.0008, 0.0000,0.0021, 0.0016, 0.0005,0.0000,-.0017, -.0054, -.0089)
   v4 <-c(-.2694,-.1918,-.2579,-.3288, -.4166,-.4460,-.3700, -.3763,-.3764, -.3349,-.3491,-.2459)
   v5 <-c(-.4395,-.4282, -.3044, -.1505,-.0363,0.091, 0.1290,0.0370,-.0543,-.1581, -.2194,-.3470)

   im <- as.integer(mnd)

# Calculate trend

   o.i$TM0<- o.i$TAM - k[im] - (v1[im]*o.i$Z) - (v2[im]*o.i$Fennomean) - (v3[im]*o.i$Fennomin) - (v4[im]*o.i$Lat) - (v5[im]*o.i$Long) 
   o.i$TM0<- o.i$TAM - k[im] - (v1[im]*o.i$Z) - (v2[im]*o.i$Fennomean) - (v3[im]*o.i$Fennomin) - (v4[im]*o.i$Lat) - (v5[im]*o.i$Long) 


# Variogram models. 
   sill<-c(7,5,1.3,.33,.7,1,.6,.19,.3,1.3,3.5,7)
   nugget<-c(0,0,0,0,0,0,0,0,0,0,0,0)
   modell<-c("Exp","Exp","Exp","Exp","Exp","Exp","Sph","Exp","Exp","Exp","Exp","Exp")
   range<-c(250,250,200,100,75,150,500,100,75,175,200,250)
   gam <-vgm(sill[im],modell[im],range[im],nugget=0)
   gam2 <-vgm(sill[im],modell[im],range[im])

   o.i$X<-o.i$X
   o.i$Y<-o.i$Y
   tam0 <- o.i$TM0

#   inter.o.i<-krige(TM0~1,~X+Y,o.i,newdata=normal.koord,model=gam)
   inter.o.i.rr<-krige(TM0~1,~X+Y,o.i,newdata=o.i,model=gam,nmax=6)

   o.i$pred<- krige.cv(TM0~1,~X+Y,o.i,model=gam,nmax=6)$var1.pred


# Put trend back on.

   output <- o.i
   output$TAMest.0<- round(inter.o.i.rr$var1.pred,digits=3) 
   output$TAMest.Z<- round(inter.o.i.rr$var1.pred + k[im] + (v1[im]*output$Z) + (v2[im]*output$Fennomean) + (v3[im]*output$Fennomin) + (v4[im]*output$Lat) + (v5[im]*output$Long),digits=2) 
   o.i$predZ<- round(o.i$pred + k[im] + (v1[im]*o.i$Z) + (v2[im]*o.i$Fennomean) + (v3[im]*o.i$Fennomin) + (v4[im]*o.i$Lat) + (v5[im]*o.i$Long),digits=2) 
   output$predZ.xv<-o.i$predZ
   o.i$Date<-date
   o.i$Daynr <- kl+1
   o.i$Jday <- Jday 


   output$Day<-  as.integer(day)
   output$Month<-  as.integer(mnd)
   output$Year<-  as.integer(year)


   write.table(output,file=paste("TG.xv/tx.xv_",date,".txt",sep=""),row.names=FALSE,col.names=TRUE)
   if (kl==0){write.table(o.i,file="tg_output_xv_db2.txt",row.names=FALSE,col.names=TRUE,append=FALSE)}

   if (kl>0){write.table(o.i,file="tg_output_xv_db2.txt",row.names=FALSE,col.names=FALSE,append=TRUE)}

   kl <- kl+1
  }
 }
}



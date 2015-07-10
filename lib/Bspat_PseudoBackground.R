# funzione per costruire il background: INIZIO
AA<-matrix(nrow=5,ncol=5,data=0.)
vv<-vector(length=5,mode="numeric")
aux_vv<-vector(length=5,mode="numeric")
id5<-matrix(data=0.,nrow=5,ncol=5)
id5[row(id5)==col(id5)]=1.

AA7<-matrix(nrow=7,ncol=7,data=0.)
vv7<-vector(length=7,mode="numeric")
aux_vv7<-vector(length=7,mode="numeric")
id7<-matrix(data=0.,nrow=7,ncol=7)
id7[row(id7)==col(id7)]=1.

XZinv<-function(b_x,b_z,b_yo,dz) {
  aux_length<-length(b_yo)
  param<-vector(length=8,mode="numeric")
  resvec<-vector(length=aux_length,mode="numeric")
  resvecz<-vector(length=aux_length,mode="numeric")
  aux_yb<-vector(length=aux_length,mod="numeric")
  aux_yb_noinv<-vector(length=aux_length,mode="numeric")
  b_yb<-vector(length=aux_length,mode="numeric")
  aux_vv1_min<-NA
  aux_vv2_min<-NA
  aux_vv3_min<-NA
  aux_vv4_min<-NA
  aux_vv5_min<-NA
  aux_length<-length(b_yo)
#  print("X  Z  Yo")
#  print(matrix(cbind(b_x,b_z,b_yo),ncol=3,byrow=FALSE))
#  print(paste("aux_length = ",aux_length,sep=""))
  if (aux_length<=0) { return(NA) }
  avx<-mean(b_x)
  avz<-mean(b_z)
  avt<-mean(b_yo)
  vv[1]<-sum(b_yo)
#  print(paste("avx=",avx," avz=",avz," avt=",avt," vv[1]=",vv[1],sep=""))
  resmin<-1e+25
  zinv<-0
# background senza inversione
  ax2<-sum((b_x-avx)**2.)
  az2<-sum((b_z-avz)**2.)
  azx<-sum((b_z-avz)*(b_x-avx))
  axt<-sum((b_yo-avt)*(b_x-avx))
  azt<-sum((b_yo-avt)*(b_z-avz))
  det<-ax2*az2-azx*azx
  alp1<-(az2*axt-azx*azt)/det
  gam1<-(ax2*azt-azx*axt)/det
  aux_yb_noinv<-avt+alp1*(b_x-avx)+gam1*(b_z-avz)
  res_noinv<-sum((b_yo-aux_yb_noinv)**2)
#  print("Yb noinv:")
#  print(aux_yb_noinv)
# background con inversione
  for (z in seq(10,2000,by=10) ) {
    if (length(b_x[b_z>z])<6|length(b_x[b_z<=z])<6) {
      next
    }
#    print(paste("=#=# z >",z,sep=""))
#    print("b_x[b_z>z] b_z[b_z>z]")
#    print(matrix(cbind(b_x[b_z>z],b_z[b_z>z]),ncol=2,byrow=FALSE))
#    print("b_x[b_z<=z] b_z[b_z<=z]")
#    print(matrix(cbind(b_x[b_z<=z],b_z[b_z<=z]),ncol=2,byrow=FALSE))
    AA[1,1]<-length(b_yo)
    AA[1,2]<-sum(b_x[b_z>z]-avx)
    AA[1,3]<-sum(b_z[b_z>z]-z)
    AA[1,4]<-sum(b_x[b_z<=z]-avx)
    AA[1,5]<-sum(b_z[b_z<=z]-z)
    AA[2,1]<-AA[1,2]
    AA[2,2]<-sum((b_x[b_z>z]-avx)**2)
    AA[2,3]<-sum((b_z[b_z>z]-z)*(b_x[b_z>z]-avx))
    AA[2,4]<-0.
    AA[2,5]<-0.
    AA[3,1]<-AA[1,3]
    AA[3,2]<-AA[2,3]
    AA[3,3]<-sum((b_z[b_z>z]-z)**2)
    AA[3,4]<-0.
    AA[3,5]<-0.
    AA[4,1]<-AA[1,4]
    AA[4,2]<-0.
    AA[4,3]<-0.
    AA[4,4]<-sum((b_x[b_z<=z]-avx)**2)
    AA[4,5]<-sum((b_z[b_z<=z]-z)*(b_x[b_z<=z]-avx))
    AA[5,1]<-AA[1,5]
    AA[5,2]<-0.
    AA[5,3]<-0.
    AA[5,4]<-sum((b_x[b_z<=z]-avx)*(b_z[b_z<=z]-z))
    AA[5,5]<-sum((b_z[b_z<=z]-z)**2)
    vv[2]<-sum(b_yo[b_z>z]*(b_x[b_z>z]-avx))
    vv[3]<-sum(b_yo[b_z>z]*(b_z[b_z>z]-z))
    vv[4]<-sum(b_yo[b_z<=z]*(b_x[b_z<=z]-avx))
    vv[5]<-sum(b_yo[b_z<=z]*(b_z[b_z<=z]-z))
#    print("AA")
#    print(AA)
#    print("id5")
#    print(id5)
    InvAA<-try(solve(AA,id5))
    if (inherits(InvAA,"try-error")) {
      quit(status=1)
    }
#    print("InvAA")
#    print(InvAA)
#    print("InvAA%*%AA")
#    print(InvAA%*%AA)
    aux_vv<-InvAA %*% vv
#    print(paste("aux_vv [1,2,3,4,5]",aux_vv[1],aux_vv[2],aux_vv[3],aux_vv[4],aux_vv[5],sep=" "))
#    print(paste("    vv [1,2,3,4,5]",vv[1],vv[2],vv[3],vv[4],vv[5],sep=" "))
    aux_yb<-NA
    zab<-z+dz
    zbe<-z-dz
    bfab<-aux_vv[1]+aux_vv[2]*(b_x-avx)+aux_vv[3]*(b_z-z)
    bfbe<-aux_vv[1]+aux_vv[4]*(b_x-avx)+aux_vv[5]*(b_z-z)
    aux_yb[b_z>zab]<-aux_vv[1]+aux_vv[2]*(b_x[b_z>zab]-avx)+aux_vv[3]*(b_z[b_z>zab]-z)
    aux_yb[b_z<=zbe]<-aux_vv[1]+aux_vv[4]*(b_x[b_z<=zbe]-avx)+aux_vv[5]*(b_z[b_z<=zbe]-z)
    aux_yb[(b_z>zbe)&(b_z<=zab)]<-( bfab[(b_z>zbe)&(b_z<=zab)]*(b_z[(b_z>zbe)&(b_z<=zab)]-zbe) +
                                      bfbe[(b_z>zbe)&(b_z<=zab)]*(zab-b_z[(b_z>zbe)&(b_z<=zab)]) ) /
                                    (zab-zbe)
    aux_yb<-matrix(aux_yb,ncol=1)
#    aux_yb[b_z>z]<-aux_vv[1]+aux_vv[2]*(b_x[b_z>z]-avx)+aux_vv[3]*(b_z[b_z>z]-z)
#    aux_yb[b_z<=z]<-aux_vv[1]+aux_vv[4]*(b_x[b_z<=z]-avx)+aux_vv[5]*(b_z[b_z<=z]-z)
#    aux_yb<-matrix(aux_yb,ncol=1)
#    print(paste("aux_vv [1,2,3,4,5]",aux_vv[1],aux_vv[2],aux_vv[3],aux_vv[4],aux_vv[5],sep=" "))
#    print("z x yo yb")
#    ii<-order(b_z)
#    print(matrix(cbind(b_z[ii],b_x[ii],b_yo[ii],aux_yb[ii]),nrow=aux_length,ncol=4,byrow=FALSE))
    res<-sum((b_yo-aux_yb)**2)
    resvec[as.integer(z/10)]<-res
    resvecz[as.integer(z/10)]<-z
#    print(paste("        res resmin",res,resmin,sep=" "))
    if ( abs(aux_vv[2]<=2)&abs(aux_vv[4])<=2&res<resmin ) {
#      print(paste(" ZINV nuovo = ",z,sep=""))
      b_yb<-aux_yb
      aux_vv1_min<-aux_vv[1]
      aux_vv2_min<-aux_vv[2]
      aux_vv3_min<-aux_vv[3]
      aux_vv4_min<-aux_vv[4]
      aux_vv5_min<-aux_vv[5]
      b_yb<-aux_yb
      zinv<-z
      resmin<-res
    }
  }
#  postscript(file="caccola.ps")
#  plot(b_yb,b_z,col="gray",xlim=c(min(c(b_yb,b_yo)),max(c(b_yb,b_yo))),ylim=c(0,3500))
#  points(b_yo,b_z,col="red")
#  abline(v=(-16))
#  par(new=T)
#  plot(resvec,resvecz,xlim=c(0,10000000),type="l",ylim=c(0,3500),axes=F)
#  axis(3)
#  dev.off()
#  print(paste(" ZINV = ",zinv,sep=""))
#  print("z x yo yb")
#  ii<-order(b_z)
#  print(matrix(cbind(b_z[ii],b_x[ii],b_yo[ii],b_yb[ii]),nrow=aux_length,ncol=4,byrow=FALSE))
  if (resmin<res_noinv) {
    param[1]<-zinv
    param[2]<-avx
    param[3]<-aux_vv1_min
    param[4]<-resmin
    param[5]<-aux_vv2_min
    param[6]<-aux_vv3_min
    param[7]<-aux_vv4_min
    param[8]<-aux_vv5_min
  } else {
    param[1]<-avz
    param[2]<-avx
    param[3]<-avt
    param[4]<-res_noinv
    param[5]<-alp1
    param[6]<-gam1
    param[7]<-(-99999.)
    param[8]<-(-99999.)
  }
  return(param)
}

XYZinv<-function(b_x,b_y,b_z,b_yo,dz) {
  numdata<-length(b_yo)
  yb.aux<-vector(length=numdata,mode="numeric")
  param<-vector(length=11,mode="numeric")
  mx<-mean(b_x)
  my<-mean(b_y)
  mz<-mean(b_z)
  myo<-mean(b_yo)
  devx<-b_x-mx
  devy<-b_y-my
  devz<-b_z-mz
  devyo<-b_yo-myo
  sumyo<-sum(b_yo)
# background NOT allowing for inversion in the vertical profile
  ide<-matrix(ncol=3,nrow=3,data=0.)
  ide[row(ide)==col(ide)]<-1
# Ax=b
  res.noinv<-NA
  x.noinv<-vector(length=3,mode="numeric")
  x.noinv<-NA
  A<-matrix(ncol=3,nrow=3,data=0.)
  b<-vector(length=3)
  A[1,1]<-sum(devx**2.)
  A[1,2]<-sum(devx*devy)
  A[1,3]<-sum(devx*devz)
  A[2,1]<-A[1,2]
  A[2,2]<-sum(devy**2.)
  A[2,3]<-sum(devy*devz)
  A[3,1]<-A[1,3]
  A[3,2]<-A[2,3]
  A[3,3]<-sum(devz**2.)
  b[1]<-sum(devyo*devx)
  b[2]<-sum(devyo*devy)
  b[3]<-sum(devyo*devz)
  InvA<-solve(A,ide)
  x.noinv<-InvA %*% b
  yb.aux[]<-myo+ x.noinv[1]*devx+ x.noinv[2]*devy+ x.noinv[3]*devz
  res.noinv<-(mean((yb.aux-b_yo)**2))**0.5
# background allowing for inversion in the vertical profile
  ide<-matrix(ncol=7,nrow=7,data=0.)
  ide[row(ide)==col(ide)]=1
# Ax=b
  res.inv<-NA
  zinv<-NA
  x<-vector(length=7,mode="numeric")
  xmin<-vector(length=7,mode="numeric")
  xmin[]<-NA
  A<-matrix(ncol=7,nrow=7,data=0.)
  A[1,1]<-numdata
  b<-vector(length=7)
  b[1]<-sumyo
  for (z in seq(2000,10,by=-10) ) {
    devz<-b_z-z
    abv<-which(b_z>z)
    blw<-which(b_z<=z)
    A[,]<-0
    b[]<-0
    if (length(abv)<=8) next
    if (length(blw)<=8) break
    if (length(abv)>0) {
      A[1,2]<-sum(devx[abv])
      A[1,3]<-sum(devy[abv])
      A[1,4]<-sum(devz[abv])
      A[2,1]<-sum(devx[abv])
      A[2,2]<-sum(devx[abv]**2.)
      A[2,3]<-sum(devx[abv]*devy[abv])
      A[2,4]<-sum(devx[abv]*devz[abv])
      A[3,1]<-A[1,3]
      A[3,2]<-A[2,3]
      A[3,3]<-sum(devy[abv]**2.)
      A[3,4]<-sum(devy[abv]*devz[abv])
      A[4,1]<-A[1,4]
      A[4,2]<-A[2,4]
      A[4,3]<-A[3,4]
      A[4,4]<-sum(devz[abv]**2.)
      b[2]<-sum(b_yo[abv]*devx[abv])
      b[3]<-sum(b_yo[abv]*devy[abv])
      b[4]<-sum(b_yo[abv]*devz[abv])
    }
    if (length(blw)>0) {
      A[1,5]<-sum(devx[blw])
      A[1,6]<-sum(devy[blw])
      A[1,7]<-sum(devz[blw])
      A[5,5]<-sum(devx[blw]**2.)
      A[5,6]<-sum(devx[blw]*devy[blw])
      A[5,7]<-sum(devx[blw]*devz[blw])
      A[6,5]<-A[5,6]
      A[6,6]<-sum(devy[blw]**2.)
      A[6,7]<-sum(devy[blw]*devz[blw])
      A[7,5]<-A[5,7]
      A[7,6]<-A[6,7]
      A[7,7]<-sum(devz[blw]**2.)
      A[5,1]<-A[1,5]
      A[6,1]<-A[1,6]
      A[7,1]<-A[1,7]
      b[5]<-sum(b_yo[blw]*devx[blw])
      b[6]<-sum(b_yo[blw]*devy[blw])
      b[7]<-sum(b_yo[blw]*devz[blw])
    }
    InvA<-solve(A,ide)
    x<-InvA %*% b
    yb.aux[]<-NA
    zab<-z+dz
    zbe<-z-dz
    abvzab<-which(b_z>zab)
    blwzbe<-which(b_z<=zbe)
    zbe_zab<-which(((b_z>zbe)&(b_z<=zab)))
    bfab<-x[1]+ x[2]*devx+ x[3]*devy+ x[4]*(b_z-z)
    bfbe<-x[1]+ x[5]*devx+ x[6]*devy+ x[7]*(b_z-z)
    if (length(abvzab)>0)  yb.aux[abvzab]<-x[1]+ x[2]*devx[abvzab]+ x[3]*devx[abvzab]+ x[4]*(b_z[abvzab]-z)
    if (length(blwzbe)>0)  yb.aux[blwzbe]<-x[1]+ x[5]*devx[blwzbe]+ x[6]*devy[blwzbe]+ x[7]*(b_z[blwzbe]-z)
    if (length(zbe_zab)>0) yb.aux[zbe_zab]<- ( bfab[zbe_zab]*(b_z[zbe_zab]-zbe) +
                                               bfbe[zbe_zab]*(zab-b_z[zbe_zab]) ) / (2.*dz)
    res<-(mean((b_yo-yb.aux)**2))**0.5
    if (is.na(res.inv)) {
      xmin<-x
      zinv<-z
      res.inv<-res
    } else {
      if ( abs(x[2])<=2 & abs(x[3])<=2 & abs(x[5])<=2 & abs(x[6])<=2 & res<res.inv ) {
        xmin<-x
        zinv<-z
        res.inv<-res
      }
    }
  } # end cycle over z
  if (is.na(res.inv)) {
    param[1] <-mz
    param[2] <-mx
    param[3] <-my
    param[4] <-myo
    param[5] <-res.noinv
    param[6] <-x.noinv[1]
    param[7] <-x.noinv[2]
    param[8] <-x.noinv[3]
    param[9] <-NA
    param[10]<-NA
    param[11]<-NA
  } else {
    if (res.inv<=res.noinv) {
      param[1] <-zinv
      param[2] <-mx
      param[3] <-my
      param[4] <-xmin[1]
      param[5] <-res.inv
      param[6] <-xmin[2]
      param[7] <-xmin[3]
      param[8] <-xmin[4]
      param[9] <-xmin[5]
      param[10]<-xmin[6]
      param[11]<-xmin[7]
    } else {
      param[1] <-mz
      param[2] <-mx
      param[3] <-my
      param[4] <-myo
      param[5] <-res.noinv
      param[6] <-x.noinv[1]
      param[7] <-x.noinv[2]
      param[8] <-x.noinv[3]
      param[9] <-NA
      param[10]<-NA
      param[11]<-NA
    }
  }
  return(param)
}

XYZinv_NOapprox<-function(b_x,b_y,b_z,b_yo,dz) {
# In case of inversion
# 1.Zinv;2.mean(X);3.mean(Y);4.Tinv;5.res;6.alphaA;7.betaA;8.gammaA;9.alphaB;10.betaB;11.gammaB
# In case of NO inversion
# 1.mean(Z);2.mean(X);3.mean(Y);4.Tinv;5.res;6.alpha;7.beta;8.gamma;9.NA;10.NA;11.NA
  numdata<-length(b_yo)
  yb.aux<-vector(length=numdata,mode="numeric")
  param<-vector(length=11,mode="numeric")
  mx<-mean(b_x)
  my<-mean(b_y)
  mz<-mean(b_z)
  myo<-mean(b_yo)
  devx<-b_x-mx
  devy<-b_y-my
  devz<-b_z-mz
  devyo<-b_yo-myo
  dz2<-dz**2
# background NOT allowing for inversion in the vertical profile
  ide<-matrix(ncol=3,nrow=3,data=0.)
  ide[row(ide)==col(ide)]<-1
# Ax=b
  res.noinv<-NA
  x.noinv<-vector(length=3,mode="numeric")
  x.noinv<-NA
  A<-matrix(ncol=3,nrow=3,data=0.)
  b<-vector(length=3)
  A[1,1]<-sum(devx**2.)
  A[1,2]<-sum(devx*devy)
  A[1,3]<-sum(devx*devz)
  A[2,1]<-A[1,2]
  A[2,2]<-sum(devy**2.)
  A[2,3]<-sum(devy*devz)
  A[3,1]<-A[1,3]
  A[3,2]<-A[2,3]
  A[3,3]<-sum(devz**2.)
  b[1]<-sum(devyo*devx)
  b[2]<-sum(devyo*devy)
  b[3]<-sum(devyo*devz)
  InvA<-solve(A,ide)
  x.noinv<-InvA %*% b
  yb.aux[]<-myo+ x.noinv[1]*devx+ x.noinv[2]*devy+ x.noinv[3]*devz
  res.noinv<-(mean((yb.aux-b_yo)**2))**0.5
# background allowing for inversion in the vertical profile
  ide<-matrix(ncol=7,nrow=7,data=0.)
  ide[row(ide)==col(ide)]=1
# Ax=b
  res.inv<-NA
  zinv<-NA
  x<-vector(length=7,mode="numeric")
  xmin<-vector(length=7,mode="numeric")
  xmin[]<-NA
  A<-matrix(ncol=7,nrow=7,data=0.)
  AA<-matrix(ncol=3,nrow=3,data=0.)
  b<-vector(length=7)
  bb<-vector(length=3)
  for (z in seq(2000,10,by=-10) ) {
    devz<-b_z-z
    abv<-which(b_z>(z+dz))
    blw<-which(b_z<=(z-dz))
    inb<-which(b_z>(z-dz) & b_z<=(z+dz))
    A[,]<-0
    b[]<-0
    bb[]<-0
    if (length(abv)<=8) next
    if (length(blw)<=8) break
    AA[1,1]<-sum(devx[inb]**2)/(4*dz2)
    AA[1,2]<-sum(devx[inb]*devy[inb])/(4*dz2)
    AA[1,3]<-sum(devx[inb]*devz[inb])/(4*dz2)
    AA[2,1]<-AA[1,2]
    AA[2,2]<-sum(devy[inb]**2)/(4*dz2)
    AA[2,3]<-sum(devy[inb]*devz[inb])/(4*dz2)
    AA[3,1]<-AA[1,3]
    AA[3,2]<-AA[2,3]
    AA[3,3]<-sum(devz[inb]**2)/(4*dz2)
    A[1,1]<-length(abv)+length(blw)+length(inb)/dz2
    A[1,2]<-sum(devx[abv])+sum(devx[inb])/(2*dz2)
    A[1,3]<-sum(devy[abv])+sum(devy[inb])/(2*dz2)
    A[1,4]<-sum(devz[abv])+sum(devz[inb])/(2*dz2)
    A[1,5]<-sum(devx[blw])+sum(devx[inb])/(2*dz2)
    A[1,6]<-sum(devy[blw])+sum(devy[inb])/(2*dz2)
    A[1,7]<-sum(devz[blw])+sum(devz[inb])/(2*dz2)
    A[2,1]<-sum(devx[abv])
    A[2,2]<-sum(devx[abv]**2.)
    A[2,3]<-sum(devx[abv]*devy[abv])
    A[2,4]<-sum(devx[abv]*devz[abv])
    A[3,1]<-A[1,3]
    A[3,2]<-A[2,3]
    A[3,3]<-sum(devy[abv]**2.)
    A[3,4]<-sum(devy[abv]*devz[abv])
    A[4,1]<-A[1,4]
    A[4,2]<-A[2,4]
    A[4,3]<-A[3,4]
    A[4,4]<-sum(devz[abv]**2.)
    A[5,5]<-sum(devx[blw]**2.)
    A[5,6]<-sum(devx[blw]*devy[blw])
    A[5,7]<-sum(devx[blw]*devz[blw])
    A[6,5]<-A[5,6]
    A[6,6]<-sum(devy[blw]**2.)
    A[6,7]<-sum(devy[blw]*devz[blw])
    A[7,5]<-A[5,7]
    A[7,6]<-A[6,7]
    A[7,7]<-sum(devz[blw]**2.)
    A[5,1]<-A[1,5]
    A[6,1]<-A[1,6]
    A[7,1]<-A[1,7]
    A[2:4,2:4]<-A[2:4,2:4]+AA
    A[2:4,5:7]<-AA
    A[5:7,2:4]<-AA
    A[5:7,5:7]<-A[5:7,5:7]+AA
    bb[1]<-sum(b_yo[inb]*devx[inb])/(2*dz)
    bb[2]<-sum(b_yo[inb]*devy[inb])/(2*dz)
    bb[3]<-sum(b_yo[inb]*devz[inb])/(2*dz)
    b[1]<-sum(b_yo[abv])+sum(b_yo[blw])+sum(b_yo[inb])/dz
    b[2]<-sum(b_yo[abv]*devx[abv])+bb[1]
    b[3]<-sum(b_yo[abv]*devy[abv])+bb[2]
    b[4]<-sum(b_yo[abv]*devz[abv])+bb[3]
    b[5]<-sum(b_yo[blw]*devx[blw])+bb[1]
    b[6]<-sum(b_yo[blw]*devy[blw])+bb[2]
    b[7]<-sum(b_yo[blw]*devz[blw])+bb[3]
#
    InvA<-solve(A,ide)
    x<-InvA %*% b
    yb.aux[]<-NA
    zab<-z+dz
    zbe<-z-dz
    abvzab<-which(b_z>zab)
    blwzbe<-which(b_z<=zbe)
    zbe_zab<-which(((b_z>zbe)&(b_z<=zab)))
    bfab<-x[1]+ x[2]*devx+ x[3]*devy+ x[4]*(b_z-z)
    bfbe<-x[1]+ x[5]*devx+ x[6]*devy+ x[7]*(b_z-z)
    if (length(abvzab)>0)  yb.aux[abvzab]<-x[1]+ x[2]*devx[abvzab]+ x[3]*devy[abvzab]+ x[4]*(b_z[abvzab]-z)
    if (length(blwzbe)>0)  yb.aux[blwzbe]<-x[1]+ x[5]*devx[blwzbe]+ x[6]*devy[blwzbe]+ x[7]*(b_z[blwzbe]-z)
    if (length(zbe_zab)>0) yb.aux[zbe_zab]<- ( bfab[zbe_zab]*(b_z[zbe_zab]-zbe) +
                                               bfbe[zbe_zab]*(zab-b_z[zbe_zab]) ) / (2.*dz)
    res<-(mean((b_yo-yb.aux)**2))**0.5
    if (is.na(res.inv)) {
      xmin<-x
      zinv<-z
      res.inv<-res
    } else {
      if ( abs(x[2])<=2 & abs(x[3])<=2 & abs(x[5])<=2 & abs(x[6])<=2 & res<res.inv ) {
        xmin<-x
        zinv<-z
        res.inv<-res
      }
    }
  } # end cycle over z
  if (is.na(res.inv)) {
    param[1] <-mz
    param[2] <-mx
    param[3] <-my
    param[4] <-myo
    param[5] <-res.noinv
    param[6] <-x.noinv[1]
    param[7] <-x.noinv[2]
    param[8] <-x.noinv[3]
    param[9] <-NA
    param[10]<-NA
    param[11]<-NA
  } else {
    if (res.inv<=res.noinv) {
      param[1] <-zinv
      param[2] <-mx
      param[3] <-my
      param[4] <-xmin[1]
      param[5] <-res.inv
      param[6] <-xmin[2]
      param[7] <-xmin[3]
      param[8] <-xmin[4]
      param[9] <-xmin[5]
      param[10]<-xmin[6]
      param[11]<-xmin[7]
    } else {
      param[1] <-mz
      param[2] <-mx
      param[3] <-my
      param[4] <-myo
      param[5] <-res.noinv
      param[6] <-x.noinv[1]
      param[7] <-x.noinv[2]
      param[8] <-x.noinv[3]
      param[9] <-NA
      param[10]<-NA
      param[11]<-NA
    }
  }
  return(param)
}

XYZnoinv<-function(b_x,b_y,b_z,b_yo,dz) {
  numdata<-length(b_yo)
  yb.aux<-vector(length=numdata,mode="numeric")
  param<-vector(length=11,mode="numeric")
  mx<-mean(b_x)
  my<-mean(b_y)
  mz<-mean(b_z)
  myo<-mean(b_yo)
  devx<-b_x-mx
  devy<-b_y-my
  devz<-b_z-mz
  devyo<-b_yo-myo
  sumyo<-sum(b_yo)
# background NOT allowing for inversion in the vertical profile
  ide<-matrix(ncol=3,nrow=3,data=0.)
  ide[row(ide)==col(ide)]<-1
# Ax=b
  res.noinv<-NA
  x.noinv<-vector(length=3,mode="numeric")
  x.noinv<-NA
  A<-matrix(ncol=3,nrow=3,data=0.)
  b<-vector(length=3)
  A[1,1]<-sum(devx**2.)
  A[1,2]<-sum(devx*devy)
  A[1,3]<-sum(devx*devz)
  A[2,1]<-A[1,2]
  A[2,2]<-sum(devy**2.)
  A[2,3]<-sum(devy*devz)
  A[3,1]<-A[1,3]
  A[3,2]<-A[2,3]
  A[3,3]<-sum(devz**2.)
  b[1]<-sum(devyo*devx)
  b[2]<-sum(devyo*devy)
  b[3]<-sum(devyo*devz)
  InvA<-solve(A,ide)
  x.noinv<-InvA %*% b
  yb.aux[]<-myo+ x.noinv[1]*devx+ x.noinv[2]*devy+ x.noinv[3]*devz
  res.noinv<-(mean((yb.aux-b_yo)**2))**0.5
  param[1] <-mz
  param[2] <-mx
  param[3] <-my
  param[4] <-myo
  param[5] <-res.noinv
  param[6] <-x.noinv[1]
  param[7] <-x.noinv[2]
  param[8] <-x.noinv[3]
  param[9] <-NA
  param[10]<-NA
  param[11]<-NA
  return(param)
}

#==============================================================================
NonLinVertProf_2<-function(par,x,y,z,yo) {
  n<-length(x)
  h0<-par[1]
  h1<-par[1]+par[2]
#  if (h1<0) h1<-100
  T0<-par[3]
  gamma<-par[4]
  a<-par[5]
  alphaa<-par[6]
  alphab<-par[7]
  betaa<-par[8]
  betab<-par[9]
  yb<-vector(length=n,mode="numeric")
  x0<-mean(x)
  y0<-mean(y)
  ab<-which(z>=h1)
  bl<-which(z<=h0)
  bw<-which((z>h0)&(z<h1))
  if (length(ab)>0) yb[ab]<-T0+gamma*z[ab]  +alphaa*(x[ab]-x0)+betaa*(y[ab]-y0)
  if (length(bl)>0) yb[bl]<-T0+gamma*z[bl]-a+alphab*(x[bl]-x0)+betab*(y[bl]-y0)
  if (length(bw)>0) {
    ybAB.bw<-vector(length=length(bw),mode="numeric")
    ybBL.bw<-vector(length=length(bw),mode="numeric")
    ybAB.bw<-alphaa*(x[bw]-x0)+betaa*(y[bw]-y0)
    ybBL.bw<-alphab*(x[bw]-x0)+betab*(y[bw]-y0)
    yb[bw]<-T0+gamma*z[bw]-a/2.*(1+cos(pi*(z[bw]-h0)/(h1-h0))) + 
            (ybAB.bw*(z[bw]-h0)+ybBL.bw*(h1-z[bw]))/(h1-h0)
  }
#          alphaa*(x[bw]-x0)+betaa*(y[bw]-y0)+alphab*(x[bw]-x0)+betab*(y[bw]-y0))/(h1-h0)
#  yb[ab]<-T0-gamma*z[ab]NonVertProf
#  yb[bl]<-T0-gamma*z[bl]-a
#  yb[bw]<-T0-gamma*z[bw]-a/2.*(1+cos(pi*(z[bw]-h0)/(h1-h0)))
  J<-log(sum((yo-yb)**2.))
  return(J)
}
 
XYZinv_step2<-function(param,b_x,b_y,b_z,b_yo) {
##  param<-vector(mode="numeric",length=9)
  param.out<-vector(mode="numeric",length=11)
# vertical profile involves two linear section, one at upper levels 
# (above elevation h1) and one at lower levels (below elevation h0).
# Both with lapse rate gamma. T0 denotes the intercept (temperature at z0)
# of the upper linear section. The two sections are shifted against each 
# other by a temperature contrast "a", and they are matched in the 
# intermediate elevation range (h0,h1) by a smooth step function, imitating
# an inversion layer. (Frei, 2013)
# 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
#------------------------------------------------------------------------------
#  lm<-lm(b_z~b_yo)
#  m<-lm$coefficients[2]
#  q<-lm$coefficients[1]
# Alpha and Beta limit C/m -> 0.00008 C/m = 8 C/100Km
  ABlim<-0.00008
  gamma.inflim<--0.012
  gamma.suplim<--0.0001
  gamma.def<--0.006
  par2.inflim<-50
  par2.suplim<-300
  a.inflim<--10
  a.suplim<-10
  n.q10<-length(b_z)/20
  z.q10<-quantile(b_z,probs=0.1)
  z.q80<-quantile(b_z,probs=0.8)
#
  if (param[1]>=z.q80) param[1]<-(z.q80-10)
  if (param[2]<=par2.inflim  | param[2]>=par2.suplim)  param[2]<-(par2.inflim+par2.suplim)/2 
  if (param[4]<=gamma.inflim | param[4]>=gamma.suplim)  param[4]<-gamma.def 
  if (abs(param[6])>=ABlim)  param[6]<-0 
  if (abs(param[7])>=ABlim)  param[7]<-0 
  if (abs(param[8])>=ABlim)  param[8]<-0 
  if (abs(param[9])>=ABlim)  param[9]<-0 
#  n.iter<-1
#  n.iter.lim<-10
#  while (n.iter<n.iter.lim) {
    ui<-matrix(ncol=9,nrow=15,data=NA)
    ci<-vector(length=15)
    ui[1,]<-c(-1,0,0,0,0,0,0,0,0)
    ui[2,]<-c(0 ,1,0,0,0,0,0,0,0)
    ui[3,]<-c(0 ,-1,0,0,0,0,0,0,0)
    ui[4,]<-c(0,0,0,0,0,-1,0,0,0)
    ui[5,]<-c(0,0,0,0,0,0,-1,0,0)
    ui[6,]<-c(0,0,0,0,0,0,0,-1,0)
    ui[7,]<-c(0,0,0,0,0,0,0,0,-1)
    ui[8,]<-c(0,0,0,0,0, 1,0,0,0)
    ui[9,]<-c(0,0,0,0,0,0, 1,0,0)
    ui[10,]<-c(0,0,0,0,0,0,0, 1,0)
    ui[11,]<-c(0,0,0,0,0,0,0,0,1)
    ui[12,]<-c(0,0,0, 1,0,0,0,0,0)
    ui[13,]<-c(0,0,0,-1,0,0,0,0,0)
    ui[14,]<-c(0,0,0,0, 1,0,0,0,0)
    ui[15,]<-c(0,0,0,0,-1,0,0,0,0)
#    if (n.iter==1) h1lim<-0
#    if (n.iter>1) {
#      h1lim<-opt$par[2]+25
#      if (h1lim>290) h1lim<-290
#      param[2]<-h1lim+1
#    }
    ci[1:15]<-c(-z.q80,par2.inflim,-par2.suplim,
                -ABlim,-ABlim,-ABlim,-ABlim,-ABlim,-ABlim,-ABlim,-ABlim,
                gamma.inflim,-gamma.suplim,a.inflim,-a.suplim)
    opt<-constrOptim(param,NonLinVertProf_2,grad=NULL, ui=ui, ci=ci, mu = 1e-04, 
                     control=list(parscale=c(1000,100,50,0.01,10,0.0001,0.0001,0.0001,0.0001),
                     ndeps=c(.01,.01,.01,0.01,0.01,0.01,0.01,0.01,0.01),
                     maxit=10000),
                     outer.iterations = 100, outer.eps = 1e-05, x=b_x, y=b_y, z=b_z, yo=b_yo,
                     hessian = FALSE)
#    n.between<-length(b_z[b_z>=opt$par[1] & b_z<=(opt$par[1]+opt$par[2])])
##    print(paste(n.between,n.q10))
#    if (n.between>=n.q10) break
#    n.iter<-n.iter+1
#  }
#  if (opt$convergence==0 & n.iter<n.iter.lim) {
  if (opt$convergence==0) {
    param.out<-c(opt$par,mean(b_x),mean(b_y))
  } else {
    param.out<-c(rep(NA,11))
  }
  return(param.out)
}

#==============================================================================
NonLinVertProf_1<-function(par,x,y,z,yo,dz) {
# 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
  n<-length(x)
  yb<-vector(length=n,mode="numeric")
  yb.ab<-vector(length=n,mode="numeric")
  yb.bl<-vector(length=n,mode="numeric")
  zinv<-par[1]
  dz<-par[2]
  tinv<-par[3]
  alphaA<-par[6]
  alphaB<-par[7]
  betaA<-par[8]
  betaB<-par[9]
  gammaA<-par[4]
  gammaB<-par[5]
#
  x0<-mean(x)
  y0<-mean(y)
  devx<-x-x0
  devy<-y-y0
  devz<-z-zinv
  zab<-zinv+dz
  zbl<-zinv-dz
  ab<-which(z>zab)
  bl<-which(z<=zbl)
  bw<-which((z>zbl)&(z<=zab))
  yb.ab<-tinv+ alphaA*devx+ betaA*devy+ gammaA*devz
  yb.bl<-tinv+ alphaB*devx+ betaB*devy+ gammaB*devz
  if (length(ab)>0) yb[ab]<-yb.ab[ab]
  if (length(bl)>0) yb[bl]<-yb.bl[bl]
  if (length(bw)>0) yb[bw]<-(yb.ab[bw]*(z[bw]-zbl)+yb.bl[bw]*(zab-z[bw]))/(2.*dz)
  J<-log(sum((yo-yb)**2.))
  return(J)
}


#+ Background 1: allowing for inversion in the vertical profile
XYZinv_step1<-function(param,b_x,b_y,b_z,b_yo) {
# 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
  param.out<-vector(mode="numeric",length=11)
#------------------------------------------------------------------------------
# Alpha and Beta limit C/m -> 0.00008 C/m = 8 C/100Km
  ABlim<-0.00008
  gamma.inflim<--0.012
  gamma.suplim<--0.0001
  gamma.def<--0.006
  gamma.bl.inflim<--0.012
  gamma.bl.suplim<-0.010
  par2.inflim<-40
  par2.suplim<-60
  par2.def<-50
#
  z.q20<-quantile(b_z,probs=0.2)
  z.q80<-quantile(b_z,probs=0.8)
  dz.aux<-10
  if ((z.q80-z.q20)<10) dz.aux<-1
#
  if (param[1]>=z.q80) param[1]<-(z.q80-dz.aux)
  if (param[1]<=z.q20) param[1]<-(z.q20+dz.aux)
  if (param[2]<=par2.inflim | param[2]>=par2.suplim)  param[2]<-par2.def 
  if (param[4]<=gamma.inflim | param[4]>=gamma.suplim)  param[4]<-gamma.def 
  if (param[5]<=gamma.inflim | param[5]>=gamma.suplim)  param[5]<-gamma.def
  if (abs(param[6])>=ABlim)  param[6]<-0 
  if (abs(param[7])>=ABlim)  param[7]<-0 
  if (abs(param[8])>=ABlim)  param[8]<-0 
  if (abs(param[9])>=ABlim)  param[9]<-0 
#  print(paste("zinv dz tinv gA gB aA aB bA bB:",
#    round(param[1],1),
#    round(param[2],0),
#    round(param[3],1),
#    round(param[4],6),
#    round(param[5],6),
#    round(param[6],6),
#    round(param[7],6),
#    round(param[8],6),
#    round(param[9],6)))
# constrOptim; package:stats
# The feasible region is defined by ‘ui %*% param - ci >= 0’. The
# starting value must be in the interior of the feasible region, but
# the minimum may be on the boundary.
# then: ( 1 0 0 ...) * par -   inflim  --> par >= inflim 
#       (-1 0 0 ...) * par - (-suplim) --> par <= suplim
  ui<-matrix(ncol=9,nrow=16,data=NA)
  ci<-vector(length=16)
  ui[1,]<-c(1,0,0,0,0,0,0,0,0)
  ui[2,]<-c(-1,0,0,0,0,0,0,0,0)
  ui[3,]<-c(0,-1,0,0,0,0,0,0,0)
  ui[4,]<-c(0 ,1,0,0,0,0,0,0,0)
  ui[5,]<-c(0,0,0,0,0,-1,0,0,0)
  ui[6,]<-c(0,0,0,0,0,0,-1,0,0)
  ui[7,]<-c(0,0,0,0,0,0,0,-1,0)
  ui[8,]<-c(0,0,0,0,0,0,0,0,-1)
  ui[9,]<-c(0,0,0,0,0, 1,0,0,0)
  ui[10,]<-c(0,0,0,0,0,0, 1,0,0)
  ui[11,]<-c(0,0,0,0,0,0,0, 1,0)
  ui[12,]<-c(0,0,0,0,0,0,0,0,1)
  ui[13,]<-c(0,0,0, 1,0,0,0,0,0)
  ui[14,]<-c(0,0,0,-1,0,0,0,0,0)
  ui[15,]<-c(0,0,0,0,1,0,0,0,0)
  ui[16,]<-c(0,0,0,0,-1,0,0,0,0)
  ci[1:16]<-c(z.q20,-z.q80,-par2.suplim,par2.inflim,
              -ABlim,-ABlim,-ABlim,-ABlim,-ABlim,-ABlim,-ABlim,-ABlim,
              gamma.inflim,-gamma.suplim,gamma.bl.inflim,-gamma.bl.suplim)
  opt<-constrOptim(param,NonLinVertProf_1,grad=NULL, ui=ui, ci=ci, mu = 1e-04, 
                   control=list(parscale=c(1000,100,50,0.01,0.01,0.0001,0.0001,0.0001,0.0001),
                   ndeps=c(.01,.01,.01,0.01,0.01,0.01,0.01,0.01,0.01),
                   maxit=10000),
                   outer.iterations = 100, outer.eps = 1e-05, x=b_x, y=b_y, z=b_z, yo=b_yo,
                   hessian = FALSE)
  if (opt$convergence==0) {
    param.out<-c(opt$par,mean(b_x),mean(b_y))
  } else {
    param.out<-c(rep(NA,11))
  }
  return(param.out)
}


#==============================================================================
#+ Background 0: NOT allowing for inversion in the vertical profile
XYZnoinv_step0<-function(b_x,b_y,b_z,b_yo) {
# Alpha and Beta limit C/m -> 0.00008 C/m = 8 C/100Km
  ABlim<-0.00008
  gamma.inflim<--0.012
  gamma.suplim<--0.0001
  gamma.def<--0.006
# 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
# 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
  param.out<-vector(mode="numeric",length=11)
  param.out[]<-NA
#
  yb<-vector(length=length(b_yo),mode="numeric")
  mx<-mean(b_x)
  my<-mean(b_y)
#  mz<-mean(b_z)
  myo<-mean(b_yo)
  devx<-b_x-mx
  devy<-b_y-my
#  devz<-b_z-mz
  devyo<-b_yo-myo
# background NOT allowing for inversion in the vertical profile
  ide<-matrix(ncol=3,nrow=3,data=0.)
  ide[row(ide)==col(ide)]<-1
# Ax=b
  res.noinv<-NA
  x.noinv<-vector(length=3,mode="numeric")
  x.noinv<-NA
  A<-matrix(ncol=3,nrow=3,data=0.)
  b<-vector(length=3)
  A[1,1]<-sum(devx**2.)
  A[1,2]<-sum(devx*devy)
  A[1,3]<-sum(devx*b_z)
  A[2,1]<-A[1,2]
  A[2,2]<-sum(devy**2.)
  A[2,3]<-sum(devy*b_z)
  A[3,1]<-A[1,3]
  A[3,2]<-A[2,3]
  A[3,3]<-sum(b_z**2.)
  A.det<-A[1,1]*(A[2,2]*A[3,3]-A[2,3]*A[3,2])-
         A[1,2]*(A[2,1]*A[3,3]-A[2,3]*A[3,1])+
         A[1,3]*(A[2,1]*A[3,2]-A[2,2]*A[3,1])
  if (abs(A.det)<(1e-08)) {
    param.out[3]<-myo
    param.out[4]<-gamma.def
    param.out[6]<-0.
    param.out[8]<-0.
    param.out[10]<-mx
    param.out[11]<-my
    return(param.out)
  }
  b[1]<-sum(devyo*devx)
  b[2]<-sum(devyo*devy)
  b[3]<-sum(devyo*b_z)
  InvA<-solve(A,ide)
  x.noinv<-InvA %*% b
#  print(paste("mnZ mxZ alpha beta gamma:",min(b_z),max(b_z),
#              round(x.noinv[1],6),round(x.noinv[2],6),round(x.noinv[3],6)))
  if (x.noinv[1]<(-ABlim)) x.noinv[1]<--ABlim
  if (x.noinv[1]>ABlim) x.noinv[1]<-ABlim
  if (x.noinv[2]<(-ABlim)) x.noinv[2]<--ABlim
  if (x.noinv[2]>ABlim) x.noinv[2]<-ABlim
  if (x.noinv[3]<gamma.inflim | x.noinv[3]>=gamma.suplim)  x.noinv[3]<-gamma.def 
  param.out[3]<-myo
  param.out[4]<-x.noinv[3]
  param.out[6]<-x.noinv[1]
  param.out[8]<-x.noinv[2]
  param.out[10]<-mx
  param.out[11]<-my
  return(param.out)
}

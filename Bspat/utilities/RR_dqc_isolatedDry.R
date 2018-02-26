RR_dqc_isolatedDry<-function(obs,
                             rrinf=1,
                             pmax=20,
                             dmax=150,
                             n.sector=16,
                             dmin.dry=150) {
#------------------------------------------------------------------------------
  require(tripack)
# [] Contiguous NO-Rain (nor) Area identification
  p<-length(obs$stid)
  pmax<-min(p,pmax)
  psub<-vector(mode="numeric",length=p)
  psub[]<-NA
  sector.angle<-360/n.sector
  ydqc.flag<-vector(mode="numeric",length=p)
  ydqc.flag[]<-0
# First Guess
  dry_ix<-which(obs$yo<rrinf)
  if (length(dry_ix)==0) return(ydqc.flag)
  # vectors:
  # +nnor.vec>number of no-rain areas
  # +Lnor.vec[i]>(i=1,..,nnor.vec) number of stns in i-th no-rain area 
  # +nor.vec[i,n]-> n-th stn (i.e. pointers to VecS) in i-th no-rain area
  #  matrix(i,j) i=1,nnor.vec j=1,Lnor.vec[i]
  #  note: an no-rain areas could contain both wet and dry stations
  nor.vec<-matrix(data=NA,ncol=p,nrow=p)
  Lnor.vec<-vector(mode="numeric",length=p)
  Lnor.vec[]<-NA
  nnor.vec<-0
  nor.aux<-matrix(data=NA,ncol=p,nrow=p)
  Lnor.aux<-vector(mode="numeric",length=p)
  Lnor.aux[]<-NA
  nnor.aux<-0
  for (b in dry_ix) {  # START: Cycle over dry stations 
    # a. identify the closest stations to the b-th dry station
    #  outputs: -close2b-> pointer to VecS of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -psub[b]-> actual number of stations used
    Disth<-sqrt((obs$y[b]-obs$y)**2.+(obs$x[b]-obs$x)**2.)/1000.
    if (any(Disth<=dmax)) {
      ix<-which(Disth<=dmax)
      psub[b]<-min(pmax,length(ix))
      c2b<-order(Disth,decreasing=F)[1:psub[b]]
    } else {
      psub[b]<-0
    }
    if (psub[b]<3) next
    close2b<-c2b[1:psub[b]]
    rm(c2b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
    aux<-cbind(obs$x[close2b],obs$y[close2b])
    nokloop<-1
    while (anyDuplicated(aux)) {
      if (nokloop%%2==0) {
        obs$x[close2b][which(duplicated(aux))]<-obs$x[close2b][which(duplicated(aux))]+1
      } else {
        obs$y[close2b][which(duplicated(aux))]<-obs$y[close2b][which(duplicated(aux))]+1
      }
      aux<-cbind(obs$x[close2b],obs$y[close2b])
      nokloop<-nokloop+1
      if (nokloop>10) {
        print("ERROR: check the input data for too many duplicated stations")
        quit(status=1)
      }
    }
    tri.rr<-tri.mesh(obs$x[close2b],obs$y[close2b])
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th dry station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th dry station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,obs$x[b],obs$y[b])
    # note: due to the fact that (obs$x[b],obs$y[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-obs$x[b]+sin(sector.angle*(cont.clock-1)*pi/180.)*1000
      y.aux<-obs$y[b]+cos(sector.angle*(cont.clock-1)*pi/180.)*1000
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,
               nodes.aux[which( (!(nodes.aux %in% nodes)) &
                                  (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify nor 
    #    merge in a (new) nor all the (old) temporary nor (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary nor merged
    aux<-vector()
    if (nnor.aux>0) {
      for (nor in 1:nnor.aux) {
        if (Lnor.aux[nor]==0) next
        if ( any(nodes %in% nor.aux[nor,1:Lnor.aux[nor]]) ) {
          aux<-c(aux[which(!(aux %in% nor.aux[nor,1:Lnor.aux[nor]]))],nor.aux[nor,1:Lnor.aux[nor]])
          nor.aux[nor,]<-NA
          Lnor.aux[nor]<-0
        }
      }
    }
    nnor.aux<-nnor.aux+1
    Lnor.aux[nnor.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    nor.aux[nnor.aux,1:Lnor.aux[nnor.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the dry stations
# Reorganise no-rain area labels
  y.nor<-vector(length=p)
  y.nor[1:p]<-NA
  nnor.vec<-0
  if (nnor.aux>0) {
    for (nor in 1:nnor.aux) {
      if (Lnor.aux[nor]==0) next
      nnor.vec<-nnor.vec+1
      Lnor.vec[nnor.vec]<-Lnor.aux[nor]
      nor.vec[nnor.vec,1:Lnor.vec[nnor.vec]]<-nor.aux[nor,1:Lnor.aux[nor]]
      y.nor[nor.vec[nnor.vec,1:Lnor.vec[nnor.vec]]]<-nnor.vec
    }
  }
  rm(nor.aux,Lnor.aux,nnor.aux)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# identify dry-stations surrounded only by wet-stations (close enough)
  if (nnor.vec>0) {
    for (j in 1:nnor.vec) {
      nor.j<-nor.vec[j,1:Lnor.vec[j]]
      indx.dry.stns<-which(nor.j %in% dry_ix)
      n.dry<-length(indx.dry.stns)
      if (n.dry==1 | (n.dry<4 & (n.dry/Lnor.vec[j])<0.2) ) {
        nor.j.dry<-nor.j[indx.dry.stns]
        aux.dist<-(outer(obs$y[nor.j.dry],obs$y[nor.j],FUN="-")**2.+
                   outer(obs$x[nor.j.dry],obs$x[nor.j],FUN="-")**2.)**0.5/1000.
        min.dist.from.dry.stn<-min(aux.dist[aux.dist>0])
        if (min.dist.from.dry.stn<dmin.dry) 
          ydqc.flag[nor.j.dry]<-200
      }
    }
    rm(nor.j,indx.dry.stns)
  }
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.dry.stn,nor.j.dry)
  ydqc.flag
}

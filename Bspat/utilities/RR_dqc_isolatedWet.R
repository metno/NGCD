RR_dqc_isolatedWet<-function(obs,
                             rrinf=1,
                             pmax=20,
                             dmax=150,
                             n.sector=16,
                             dmin.wet=100) {
#------------------------------------------------------------------------------
  require(tripack)
# [] Precipitation events identification (eve) (contiguous rain areas)
  p<-length(obs$stid)
  pmax<-min(p,pmax)
  psub<-vector(mode="numeric",length=p)
  psub[]<-NA
  sector.angle<-360/n.sector
  ydqc.flag<-vector(mode="numeric",length=p)
  ydqc.flag[]<-0
# First Guess
  wet_ix<-which(obs$yo>=rrinf)
  if (length(wet_ix)==0) return(ydqc.flag) 
  # vectors:
  # +neve.vec>number of events 
  # +Leve.vec[i]>(i=1,..,neve.vec) number of stns in i-th event
  # +eve.vec[i,n]-> n-th stn (i.e. pointers to VecS) in i-th event
  #  matrix(i,j) i=1,neve.vec j=1,Leve.vec[i]
  #  note: an event could contain both wet and dry stations
  eve.vec<-matrix(data=NA,ncol=p,nrow=p)
  Leve.vec<-vector(mode="numeric",length=p)
  Leve.vec[]<-NA
  neve.vec<-0
  eve.aux<-matrix(data=NA,ncol=p,nrow=p)
  Leve.aux<-vector(mode="numeric",length=p)
  Leve.aux[]<-NA
  neve.aux<-0
  for (b in wet_ix) {  # START: Cycle over wet stations 
    # a. identify the closest stations to the b-th wet station
    #  outputs: -close2b-> position (i.e. pointer to VecS) of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -Lsubsample.vec[b]-> actual number of stations used
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
    #    to adjacent nodes (respect to the b-th wet station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th wet station
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
      nodes<-c(nodes,nodes.aux[which( (!(nodes.aux %in% nodes)) & (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify eve 
    #    merge in a (new) eve all the (old) temporary eve (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary eve merged
    aux<-vector()
    if (neve.aux>0) {
      for (eve in 1:neve.aux) {
        if (Leve.aux[eve]==0) next
        if ( any(nodes %in% eve.aux[eve,1:Leve.aux[eve]]) ) {
          aux<-c(aux[which(!(aux %in% eve.aux[eve,1:Leve.aux[eve]]))],eve.aux[eve,1:Leve.aux[eve]])
          eve.aux[eve,]<-NA
          Leve.aux[eve]<-0
        }
      }
    }
    neve.aux<-neve.aux+1
    Leve.aux[neve.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    eve.aux[neve.aux,1:Leve.aux[neve.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the wet stations
# reorganise event labels
  y.eve<-vector(length=p)
  y.eve[1:p]<-NA
  neve.vec<-0
  if (neve.aux>0) {
    for (eve in 1:neve.aux) {
      if (Leve.aux[eve]==0) next
      neve.vec<-neve.vec+1
      Leve.vec[neve.vec]<-Leve.aux[eve]
      eve.vec[neve.vec,1:Leve.vec[neve.vec]]<-eve.aux[eve,1:Leve.aux[eve]]
      y.eve[eve.vec[neve.vec,1:Leve.vec[neve.vec]]]<-neve.vec
    }
  }
  rm(eve.aux,Leve.aux,neve.aux)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# wet-stations surrounded only by dry-stations (close enough)
  if (neve.vec>0) {
    for (j in 1:neve.vec) {
      eve.j<-eve.vec[j,1:Leve.vec[j]]
      indx.wet.stns<-which(eve.j %in% wet_ix)
      n.wet<-length(indx.wet.stns)
      if (n.wet==1) {
        eve.j.wet<-eve.j[indx.wet.stns]
        aux.dist<-(outer(obs$y[eve.j.wet],obs$y[eve.j],FUN="-")**2.+
                   outer(obs$x[eve.j.wet],obs$x[eve.j],FUN="-")**2.)**0.5/1000.
        min.dist.from.wet.stn<-min(aux.dist[aux.dist>0])
        if (min.dist.from.wet.stn<dmin.wet) 
          ydqc.flag[eve.j.wet]<-300
      }
    }
  }
#  rm(eve.j,indx.wet.stns)
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.wet.stn,eve.j.wet)
  ydqc.flag
}

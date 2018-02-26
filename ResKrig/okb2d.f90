!***********************************************************************
	
      subroutine okb2d(nd,x,y,vr,nx,ny,xmn,ymn,xsiz,ysiz,nxdis,nydis, &
                      ndmin,ndmax,radius,nst,c0,it,aa,cc,ang,anis, &
                      idbg,ldbg,est,estv)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1992 Stanford Center for Reservoir Forecasting.  All   %
! rights reserved.  Distributed with: C.V. Deutsch and A.G. Journel.   %
! ``GSLIB: Geostatistical Software Library and User's Guide,'' Oxford  %
! University Press, New York, 1992.                                    %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!c-----------------------------------------------------------------------
!c
!c             Ordinary Kriging of a 2-D Rectangular Grid
!c             ******************************************
!c
!c This subroutine estimates point or block values of one variable by
!c ordinary kriging.  All of the samples are rescanned for each block
!c estimate; this makes the program simple but inefficient.  The data
!c should NOT contain any missing values.  Unestimated points are
!c returned as -1.0e21
!c
!c
!c
!c INPUT VARIABLES:
!c
!c   nd               Number of data (no missing values)
!c   x(nd),y(nd)      X and Y coordinates of the data
!c   vr(nd)           Data values
!c   nx,ny            Number of blocks in X and Y
!c   xmn              Coordinate at the center of the first X Block
!c   ymn              Coordinate at the center of the first Y Block
!c                      all blocks or grid nodes will have X and Y
!c                      coordinates equal to or greater than xmn and ymn
!c   xsiz,ysiz        X and Y spacing of the grid nodes (block size)
!c   nxdis            Number of discretization points/block in X
!c   nydis            Number of discretization points/block in Y
!c   ndmin            Minimum number of samples required for kriging
!c   ndmax            Maximum number of samples to use in kriging.
!c   radius           Maximum search radius (isotropic)
!c
!c   nst              Number of variogram structures
!c   c0               Nugget effect
!c   it(nst)          Variogram type (1=sph,2=exp,3=gaus,4=powr)
!c   aa(nst)          Range parameter except for the power model ("aa" is
!c                      the power parameter)
!c   cc(nst)          Contribution of each nested structure except for
!c                      the power model where "cc" is the slope
!c   ang(nst)         Azimuth angle for anisotropy (measured clockwise
!c                      from the Y axis) for each variogram structure.
!c   anis(nst)        Anisotropy ratio, i.e., the range in the direction
!c                      90 degrees from "ang" divided by the range given
!c                      by "aa".  For example, if the variogram model has
!c                      a range of 20 at an azimuth of N30E and a range 
!c                      of 10 in the direction N120E the parameters 
!c                      should be entered as:  aa=20, ang=30, anis=0.5
!c   idbg             The debug level:
!c                      0. None
!c                      1. Input parameters and summary statistics
!c                      2. The data used for each kriging along with the
!c                         kriging weights are printed.
!c                      3. Kriging weights, values and variances / block
!c                         kriging matrices also printed (serious)
!c   ldbg             The output unit for the debugging output
!c
!c
!c
!c OUTPUT VARIABLES:
!c
!c   est(nx,ny)        the kriged values
!c   estv(nx,ny)       the kriging variances
!c
!c
!c
!c EXTERNAL REFERENCES:
!c
!c   cova2            Calculates the covariance given a variogram model
!c   ksol             Linear system solver using Gaussian Elimination
!c                    (double precision)
!c
!c
!c
!c Original:  A.G. Journel                                           1978
!c Revisions: B.E. Buxton                                       Apr. 1983
!c            C.V. Deutsch                                      July 1990
!c-----------------------------------------------------------------------
      include  'OKB2D.INC'
      implicit real (a-h,o-z)
      real      x(MAXDAT),y(MAXDAT),vr(MAXDAT),aa(MAXNST),cc(MAXNST)
      real      ang(MAXNST),anis(MAXNST),est(MAXX,MAXY),estv(MAXX,MAXY)
      real      xdb(MAXDIS),ydb(MAXDIS),xa(MAXSAM),ya(MAXSAM),
      real      vra(MAXSAM),dist(MAXSAM)
      real*8    r(MAXSAM+1),rr(MAXSAM+1),s(MAXSAM+1),a(MAXKRG)
      integer   it(MAXNST),nums(MAXSAM)
      logical   first
      data      first/.true./,PMX/9999.0/
      
!c
!c Echo the input parameters if debugging flag is >2:
!c
      if(idbg.gt.2) then
            write(ldbg,*) 'OKB2D Parameters'
            write(ldbg,*)
            write(ldbg,*) 'Variogram Parameters for ',nst,' structures:'
            write(ldbg,*) '  Nugget effect:         ',c0
            write(ldbg,*) '  Types of variograms:   ',(it(i),i=1,nst)
            write(ldbg,*) '  Contribution cc        ',(cc(i),i=1,nst)
            write(ldbg,*) '  Ranges:                ',(aa(i),i=1,nst)
            write(ldbg,*) '  Angle for Continuity:  ',(ang(i),i=1,nst)
            write(ldbg,*) '  Anisotropy Factors:    ',(anis(i),i=1,nst)
            write(ldbg,*) ' '
            write(ldbg,*) 'Grid for Kriging:'
            write(ldbg,*) '  Number of X and Y Blocks:',nx,ny
            write(ldbg,*) '  Origin of X and Y Blocks:',xmn,ymn
            write(ldbg,*) '  Size   of X and Y Blocks:',xsiz,ysiz
            write(ldbg,*) ' '
            write(ldbg,*) 'Discretization of blocks:  ',nxdis,nydis
            write(ldbg,*) 'Search Radius:             ',radius
            write(ldbg,*) 'Minimum number of samples: ',ndmin
            write(ldbg,*) 'Maximum number of samples: ',ndmax
            write(ldbg,*) ' '
      endif
!c
!c Echo the input data if debugging flag >1:
!c
      if(idbg.ge.4) then
            do 1 id=1,nd
                  write(ldbg,99) id,x(id),y(id),vr(id)
 99               format('Data: ',i5,' at ',2f12.3,' value: ',f12.5)
 1          continue
      endif
!c
!c Set up the discretization points per block.  Figure out how many
!c are needed, the spacing, and fill the xdb and ydb arrays with the
!c offsets relative to the block center (this only gets done once):
!c
      ndb  = nxdis * nydis
      if(ndb.gt.MAXDIS) then
            write(*,*) 'ERROR OKB2D: Too many discretization points '
            write(*,*) '             Increase MAXDIS or lower n[xy]dis'
            stop
      endif
      xdis = xsiz  / amax1(real(nxdis),1.0)
      ydis = ysiz  / amax1(real(nydis),1.0)
      xloc = -0.5*(xsiz+xdis)
      i    = 0
      do 2 ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do 2 iy=1,nydis
                  yloc = yloc + ydis
                  i = i+1
                  xdb(i) = xloc
                  ydb(i) = yloc
 2    continue
!c
!c Initialize accumulators:
!c
      uk   = 0.0
      vk   = 0.0
      nk   = 0
      cbb  = 0.0
      rad2 = radius*radius
!c
!c Calculate Block Covariance. Check for point kriging.
!c
      cov   = cova2(xdb(1),ydb(1),xdb(1),ydb(1),nst,c0,PMX,cc, &
                   aa,it,ang,anis,first)
!c
!c Keep this value to use for the unbiasedness constraint:
!c
      unbias = cov
      first  = .false.
      if (ndb.le.1) then
            cbb = cov
      else
            do 3 i=1,ndb
            do 3 j=1,ndb
                  cov = cova2(xdb(i),ydb(i),xdb(j),ydb(j),nst,c0, &
                             PMX,cc,aa,it,ang,anis,first)
                  if(i.eq.j) cov = cov - c0
                  cbb = cbb + cov
 3          continue
            cbb = cbb/real(ndb*ndb)
      endif
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      endif
!c
!c MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
!c
      do 4 iy=1,ny
      yloc = ymn + (iy-1)*ysiz
      do 4 ix=1,nx
      
     
            xloc = xmn + (ix-1)*xsiz
!c
!c Find the nearest samples within each octant: First initialize
!c the counter arrays:
!c
            na = 0
            do 5 isam=1,ndmax
                  dist(isam) = 1.0e+20
                  nums(isam) = 0
 5          continue
!c
!c Scan all the samples (this is inefficient and the user with lots of
!c data should move to ktb3d):
!c
            do 6 id=1,nd
                  dx = x(id) - xloc
                  dy = y(id) - yloc
                  h2 = dx*dx + dy*dy
                  if(h2.gt.rad2) goto 6
!c
!c Do not consider this sample if there are enough close ones:
!c
                  if(na.eq.ndmax.and.h2.gt.dist(na)) goto 6
!c
!c Consider this sample (it will be added in the correct location):
!c
                  if(na.lt.ndmax) na = na + 1
                  nums(na)           = id
                  dist(na)           = h2
                  if(na.eq.1) goto 6
!c
!c Sort samples found thus far in increasing order of distance:
!c
                  n1 = na-1
                  do 7 ii=1,n1
                        k=ii
                        if(h2.lt.dist(ii)) then
                              jk = 0
                              do 8 jj=k,n1
                                    j  = n1-jk
                                    jk = jk+1
                                    j1 = j+1
                                    dist(j1) = dist(j)
                                    nums(j1) = nums(j)
 8                            continue
                              dist(k) = h2
                              nums(k) = id
                              goto 6
                        endif
 7                continue
 6          continue
!c
!c Is there enough samples?
!c
            if(na.lt.ndmin) then
                  if(idbg.ge.2) &
                 write(ldbg,*) 'Block ',ix,iy, 'not estimated'
                  est(ix,iy)  = UNEST
                  estv(ix,iy) = UNEST
                  goto 4
            endif
!c
!c Put coordinates and values of neighborhood samples into xa,ya,vra:
!c
            do 10 ia=1,na
                  jj      = nums(ia)
                  xa(ia)  = x(jj)
                  ya(ia)  = y(jj)
                  vra(ia) = vr(jj)
 10         continue
!c
!c Handle the situation of only one sample:
!c
            if(na.eq.1) then
                  cb1 = cova2(xa(1),ya(1),xa(1),ya(1),nst,c0, &
                             PMX,cc,aa,it,ang,anis,first)
                  xx  = xa(1) - xloc
                  yy  = ya(1) - yloc
!c
!c Establish Right Hand Side Covariance:
!c
                  if(ndb.le.1) then
                        cb = cova2(xx,yy,xdb(1),ydb(1),nst,c0, &
                                  PMX,cc,aa,it,ang,anis,first)
                  else
                        cb  = 0.0
                        do 12 i=1,ndb
                              cb = cb + cova2(xx,yy,xdb(i),ydb(i),nst, &
                                       c0,PMX,cc,aa,it,ang,anis,first)
                              dx = xx - xdb(i)
                              dy = yy - ydb(i)
                              if((dx*dx+dy*dy).lt.EPSLON) &
                             cb = cb - c0
 12                     continue
                        cb = cb / real(ndb)
                  end if
                  est(ix,iy)  = vra(1)
                  estv(ix,iy) = cbb - 2.0*cb + cb1
            else
!c
!c Solve the Kriging System with more than one sample:
!c
                  neq = na + 1
                  nn  = (neq + 1)*neq/2
!c
!c Set up kriging matrices:
!c
                  in=0
                  do 13 j=1,na
!c
!c Establish Left Hand Side Covariance Matrix:
!c
                        do 14 i=1,j
                              in = in + 1
                              a(in) = dble( cova2(xa(i),ya(i),xa(j), &
                                           ya(j),nst,c0,PMX,cc,aa,it, &
                                           ang,anis,first) )
 14                     continue
                        xx = xa(j) - xloc
                        yy = ya(j) - yloc
!c
!c Establish Right Hand Side Covariance:
!c
                        if(ndb.le.1) then
                              cb = cova2(xx,yy,xdb(1),ydb(1),nst,c0, &
                                        PMX,cc,aa,it,ang,anis,first)
                        else
                              cb  = 0.0
                              do 15 j1=1,ndb
                              cb = cb + cova2(xx,yy,xdb(j1),ydb(j1),nst, &
                                       c0,PMX,cc,aa,it,ang,anis,first)
                              dx = xx - xdb(j1)
                              dy = yy - ydb(j1)
                              if((dx*dx+dy*dy).lt.EPSLON) &
                             cb = cb - c0
 15                           continue
                              cb = cb / real(ndb)
                        end if
                        r(j)  = dble(cb)
                        rr(j) = r(j)
 13               continue
!c
!c Set the unbiasedness constraint:
!c
                  do 16 i=1,na
                        in    = in + 1
                        a(in) = dble(unbias)
 16               continue
                  in      = in + 1
                  a(in)   = 0.0
                  r(neq)  = dble(unbias)
                  rr(neq) = r(neq)
!c
!c Write out the kriging Matrix if Seriously Debugging:
!c
                  if(idbg.ge.3) then
                        write(ldbg,101) ix,iy
                        is = 1
                        do 100 i=1,neq
                              ie = is + i - 1
                              write(ldbg,102) i,r(i),(a(j),j=is,ie)
                              is = is + i
 100                    continue
 101                    format(/,'Kriging Matrices for Node: ',2i4, &
                                ' RHS first')
 102                    format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
                  endif
!c
!c Solve the Kriging System:
!c
                  call ksol(1,neq,1,a,r,s,ising)
!c
!c Write a warning if the matrix is singular:
!c
                  if(ising.ne.0) then
                        write(*,*) 'WARNING OKB2D: singular matrix'
                        write(*,*) '               for block',ix,iy
                        est(ix,iy)  = UNEST
                        estv(ix,iy) = UNEST
                        goto 4
                  endif
!c
!c Write the kriging weights and data if requested:
!c
                  if(idbg.ge.2) then
                        write(ldbg,*) '       '
                        write(ldbg,*) 'BLOCK: ',ix,iy
                        write(ldbg,*) '       '
                        write(ldbg,*) '  Lagrange multiplier: ',s(neq)
                        write(ldbg,*) '  BLOCK EST: x,y,vr,wt '
                        do 18 i=1,na
 18                     write(ldbg,'(4f12.3)') xa(i),ya(i),vra(i),s(i)
                  endif
!c
!c Compute the estimate and the kriging variance:
!c
                  ook  = 0.0
                  ookv = cbb-real(s(na+1))
                  do 19 i=1,na
                        ook  = ook  + real(s(i))*vra(i)
                        ookv = ookv - real(s(i))*real(rr(i))
 19               continue
                  est(ix,iy)  = ook
                  estv(ix,iy) = ookv
            endif
!c
!c Accumulate statistics of kriged blocks:
!c
            nk = nk + 1
            uk = uk + est(ix,iy)
            vk = vk + est(ix,iy)*est(ix,iy)
!c
!c Write more details if requested:
!c
            if(idbg.ge.2) then
                  write(ldbg,*) '  est  ',est(ix,iy)
                  write(ldbg,*) '  estv ',estv(ix,iy)
                  write(ldbg,*) ' '
            endif
!c
!c END OF MAIN LOOP OVER ALL THE BLOCKS:
!c

      
      if(ix==500.and.mod(iy,100)==0)write(99,*) ix,iy,est(ix,iy)
      if(est(ix,iy)<-50. .or. est(ix,iy)>50.) then
	    write(99,*)"Something's wrong here: ", ix,iy,est(ix,iy)	  
      endif	  

 4    continue
!c
!c Write statistics of kriged values:
!c
      if(nk.gt.0.and.idbg.gt.0) then
            vk = (vk-uk*uk/real(nk))/real(nk)
            uk = uk/real(nk)
            write(ldbg,*) ' '
            write(ldbg,*) 'Estimated  ',nk,' blocks '
            write(ldbg,*) '  average  ',uk
            write(ldbg,*) '  variance ',vk
      endif


      return
      end
 
 
 
      real function cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it, &
                         ang,anis,first)
!c-----------------------------------------------------------------------
!c
!c              Covariance Between Two Points (2-D Version)
!c              *******************************************
!c
!c This function returns the covariance associated with a variogram model
!c that is specified by a nugget effect and possibly four different
!c nested varigoram structures.  The anisotropy definition can be
!c different for each of the nested structures (spherical, exponential,
!c gaussian, or power).
!c
!c
!c
!c INPUT VARIABLES:
!c
!c   x1,y1            Coordinates of first point
!c   x2,y2            Coordinates of second point
!c   nst              Number of nested structures (max. 4).
!c   c0               Nugget constant (isotropic).
!c   PMX              Maximum variogram value needed for kriging when
!c                      using power model.  A unique value of PMX is
!c                      used for all nested structures which use the
!c                      power model.  therefore, PMX should be chosen
!c                      large enough to account for the largest single
!c                      structure which uses the power model.
!c   cc(nst)          Multiplicative factor of each nested structure.
!c   aa(nst)          Parameter "a" of each nested structure.
!c   it(nst)          Type of each nested structure:
!c                      1. spherical model of range a;
!c                      2. exponential model of parameter a;
!c                           i.e. practical range is 3a
!c                      3. gaussian model of parameter a;
!c                           i.e. practical range is a*sqrt(3)
!c                      4. power model of power a (a must be gt. 0  and
!c                           lt. 2).  if linear model, a=1,c=slope.
!c   ang(nst)         Azimuth angle for the principal direction of
!c                      continuity (measured clockwise in degrees from Y)
!c   anis(nst)        Anisotropy (radius in minor direction at 90 degrees
!c                      from "ang" divided by the principal radius in 
!c                      direction "ang")
!c   first            A logical variable which is set to true if the
!c                      direction specifications have changed - causes
!c                      the rotation matrices to be recomputed.
!c
!c
!c
!c OUTPUT VARIABLES: returns "cova2" the covariance obtained from the
!c                   variogram model.
!c
!c
!c
!c AUTHOR: C.V. Deutsch                                   DATE: July 1990
!c-----------------------------------------------------------------------
      parameter(DTOR=3.14159265/180.0,EPSLON=0.0000001)
      real      aa(*),cc(*),ang(*),anis(*),rotmat(4,4),maxcov
      integer   it(*)
      logical   first
      save      rotmat,maxcov
!c
!c The first time around, re-initialize the cosine matrix for the
!c variogram structures:
!c
      if(first) then
            maxcov = c0
            do 1 is=1,nst
                  azmuth       = (90.0-ang(is))*DTOR
                  rotmat(1,is) =  cos(azmuth)
                  rotmat(2,is) =  sin(azmuth)
                  rotmat(3,is) = -sin(azmuth)
                  rotmat(4,is) =  cos(azmuth)
                  if(it(is).eq.4) then
                        maxcov = maxcov + PMX
                  else
                        maxcov = maxcov + cc(is)
                  endif
 1          continue
      endif
!c
!c Check for very small distance:
!c
      dx = x2-x1
      dy = y2-y1
      if((dx*dx+dy*dy).lt.EPSLON) then
            cova2 = maxcov
            return
      endif
!c
!c Non-zero distance, loop over all the structures:
!c
      cova2 = 0.0
      do 2 is=1,nst
!c
!c Compute the appropriate structural distance:
!c
            dx1 = (dx*rotmat(1,is) + dy*rotmat(2,is))
            dy1 = (dx*rotmat(3,is) + dy*rotmat(4,is))/anis(is)
            h   = sqrt(dx1*dx1+dy1*dy1)
            if(it(is).eq.1) then
!c
!c Spherical model:
!c
                  hr = h/aa(is)
                  if(hr.ge.1.0) goto 2
                  cova2 = cova2 + cc(is)*(1.-hr*(1.5-.5*hr*hr))
            else if(it(is).eq.2) then
!c
!c Exponential model:
!c
                  cova2 = cova2 +cc(is)*exp(-h/aa(is))
            else if(it(is).eq. 3) then
!c
!c Gaussian model:
!c
                  hh=-(h*h)/(aa(is)*aa(is))
                  cova2 = cova2 +cc(is)*exp(hh)
            else
!c
!c Power model:
!c
                  cov1  = PMX - cc(is)*(h**aa(is))
                  cova2 = cova2 + cov1
            endif
 2    continue
      return
      end
 
 
 
      subroutine ksol(nright,neq,nsb,a,r,s,ising)
!c-----------------------------------------------------------------------
!c
!c                Solution of a System of Linear Equations
!c                ****************************************
!c
!c
!c
!c INPUT VARIABLES:
!c
!c   nright,nsb       number of columns in right hand side matrix.
!c                      for OKB2D: nright=1, nsb=1
!c   neq              number of equations
!c   a()              upper triangular left hand side matrix (stored 
!c                      columnwise)
!c   r()              right hand side matrix (stored columnwise)
!c                      for okb2d, one column per variable
!c
!c
!c
!c OUTPUT VARIABLES:
!c
!c   s()              solution array, same dimension as  r  above.
!c   ising            singularity indicator
!c                      0,  no singularity problem
!c                     -1,  neq .le. 1
!c                      k,  a null pivot appeared at the kth iteration
!c
!c
!c
!c PROGRAM NOTES:
!c
!c   1. Requires the upper triangular left hand side matrix.
!c   2. Pivots are on the diagonal.
!c   3. Does not search for max. element for pivot.
!c   4. Several right hand side matrices possible.
!c   5. USE for ok and sk only, NOT for UK.
!c
!c
!c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
!c
!c If there is only one equation then set ising and return:
!c
      if(neq.le.1) then
            ising = -1
            return
      endif
!c
!c Initialize:
!c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
!c
!c Start triangulation:
!c
      do 70 k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do 60 iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do 50 i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do 30 j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
 30                     continue
                        do 40 llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
 40                     continue
 50               continue
 60         continue
 70   continue
!c
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
!c
!c Finished triangulation, start solving back:
!c
      do 140 iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do 100 llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
 100        continue
            i=neq
            kk=ij
            do 130 ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do 120 llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do 110 j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
 110                    continue
                        s(ll1)=ap*piv
 120              continue
 130        continue
 140  continue
!c
!c Finished solving back, return:
!c
      return
      end




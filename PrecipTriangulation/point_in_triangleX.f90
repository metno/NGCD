      program point_in_triangle
      
      implicit none
	  
      integer, parameter :: node_num_max = 2500, triangle_num_max = 5000
      integer, parameter :: dim_num = 2
  
      integer :: i,j,ii,io,jj,n,tri,point_within_how_many_triangles,triangle_num,node_num,it,stnr,inint,stnrx
	  
      integer (kind=4 ), dimension(1,2) :: x_grid_indices_norge
      integer (kind=4 ), dimension(1,1) :: point_within_triangle,buffer
      integer (kind=4 ), dimension (3,triangle_num_max) :: triangle_nodes_index 	  
      integer (kind=4 ) :: istatus,itmp
      real ( kind = 8 ), dimension (dim_num,node_num_max) :: node_xy 
      real ( kind = 8 ), dimension (node_num_max) :: rr,z_amsl 
      real ( kind = 8 ), dimension (1) :: x
      real ( kind = 8 ), dimension (1) :: y
      real ( kind = 8 ), dimension (3) :: x_triangle,y_triangle
      real ( kind = 8 ) :: invDenom,u,v,dot00,dot01,dot02,dot11,dot12,diffz,rr2_grd,z_grd,distance,distance_to_nearest_neighbour
      real ( kind = 4 ) :: unew,vnew,u_plus_v,yyyyy
      real ( kind = 8 ) :: mean_rr1km_out,std_rr1km_out,zfac
      real ( kind = 8 ), dimension(1,1)  :: rr_grd,z_amsl_grd,dem1km,rr1km
      real ( kind = 8 ), dimension(1,1)  :: z_amsl_grd_new,rr_grd_new
      real ( kind = 8 ) :: x1,x2,x3,h1,h2,h3,Denom,y1,y2,y3,d_rr_over_d_x,d_rr_over_d_y,d_z_over_d_x,d_z_over_d_y,rr1,rr2,rr3,rrx
      character cdum*1	  
          
!---- Read daily precipitation data from file "rr_x.par" ----	  	  
      zfac=0.001

      open(1,file='rr_x.par',status='old')
      do j=1,10000
        read(1,*,iostat=io)(node_xy(i,j),i=1,2),rr(j),it,stnr,z_amsl(j)
        if(io<0) exit		
        write(79,*)(node_xy(i,j),i=1,2),rr(j)
      enddo
      close(1)
      node_num=j-1
!      write(*,*)'Number of nodes/stations =', node_num

!---- Read triangle nodes computed by Geompack ----	  

      open(1,file='triangle_nodes_X.dat',status='old')
      do j=1,10000
        read(1,*,iostat=io)(triangle_nodes_index(i,j),i=1,3)
        if(io<0) exit			
        write(79,*)(triangle_nodes_index(i,j),i=1,3)      
      enddo
      close(1)
      triangle_num=j-1
!      write(*,*)'Number of triangles = ', triangle_num

!---- Read candidate station ----	  
      open(1,file='x_station.par',status='old')
      read(1,*,iostat=io)x(1),y(1),rrx,it,stnrx,dem1km(1,1)
!      write(*,*)x(1),y(1),rrx,it,stnrx,dem1km(1,1)
!----
!---- Loop over all grid points of seNorge.no ----
!----

      do i=1,1
        do j=1,1
!          if (buffer(i,j)<0) cycle	
          point_within_how_many_triangles=0

          point_within_triangle(i,j)=-999
!          point_within_triangle(i,j)=10000

!---- Loop over all triangles ----

          do tri=1,triangle_num 

!---- Set coordinates of nodes 1, 2 and 3 of current triangle ----

            x_triangle(1)=node_xy(1,triangle_nodes_index(1,tri))
            y_triangle(1)=node_xy(2,triangle_nodes_index(1,tri))
            x_triangle(2)=node_xy(1,triangle_nodes_index(2,tri))
            y_triangle(2)=node_xy(2,triangle_nodes_index(2,tri))
            x_triangle(3)=node_xy(1,triangle_nodes_index(3,tri))
            y_triangle(3)=node_xy(2,triangle_nodes_index(3,tri))

	
!---- If grid point not within square comprised of current triangle points skip some operations to save time ----
            if(x(i) < min(x_triangle(1),x_triangle(2),x_triangle(3))) goto 33333
            if(x(i) > max(x_triangle(1),x_triangle(2),x_triangle(3))) goto 33333		 
            if(y(j) < min(y_triangle(1),y_triangle(2),y_triangle(3))) goto 33333		 
            if(y(j) > max(y_triangle(1),y_triangle(2),y_triangle(3))) goto 33333		 		 

!---- 
!---- Perform "Point in triangle test" according to Barycentric Technique (see http://www.blackpawn.com/texts/pointinpoly/default.html) ----
!----

            dot00=(x_triangle(3)-x_triangle(1))*(x_triangle(3)-x_triangle(1))&
                 & + (y_triangle(3)-y_triangle(1))*(y_triangle(3)-y_triangle(1))
            dot01=(x_triangle(3)-x_triangle(1))*(x_triangle(2)-x_triangle(1))&
                 & + (y_triangle(3)-y_triangle(1))*(y_triangle(2)-y_triangle(1))
            dot02=(x_triangle(3)-x_triangle(1))*(x(i)-x_triangle(1))&
                 & + (y_triangle(3)-y_triangle(1))*(y(j)-y_triangle(1))
            dot11=(x_triangle(2)-x_triangle(1))*(x_triangle(2)-x_triangle(1))&
                 & + (y_triangle(2)-y_triangle(1))*(y_triangle(2)-y_triangle(1))
            dot12=(x_triangle(2)-x_triangle(1))*(x(i)-x_triangle(1))&
                 & + (y_triangle(2)-y_triangle(1))*(y(j)-y_triangle(1))
  	
  
            invDenom = 1.d0 / (dot00 * dot11 - dot01 * dot01)
            u = (dot11 * dot02 - dot01 * dot12) * invDenom
            v = (dot00 * dot12 - dot01 * dot02) * invDenom 
  
            !unew = int(u*100.)/100. ! Truncate/Eliminate some rounding errors in last decimals (might not be required/give wrong results)
            !vnew = int(v*100.)/100. ! Truncate/Eliminate some rounding errors in last decimals (might not be required/give wrong results)
            !u_plus_v=unew+vnew
   
            if(u >= 0. .and. v >= 0. .and. u + v <= 1.) then
              point_within_triangle(i,j)=tri
              point_within_how_many_triangles=point_within_how_many_triangles+1
            endif
!            write(41,*)dot00,dot01,dot02,dot11,dot12
!            write(41,'(a3,3i5)')'PIN',i,j,point_within_triangle(i,j)
33333       continue

          enddo

          if(mod(j,100)==0 .and. mod(i,100)==0) then
            write(79,*)i,j,point_within_triangle(i,j),point_within_how_many_triangles
          endif
          if(point_within_how_many_triangles > 1) then
            write(79,*)i,j,point_within_how_many_triangles
            write(79,*)'>>>>>>>>> SOMETHING IS WRONG! CHECK PROGRAM! <<<<<<<<<'
            !stop
          endif
 
        enddo
      enddo	
	  
!----
!---- Loop over all grid points  ----
!----
     

      do i=1,1
        do j=1,1
!          if (buffer(i,j) < 0)cycle
          if (point_within_triangle(i,j) > 0) then
             tri=point_within_triangle(i,j)
             x1=node_xy(1,triangle_nodes_index(1,tri))
             y1=node_xy(2,triangle_nodes_index(1,tri))
             x2=node_xy(1,triangle_nodes_index(2,tri))
             y2=node_xy(2,triangle_nodes_index(2,tri))
             x3=node_xy(1,triangle_nodes_index(3,tri))
             y3=node_xy(2,triangle_nodes_index(3,tri))
             rr1=rr(triangle_nodes_index(1,tri))
             rr2=rr(triangle_nodes_index(2,tri))
             rr3=rr(triangle_nodes_index(3,tri))
             !---- Compute gradients and interpolated value at grid point using triangle ----		 
             Denom=(x3-x1)*(y2-y1)-(x2-x1)*(y3-y1)
             d_rr_over_d_x=((y2-y1)*(rr3-rr1)-(y3-y1)*(rr2-rr1))/Denom
             d_rr_over_d_y=((x3-x1)*(rr2-rr1)-(x2-x1)*(rr3-rr1))/Denom
             itmp=nint(d_rr_over_d_x*10000)
             d_rr_over_d_x=itmp/10000.
             itmp=nint(d_rr_over_d_y*10000)
             d_rr_over_d_y=itmp/10000.

             rr_grd(i,j)=rr1+(x(i)-x1)*d_rr_over_d_x+(y(j)-y1)*d_rr_over_d_y
             itmp=nint(rr_grd(i,j)*100.)
             rr_grd(i,j)=itmp/100.
             write(41,*)i,j,rr_grd(i,j),(x(i)-x1),d_rr_over_d_x,(y(j)-y1),d_rr_over_d_y
             h1=z_amsl(triangle_nodes_index(1,tri))
             h2=z_amsl(triangle_nodes_index(2,tri))
             h3=z_amsl(triangle_nodes_index(3,tri))		 
             !---- Compute gradients and interpolated value at grid point using triangle ----		 			
             Denom=(x3-x1)*(y2-y1)-(x2-x1)*(y3-y1)
             d_z_over_d_x=((y2-y1)*(h3-h1)-(y3-y1)*(h2-h1))/Denom
             d_z_over_d_y=((x3-x1)*(h2-h1)-(x2-x1)*(h3-h1))/Denom
             z_amsl_grd(i,j)=h1+(x(i)-x1)*d_z_over_d_x+(y(j)-y1)*d_z_over_d_y
             itmp=nint(z_amsl_grd(i,j)*100)
             z_amsl_grd(i,j)=itmp/100.
          else
             rr_grd(i,j)=-999	
             z_amsl_grd(i,j)=-999
          endif

       enddo
    enddo
!    print*,'Point in triangle complete'
    
	  
      open(1,file='rr_x_grd.txt')
      write(1,'(i8,4f7.1)')stnrx,rrx,rr_grd(1,1),dem1km(1,1),z_amsl_grd(1,1)
      close(1)

	  
	  
	  
	  
!---- Use nearest neighbour for points of Norway that are not covered by any triangles ----
      write(79,*)'Fill points of Norway that are not covered by any triangles with nearest neighbour:' 
      do i=1,1
        do j=1,1
          rr_grd_new(i,j)=-999			
          z_amsl_grd_new(i,j)=-999				
!          if (buffer(i,j)<0) cycle     

          if(rr_grd(i,j)==-999 .and. dem1km(i,j)>=0.) then		  
            distance_to_nearest_neighbour=10.e6
            do ii=max(i-100,1),min(i+100,1985)       !  --> Find closest neighbour from triangulated precipitation field (by going 100 km to west, east, north and south)
            do jj=max(j-100,1),min(j+100,2095)	
              if(rr_grd(ii,jj) < 0.) cycle         !  --> If point has no triangulated precipitation, cycle...
              if((ii-i)**2+(jj-j)**2 > 7225) cycle !  --> If point more than 100 km away, cycle... (Skip square root to save some time.)
              distance=(ii-i)**2+(jj-j)**2         !  --> Calculate distance to point (Skip square root to save some time.)
              if(distance < distance_to_nearest_neighbour) then
                distance_to_nearest_neighbour = distance
                rr_grd_new(i,j) = rr_grd(ii,jj)			
                z_amsl_grd_new(i,j) = z_amsl_grd(ii,jj)				
              endif			  
            enddo		
            enddo
          endif		  
        enddo
      enddo
      n=0	  
      do i=1,1
        do j=1,1
           if (dem1km(i,j)<-99) cycle
           if(rr_grd_new(i,j)>=0.) then
              n=n+1
              rr_grd(i,j)=rr_grd_new(i,j)
              z_amsl_grd(i,j)=z_amsl_grd_new(i,j)
           endif
 
          if (rr_grd(i,j)>=0) write(42,'(2i10,f10.2,f8.1)')(i*1000)+3863171,(j*1000)+3483615,rr_grd(i,j),z_amsl_grd_new(i,j)
        enddo
     enddo
      write(79,'(i10,a)')n,' points were outside triangles. These are filled with values from nearest triangle!'
	  

      do i=1,1
        do j=1,1     

          rr1km(i,j) = -999	   		
		
          z_grd = z_amsl_grd(i,j)
          rr2_grd = rr_grd(i,j)

          if(z_amsl_grd(i,j) < 0) cycle ! --> "Cycle" means jump to next execution of loop ("do" command) by increasing i or j.
		  
          diffz = dem1km(i,j) - z_grd
		  
          if (rr2_grd > .01) then
            if (z_grd < 1000.) then
              if (dem1km(i,j) < 1000.)  then
                rr1km(i,j) = rr2_grd + (rr2_grd * zfac * diffz)
              else
                rr1km(i,j) = rr2_grd + (rr2_grd * zfac * (1000. - z_grd)) + (rr2_grd * (zfac/2.) * (dem1km(i,j) - 1000.))
              endif				
            else
              if (dem1km(i,j) > 1000.) then
                rr1km(i,j) = rr2_grd + (rr2_grd * (zfac/2.) * diffz)
              else
                rr1km(i,j) = rr2_grd + (rr2_grd * (zfac/2.) * (1000. - z_grd)) + (rr2_grd * zfac * (dem1km(i,j) - 1000.))
              endif				
            endif			  
          else
            rr1km(i,j) = 0.0
          endif

          if ( rr1km(i,j) < 0.0) rr1km(i,j) = 0.0		  

        enddo
      enddo	

	  
!---- Write final gridded and terrain-corrected precipitation file ----
	  
      open(1,file='rr1km_out_X.txt')
      write(1,'(i8,5f7.1)')stnrx,rrx,rr_grd(1,1),dem1km(1,1),z_amsl_grd(1,1),rr1km(1,1)
      close(1)



	  
	  
      end program point_in_triangle
      

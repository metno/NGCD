!C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!C                                                                      %
!C Copyright (C) 1992 Stanford Center for Reservoir Forecasting.  All   %
!C rights reserved.  Distributed with: C.V. Deutsch and A.G. Journel.   %
!C ``GSLIB: Geostatistical Software Library and User's Guide,'' Oxford  %
!C University Press, New York, 1992.                                    %
!C                                                                      %
!C The programs in GSLIB are distributed in the hope that they will be  %
!C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
!C responsibility to anyone for the consequences of using them or for   %
!C whether they serve any particular purpose or work at all, unless he  %
!C says so in writing.  Everyone is granted permission to copy, modify  %
!C and redistribute the programs in GSLIB, but only under the condition %
!C that this notice and the above copyright notice remain intact.       %
!C                                                                      %
!C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!c-----------------------------------------------------------------------
!c
!c             Ordinary Kriging of a 2-D Rectangular Grid
!c             ******************************************
!c
!c The following Parameters control static dimensioning within okb2d:
!c
!c   MAXX      maximum nodes in X
!c   MAXY      maximum nodes in Y
!c   MAXDAT    maximum number of data points
!c   MAXSAM    maximum number of data points to use in one kriging system
!c   MAXDIS    maximum number of discretization points per block
!c   MAXNST    maximum number of nested structures
!c
!c   MAXKRG    MAXSAM*MAXSAM - used for dimensioning
!c
!c   UNEST     Assigned to unestimated blocks
!c
!c-----------------------------------------------------------------------
!c
!c User Adjustable:
!c
      parameter(MAXX  =  1550, &
               MAXY   =  2020, &
               MAXDAT =  1200,  &! NB! Was 1000 in original file!
               MAXSAM =  120, &
               MAXDIS =   64,  &! originally 64
               MAXNST =    1)  ! originally 4
!c
!c Fixed
!c
      parameter(MAXKD=MAXSAM+1,MAXKRG=MAXKD*MAXKD,UNEST=-999., &
               EPSLON=0.0000001,VERSION=1.202)


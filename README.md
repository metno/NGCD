NGCD
====
The NGCD is an high-resolution observational gridded dataset for daily temperature and precipitation over Fennoscandia (Finland, Sweden and Norway). 

Direct data sources are: 
1. European Climate Assessment & Dataset project [ECA&D](ecad.eu)
1. [Norwegian Meteorological Institute Climate Database](frost.met.no)

Note that most of the ECA&D data used outside Norway has been provided to ECA&D by The Swedish Meteorological and Hydrological Institute and The Finnish Meteorological Institute.

See the [wiki](https://github.com/metno/NGCD/wiki) for the software description.

Requirements
------------
Required libraries are:

* gdal


Required programs are:
* R
* Fortran compiler (gfortran or similar)
* C compiler (gcc or similar)

Required R-packages are:

* sp
* raster
* rgdal
* gstat
* foreign
* argparser
* ncdf4
* igraph
* tripack
* cluster

External Fortran libraries:
* GSLIB (www.gslib.com)
* Geompack (https://people.sc.fsu.edu/~jburkardt/f77_src/geompack/geompack.html)

Copyright and license
---------------------
Copyright (C) 2018 MET Norway. NGCD software is licensed under [GPLversion 2](https://github.com/metno/NGCD/blob/master/LICENSE) or (at your option) any later version.

External libraries such as GSLIB and Geompack might have different licences.

Contact
-------
E-mails: cristianl@met.no and oleet@met.no

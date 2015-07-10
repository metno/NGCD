ncout.spec.list <- alist(

    # precipitation sums
    # -----------------------------------------------------------------
#    list( pname="RprelimD",  
#          opts =  alist(
#                       fname.nc = "AUTO",
#                       t.ref = "190001010000",
#                       t.unit = "D",
#                       var.name = "RprelimD",
#                       var.longname = "daily precipitation sum (preliminary)",
#                       var.unit = "millimeter",
#                       var.mv = -999.99,
#                       reference = paste("",sep=""))
#        ),

    list( pname="elevation",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "elevation",
                       var.longname = "elevation",
                       var.unit = "meter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="IDI",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "IDI",
                       var.longname = "Integral Data Influence",
                       var.unit = "adimensional",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="IDIms",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "IDIms",
                       var.longname = "Integral Data Influence (multiscale modeling)",
                       var.unit = "adimensional",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    

    list( pname="PREC1d",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "precipitation_amount",
                       var.longname = "daily precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="PREC1h",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "precipitation_amount",
                       var.longname = "hourly precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="PREC1hRT",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "precipitation_amount",
                       var.longname = "hourly precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="PREC3h",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "precipitation_amount_3h",
                       var.longname = "3 hour precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="PREC3hRT",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "precipitation_amount_3h",
                       var.longname = "3 hour precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="RR_1",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "precipitation_amount",
                       var.longname = "hourly precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="RR_3",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "precipitation_amount_3h",
                       var.longname = "3 hour precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="RR",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "precipitation_amount",
                       var.longname = "daily precipitation sum",
                       var.unit = "millimeter",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="TA",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "temperature",
                       var.longname = "air temperature",
                       var.unit = "Celsius degrees",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="TEMP1h",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "H",
                       var.name = "temperature",
                       var.longname = "air temperature",
                       var.unit = "Celsius degrees",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),
    
    list( pname="TAM",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "mean_temperature",
                       var.longname = "daily mean temperature",
                       var.unit = "Celsius degrees",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="TEMP1d",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "mean_temperature",
                       var.longname = "daily mean temperature",
                       var.unit = "Celsius degrees",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="TAMRR",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "mean_temperature",
                       var.longname = "daily mean temperature",
                       var.unit = "Celsius degrees",
#                       var.mv = -999.99,
                       var.mv = "-999.99",
#                       reference = paste("",sep=""))
                       reference = "")
        ),

    list( pname="RhiresD",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "RhiresD",
                       var.longname = "daily precipitation sum",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RhiresM", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "RhiresM",
                       var.longname = "monthly precipitation sum",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RhiresY",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "RhiresY",
                       var.longname = "yearly precipitation sum",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RapdD",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "RapdD",
                       var.longname = "daily precipitation sum",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("Isotta, F.A. et al. 2013: ",
                            "The climate of daily precipitation in the Alps: development ",
                            "and analysis of a high-resolution grid dataset from pan-Alpine ",
                            "rain-gauge data. Int. J. Climatol., in press. (Please check ",
                            "publication status!)",sep=""))
        ),

    list( pname="Ral20D",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "Ral20D",
                       var.longname = "daily precipitation sum",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("Isotta, F.A. et al. 2013: ",
                            "The climate of daily precipitation in the Alps: development ",
                            "and analysis of a high-resolution grid dataset from pan-Alpine ",
                            "rain-gauge data. Int. J. Climatol., in press. (Please check ",
                            "publication status!)",sep=""))
        ),

    list( pname="RkedradD",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "RkedradD",
                       var.longname = "daily precipitation sum",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),


    # precipitation anomalies
    # -----------------------------------------------------------------
    list( pname="RanomM6190",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "RanomM6190",
                       var.longname = "monthly precipitation anomaly wrt 1961-1990",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RanomY6190", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "RanomY6190",
                       var.longname = "yearly precipitation anomaly wrt 1961-1990",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RanomM8110",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "RanomM8110",
                       var.longname = "monthly precipitation anomaly wrt 1981-2010",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RanomY8110", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "RanomY8110",
                       var.longname = "yearly precipitation anomaly wrt 1981-2010",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    # precipitation norm values
    # -----------------------------------------------------------------
    list( pname="RnormM6190",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "RnormM6190",
                       var.longname = "mean monthly precipitation 1961-1990",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RnormY6190", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "RnormY6190",
                       var.longname = "mean yearly precipitation 1961-1990",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RnormM8110",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "RnormM8110",
                       var.longname = "mean monthly precipitation 1981-2010",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="RnormY8110", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "RnormY8110",
                       var.longname = "mean yearly precipitation 1981-2010",
                       var.unit = "millimeter",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="R8110m6190M",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "R8110m6190M",
                       var.longname = "ratio of mean monthly precipitation between norm periods 1981-2010 / 1961-1990",
                       var.unit = "-",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="R8110m6190Y", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "R8110m6190Y",
                       var.longname = "ratio of mean yearly precipitation between norm periods 1981-2010 / 1961-1990",
                       var.unit = "-",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),


    # relative sunshine duration
    # ------------------------------------------------------------------
    list( pname="SrelD", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "SrelD",
                       var.longname = "daily sunshine duration relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),
                         
    list( pname="SrelM", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "SrelM",
                       var.longname = "monthly sunshine duration relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SrelY", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "SrelY",
                       var.longname = "yearly sunshine duration relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    # sunshine anomalies
    # ------------------------------------------------------------------
    list( pname="SanomM6190", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "SanomM6190",
                       var.longname = "anomaly of monthly sunshine duration wrt 1961-1990",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SanomY6190", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "SanomY6190",
                       var.longname = "anomaly of yearly sunshine duration wrt 1961-1990",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SanomM8110", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "SanomM8110",
                       var.longname = "anomaly of monthly sunshine duration wrt 1981-2010",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SanomY8110", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "SanomY8110",
                       var.longname = "anomaly of yearly sunshine duration wrt 1981-2010",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    # sunshine duration norm values
    # -----------------------------------------------------------------
    list( pname="SnormM6190",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "SnormM6190",
                       var.longname = "mean monthly sunshine duration 1961-1990 relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SnormY6190", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "SnormY6190",
                       var.longname = "mean yearly sunshine duration 1961-1990 relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SnormM8110",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "SnormM8110",
                       var.longname = "mean monthly sunshine duration 1981-2010 relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="SnormY8110", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "SnormY8110",
                       var.longname = "mean yearly sunshine duration 1981-2010 relative to max possible",
                       var.unit = "percent",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="S8110m6190M",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "S8110m6190M",
                       var.longname = "ratio of mean monthly sunshine duration between norm periods 1981-2010 / 1961-1990",
                       var.unit = "-",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),

    list( pname="S8110m6190Y", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "S8110m6190Y",
                       var.longname = "ratio of mean yearly sunshine duration between norm periods 1981-2010 / 1961-1990",
                       var.unit = "-",
                       var.mv = -999.99,
                       reference = paste("",sep=""))
        ),


    # absolute temperature
    # -----------------------------------------------------------------
    list( pname="TabsD",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TabsD",
                       var.longname = "daily mean temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786", sep=""))
        ),

    list( pname="TminD",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TminD",
                       var.longname = "daily minimum temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TmaxD",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TmaxD",
                       var.longname = "daily maximum temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TabsM",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "TabsM",
                       var.longname = "monthly mean temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TminM",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "TminM",
                       var.longname = "monthly mean minimum temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TmaxM",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "TmaxM",
                       var.longname = "monthly mean maximum temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TabsY",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "TabsY",
                       var.longname = "yearly mean temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TminY",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "TminY",
                       var.longname = "yearly mean minimum temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TmaxY",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "TmaxY",
                       var.longname = "yearly mean maximum temperature",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TmosxD",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TmosxD",
                       var.longname = "daily maximum temperature forecast",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TmosnD",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TmosnD",
                       var.longname = "daily minimum temperature forecast",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    # temperature norm values
    # -----------------------------------------------------------------
    list( pname="TnormD6190",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "D",
                       var.name = "TnormD6190",
                       var.longname = "calendar day mean temperature 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TnormM6190",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "TnormM6190",
                       var.longname = "mean monthly temperature 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TnormY6190", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "TnormY6190",
                       var.longname = "mean yearly temperature 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TnormD8110",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "D",
                       var.name = "TnormD8110",
                       var.longname = "calendar day mean temperature 1981-2010",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TnormM8110",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "TnormM8110",
                       var.longname = "mean monthly temperature 1981-2010",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TnormY8110", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "TnormY8110",
                       var.longname = "mean yearly temperature 1981-2010",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="T8110m6190M",  
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "M",
                       var.name = "T8110m6190M",
                       var.longname = "difference in mean monthly temperature between norm periods 1981-2010 / 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="T8110m6190Y", 
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "000001010000",
                       t.unit = "Y",
                       var.name = "T8110m6190Y",
                       var.longname = "difference in mean yearly temperature between norm periods 1981-2010 / 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    # temperature anomalies
    # -----------------------------------------------------------------
    list( pname="TanomD6190",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TanomD6190",
                       var.longname = "daily temperature anomaly wrt 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TanomM6190",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "TanomM6190",
                       var.longname = "monthly temperature anomaly wrt 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TanomY6190",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "TanomY6190",
                       var.longname = "yearly temperature anomaly wrt 1961-1990",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TanomD8110",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "D",
                       var.name = "TanomD8110",
                       var.longname = "daily temperature anomaly wrt 1981-2010",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TanomM8110",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "M",
                       var.name = "TanomM8110",
                       var.longname = "monthly temperature anomaly wrt 1981-2010",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        ),

    list( pname="TanomY8110",   
          opts =  alist(
                       fname.nc = "AUTO",
                       t.ref = "190001010000",
                       t.unit = "Y",
                       var.name = "TanomY8110",
                       var.longname = "yearly temperature anomaly wrt 1981-2010",
                       var.unit = "degree",
                       var.mv = -999.99,
                       reference = paste("Frei C., 2013: Interpolation of temperatures ",
                                         "in a mountainous region using nonlinear ",
                                         "profiles and non-Euclidean distances. ",
                                         "Int. J. Climatol., DOI: 10.1002/joc.3786",sep=""))
        )

  )


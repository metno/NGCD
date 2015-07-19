#!/bin/bash
#===============================================================================
# <seNorge2_PREC3hRT.sh>
#
# DESCRIPTION:
# ===========
# Spatial Interpolation of 3-hourly accumulated precipitation (Real-Time mode).
#
# COMMAND LINE:
# =============
#  >seNorge2_PREC3hRT.sh -s yyyy.mm.dd  (date start)
#                        -e yyyy.mm.dd  (date end)
#                        -c config.file (configuration file)
#                        -p config.par  (parameter name in configuration file)
#                        -l log.directory (log directory for Bspat)
#===============================================================================
function trim()
{
    local var=$1;
    var="${var#"${var%%[![:space:]]*}"}";   # remove leading whitespace characters
    var="${var%"${var##*[![:space:]]}"}";   # remove trailing whitespace characters
    var="${var%%,}";   # remove trailing comma characters
    echo -n "$var";
}
# whereis R?
  R=/usr/bin/R
#----------------------------
# Read command line arguments
#----------------------------
  flag_s=0
  flag_e=0
  flag_c=0
  flag_p=0
  flag_l=0
  while getopts "s:e:c:p:l:" Option
  do
    case $Option in
    s ) flag_s=1
    DATESTART=$OPTARG
    ;;
    e ) flag_e=1
    DATEEND=$OPTARG
    ;;
    c ) flag_c=1
    CONFIG_FILE=$OPTARG
    ;;
    p ) flag_p=1
    CONFIG_PAR=$OPTARG
    ;;
    l ) flag_l=1
    LOGDIR=$OPTARG
    ;;
    * ) echo " not recognized Option ";;
    esac
  done
#------------------------------------------------------------------------------
  # time-related operations
  YYYYbeg=${DATESTART:0:4}
  MMbeg=${DATESTART:5:2}
  DDbeg=${DATESTART:8:2}
  HHbeg=01
  YYYYend=${DATEEND:0:4}
  MMend=${DATEEND:5:2}
  DDend=${DATEEND:8:2}
  HHend=23
  SECbeg=`date +%s -d "$YYYYbeg-$MMbeg-$DDbeg $HHbeg:00:00"`
  SECend=`date +%s -d "$YYYYend-$MMend-$DDend $HHend:00:00"`
#------------------------------------------------------------------------------
# checks
  if [ "$flag_s" -eq 0 ] || [ "$flag_e" -eq 0 ] || [ "$SECbeg" -gt "$SECend" ]
  then
    echo "error in dates"
    exit 1
  fi
  if [ "$flag_c" -eq 0 ] || [ "$flag_p" -eq 0 ]
  then
    echo "error in configuration file specification"
    exit 1
  fi
  if [ "$flag_l" -eq 0 ]
  then
    echo "error in log file specification"
    exit 1
  fi
#------------------------------------------------------------------------------
# read config file and setup path variable
  flag_var=0
  while read line; do
    name=$(trim "${line%%=*}");
    value=$(trim "${line#*=}");
    aux="\"$CONFIG_PAR\""
    aux1="main.path"
    if [ "$value" == $aux ]; then
      flag_var=1
    fi
    if [ "$flag_var" -eq 1 ] && [ "$name" == $aux1 ]; then
        MAINDIR=$value
        break
      fi
  done < $CONFIG_FILE
  MAINDIR=${MAINDIR:1:-1}
#------------------------------------------------------------------------------
# log
  echo "seNorge2_PREC3hRT.sh "`date +%Y-%m-%d" "%H:%M`" > elaborations from "$DATESTART" UTC to "$DATEEND" UTC"
  echo "configuration file: "$CONFIG_FILE" configuration parameter:"$CONFIG_PAR
  echo "main directory:"$MAINDIR
  vis="vis-m1"
  if [ "$CONFIG_PAR" == "$vis" ]; then
    export R_LIBS=/home/cristianl/programs/spatial_interpolation/lib/rpacks
    echo $R_LIBS
  fi
#------------------------------------------------------------------------------
# Variables
  Bspat=$MAINDIR/Bspat_PREC3hRT/Bspat_PREC3hRT.R
  BLACKL=$MAINDIR/etc/blacklists/seNorge2_PREC1h_blacklist.txt
  ERROBS=$MAINDIR/etc/suspect_observations/seNorge2_PREC1h_suspect_observations.txt
#--------------------------------------------------
#  Clean temporary directories, if needed 
#--------------------------------------------------
#  echo "seNorge2_PREC3hRT.sh "`date +%Y-%m-%d" "%H:%M`" > clean directories"
#  find $OUTDIR/ -mtime +1 -exec rm -vf {} \;
#  rm -v $TOFTP_DIR/*.bz2
#  rm -v $TOFTP_DIR/ok
#  rm -v $FROMFTP_DIR/*.bz2
#  rm -v $FROMFTP_DIR/ok
#-----------------------------------------------------
# Statistical Interpolation 
#-----------------------------------------------------
  SECcur=$SECbeg
  SECnxt=$(( SECcur+86400 ))
  while (( "$SECcur" <= "$SECend" )) 
  do
    DATEcur=`date --date="1970-01-01 $SECcur sec UTC" +%Y.%m.%d`
    DATEnxt=`date --date="1970-01-01 $SECnxt sec UTC" +%Y.%m.%d`
    echo "================================================================================"
    echo "$R --vanilla $DATEcur.01 $DATEcur.03 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.01_$DATEcur.03.log 2>&1"
    $R --vanilla $DATEcur.01 $DATEcur.03 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.01_$DATEcur.03.log 2>&1
    echo "$R --vanilla $DATEcur.04 $DATEcur.06 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.04_$DATEcur.06.log 2>&1"
    $R --vanilla $DATEcur.04 $DATEcur.06 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.04_$DATEcur.06.log 2>&1
    echo "$R --vanilla $DATEcur.07 $DATEcur.09 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.07_$DATEcur.09.log 2>&1"
    $R --vanilla $DATEcur.07 $DATEcur.09 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.07_$DATEcur.09.log 2>&1
    echo "$R --vanilla $DATEcur.10 $DATEcur.12 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.10_$DATEcur.12.log 2>&1"
    $R --vanilla $DATEcur.10 $DATEcur.12 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.10_$DATEcur.12.log 2>&1
    echo "$R --vanilla $DATEcur.13 $DATEcur.15 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.13_$DATEcur.15.log 2>&1"
    $R --vanilla $DATEcur.13 $DATEcur.15 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.13_$DATEcur.15.log 2>&1
    echo "$R --vanilla $DATEcur.16 $DATEcur.18 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.16_$DATEcur.18.log 2>&1"
    $R --vanilla $DATEcur.16 $DATEcur.18 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.16_$DATEcur.18.log 2>&1
    echo "$R --vanilla $DATEcur.19 $DATEcur.21 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.19_$DATEcur.21.log 2>&1"
    $R --vanilla $DATEcur.19 $DATEcur.21 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.19_$DATEcur.21.log 2>&1
    echo "$R --vanilla $DATEcur.22 $DATEnxt.00 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.22_$DATEnxt.00.log 2>&1"
    $R --vanilla $DATEcur.22 $DATEnxt.00 $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC3hRT_$DATEcur.22_$DATEnxt.00.log 2>&1
    SECcur=$(( SECcur+86400 ))
    SECnxt=$(( SECcur+86400 ))
  done
  echo "seNorge2_PREC3hRT.sh "`date +%Y-%m-%d" "%H:%M`" > end main"
#--------------------------
# Exit 
#--------------------------
  echo "seNorge2_PREC3hRT.sh "`date +%Y-%m-%d" "%H:%M`" > success!"
  exit 0

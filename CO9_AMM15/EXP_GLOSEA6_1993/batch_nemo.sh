
iy1=1996 # initial year
iy2=1996 # end year
nit0=1215001 # namelist iteration start
ms=1 # months start
im2=12

tpd=1080

SRUN_CMD=$1

## year and month are the formatted version of running month and year to be used for script launching
## y and m are the numeric one to be used in counting
if [ $im2 -ge 13 ]; then im2=12; fi   # if out of step at end of year just run to end of year

# loop for years
yy=$iy1
while [  $yy -le $iy2 ]; do

   # loop for years
   year=${yy}
   
   #loop for months
   if [ $yy -eq $iy1 ]; then mm=$ms; else mm=1; fi
   while [ $mm -le $im2 ]; do
      
      if [ $mm -lt 10 ] ; then
         month=0$mm
      else
         month=$mm
      fi
      
      date
      
      echo "Running year $year, month $month ..."
      #get run length in iterations (nit)
      case $month in
              04|06|09|11) nit=$((30*tpd)) ;;
              #february considering leap years:
              02) if [ `expr $year % 4` -ne 0 -o `expr $year % 100` -eq 0 -a `expr $year % 400` -ne 0 ]; then nit=$((28*tpd)); else nit=$((29*tpd)); fi ;;
              #february not considering leap years:
              #02) nit=$((28*tpd)) ;;
              *) nit=$((31*tpd)) ;;
      esac
      
      nitend=`expr $nit0 + $nit`
      nitend=`expr $nitend - 1`
      
      if [ $nit0 -eq 1 ] ; then
         rstart=.false.  # start from rest
         rstctl=0
         inits=.true.    # start from rest
         nitrst=00000000
      else
         rstart=.true.   # start from restart
         rstctl=2
         inits=.false.   # start from restart
         nitrst=$year$month'01'
      fi
      
      cat namelist.template \
          | sed "s,__IT000__,$nit0,g" \
          | sed "s,__ITEND__,$nitend,g" \
          | sed "s,__RSTART__,$rstart,g" \
          | sed "s,__INITS__,$inits,g" \
          | sed "s,__RSTNIT__,$nitrst,g" \
          | sed "s,__RSTCTL__,$rstctl,g" \
          > namelist_cfg

      # move previous restarts
      mv RESTART_OUT/CHAMFER_lat_north_${nitrst}_restart* RESTART_IN/
      
      # this is the RUN bit ------>>>>>
      
      echo $SRUN_CMD
      eval $SRUN_CMD
      
      ## end of run bit <-----
      date
      ##
      ### move output files
      OUTDIR=$(pwd)/$year/$month
      mkdir -p $OUTDIR
      mkdir -p ${OUTDIR}/RESTARTS
      mv ocean.output $OUTDIR/ocean.output
      mv namelist_cfg $OUTDIR/namelist_cfg
      mv *25hour*  $OUTDIR/.
      mv CHAMFER_north_lat_1d_*  $OUTDIR/.
      mv mon*  $OUTDIR/.
      mv shelf*  $OUTDIR/.
      mv run.stat $OUTDIR/.
      mv RESTART_IN/CHAMFER_lat_north_${nitrst}_restart* RESTARTS/.
      
      current_stp=`sed -n 1,1p time.step`
      if [ ! $current_stp -eq $nitend ]
      then
         exit
      fi
      
      nit0=`expr $nitend + 1`
      
      mm=$((mm+1))
   done # 1..12 months
   yy=$((yy+1))
done # year and year+1

#if [ $y1 -lt 2006 ]; then
#echo qsub -v I=$I,J=$J,IJ=$IJ,m=$m1,y=$y1,nit0=$nit0  -o $RUNDIR/AMM7-$IJ-$y1-$month1 -N AMM7$y1$month1 monthlyrun.pbs
#qsub -v I=$I,J=$J,IJ=$IJ,m=$m1,y=$y1,nit0=$nit0 -o $RUNDIR/NE-AMM7-$IJ-$y1-$month1 -N NE$y1$month1  monthlyrun.pbs
#fi

exit



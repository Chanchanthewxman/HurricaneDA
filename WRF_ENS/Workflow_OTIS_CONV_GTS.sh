#!/bin/bash 


source ~/.bashrc
module restore intel  # Intel Compiler Version 19

#!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!!

#load configuration files, functions, parameters
export CONFIG_FILE=/expanse/lustre/projects/pen116/cpruett/WRF_ENS/CONV_GTS.OTIS
#!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!!

source "$CONFIG_FILE"
source "$SCRIPT_DIR"/util.sh

if [[ ! -d "$ENS_DIR" ]]; then mkdir -p "$ENS_DIR"; fi
cd "$ENS_DIR" || exit

####total_ntasks####
if [ "$JOB_SUBMIT_MODE" == 1 ]; then
  if [[ $HOSTTYPE == "stampede" || $HOSTTYPE == "Expanse" ]]; then
    export total_ntasks=$SLURM_NTASKS
  fi
  if [[ $HOSTTYPE == "jet" ]]; then
    export total_ntasks=$PBS_NP
  fi
else
  export total_ntasks=9999999
fi

#################################
#date

export DATE=$DATE_CYCLE_END
export PREVDATE=$DATE_CYCLE_END
export NEXTDATE=$DATE_CYCLE_END

while [[ $NEXTDATE == $DATE_CYCLE_END ]]; do

  export OWMAX=$(min ${OBS_WIN_MAX[@]})
  export OWMIN=$(min ${OBS_WIN_MIN[@]})
  export DT=$TIME_STEP
  export MPS=$(min ${MINUTE_PER_SLOT[@]})
  export FCSTM=$(min ${FORECAST_MINUTES[@]})
  export CP=$(diff_time "$DATE" "$DATE_END")

  # -----------------------------------------------------------------
  # calculate start_date and run_minutes, used by namelist_wrf.sh to &
  # generate correct time in namelist.input
  # -----------------------------------------------------------------
  if $RUN_4DVAR; then
    export start_date_cycle=$DATE
    export run_minutes_cycle=$(echo "$CP"+"$OWMAX" |bc)
    if [[ $DATE -ge $DATE_CYCLE_START ]]; then
      export start_date_cycle=$(advance_time "$start_date_cycle" "$OWMIN")  
      export run_minutes_cycle=$(echo "$run_minutes"+"$OWMIN" |bc)  
    fi
  else
    export start_date_cycle=$DATE
    export run_minutes_cycle=$CP
  fi
  if $RUN_DETERMINISTIC; then
    export run_minutes_forecast=$(diff_time "$DATE" "$DATE_END")
  else
    export run_minutes_forecast=$(max "$CP" "$FCSTM")
  fi
 
  # ------------------
  # LBDATE
  # ------------------
  #export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$LBC_INTERVAL" |bc`
  export minute_off=$(echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$BC_INTERVAL" |bc)
  if [[ $DATE == $(advance_time "$start_date_cycle" -"$minute_off") ]]; then
    export LBDATE=$(advance_time "$start_date_cycle" -"$minute_off")
    #export LBDATE=$DATE_START
  fi

  # ------------------
  # FCSTDATE
  # ------------------
  export minute_off=$(echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$FCST_INTERVAL" |bc)  
  export FCSTDATE=$(advance_time "$start_date_cycle" -"$minute_off")

  export NEXTDATE=$(advance_time "$DATE" "$CP")
  echo "----------------------------------------------------------------------"
  echo "FREE ENSEMBLE FORECAST: $(wrf_time_string "$DATE") => $(wrf_time_string "$NEXTDATE")"

  # ------------------
  # run components
  # ------------------
  "$ENS_SCRIPT_DIR"/wrf_ens.sh

  wait
  date

  # ------------------
  # check errors  
  # ------------------
  if [[ $(cat stat) == "error" ]]; then
    echo CYCLING STOP DUE TO FAILED COMPONENT
    exit 1
  fi

  # ------------------
  #advance to next cycle
  # ------------------
  export PREVDATE=$DATE
  export DATE=$NEXTDATE

done
echo CYCLING COMPLETE
echo bottom "$MODULEPATH_ROOT"

#!/bin/bash --login
. "$CONFIG_FILE"

rundir="$ENS_DIR"
if [[ ! -d $rundir ]]; then mkdir -p "$rundir"; echo running > "$rundir"/stat; fi

cd "$rundir" || exit
if [[ $(cat stat) == "complete" ]]; then exit; fi

#Setup for wrf run
echo "  Running Long-Range WRF ensemble..."

tid=0  #does not start from 0, because the wrf forecast runs with ens at the same time.
nt=$((total_ntasks/$wrf_ntasks))

if $RECENTER; then NUM_ENS=$(expr "$NUM_ENS" + 1); fi

for NE in $(seq 1 "$NUM_ENS"); do
#for NE in `seq 1 1`; do
  id=$(expr "$NE" + 1000 |cut -c2-)
  if [[ ! -d $id ]]; then mkdir "$id"; fi
  touch "$id"/rsl.error.0000
  if [[ $(tail -n1 "$id"/rsl.error.0000 |grep SUCCESS) ]]; then continue; fi

  cd "$id" || exit
  #ln -fs $WRF_DIR/run/* .
  ln -fs "$WRF_PRESET_DIR"/run/* .
  rm -f namelist.*

  for n in $(seq 1 "$MAX_DOM"); do
    dm=d$(expr "$n" + 100 |cut -c2-)
    #if [[ $NE -eq $NUM_ENS ]]; then
    #  ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm} wrfinput_$dm
    #  ln -fs $WORK_DIR/fc/$DATE/wrfbdy_d01 wrfbdy_d01
    #else
    ln -fs "$WORK_DIR"/fc/"$DATE"/wrfinput_"${dm}"_"$id" wrfinput_"$dm"
    ln -fs "$WORK_DIR"/fc/"$DATE"/wrfbdy_d01_"$id" wrfbdy_d01
    #fi
    #ncl $SCRIPT_DIR/util_change_nc_time.ncl 'ncfile="wrfinput_d01"' 'time="'`wrf_time_string $DATE`'"'
  done


  export start_date=$start_date_cycle
  export run_minutes=$run_minutes_cycle
  export inputout_interval=$run_minutes
  export inputout_begin=0
  export inputout_end=$run_minutes
  export GET_PHYS_FROM_FILE=false

  export wrf_for=forecast
  "$SCRIPT_DIR"/namelist_wrf_realtime.sh wrfw "$RUN_DOMAIN" > namelist.input
  
  echo Test1 
  
cat > run_wrf_ens.sh << EOF
#!/bin/bash -x
#SBATCH -A pen116
#SBATCH -J wrfens
#SBATCH -p compute
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH -t 8:00:00
#SBATCH --mem=249208M
#SBATCH -o run_wrf_ens.batch
#SBATCH --mail-type=FAIL # only email on job failure
#SBATCH --mail-user=cmp7089@psu.edu

ibrun ./wrf.exe >& wrf.log
	
EOF
  
  echo Test2
  sbatch run_wrf_ens.sh &> job_submit.log
  cd ..
done

echo complete > stat

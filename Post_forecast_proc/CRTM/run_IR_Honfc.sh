#!/bin/bash
# Script to forward calculate model state to IR Tb

#SBATCH -A pen116
#SBATCH -J toIR
#SBATCH -p compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH -e IR.e%j       # Name of stderr error file
#SBATCH -o IR.o%j       # Name of stdout output file
##SBATCH --mail-user=cmp7089@psu.edu
##SBATCH --mail-type=all    # Send email at begin and end of job

module restore intel
source util.sh

# Parent paths
fc_dir=/expanse/lustre/scratch/cpruett/temp_project/BERYL/NoDA/WRF_FREE # the location of forecasts; Usuallly under scratch directory
toHx_dir=/expanse/lustre/scratch/cpruett/temp_project/BERYL/CRTM_OUT/NoDA # the location of model-equivalent IR TBs

############ User control parameters
max_num_of_crtm=2				# Max number of CRTM.exe to run concurrently 
								# (make same as # of nodes requested)
cores_per_crtm=128				# Number of cores given to each crtm.exe 
								# (make same as # of cores per node)
date_st=202406260000			# Start date  
date_ed=202407091200			# End date
time_int=60						# Time interval btwn cycles in minutes
dom=1                           # Domain you are running it on 

##### Initialize counting variable to keep track of number of active CRTM.exe 
num_crtm=0

# Initialize date variable
DAtime=$date_st

# Iterate thru dates
while [[ $DAtime -le $date_ed ]]; do

  # figure out the time for this ensemble
  year=${DAtime:0:4}
  month=${DAtime:4:2}
  day=${DAtime:6:2} 
  hour=${DAtime:8:2}
  minute=${DAtime:10:2}  

  wrffile=${fc_dir}/wrfout_d0${dom}_${year}-${month}-${day}_${hour}:${minute}:00
  outfile=${toHx_dir}/GOES16_Ch7_wrfout_d0${dom}_${year}-${month}-${day}_${hour}:${minute}:00.bin

  echo $wrffile $outfile
      
  # Run the CRTM.exe

  ibrun -n $cores_per_crtm -o $(($num_crtm*$cores_per_crtm)) XbtoIR_crtm.exe $wrffile $outfile >& log_XbtoIR &


  # Increment number of active CRTM programs
  num_crtm=$(($num_crtm + 1))

  # If already has max_num_of_crtm CRTM programs active, wait for all processes to clear
  if [[ $num_crtm -ge $max_num_of_crtm ]]; then
    wait
    num_crtm=0
  fi

  # Increment date
  DAtime=`advance_time $DAtime $time_int`

done # End looping over dates

wait # Trick: the background process will leave the compute node once the process goes over the shell script. If at that moment the crtm processes are not finished yet, these crtm processes will crash.












#!/bin/bash
# Script to download the GOES data from AWS S3.
# Author: Zhu (Judy) Yao. Nov, 2024
# note: only tested for goes16 satellite.


#--------- User Defined --------------------

# Specify rclone locaiton
export rclone='/home/cpruett/rclone/rclone-v1.67.0-linux-amd64/rclone'
# Specify the destination directory
DownloadTo=/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/IR_DATA
# AMS S3 bucket name
Satellite=noaa-goes16
# Year
Year=2024
# Product name
Product=ABI-L2-CMIPF
# Channel name (e.g.,C01,C16)
CH=C08
# Specify date range (YYYYMMDD format)
StartDate=20240626
EndDate=20240709
# -------------------------------------------

# Convert date range to Day of Year (DOY)
StartDOY=$(date -d ${StartDate} +%j)
EndDOY=$(date -d ${EndDate} +%j)

# List available dates on AWS
Dates=$(${rclone} lsf --dirs-only remote:${Satellite}/${Product}/${Year})
DatesArray=($Dates)

# Filter dates by DOY
FilteredDates=()
for date in ${DatesArray[@]}; do
  CleanDate=${date%/}
  if [[ $CleanDate -ge $StartDOY && $CleanDate -le $EndDOY ]]; then
    FilteredDates+=($CleanDate)
  fi
done

# Loop through filtered dates
for time in ${FilteredDates[@]}; do
  # List files that have channel number as the keyword under that date
  files=$(${rclone} lsf remote:${Satellite}/${Product}/${Year}/${time}/ --recursive | grep ${CH})
  filesArray=($files)

  # Make directory on the server
  outdir=${DownloadTo}/${CH}/${time}/
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi

  # Download files if not already present
  for ff in "${filesArray[@]}"; do
    localfile="${outdir}/$(basename $ff)"
    if [[ -f "$localfile" ]]; then
      echo "File already exists: $localfile, skipping."
    else
      echo "Downloading: ${Satellite}/${Product}/${Year}/${time}/$ff"
      ${rclone} copy remote:${Satellite}/${Product}/${Year}/${time}/${ff} "$outdir" --progress
    fi
  
  done
done


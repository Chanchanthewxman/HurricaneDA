#!/bin/bash

# Set start and end datetimes (YYYYMMDDHH format)
START_DATE="2024062512"
END_DATE="2024070912"

DA_CYCLE_START_DATE="2024062600"
DA_CYCLE_END_DATE="2024062612"

# Time step in hours (e.g., 6 for GFS)
TIME_STEP=6
VALID_HOURS=("00" "06" "12" "18")

# Base URL
BASE_URL="https://data.rda.ucar.edu/d084001"

# Base Directory
BASE_DIR="/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF"

# Directory to store downloaded files
STORM_NAME_DIR="$BASE_DIR/BERYL"
mkdir -p "$STORM_NAME_DIR"

DATA_DIR="$STORM_NAME_DIR/Data"
mkdir -p "$DATA_DIR"

MODEL_DATA_DIR="$DATA_DIR/GFS_DATA"
mkdir -p "$MODEL_DATA_DIR"

OUTPUT_DIR="$MODEL_DATA_DIR/${START_DATE}_${END_DATE}"
mkdir -p "$OUTPUT_DIR"

# -------- Helper Functions -------- #
function is_valid_hour() {
    local hour="$1"
    for valid in "${VALID_HOURS[@]}"; do
        if [[ "$hour" == "$valid" ]]; then
            return 0
        fi
    done
    return 1
}

function download_forecast() {
    local init_datetime="$1"
    local end_datetime="$2"
    local init_epoch=$(date -d "${init_datetime:0:8} ${init_datetime:8:2}:00" +%s)
    local end_epoch=$(date -d "${end_datetime:0:8} ${end_datetime:8:2}:00" +%s)
    local total_hours=$(( (end_epoch - init_epoch) / 3600 ))

    echo "Processing initialization: $init_datetime -> $end_datetime (forecast range: ${total_hours}h)"

    for (( fhour=0; fhour<=total_hours; fhour+=TIME_STEP )); do
        f_str=$(printf "f%03d" "$fhour")
        yyyymmdd="${init_datetime:0:8}"
        hh="${init_datetime:8:2}"
        YEAR="${yyyymmdd:0:4}"

        FILE_NAME="gfs.0p25.${yyyymmdd}${hh}.${f_str}.grib2"
        URL="$BASE_URL/$YEAR/$yyyymmdd/$FILE_NAME"
        FILE_PATH="$OUTPUT_DIR/$FILE_NAME"

        if [[ -f "$FILE_PATH" ]]; then
            echo "Already exists: $FILE_NAME"
        else
            echo "Downloading: $URL"
            wget -q --show-progress -P "$OUTPUT_DIR" "$URL"
            if [[ $? -ne 0 ]]; then
                echo "Failed to download: $URL"
            fi
        fi
    done
}

# -------- Validity Checks -------- #
for dt in "$START_DATE" "$END_DATE" "$DA_CYCLE_START_DATE" "$DA_CYCLE_END_DATE"; do
    hour="${dt:8:2}"
    if ! is_valid_hour "$hour"; then
        echo "ERROR: Invalid hour in datetime: $dt. Must be one of: ${VALID_HOURS[*]}"
        exit 1
    fi
done

# -------- Step 1: Forecast from START_DATE to DA_CYCLE_START_DATE -------- #
download_forecast "$START_DATE" "$DA_CYCLE_START_DATE"

# -------- Step 2: Forecasts initialized at DA_CYCLE_START_DATE to END_DATE -------- #
da_start_epoch=$(date -d "${DA_CYCLE_START_DATE:0:8} ${DA_CYCLE_START_DATE:8:2}:00" +%s)
da_end_epoch=$(date -d "${DA_CYCLE_END_DATE:0:8} ${DA_CYCLE_END_DATE:8:2}:00" +%s)

# Step through DA cycle interval in TIME_STEP increments
current_epoch="$da_start_epoch"
while (( current_epoch <= da_end_epoch )); do
    current_datetime=$(date -d "@$current_epoch" +"%Y%m%d%H")
    current_hour="${current_datetime:8:2}"

    if is_valid_hour "$current_hour"; then
        download_forecast "$current_datetime" "$END_DATE"
    fi

    current_epoch=$((current_epoch + TIME_STEP * 3600))
done

echo "Download completed!"


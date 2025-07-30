import os
import requests
import pandas as pd
from datetime import datetime

def parse_hurdat_best_track(url, storm_name, year, output_dir):
    storm_name = storm_name.upper()
    response = requests.get(url)
    response.raise_for_status()
    lines = response.text.splitlines()

    storm_data = []
    for i, line in enumerate(lines):
        if ',' in line:
            parts = line.split(',')
            header_year = parts[0][4:8]
            header_name = parts[1].strip().upper()
            if header_name == storm_name and header_year == str(year):
                num_entries = int(parts[2])
                storm_data = lines[i + 1 : i + 1 + num_entries]
                break

    if not storm_data:
        raise ValueError(f"Storm '{storm_name}' in year {year} not found.")
    
    # Parse data lines into a DataFrame
    records = []
    for entry in storm_data:
        parts = [p.strip() for p in entry.split(',')]
        # Format datetime
        dt_raw = parts[0] + parts[1]  # yyyymmddHHMM
        dt = datetime.strptime(dt_raw, "%Y%m%d%H%M").strftime("%Y-%m-%d %H:%M:%S")

        lat = float(parts[4][:-1]) * (-1 if parts[4].endswith('S') else 1)
        lon = float(parts[5][:-1]) * (-1 if parts[5].endswith('W') else 1)
        wind_kt = int(parts[6])
        pressure_mb = int(parts[7])

        records.append({
            "Datetime": dt,
            "Center_Lon": lon,
            "Center_Lat": lat,
            "Max_Wind": round(wind_kt * 0.514444, 6),
            "Min_Press": pressure_mb
        })

    df = pd.DataFrame(records)

    # Save DataFrame to user-specified output directory
    output_file = f"{output_dir}/{storm_name}_{year}_best_track.csv"
    df.to_csv(output_file, index=False)
    print(f"Saved best track to: {output_file}")

    return df

# Example usage
if __name__ == "__main__":
    hurdat_url = "https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2024-040425.txt"
    storm_name = "BERYL"
    year = 2024
    output_dir = f"/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/{storm_name}/Data/HURDAT"

    os.makedirs(output_dir, exist_ok=True)

    df = parse_hurdat_best_track(hurdat_url, storm_name, year, output_dir)


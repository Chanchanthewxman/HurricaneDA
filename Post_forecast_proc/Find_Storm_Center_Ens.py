import numpy as np
import wrf
import scipy.ndimage
from netCDF4 import Dataset
import os
import glob
import csv
from datetime import datetime, timedelta

"""This script is based on the tracking algorithm described in Hartman et al. (2023)."""

def process_wrf_file(wrf_file, resolution):
    """Extracts SLP, wind speed, 700–850 hPa circulation, and 200–850 hPa thickness anomaly."""
    ncfile = Dataset(wrf_file)
    size = round(100/resolution)

    # --- Basic fields ---
    slp = wrf.getvar(ncfile, "slp")
    u10 = wrf.getvar(ncfile, "U10")
    v10 = wrf.getvar(ncfile, "V10")
    wspd = np.sqrt(u10**2 + v10**2)

    # --- Coordinates ---
    lats, lons = wrf.latlon_coords(slp)

    # Smooth SLP with 100 km boxcar (~size=5 for 3 km grid)
    slp_smooth = scipy.ndimage.uniform_filter(slp, size=size)

    # --- Relative vorticity and pressure levels ---
    avo = wrf.getvar(ncfile, "avo")  # Absolute vorticity
    pressure = wrf.getvar(ncfile, "pressure")

    # Coriolis parameter f = 2 * Omega * sin(lat)
    Omega = 7.2921e-5  # rad/s
    lat_radians = np.deg2rad(lats)
    f = 2 * Omega * np.sin(lat_radians)

    # Interpolate vorticity to 700 and 850 hPa
    avo_700 = wrf.interplevel(avo, pressure, 700)
    avo_850 = wrf.interplevel(avo, pressure, 850)

    # Compute relative vorticity = absolute vorticity - f
    rv_700 = avo_700 - f
    rv_850 = avo_850 - f

    # Apply 100 km × 100 km boxcar smoothing to each vorticity layer
    rv_700_smooth = scipy.ndimage.uniform_filter(rv_700, size=size)
    rv_850_smooth = scipy.ndimage.uniform_filter(rv_850, size=size)

    # Layer-averaged vorticity (approximation of circulation)
    circulation_smooth = (rv_700_smooth + rv_850_smooth) / 2

    # --- Geopotential height for thickness calculation ---
    g = 9.81  # m/s^2
    geopt = wrf.getvar(ncfile, "geopt")  # geopotential (m^2/s^2)
    geo_height = geopt / g  # geopotential height in meters
    geo_height_200 = wrf.interplevel(geo_height, pressure, 200)
    geo_height_850 = wrf.interplevel(geo_height, pressure, 850)
    thickness = geo_height_200 - geo_height_850

    # Subtract domain-mean thickness to get anomaly
    thickness_anom = thickness - np.nanmean(thickness)

    # Smooth thickness anomaly
    thickness_anom_smooth = scipy.ndimage.uniform_filter(thickness_anom, size=size)
    
    ncfile.close()

    return lons, lats, wspd, slp_smooth, circulation_smooth, thickness_anom_smooth

def distance(lats, lons, ref_lat, ref_lon):
    """
    Calculates great-circle distance in km between a reference point and one or more lat/lon points.
    
    Parameters:
        lats (float or np.ndarray): Latitude(s) of target point(s)
        lons (float or np.ndarray): Longitude(s) of target point(s)
        ref_lat (float): Reference latitude
        ref_lon (float): Reference longitude

    Returns:
        np.ndarray or float: Distance(s) in kilometers
    """
    return np.sqrt((lats - ref_lat)**2 + (lons - ref_lon)**2) * 111.0

def extrema(var, lats, lons, mode='minima'):
    """
    Finds the latitude and longitude of the min or max of a 2D variable.

    Parameters:
        var (np.ndarray): 2D variable field (e.g., SLP or thickness anomaly)
        lats (np.ndarray): 2D array of latitudes (same shape as var)
        lons (np.ndarray): 2D array of longitudes (same shape as var)
        mode (str): 'minima' or 'maxima'

    Returns:
        tuple: (lat, lon) of the min or max location
    """
    flat_var = var.values.flatten() if hasattr(var, 'values') else var.flatten()
    flat_lats = lats.values.flatten() if hasattr(lats, 'values') else lats.flatten()
    flat_lons = lons.values.flatten() if hasattr(lons, 'values') else lons.flatten()

    if mode == 'minima':
        idx = np.nanargmin(var)
    elif mode == 'maxima':
        idx = np.nanargmax(var)
    else:
        raise ValueError("mode must be 'minima' or 'maxima'")

    return flat_lats[idx], flat_lons[idx], flat_var[idx]

def centroid(lat1, lon1, lat2, lon2, lat3, lon3):
    """
    Computes the centroid (average lat/lon) of three coordinate pairs.

    Parameters:
        lat1, lon1, lat2, lon2, lat3, lon3 (float): Coordinates of the three points

    Returns:
        tuple: (centroid_lat, centroid_lon)
    """
    centroid_lat = (lat1 + lat2 + lat3) / 3.0
    centroid_lon = (lon1 + lon2 + lon3) / 3.0
    return centroid_lat, centroid_lon

def find_storm_center(wrf_file, resolution, best_lon, best_lat, max_dist, max_dist_centroid):
    """
    Finds the tropical cyclone (TC) center using SLP, vorticity at 850 hPa, and surface wind speed.
    """
    # Process wrf file
    lons, lats, wspd, slp_smooth, circulation_smooth, thickness_anom_smooth = process_wrf_file(wrf_file, resolution)

    # Compute distances from the best track estimate
    distances = distance(lats, lons, best_lat, best_lon)

    # Mask grid points beyond max_dist
    mask = distances <= max_dist
    masked_circulation = np.where(mask, circulation_smooth, np.nan)
    masked_thickness_anom = np.where(mask, thickness_anom_smooth, np.nan)
    masked_slp = np.where(mask, slp_smooth, np.nan)

    # Find maxima and minima
    circ_max_lat, circ_max_lon, _ = extrema(masked_circulation, lats, lons, 'maxima')
    thk_anom_max_lat, thk_anom_max_lon, _ = extrema(masked_thickness_anom, lats, lons, 'maxima')
    slp_min_lat, slp_min_lon, _ = extrema(masked_slp, lats, lons, 'minima')

    circ_thk_dist = distance(thk_anom_max_lat, thk_anom_max_lon, circ_max_lat, circ_max_lon)
    if circ_thk_dist <= max_dist:
        circ_slp_dist = distance(slp_min_lat, slp_min_lon, circ_max_lat, circ_max_lon)
        thk_slp_dist = distance(slp_min_lat, slp_min_lon, thk_anom_max_lat, thk_anom_max_lon)
        if (circ_slp_dist <= max_dist) and (thk_slp_dist <= max_dist):
            print(f"All three extrema are within {max_dist} km of one another!")
            center_lat, center_lon = centroid(circ_max_lat, circ_max_lon, thk_anom_max_lat, thk_anom_max_lon, slp_min_lat, slp_min_lon)
            distances_centroid = distance(lats, lons, center_lat, center_lon)
            mask_centroid = distances_centroid <= max_dist_centroid
            masked_slp_centroid = np.where(mask_centroid, slp_smooth, np.nan)
            _, _, min_slp = extrema(masked_slp_centroid, lats, lons, 'minima')
            masked_wspd_centroid = np.where(mask_centroid, wspd, np.nan)
            _, _, max_wind_speed = extrema(masked_wspd_centroid, lats, lons, 'maxima')
        else:
            print(f"Remasking slp for values within {max_dist} km of other extrema!")
            dist_to_circ = distance(lats, lons, circ_max_lat, circ_max_lon)
            dist_to_thk = distance(lats, lons, thk_anom_max_lat, thk_anom_max_lon)
            mask_for_slp = (dist_to_circ <= max_dist) & (dist_to_thk <= max_dist)
            masked_slp = np.where(mask_for_slp, slp_smooth, np.nan)
            slp_min_lat, slp_min_lon, _ = extrema(masked_slp, lats, lons, 'minima')
            center_lat, center_lon = centroid(circ_max_lat, circ_max_lon, thk_anom_max_lat, thk_anom_max_lon, slp_min_lat, slp_min_lon)
            distances_centroid = distance(lats, lons, center_lat, center_lon)
            mask_centroid = distances_centroid <= max_dist_centroid
            masked_slp_centroid = np.where(mask_centroid, slp_smooth, np.nan)
            _, _, min_slp = extrema(masked_slp_centroid, lats, lons, 'minima')
            masked_wspd_centroid = np.where(mask_centroid, wspd, np.nan)
            _, _, max_wind_speed = extrema(masked_wspd_centroid, lats, lons, 'maxima')
    else:
        print(f"Distance between circulation and thickness is {circ_thk_dist} km!")
        dist_to_circ = distance(lats, lons, circ_max_lat, circ_max_lon)
        mask_around_circ = dist_to_circ <= max_dist
        masked_thickness_anom = np.where(mask_around_circ, thickness_anom_smooth, np.nan)
        masked_slp = np.where(mask_around_circ, slp_smooth, np.nan)
        thk_anom_max_lat, thk_anom_max_lon, _ = extrema(masked_thickness_anom, lats, lons, 'maxima')
        slp_min_lat, slp_min_lon, _ = extrema(masked_slp, lats, lons, 'minima')
        center_lat, center_lon = centroid(circ_max_lat, circ_max_lon, thk_anom_max_lat, thk_anom_max_lon, slp_min_lat, slp_min_lon)
        distances_centroid = distance(lats, lons, center_lat, center_lon)
        mask_centroid = distances_centroid <= max_dist_centroid
        masked_slp_centroid = np.where(mask_centroid, slp_smooth, np.nan)
        _, _, min_slp = extrema(masked_slp_centroid, lats, lons, 'minima')
        masked_wspd_centroid = np.where(mask_centroid, wspd, np.nan)
        _, _, max_wind_speed = extrema(masked_wspd_centroid, lats, lons, 'maxima')

    return center_lon, center_lat, max_wind_speed, min_slp

def pad_member_id(member_num):
    return f"{member_num:03d}"

def calculate_expected_rows(start, end, step_hours):
    return ((end - start).total_seconds() // 3600) // step_hours + 1

if __name__ == "__main__":
    # ---------- Configuration -------------------------
    Storm = 'Otis'
    Exper_name = 'CONV_GTS_H23'
    resolution = 9
    wrf_dir = '/expanse/lustre/scratch/cpruett/temp_project/OTIS/WRF_ENS/CONV_GTS'
    start_time_str = '2023-10-16_06:00:00'
    end_time_str = '2023-10-21_06:00:00'
    timestep = 1 # 1 means no skip, 2 means skip every other hour
    best_lon = -89.0 # Ensure it is a float by using decimal
    best_lat = 9.0 # Ensure it is a float by using decimal
    max_dist, max_dist_centroid = 300, 100
    num_ensembles = 40 # Ensure integer
    track_dir = f'./{Storm}_tracks'
    exp_dir = f'{track_dir}/{Exper_name}_tracks'
    parent_csv_dir = f'{exp_dir}/{start_time_str}_to_{end_time_str}'
    # -------------------------------------------------------

    start_date = datetime.strptime(start_time_str, "%Y-%m-%d_%H:%M:%S")
    end_date = datetime.strptime(end_time_str, "%Y-%m-%d_%H:%M:%S")
    expected_rows = int(calculate_expected_rows(start_date, end_date, timestep))

    os.makedirs(track_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)
    os.makedirs(parent_csv_dir, exist_ok=True)

    for ens_num in range(1, num_ensembles + 1):
        member_id = pad_member_id(ens_num)
        member_dir = os.path.join(wrf_dir, member_id)
        if not os.path.isdir(member_dir):
            print(f"Warning: Directory not found for ensemble {member_id}: {member_dir}")
            continue

        csv_filename = os.path.join(parent_csv_dir, f"storm_center_{member_id}.csv")
        existing_dates = set()
        resume_date = start_date

        # Check if the CSV file already exists and determine the resume point
        
        mode = 'w'
        if os.path.exists(csv_filename):
            with open(csv_filename, mode='r') as file:
                reader = csv.reader(file)
                next(reader, None)  # skip header
                for row in reader:
                    if row:
                        try:
                            existing_dt = datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S')
                            existing_dates.add(existing_dt)
                        except:
                            continue

            if len(existing_dates) >= expected_rows:
                print(f"[{member_id}] CSV already complete. Skipping.")
                continue
            elif existing_dates:
                latest_dt = max(existing_dates)
                resume_date = latest_dt + timedelta(hours=timestep)
                print(f"[{member_id}] Resuming from {resume_date}")
                mode = 'a'
            else:
                print(f"[{member_id}] Empty or malformed CSV found. Restarting.")

        with open(csv_filename, mode=mode, newline='') as file:
            writer = csv.writer(file)
            if os.stat(csv_filename).st_size == 0:
                writer.writerow(["Datetime", "Center_Lon", "Center_Lat", "Max_Wind", "Min_Press"])

            current_date = resume_date 
            if resume_date == start_date:
                local_best_lon, local_best_lat = best_lon, best_lat
            else:
                # Get the last known center from the CSV
                last_row = sorted(existing_dates)[-1]  # latest datetime
                with open(csv_filename, mode='r') as file:
                    reader = csv.reader(file)
                    next(reader, None)  # skip header
                    found_coords = False
                    for row in reader:
                        if row[0] == last_row.strftime('%Y-%m-%d %H:%M:%S'):
                            local_best_lon = float(row[1])
                            local_best_lat = float(row[2])
                            found_coords = True
                            break
                    if not found_coords:
                        print(f"[{member_id}] Warning: Could not find coordinates for {last_row}. Using default starting point.")
                        local_best_lon, local_best_lat = best_lon, best_lat

            local_max_dist = max_dist if resume_date == start_date else 150
            
            while current_date <= end_date:
                wrf_file = os.path.join(member_dir, f"wrfout_d01_{current_date.strftime('%Y-%m-%d_%H:%M:%S')}")
                if os.path.exists(wrf_file):
                    center_lon, center_lat, max_wind_speed, min_slp = find_storm_center(
                            wrf_file, resolution, local_best_lon, local_best_lat,
                            local_max_dist, max_dist_centroid)
                    writer.writerow([
                        current_date.strftime('%Y-%m-%d %H:%M:%S'),
                        center_lon, center_lat, max_wind_speed, min_slp
                    ])
                    print(f"[{member_id}] File processed: {wrf_file}")
                    local_best_lon, local_best_lat = center_lon, center_lat
                    local_max_dist = 100  # shrink search radius after first frame
                else:
                    print(f"[{member_id}] File not found: {wrf_file}")

                current_date += timedelta(hours=timestep)


import numpy as np
import wrf
from scipy.ndimage import uniform_filter, uniform_filter1d
from netCDF4 import Dataset, num2date
import os
import glob
import csv
from datetime import datetime, timedelta

"""This script is based on the tracking algorithm described in Hartman et al. (2023)."""

def nanuniform_filter(data, size):
    """Apply a uniform filter that ignores NaNs."""
    data_nan = np.isnan(data)
    data_filled = np.where(data_nan, 0, data)

    kernel_sum = uniform_filter((~data_nan).astype(float), size=size, mode='constant', cval=0.0)
    data_sum = uniform_filter(data_filled, size=size, mode='constant', cval=0.0)

    with np.errstate(invalid='ignore', divide='ignore'):
        result = data_sum / kernel_sum
        result[kernel_sum == 0] = np.nan  # if kernel area was all-NaN, set result to NaN

    return result

# --- Debugging checks ---
def check_array(name, arr):
    print(f"\n--- {name} ---")
    print(f"Shape: {arr.shape}")
    print(f"Min: {np.nanmin(arr)}, Max: {np.nanmax(arr)}, Mean: {np.nanmean(arr)}")
    if np.all(np.isnan(arr)):
        print(f"Warning: {name} is entirely NaN!")
    elif arr.size == 0:
        print(f"Warning: {name} is empty!")

def sanity_check(variable, label, cmap):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.pcolormesh(variable, cmap=cmap)
    plt.colorbar(label=f'{label}')
    plt.title(f'{label}')
    plt.xlabel('X grid index')
    plt.ylabel('Y grid index')
    plt.savefig(f"./sanity_check_{label}.png", dpi=400)

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
    slp_smooth = nanuniform_filter(slp, size=size)

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
    rv_700_smooth = nanuniform_filter(rv_700, size=size)
    rv_850_smooth = nanuniform_filter(rv_850, size=size)

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
    thickness_anom_smooth = nanuniform_filter(thickness_anom, size=size)

    ncfile.close()
    
    return lons, lats, wspd, slp_smooth, circulation_smooth, thickness_anom_smooth

def process_era5_file(surface_file, upper_file, time):

    # Open ERA5 files
    surf = Dataset(surface_file)
    upper = Dataset(upper_file)

    # --- Time Index Matching ---
    time_var_s = surf.variables['valid_time']
    time_var_u = upper.variables['valid_time']

    times_s = num2date(time_var_s[:], units=time_var_s.units)
    times_u = num2date(time_var_u[:], units=time_var_u.units)

    t_idx_s = np.where(times_s == time)[0]
    t_idx_u = np.where(times_u == time)[0]

    if len(t_idx_s) == 0 or len(t_idx_u) == 0:
        raise ValueError(f"Time {time} not found in one or both files.")

    t_idx_s = t_idx_s[0]
    t_idx_u = t_idx_u[0]

    # Extract coordinates (assume same on both files)
    lats = surf.variables['latitude'][:]
    lons = surf.variables['longitude'][:]
    lons, lats = np.meshgrid(lons, lats)

    # Extract variables
    slp = surf.variables['msl'][t_idx_s] / 100  # Pa to hPa
    u10 = surf.variables['u10'][t_idx_s]
    v10 = surf.variables['v10'][t_idx_s]
    wspd = np.sqrt(u10**2 + v10**2)

    # Get pressure and geopotential levels
    levels = upper.variables['pressure_level'][:]  # should include 200, 700, 850
    z = upper.variables['z'][t_idx_u] / 9.81  # convert geopotential to height
    rv = upper.variables['vo'][t_idx_u]  # relative vorticity in s^-1

    # Extract fields at required levels
    level_idx = {}
    for lev in [200, 700, 850]:
        idxs = np.where(levels == lev)[0]
        if len(idxs) == 0:
            raise ValueError(f"Level {lev} not found in upper file")
        level_idx[lev] = idxs[0]

    z200 = z[level_idx[200], :, :]
    z850 = z[level_idx[850], :, :]
    rv700 = rv[level_idx[700], :, :]
    rv850 = rv[level_idx[850], :, :]

    circulation = (rv700 + rv850) / 2
    thickness = z200 - z850
    thickness_anom = thickness - np.nanmean(thickness)

    surf.close()
    upper.close()

    return lons, lats, wspd, slp, circulation, thickness_anom

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

def find_storm_center_WRF(wrf_file, resolution, best_lon, best_lat, max_dist, max_dist_centroid):
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

def find_storm_center_ERA5(surface_file, upper_file, time, best_lon, best_lat, max_dist, max_dist_centroid):
    """
    Finds the tropical cyclone (TC) center using SLP, vorticity at 850 hPa, and surface wind speed.
    """
    # Process wrf file
    lons, lats, wspd, slp_smooth, circulation_smooth, thickness_anom_smooth = process_era5_file(surface_file, upper_file, time)

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


if __name__ == "__main__":
    # ---------- Configuration -------------------------
    Storm = 'Beryl'
    Exper_name = 'NoDA'
    WRF = True
    resolution = 6
    wrf_dir = '/expanse/lustre/scratch/cpruett/temp_project/BERYL/NoDA/WRF_FREE/'
    surface_file = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/ERA5/era5_sfc_wind_press_2024.nc'
    upper_file = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/ERA5/era5_upper_geog_vort_2024.nc'
    start_time_str = '2024-06-26_00:00:00'
    end_time_str = '2024-07-09_12:00:00'
    timestep = 1 # 1 means no skip, 2 means skip every other hour
    best_lon = -27.83 #-89.0 # Ensure it is a float by using decimal
    best_lat = 10.58 #9.0 # Ensure it is a float by using decimal
    max_dist, max_dist_centroid = 300, 100
    track_dir = f'./{Storm}_tracks'
    exp_dir = f'{track_dir}/{Exper_name}_tracks'
    # -------------------------------------------------------
    
    start_date = datetime.strptime(start_time_str, "%Y-%m-%d_%H:%M:%S")
    end_date = datetime.strptime(end_time_str, "%Y-%m-%d_%H:%M:%S")
    csv_filename = f"{start_time_str}_to_{end_time_str}.csv"
    
    os.makedirs(track_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)
    
    csv_file_path = os.path.join(exp_dir, csv_filename)
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Datetime", "Center_Lon", "Center_Lat", "Max_Wind", "Min_Press"])
    
        current_date = start_date
        local_best_lon, local_best_lat = best_lon, best_lat  # start location for this member
        local_max_dist = max_dist
        if WRF:
            while current_date <= end_date:
                wrf_file = os.path.join(wrf_dir, f"wrfout_d01_{current_date.strftime('%Y-%m-%d_%H:%M:%S')}")
                if os.path.exists(wrf_file):
                    center_lon, center_lat, max_wind_speed, min_slp = find_storm_center_WRF( 
                            wrf_file, resolution, local_best_lon, local_best_lat, 
                            local_max_dist, max_dist_centroid)
                    writer.writerow([current_date.strftime('%Y-%m-%d %H:%M:%S'), center_lon, center_lat, max_wind_speed, min_slp])
                    print(f"File processed: {wrf_file}")
                    local_best_lon, local_best_lat = center_lon, center_lat
                    local_max_dist = 100  # shrink search radius after first frame
                else:
                    print(f"File not found: {wrf_file}")
            
                current_date += timedelta(hours=timestep)  # Adjust timestep as needed
        else:
            if (os.path.exists(surface_file) and os.path.exists(upper_file)):
                while current_date <= end_date:
                    center_lon, center_lat, max_wind_speed, min_slp = find_storm_center_ERA5(
                            surface_file, upper_file, current_date, local_best_lon, local_best_lat,
                            local_max_dist, max_dist_centroid)
                    writer.writerow([current_date.strftime('%Y-%m-%d %H:%M:%S'), center_lon, center_lat, max_wind_speed, min_slp])
                    print(f"ERA5 files processed for: {current_date}")
                    print(f"{center_lat},{center_lon}")
                    local_best_lon, local_best_lat = center_lon, center_lat
                    local_max_dist = 100  # shrink search radius after first frame
                    current_date += timedelta(hours=timestep)  # Adjust timestep as needed
            else:
                print(f"ERA5 files not found!")


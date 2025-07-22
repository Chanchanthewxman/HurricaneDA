import numpy as np
import wrf
import scipy.ndimage
from netCDF4 import Dataset
import os
import glob
import csv
from datetime import datetime, timedelta

def find_storm_center(lons, lats, slp, vorticity, wspd, max_dist, best_lon, best_lat, max_dist_pv, max_dist_press):
    """
    Finds the tropical cyclone (TC) center using SLP, vorticity at 850 hPa, and surface wind speed.
    """
    # Smooth SLP and vorticity fields
    slp_smooth = scipy.ndimage.gaussian_filter(slp, sigma=1)
    vorticity_smooth = scipy.ndimage.gaussian_filter(vorticity, sigma=1)
     
    # Compute distances from the best track estimate
    distances = np.sqrt((lats - best_lat) ** 2 + (lons - best_lon) ** 2) * 111  # Approx km
    
    # Mask grid points beyond max_dist
    mask = distances < max_dist
    masked_vort = np.where(mask, vorticity_smooth, np.nan)
    masked_wspd = np.where(mask, wspd, np.nan)
    
    # Find max vorticity point meeting wind speed condition
    pv_max_idx = np.nanargmax(masked_vort)
    pv_max_lat = lats.values.flatten()[pv_max_idx]
    pv_max_lon = lons.values.flatten()[pv_max_idx]
    
    # Find min SLP near max vorticity location
    distances_pv = np.sqrt((lats - pv_max_lat) ** 2 + (lons - pv_max_lon) ** 2) * 111
    mask_pv = distances_pv < max_dist_pv
    masked_slp = np.where(mask_pv, slp_smooth, np.nan)
    slp_min_idx = np.nanargmin(masked_slp)
    
    min_slp = slp.values.flatten()[slp_min_idx]
    center_lat = lats.values.flatten()[slp_min_idx]
    center_lon = lons.values.flatten()[slp_min_idx]

    # Find max wind speed near the low-pressure center
    distances_press = np.sqrt((lats - center_lat) ** 2 + (lons - center_lon) ** 2) * 111
    mask_press = distances_press < max_dist_press
    masked_wspd_press = np.where(mask_press, wspd, np.nan)
    max_wind_idx = np.nanargmax(masked_wspd_press)
    max_wind_speed = wspd.values.flatten()[max_wind_idx]
    
    return center_lon, center_lat, max_wind_speed, min_slp

def process_wrf_file(wrf_file):
    """Extracts SLP, vorticity (interpolated to 850 hPa), and wind speed from a WRF output file."""
    ncfile = Dataset(wrf_file)
    
    slp = wrf.getvar(ncfile, "slp")
    u10 = wrf.getvar(ncfile, "U10")
    v10 = wrf.getvar(ncfile, "V10")
    wspd = np.sqrt(u10**2 + v10**2)
    
    pressure = wrf.getvar(ncfile, "pressure")
    vorticity = wrf.getvar(ncfile, "avo")
    vorticity_850 = wrf.interplevel(vorticity, pressure, 850)
    
    lats, lons = wrf.latlon_coords(slp)
    
    return lons, lats, slp, vorticity_850, wspd

if __name__ == "__main__":
    # ---------- Configuration -------------------------
    Storm = 'Otis'
    Exper_name = 'CONV_GTS'
    wrf_dir = '/expanse/lustre/scratch/cpruett/temp_project/OTIS/WRF_DET/CONV_GTS'
    start_time_str = '2023-10-16_06:00:00'
    end_time_str = '2023-10-21_06:00:00'
    timestep = 1 # 1 means no skip, 2 means skip every other hour
    best_lon = -89.0 #-92.07 # Ensure it is a float by using decimal
    best_lat = 9.0 #8.628 # Ensure it is a float by using decimal
    max_dist, max_dist_pv, max_dist_press = 300, 150, 100
    # -------------------------------------------------------
    
    start_date = datetime.strptime(start_time_str, "%Y-%m-%d_%H:%M:%S")
    end_date = datetime.strptime(end_time_str, "%Y-%m-%d_%H:%M:%S")
    csv_filename = f"{Storm}_{Exper_name}_{start_time_str}_to_{end_time_str}.csv"
    
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Datetime", "Center_Lon", "Center_Lat", "Max_Wind", "Min_Press"])
    
        current_date = start_date
        local_best_lon, local_best_lat = best_lon, best_lat  # start location for this member
        local_max_dist = max_dist

        while current_date <= end_date:
            wrf_file = os.path.join(wrf_dir, f"wrfout_d01_{current_date.strftime('%Y-%m-%d_%H:%M:%S')}")
            if os.path.exists(wrf_file):
                lons, lats, slp, vorticity_850, wspd = process_wrf_file(wrf_file)
                center_lon, center_lat, max_wind_speed, min_slp = find_storm_center(
                        lons, lats, slp, vorticity_850, wspd, 
                        local_max_dist, local_best_lon, local_best_lat, 
                        max_dist_pv, max_dist_press)
                writer.writerow([current_date.strftime('%Y-%m-%d %H:%M:%S'), center_lon, center_lat, max_wind_speed, min_slp])
                print(f"File processed: {wrf_file}")
                local_best_lon, local_best_lat = center_lon, center_lat
                local_max_dist = 150  # shrink search radius after first frame
            else:
                print(f"File not found: {wrf_file}")
            
            current_date += timedelta(hours=timestep)  # Adjust timestep as needed


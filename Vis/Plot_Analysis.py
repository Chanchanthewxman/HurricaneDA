import os
import glob
import numpy as np
import netCDF4 as nc
import cmweather
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('agg')
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from datetime import datetime, timedelta
import time
import imageio.v2 as imageio  # Use v2 for compatibility
import zipfile
import pandas as pd
import Util_Vis
import wrf
from scipy.ndimage import uniform_filter, uniform_filter1d
import multiprocessing as mp

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

def plot_wrf_variable(wrf_file, resolution, varname, level, wind, Storm, Exper_name, output_path, timestamp):
    """
    Plots a 2D WRF variable.

    Parameters:
    - varname: str, the WRF variable name (e.g., "avo", "slp", "tk", etc.)
    - wrf_file: str, path to the WRF output file
    - level: int or float, pressure level in hPa (only needed for 3D vars like "avo")
    """
    ncfile = nc.Dataset(wrf_file)
    size = round(100/resolution)

    if varname == "slp":
        var = wrf.getvar(ncfile,varname)
        lats, lons = wrf.latlon_coords(var)
        min_threshold = 995
        max_threshold = 1030
        cmap = "viridis"
        label = "Sea Level Pressure (hPa)"
    elif varname == "avo":
        avo = wrf.getvar(ncfile,varname)
        lats, lons = wrf.latlon_coords(avo)
        # Coriolis parameter f = 2 * Omega * sin(lat)
        Omega = 7.2921e-5  # rad/s
        lat_radians = np.deg2rad(lats)
        f = 2 * Omega * np.sin(lat_radians)
        var = avo - f
        min_threshold = -30
        max_threshold = 110
        cmap = "plasma"
        label = r"Relative Voriticity ($10^{-5} s^{-1}$)"
    else:
        var = wrf.getvar(ncfile,varname)
        lats, lons = wrf.latlon_coords(var)
        min_threshold = -32
        max_threshold = 50
        cmap = "NWSRef"
        label = "Maximum Reflectivity (dBZ)"

    if level:
        pressure = wrf.getvar(ncfile, "pressure")
        var_lvl = wrf.interplevel(var, pressure, level)
        var_smooth = nanuniform_filter(var_lvl, size=size)
        if wind:
            u = wrf.getvar(ncfile, "ua")
            v = wrf.getvar(ncfile, "va")
            u = wrf.interplevel(u, pressure, level)
            v = wrf.interplevel(v, pressure, level)
    else:
        var_smooth = nanuniform_filter(var, size=size)
        if wind:
            u = wrf.getvar(ncfile, "U10")
            v = wrf.getvar(ncfile, "V10")

    levels=np.linspace(min_threshold, max_threshold, 20)
    lat_min = np.amin(lats)
    lat_max = np.amax(lats)
    lon_min = np.amin(lons)
    lon_max = np.amax(lons)

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 5), dpi=400)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)

    contour = ax.contourf(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(var_smooth), levels=levels, cmap=cmap, transform=ccrs.PlateCarree())
    
    if wind:
        strm = ax.streamplot(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(u), wrf.to_np(v), density=1.0, linewidth=1, arrowsize=1, color='k')

    # Gridlines
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2, 5))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2, 5))
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    cbar = plt.colorbar(contour, orientation='horizontal', pad=0.05, label=label)
    cbar.ax.tick_params(labelsize=9)
    cbar.ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))

    fig.suptitle(f"{Storm} (DA Exp. {Exper_name}): WRF Output at {level} hPa\n{timestamp}" if level else f"{Storm} (DA Exp. {Exper_name}): WRF Output\n{timestamp}", fontsize=9, fontweight='bold')
    plt.savefig(output_path, dpi=400)
    plt.close(fig)

def plot_tracking_wrf_variable(wrf_file, resolution, varname, level, wind, Storm, Exper_name, output_path, timestamp, center_lon, center_lat):
    """
    Plots a 2D WRF variable.

    Parameters:
    - varname: str, the WRF variable name (e.g., "avo", "slp", "tk", etc.)
    - wrf_file: str, path to the WRF output file
    - level: int or float, pressure level in hPa (only needed for 3D vars like "avo")
    """
    ncfile = nc.Dataset(wrf_file)
    size = round(100/resolution)

    domain_half_size = 3  # Degrees (half of 6-degree box)
    lat_min = center_lat - domain_half_size
    lat_max = center_lat + domain_half_size
    lon_min = center_lon - domain_half_size
    lon_max = center_lon + domain_half_size

    if varname == "slp":
        var = wrf.getvar(ncfile,varname)
        lats, lons = wrf.latlon_coords(var)
        min_threshold = 995
        max_threshold = 1030
        cmap = "viridis"
        label = "Sea Level Pressure (hPa)"
    elif varname == "avo":
        avo = wrf.getvar(ncfile,varname)
        lats, lons = wrf.latlon_coords(avo)
        # Coriolis parameter f = 2 * Omega * sin(lat)
        Omega = 7.2921e-5  # rad/s
        lat_radians = np.deg2rad(lats)
        f = 2 * Omega * np.sin(lat_radians)
        var = avo - f
        min_threshold = -30
        max_threshold = 110
        cmap = "plasma"
        label = r"Relative Voriticity ($10^{-5} s^{-1}$)"
    else:
        var = wrf.getvar(ncfile,varname)
        lats, lons = wrf.latlon_coords(var)
        min_threshold = -32
        max_threshold = 50
        cmap = "NWSRef"
        label = "Maximum Reflectivity (dBZ)"

    if level:
        pressure = wrf.getvar(ncfile, "pressure")
        var_lvl = wrf.interplevel(var, pressure, level)
        var_smooth = nanuniform_filter(var_lvl, size=size)
        if wind:
            u = wrf.getvar(ncfile, "ua")
            v = wrf.getvar(ncfile, "va")
            u = wrf.interplevel(u, pressure, level)
            v = wrf.interplevel(v, pressure, level)
    else:
        var_smooth = nanuniform_filter(var, size=size)
        if wind:
            u = wrf.getvar(ncfile, "U10")
            v = wrf.getvar(ncfile, "V10")

    levels=np.linspace(min_threshold, max_threshold, 20)

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5, 5), dpi=400)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)

    contour = ax.contourf(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(var_smooth), levels=levels, cmap=cmap, transform=ccrs.PlateCarree())

    if wind:
        strm = ax.streamplot(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(u), wrf.to_np(v), density=2.0, linewidth=1, arrowsize=1, color='k')

    # Gridlines
    lon_ticks = list(range(math.floor(lon_min), math.ceil(lon_max)+1, 1))
    lat_ticks = list(range(math.floor(lat_min), math.ceil(lat_max)+1, 1))
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    cbar = plt.colorbar(contour, orientation='horizontal', pad=0.05, label=label)
    cbar.ax.tick_params(labelsize=9)
    cbar.ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))

    fig.suptitle(f"{Storm} (DA Exp. {Exper_name}): WRF Output at {level} hPa\n{timestamp}" if level else f"{Storm} (DA Exp. {Exper_name}): WRF Output\n{timestamp}", fontsize=9, fontweight='bold')
    plt.savefig(output_path, dpi=400)
    plt.close(fig)

def process_target_time(args):
    target_time, wrf_dir, out_frame_dir, tracking_frame_dir, resolution, var, level, wind, Storm, Exper_name, storm_center_df, i, Full_domain, Tracking_domain = args
    dt_obj = datetime.strptime(target_time, "%Y%m%d%H%M")
    dt_str = dt_obj.strftime("%Y-%m-%d_%H:%M:%S")
    wrf_file = os.path.join(wrf_dir, f"wrfout_d01_{dt_str}")

    full_png_path = os.path.join(out_frame_dir, f"{target_time}.png")
    tracking_png_path = os.path.join(tracking_frame_dir, f"{target_time}.png")

    if not os.path.exists(wrf_file):
        print(f"Missing: {wrf_file} — Skipping...")
        return

    if Full_domain and os.path.exists(full_png_path):
        print(f"Skipping {target_time} — Full {var} PNG already exists")
    elif Full_domain:
        plot_wrf_variable(wrf_file, resolution, var, level[i], wind[i], Storm, Exper_name, full_png_path, target_time)
        print(f"Saved {var} full frame: {target_time}")

    if Tracking_domain and os.path.exists(tracking_png_path):
        print(f"Skipping {target_time} — Tracking {var} PNG already exists")
    elif Tracking_domain and dt_obj in storm_center_df.index:
        center_lon = storm_center_df.loc[dt_obj]["Center_Lon"]
        center_lat = storm_center_df.loc[dt_obj]["Center_Lat"]
        plot_tracking_wrf_variable(wrf_file, resolution, var, level[i], wind[i], Storm, Exper_name, tracking_png_path, target_time, center_lon, center_lat)
        print(f"Saved {var} tracking frame: {target_time}")

if __name__ == '__main__':
    # ---------- Configuration -------------------------
    Storm = 'Beryl'
    Exper_name = 'NoDA'
    resolution = 6
    variable = ['slp','avo','mdbz']
    level = [None,700,None]
    wind = [True,True,False]
    Full_domain = True
    Tracking_domain = True
    Storm_track_csv = '/expanse/lustre/projects/pen116/cpruett/Post_forecast_proc/Beryl_tracks/NoDA_tracks/2024-06-26_00:00:00_to_2024-07-09_12:00:00.csv'
    start_time_str = '202406260000'
    end_time_str = '202407091200'
    Consecutive_times = True
    vis_dir = f'./{Storm}_visuals'
    exp_dir = f'{vis_dir}/{Exper_name}_visuals'
    wrf_dir = '/expanse/lustre/scratch/cpruett/temp_project/BERYL/NoDA/WRF_FREE'
    # -------------------------------------------------------

    if not Consecutive_times:
        all_times = ['2017090500','201709050600','201709051200','201709051800','201709060000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in range(0, int(time_diff_hour)+1)]
        all_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    
    storm_center_df = Util_Vis.get_storm_center(Storm_track_csv) if Tracking_domain else None
    
    os.makedirs(vis_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)

    plot_num = len(variable)
    for i in range(0, plot_num):
        var = variable[i]
        out_frame_dir = f'{exp_dir}/{var}_frames'
        tracking_frame_dir = f'{out_frame_dir}/tracking_frames'
        gif_output_path = f'{out_frame_dir}/{Storm}_{Exper_name}_{var}.gif'
        zip_output_path = f'{out_frame_dir}/{Storm}_{Exper_name}_{var}_frames.zip'
        tracking_gif_output = f'{tracking_frame_dir}/{Storm}_{Exper_name}_{var}_tracking.gif'
        tracking_zip_output = f'{tracking_frame_dir}/{Storm}_{Exper_name}_{var}_tracking_frames.zip'

        os.makedirs(out_frame_dir, exist_ok=True)
        os.makedirs(tracking_frame_dir, exist_ok=True)

        print(f'------------ Creating {var} frames ----------------------')

        task_args = [(target_time, wrf_dir, out_frame_dir, tracking_frame_dir, resolution, var, level, wind, Storm, Exper_name, storm_center_df, i, Full_domain, Tracking_domain) for target_time in all_times]

        num_workers = min(32, mp.cpu_count())
        with mp.Pool(processes=num_workers) as pool:
            pool.map(process_target_time, task_args)

#        for target_time in all_times:
#            dt_obj = datetime.strptime(target_time, "%Y%m%d%H%M")
#            dt_str = dt_obj.strftime("%Y-%m-%d_%H:%M:%S")
#            wrf_file = os.path.join(wrf_dir, f"wrfout_d01_{dt_str}")

#            full_png_path = os.path.join(out_frame_dir, f"{target_time}.png")
#            tracking_png_path = os.path.join(tracking_frame_dir, f"{target_time}.png")

#            if not os.path.exists(wrf_file):
#                print(f"Missing: {wrf_file} — Skipping...")
#                continue

#            if (Full_domain and os.path.exists(full_png_path)):
#                print(f"Skipping {target_time} — Full {var} PNG already exists")
#            elif Full_domain:
#                plot_wrf_variable(wrf_file, resolution, var, level[i], wind[i], Storm, Exper_name, full_png_path, target_time)
#                print(f"Saved {var} full frame: {target_time}")
#            if (Tracking_domain and os.path.exists(tracking_png_path)):
#                print(f"Skipping {target_time} — Tracking {var} PNG already exists")
#            elif Tracking_domain and dt_obj in storm_center_df.index:
#                center_lon = storm_center_df.loc[dt_obj]["Center_Lon"]
#                center_lat = storm_center_df.loc[dt_obj]["Center_Lat"]
#                plot_tracking_wrf_variable(wrf_file, resolution, var, level[i], wind[i], Storm, Exper_name, tracking_png_path, target_time, center_lon, center_lat)
#                print(f"Saved {var} tracking frame: {target_time}")    

        if (os.path.exists(gif_output_path) and os.path.exists(zip_output_path)):
            print(f"Skipping {var} Full Domain GIF and zip as they already exist")
        elif Full_domain:
            print(f'------------ Creating {var} Full Domain GIF ----------------------')
            Util_Vis.create_gif(out_frame_dir, gif_output_path, fps=2)
            Util_Vis.zip_frames(out_frame_dir, zip_output_path)

        if (os.path.exists(tracking_gif_output) and os.path.exists(tracking_zip_output)):
            print(f"Skipping {var} Tracking Domain GIF and zip as they already exist")
        elif Tracking_domain:
            print(f'------------ Creating {var} Tracking Domain GIF ----------------------')
            Util_Vis.create_gif(tracking_frame_dir, tracking_gif_output, fps=0.003)
            Util_Vis.zip_frames(tracking_frame_dir, tracking_zip_output)


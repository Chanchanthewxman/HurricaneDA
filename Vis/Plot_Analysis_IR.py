import os
import glob
import numpy as np
import netCDF4 as nc
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

# Read crtm calculated IR data from one binary file
def read_simu_IR_one(Hxb_file):
    xmax = 1395  # model domain: west-east 
    ymax = 743  # model domain: south-north
    print('Hxb_file:' + Hxb_file)
    print('xmax, ymax: ' + str(xmax) + ' ' + str(ymax))

    Hxb_data = np.fromfile(Hxb_file, dtype='<f4')  # <: little endian; f: float; 4: 4 bytes
    n_ch = len(Hxb_data) / (xmax * ymax) - 2
    n_ch = int(n_ch)
    Hxb_sim = Hxb_data[:].reshape(n_ch + 2, ymax, xmax)

    dict_simu_Tb = {}
    dict_simu_Tb['Lon_x'] = Hxb_sim[0, :, :]
    dict_simu_Tb['Lat_x'] = Hxb_sim[1, :, :]
    dict_simu_Tb['Yb_x'] = Hxb_sim[2, :, :]

    return dict_simu_Tb

def plot_Tb(Storm, Exper_name, Hxb, FCtime, out_dir, channel):
    # Read simulated data
    d_simu = read_simu_IR_one(Hxb)

    # ------------------ Plot -----------------------
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()},
                           figsize=(10, 5), dpi=400)

    # Define the domain
    lat_min = np.amin(d_simu['Lat_x'])
    lat_max = np.amax(d_simu['Lat_x'])
    lon_min = np.amin(d_simu['Lon_x'])
    lon_max = np.amax(d_simu['Lon_x'])

    # Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap(0.5)

    # Plot data
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)
    cs = ax.scatter(d_simu['Lon_x'], d_simu['Lat_x'], 1, c=d_simu['Yb_x'],
                    edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,
                    transform=ccrs.PlateCarree())

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

    # Title and colorbar
    fig.suptitle(f"{Storm} (DA Exp. {Exper_name}): Simulated {channel} Tb\n{FCtime}", fontsize=9, fontweight='bold')
    cbar = fig.colorbar(cs, orientation="horizontal", pad=0.05, label='Brightness Temperature (K)')
    cbar.ax.tick_params(labelsize=9)

    # Save
    os.makedirs(out_dir, exist_ok=True)
    img_name = os.path.join(out_dir, f"{FCtime}.png")
    plt.savefig(img_name, dpi=400)
    plt.close(fig)
    print('Saved frame: ', img_name)

def plot_tracking_Tb(Storm, Exper_name, Hxb, FCtime, out_dir, channel, center_lon, center_lat):
    d_simu = read_simu_IR_one(Hxb)

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()},
                           figsize=(5, 5), dpi=400)

    domain_half_size = 3  # Degrees (half of 5-degree box)
    lat_min = center_lat - domain_half_size
    lat_max = center_lat + domain_half_size
    lon_min = center_lon - domain_half_size
    lon_max = center_lon + domain_half_size

    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap(0.5)

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)
    cs = ax.scatter(d_simu['Lon_x'], d_simu['Lat_x'], 12, c=d_simu['Yb_x'],
                    edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,
                    marker='s', transform=ccrs.PlateCarree())

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

    fig.suptitle(f"{Storm} (DA Exp. {Exper_name}): Simulated {channel} Tb\n{FCtime}", fontsize=9, fontweight='bold')
    cbar = fig.colorbar(cs, orientation="horizontal", pad=0.05, label='Brightness Temperature (K)')
    cbar.ax.tick_params(labelsize=9)

    os.makedirs(out_dir, exist_ok=True)
    img_name = os.path.join(out_dir, f"{FCtime}.png")
    plt.savefig(img_name, dpi=400)
    plt.close(fig)
    print('Saved TRACKING frame:', img_name)

if __name__ == '__main__':
    # ---------- Configuration -------------------------
    Storm = 'Beryl'
    Exper_name = 'NoDA'
    channel = 'Channel-7'
    Full_domain = True
    Tracking_domain = True
    Storm_track_csv = '/expanse/lustre/projects/pen116/cpruett/Post_forecast_proc/Beryl_tracks/NoDA_tracks/2024-06-26_00:00:00_to_2024-07-09_12:00:00.csv'
    start_time_str = '202406260000'
    end_time_str = '202407091200'
    Consecutive_times = True
    vis_dir = f'./{Storm}_visuals'
    exp_dir = f'{vis_dir}/{Exper_name}_visuals'
    out_frame_dir = f'{exp_dir}/IR_frames'
    tracking_frame_dir = f'{out_frame_dir}/tracking_frames'
    gif_output_path = f'{out_frame_dir}/{Storm}_{Exper_name}_IR.gif'
    tracking_gif_output = f'{tracking_frame_dir}/{Storm}_{Exper_name}_IR_tracking.gif'
    Analysis_dir = '/expanse/lustre/scratch/cpruett/temp_project/BERYL/CRTM_OUT/NoDA'
    # -------------------------------------------------------

    if not Consecutive_times:
        IR_times = ['2017090500','201709050600','201709051200','201709051800','201709060000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in range(0, int(time_diff_hour)+1)]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    
    storm_center_df = Util_Vis.get_storm_center(Storm_track_csv) if Tracking_domain else None
    
    os.makedirs(vis_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)
    os.makedirs(out_frame_dir, exist_ok=True)
    os.makedirs(tracking_frame_dir, exist_ok=True)
    
    print('------------ Creating frames ----------------------')
    for FCtime in IR_times:
        dt_obj = datetime.strptime(FCtime, "%Y%m%d%H%M")
        dt_str = dt_obj.strftime("%Y-%m-%d_%H:%M:%S")
        bin_path = os.path.join(Analysis_dir, f"GOES16_Ch7_wrfout_d01_{dt_str}.bin")

        full_png_path = os.path.join(out_frame_dir, f"{FCtime}.png")
        tracking_png_path = os.path.join(tracking_frame_dir, f"{FCtime}.png")

        if not os.path.exists(bin_path):
            print(f"Missing: {bin_path} — Skipping...")
            continue

        if (Full_domain and os.path.exists(full_png_path)):
            print(f"Skipping {FCtime} — Full PNG already exists")
        elif Full_domain:
            plot_Tb(Storm, Exper_name, bin_path, FCtime, out_frame_dir, channel)

        if (Tracking_domain and os.path.exists(tracking_png_path)):
            print(f"Skipping {FCtime} — Tracking PNG already exists")
        elif Tracking_domain and dt_obj in storm_center_df.index:
            center_lon = storm_center_df.loc[dt_obj]["Center_Lon"]
            center_lat = storm_center_df.loc[dt_obj]["Center_Lat"]
            plot_tracking_Tb(Storm, Exper_name, bin_path, FCtime, tracking_frame_dir, channel, center_lon, center_lat)
        

    if Full_domain:
        print('------------ Creating Full Domain GIF ----------------------')
        Util_Vis.create_gif(out_frame_dir, gif_output_path, fps=2)
        zip_output_path = f"{out_frame_dir}/{Storm}_{Exper_name}_IR_frames.zip"
        Util_Vis.zip_frames(out_frame_dir, zip_output_path)

    if Tracking_domain:
        print('------------ Creating Tracking Domain GIF ----------------------')
        Util_Vis.create_gif(tracking_frame_dir, tracking_gif_output, fps=0.003)
        zip_output_path = f"{tracking_frame_dir}/{Storm}_{Exper_name}_IR_tracking_frames.zip"
        Util_Vis.zip_frames(tracking_frame_dir, zip_output_path)

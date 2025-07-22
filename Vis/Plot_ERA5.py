import xarray as xr
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from datetime import datetime
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import pandas as pd
import math
import os
import glob
import sys
import gc
import Util_Vis

# --- Plotting Full ERA5 Function ---
def plot_era5(u, v, q, lons, lats, timestamp, level, output_path, Storm):
    lat_min = np.amin(lats)
    lat_max = np.amax(lats)
    lon_min = np.amin(lons)
    lon_max = np.amax(lons)
    step = 8  # subsample arrows

    u_plot = u[::step, ::step]
    v_plot = v[::step, ::step]
    q_plot = q[:, :]
    lon_plot = lons[::step]
    lat_plot = lats[::step]

    # Define experiment domain box
    exp_lat_min, exp_lat_max = 5.0, 35.0
    exp_lon_min, exp_lon_max = -100.0, -16.0
    box_lats = [exp_lat_min, exp_lat_min, exp_lat_max, exp_lat_max, exp_lat_min]
    box_lons = [exp_lon_min, exp_lon_max, exp_lon_max, exp_lon_min, exp_lon_min]

    min_q = 0
    max_q = 15 # (g/kg)

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 5), dpi=400)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)

    # Domain box
    ax.plot(box_lons, box_lats, transform=ccrs.PlateCarree(),
            color='red', linestyle='--', linewidth=2, label='Experiment Domain')

    # Specific humidity background
    humidity_plot = ax.contourf(lons, lats, q_plot*1000, levels=np.linspace(min_q, max_q, 20), cmap='YlGnBu', vmin=min_q, vmax=max_q, transform=ccrs.PlateCarree())

    # Wind vectors
    qv = ax.quiver(lon_plot, lat_plot, u_plot, v_plot, transform=ccrs.PlateCarree(), scale=400)
    ax.quiverkey(qv, 0.9, 0.1, 25, '25 m/s', labelpos='S', coordinates='axes',
                 color='white', labelsep=0.05, labelcolor='white', fontproperties={'size': 9})

    # Gridlines
    lon_ticks = list(range(math.ceil(lon_min)-2, math.floor(lon_max)+2, 8))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.floor(lat_max)+2, 8))
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

    cbar = plt.colorbar(humidity_plot, orientation='horizontal', pad=0.05, label='Specific Humidity (g/kg)')
    cbar.ax.tick_params(labelsize=9)

    plt.legend(loc='upper right', fontsize=9)
    fig.suptitle(f'{Storm}: ERA5 Reanalysis at {level} hPa\n{timestamp}', fontsize=9, fontweight='bold')
    plt.savefig(output_path, dpi=400)
    plt.close()

# --- Plotting Track-Domain ERA5 Function ---
def plot_tracking_era5(u, v, q, lons, lats, timestamp, level, output_path, Storm, center_lon, center_lat):
    
    domain_half_size = 2.5  # Degrees (half of 5-degree box)
    lat_min = center_lat - domain_half_size
    lat_max = center_lat + domain_half_size
    lon_min = center_lon - domain_half_size
    lon_max = center_lon + domain_half_size
    step = 5  # subsample arrows

    u_plot = u[::step, ::step]
    v_plot = v[::step, ::step]
    q_plot = q[:, :]
    lon_plot = lons[::step]
    lat_plot = lats[::step]

    min_q = 0
    max_q = 15 # (g/kg)

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5, 5), dpi=400)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)

    # Specific humidity background
    humidity_plot = ax.contourf(lons, lats, q_plot*1000, levels=np.linspace(min_q, max_q, 20), cmap='YlGnBu', vmin=min_q, vmax=max_q, transform=ccrs.PlateCarree())

    # Wind vectors
    qv = ax.quiver(lon_plot, lat_plot, u_plot, v_plot, transform=ccrs.PlateCarree(), scale=400)
    ax.quiverkey(qv, 0.9, 0.1, 25, '25 m/s', labelpos='S', coordinates='axes',
                 color='white', labelsep=0.05, labelcolor='white', fontproperties={'size': 9})

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

    cbar = plt.colorbar(humidity_plot, orientation='horizontal', pad=0.05, label='Specific Humidity (g/kg)')
    cbar.ax.tick_params(labelsize=9)

    fig.suptitle(f'{Storm}: ERA5 Reanalysis at {level} hPa\n{timestamp}', fontsize=9, fontweight='bold')
    plt.savefig(output_path, dpi=400)
    plt.close()


if __name__ == '__main__':
    # ---------- Configuration -------------------------
    Storm = 'Beryl'
    Exper_name = 'ERA5_SCOUT'
    level = 700
    Full_domain = False
    Tracking_domain = True
    Storm_track_csv = '/expanse/lustre/projects/pen116/cpruett/Post_forecast_proc/Beryl_tracks/ERA5_SCOUT_tracks/2024-06-25_12:00:00_to_2024-07-09_12:00:00.csv'
    start_time_str = '202406251200'
    end_time_str   = '202407091200'
    frequency = '1h'
    Consecutive_times = True
    vis_dir = f'./{Storm}_visuals'
    out_frame_dir = f"{vis_dir}/ERA5_{level}hPa_frames"
    tracking_frame_dir = f"{out_frame_dir}/{Exper_name}_tracking_frames"
    gif_indices = [84, 229]
    gif_basename = f"{out_frame_dir}/ERA5_{Storm}_{level}hPa"
    gif_tracking_basename = f"{tracking_frame_dir}/ERA5_{Storm}_{level}hPa_tracking"
    gif_output_path = f"{out_frame_dir}/ERA5_{Storm}_{level}hPa.gif"
    zip_output_path = f"{out_frame_dir}/ERA5_{Storm}_{level}hPa_frames.zip"
    tracking_gif_output = f"{tracking_frame_dir}/ERA5_{Storm}_{level}hPa_tracking.gif"
    tracking_zip_output = f"{tracking_frame_dir}/ERA5_{Storm}_{level}hPa_tracking_frames.zip"
    input_file = "/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/ERA5/era5_700hPa_junjuly2023.nc"
    # -------------------------------------------------------    
    
    if not Consecutive_times:
        all_times = pd.to_datetime([
            datetime(2023, 10, 16, 0, 0),
            datetime(2023, 10, 17, 5, 0),
            datetime(2023, 10, 19, 6, 0),
            datetime(2023, 10, 20, 9, 0),
            datetime(2023, 10, 20, 23, 0),
        ])
    else:
        start_time = pd.to_datetime(start_time_str, format="%Y%m%d%H%M")
        end_time = pd.to_datetime(end_time_str, format="%Y%m%d%H%M")
        all_times = pd.date_range(start=start_time, end=end_time, freq=frequency)
    
    storm_center_df = Util_Vis.get_storm_center(Storm_track_csv) if Tracking_domain else None
    
    os.makedirs(vis_dir, exist_ok=True)
    os.makedirs(out_frame_dir, exist_ok=True)
    os.makedirs(tracking_frame_dir, exist_ok=True)

    ds = xr.open_dataset(input_file)
    valid_times = pd.to_datetime(ds.valid_time.values)

    u_all = ds['u'].sel(pressure_level=level)
    v_all = ds['v'].sel(pressure_level=level)
    q_all = ds['q'].sel(pressure_level=level)
    lats = ds.latitude.values
    lons = ds.longitude.values

    print('------------ Creating frames ----------------------')
    for target_time in all_times:
        timestamp = target_time.strftime("%Y%m%d%H%M")
        dt_obj = datetime.strptime(timestamp, "%Y%m%d%H%M")
        full_png_path = os.path.join(out_frame_dir, f"{timestamp}.png")
        tracking_png_path = os.path.join(tracking_frame_dir, f"{timestamp}.png")

        if target_time not in valid_times.values:
            print(f"Skipping {target_time} (not in dataset)")
            continue

        idx = np.where(valid_times == target_time)[0][0]
        u = u_all[idx]
        v = v_all[idx]
        q = q_all[idx]

        
        if (Full_domain and os.path.exists(full_png_path)):
            print(f"Skipping {target_time} — Full PNG already exists")
        elif Full_domain:
            plot_era5(u, v, q, lons, lats, timestamp, level, full_png_path, Storm)
            print(f"Saved full frame: {target_time}")
        if (Tracking_domain and os.path.exists(tracking_png_path)):
            print(f"Skipping {target_time} — Tracking PNG already exists")
        elif Tracking_domain and dt_obj in storm_center_df.index:
            center_lon = storm_center_df.loc[dt_obj]["Center_Lon"]
            center_lat = storm_center_df.loc[dt_obj]["Center_Lat"]
            plot_tracking_era5(u, v, q, lons, lats, timestamp, level, tracking_png_path, Storm, center_lon, center_lat)
            print(f"Saved tracking frame: {target_time}")

    gif_pattern = f"{gif_basename}_part*.gif"
    gif_tracking_pattern = f"{gif_tracking_basename}_part*.gif"

    if (os.path.exists(gif_output_path) and os.path.exists(zip_output_path)):
        print(f"Skipping Full Domain GIF and zip as they already exist")
    elif (glob.glob(gif_pattern) and os.path.exists(zip_output_path)):
        print(f"Skipping Full Domain GIF and zip as they already exist")
    elif Full_domain:
        print('------------ Creating Full Domain GIF ----------------------')
       # try:
       #     Util_Vis.create_gif(out_frame_dir, gif_output_path, fps=2)
       # except np.core._exceptions._ArrayMemoryError as e:
       #     print("MemoryError encountered while creating GIF. Switching to chunked GIF creation...")
       #     Util_Vis.create_gifs_by_indices(out_frame_dir, gif_basename, gif_indices, fps=2)
        Util_Vis.create_gifs_by_indices(out_frame_dir, gif_basename, gif_indices, fps=2)
        Util_Vis.zip_frames(out_frame_dir, zip_output_path)

    if (os.path.exists(tracking_gif_output) and os.path.exists(tracking_zip_output)):
        print(f"Skipping Tracking Domain GIF and zip as they already exist")
    elif (glob.glob(gif_tracking_pattern) and os.path.exists(tracking_zip_output)):
        print(f"Skipping Tracking Domain GIF and zip as they already exist")
    elif Tracking_domain:
        print('------------ Creating Tracking Domain GIF ----------------------')
        #try:
         #   Util_Vis.create_gif(tracking_frame_dir, tracking_gif_output, fps=0.003)
        #except np.core._exceptions._ArrayMemoryError as e:
            #print("MemoryError encountered while creating GIF. Switching to chunked GIF creation...")
            
            # --- Clean up memory ---
           # for name in dir():
                #if not name.startswith('_'):
                #    obj = globals().get(name)
               #     try:
              #          if hasattr(obj, 'close'): obj.close()
             #       except: pass
            #        try:
           #             del obj
          #          except: pass

         #   gc.collect()  # Force garbage collection

        #    Util_Vis.create_gifs_by_indices(tracking_frame_dir, gif_tracking_basename, gif_indices, fps=0.003)
        Util_Vis.create_gifs_by_indices(tracking_frame_dir, gif_tracking_basename, gif_indices, fps=0.003)
        Util_Vis.zip_frames(tracking_frame_dir, tracking_zip_output)

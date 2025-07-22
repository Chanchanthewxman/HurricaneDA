import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
import os
import glob
import Util_Vis
import pandas as pd
from scipy.interpolate import griddata

def read_TSV_radiance(tsv_file):
    """
    Reads TSV-formatted radiance files and returns a dictionary of raw point data.
    """
    df = pd.read_csv(tsv_file, delim_whitespace=True, header=None)
    return {
        'lat_obs': df[3].values,
        'lon_obs': df[4].values,
        'Tb_obs': df[5].values
    }


def read_Tb(Tb_file):
    dataset = Dataset(Tb_file, 'r')
    n = 2  # Downsampling rate
    Tbs = dataset.variables['CMI'][::n, ::n]
    x = dataset.variables['x'][::n]
    y = dataset.variables['y'][::n] 
    l0 = dataset.variables['nominal_satellite_subpoint_lon'][0] * np.pi / 180

    req = 6378137.0
    rpol = 6356752.31414
    H = 42164160.0

    yy, xx = np.meshgrid(y, x)
    a = (np.sin(xx))**2 + (np.cos(xx))**2 * ((np.cos(yy))**2 + (req / rpol * np.sin(yy))**2)
    b = -2 * H * np.cos(xx) * np.cos(yy)
    c = H**2 - req**2

    r = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    sx = r * np.cos(xx) * np.cos(yy)
    sy = -r * np.sin(xx)
    sz = r * np.cos(xx) * np.sin(yy)

    lat = np.arctan(req**2 / rpol**2 * sz / np.sqrt((H - sx)**2 + sy**2))
    lon = l0 - np.arctan(sy / (H - sx))
    lat = np.degrees(lat)
    lon = np.degrees(lon)
    lat = np.transpose(lat)
    lon = np.transpose(lon)

    lat[np.iscomplex(lat)] = np.nan
    lon[np.iscomplex(lon)] = np.nan
    idx_noNan = ~np.isnan(lat)

    dataset.close()

    return {
        'lat_obs': lat[idx_noNan],
        'lon_obs': lon[idx_noNan],
        'Tb_obs': Tbs[idx_noNan],
    }

def plot_Tb(d_all, lat_min, lat_max, lon_min, lon_max, frame_path, timestamp, Satellite, channel):
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10, 5), dpi=400)
    cmap = Util_Vis.IRcmap(0.5)
    min_T = 185
    max_T = 325

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)
    cs = ax.scatter(d_all['lon_obs'], d_all['lat_obs'], 1.5, c=d_all['Tb_obs'], edgecolors='none',
                    cmap=cmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
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

    cbar = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.05, label='Brightness Temperature (K)')
    cbar.ax.tick_params(labelsize=9)
    fig.suptitle(f'{Satellite}: {channel} Tb\n{timestamp}', fontsize=9, fontweight='bold')
    plt.savefig(os.path.join(frame_path, f'{timestamp}.png'), dpi=400)
    plt.close()

def parse_datetime_from_filename(fname):
    try:
        parts = fname.split('_')
        timestamp = parts[-3][1:]  # 's20171671145342' → remove leading 's'
        dt = datetime.strptime(timestamp[:13], '%Y%j%H%M%S')  # '2017167114534' (first 13 chars)
        return dt
    except Exception as e:
        print(f"Error parsing {fname}: {e}")
        return None

def find_closest_file(target_dt, files):
    closest_file = None
    min_diff = timedelta.max
    for f in files:
        dt = parse_datetime_from_filename(os.path.basename(f))
        if dt:
            diff = abs(dt - target_dt)
            if diff < min_diff:
                min_diff = diff
                closest_file = f
    return closest_file

def plot_tracking_Tb(d_all, frame_path, timestamp, Satellite, channel, center_lat, center_lon):
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5, 5), dpi=400)
    
    domain_half_size = 2.5  # Degrees (half of 5-degree box)
    lat_min = center_lat - domain_half_size
    lat_max = center_lat + domain_half_size
    lon_min = center_lon - domain_half_size
    lon_max = center_lon + domain_half_size

    cmap = Util_Vis.IRcmap(0.5)
    min_T = 185
    max_T = 325

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)
    cs = ax.scatter(d_all['lon_obs'], d_all['lat_obs'], 10, c=d_all['Tb_obs'], edgecolors='none',
                    cmap=cmap, vmin=min_T, vmax=max_T, marker='s', transform=ccrs.PlateCarree())
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

    cbar = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.05, label='Brightness Temperature (K)')
    cbar.ax.tick_params(labelsize=9)
    fig.suptitle(f'{Satellite}: {channel} Tb\n{timestamp}', fontsize=9, fontweight='bold')
    plt.savefig(os.path.join(frame_path, f'{timestamp}.png'), dpi=400)
    plt.close()

if __name__ == '__main__':
    # ---------- Configuration -------------------------
    Storm = 'Beryl'
    Exper_name = 'NoDA'
    Satellite = 'GOES-18'
    channel = 'Channel-7'
    Full_domain = True
    Tracking_domain = False
    Storm_track_csv = '/expanse/lustre/projects/pen116/cpruett/Post_forecast_proc/Beryl_tracks/NoDA_track/2024-06-26_00:00:00_to_2024-07-09_12:00:00.csv'
    start_time_str = '202406260000'
    end_time_str = '202407091200'
    Consecutive_times = True
    interval_hours = 1
    vis_dir = f'./{Storm}_visuals'
    Observation_dir = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/OTIS/Data/IR_DATA/C07'
    Radiance_dir = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/Preprocess_Obs/toEnKFobs/GOES_IR'
    Processed = True
    # -------------------------------------------------------
    
    if not Consecutive_times:
        target_dts = [
            datetime(2023, 10, 16, 0, 0),
            datetime(2023, 10, 17, 5, 0),
            datetime(2023, 10, 19, 6, 0),
            datetime(2023, 10, 20, 9, 0),
            datetime(2023, 10, 20, 23, 0),
        ]
    else:
        start_dt = datetime.strptime(start_time_str, '%Y%m%d%H%M')
        end_dt = datetime.strptime(end_time_str, '%Y%m%d%H%M')
        target_dts = [start_dt + timedelta(hours=i) for i in range(0, int((end_dt - start_dt).total_seconds() // 3600) + 1, interval_hours)]
    
    storm_center_df = Util_Vis.get_storm_center(Storm_track_csv) if Tracking_domain else None

    os.makedirs(vis_dir, exist_ok=True)

    if Processed:
        processed_frame_dir = f'{vis_dir}/{Satellite}_{channel}_processed_frames'
        processed_tracking_frame_dir = f'{processed_frame_dir}/{Exper_name}_processed_tracking_frames'
        processed_gif_output_path = f'{processed_frame_dir}/{Satellite}_{Storm}_{channel}_processed_IR.gif'
        processed_tracking_gif_output = f'{processed_tracking_frame_dir}/{Satellite}_(Storm)_{channel}_processed_IR_tracking.gif'
    
        os.makedirs(processed_frame_dir, exist_ok=True)
        os.makedirs(processed_tracking_frame_dir, exist_ok=True)

        print('------------ Creating frames ----------------------')
        for dt in target_dts:
            timestamp_str = dt.strftime('%Y%m%d%H%M')
            radiance_filename = f"radiance_d01_{timestamp_str}_so"
            radiance_path = os.path.join(Radiance_dir, radiance_filename)

            full_png_path = os.path.join(processed_frame_dir, f"{timestamp_str}.png")
            tracking_png_path = os.path.join(processed_tracking_frame_dir, f"{timestamp_str}.png")
            
            if not os.path.exists(radiance_path):
                print(f"Missing TSV file: {radiance_path}")
                continue

            d_all = read_TSV_radiance(radiance_path)
            dt_obj = datetime.strptime(timestamp_str, "%Y%m%d%H%M")

            if (Full_domain and os.path.exists(full_png_path)):
                print(f"Skipping {timestamp_str} — Full PNG already exists")
            elif Full_domain:
                lat_min = np.nanmin(d_all['lat_obs'])
                lat_max = np.nanmax(d_all['lat_obs'])
                lon_min = np.nanmin(d_all['lon_obs'])
                lon_max = np.nanmax(d_all['lon_obs'])
                plot_Tb(d_all, lat_min, lat_max, lon_min, lon_max, processed_frame_dir, timestamp_str, Satellite, channel)
            
            if (Tracking_domain and os.path.exists(tracking_png_path)):
                print(f"Skipping {timestamp_str} — Tracking PNG already exists")
            elif Tracking_domain and dt_obj in storm_center_df.index:
                center_lon = storm_center_df.loc[dt_obj]["Center_Lon"]
                center_lat = storm_center_df.loc[dt_obj]["Center_Lat"]
                plot_tracking_Tb(d_all, tracking_frame_dir, timestamp_str, Satellite, channel, center_lat, center_lon)

            print(f"Saved frame: {timestamp_str}")

        if Full_domain:
            print('------------ Creating Full Domain GIF ----------------------')
            Util_Vis.create_gif(processed_frame_dir, processed_gif_output_path, fps=2)
            zip_output_path = f"{processed_frame_dir}/{Satellite}_{Storm}_{channel}_processed_IR_frames.zip"
            Util_Vis.zip_frames(processed_frame_dir, zip_output_path)

        if Tracking_domain:
            print('------------ Creating Tracking Domain GIF ----------------------')
            Util_Vis.create_gif(processed_tracking_frame_dir, processed_tracking_gif_output, fps=0.003)
            zip_output_path = f"{processed_tracking_frame_dir}/{Satellite}_{Storm}_{channel}_processed_IR_tracking_frames.zip"
            Util_Vis.zip_frames(processed_tracking_frame_dir, zip_output_path)
    
    else:
        julian_days = []
        julian_days = sorted(set(dt.strftime('%j') for dt in target_dts))
    
        out_frame_dir = f'{vis_dir}/{Satellite}_{channel}_frames'
        tracking_frame_dir = f'{out_frame_dir}/{Exper_name}_tracking_frames'
        gif_output_path = f'{out_frame_dir}/{Satellite}_{Storm}_{channel}_IR.gif'
        tracking_gif_output = f'{tracking_frame_dir}/{Satellite}_(Storm)_{channel}_IR_tracking.gif'

        os.makedirs(out_frame_dir, exist_ok=True)
        os.makedirs(tracking_frame_dir, exist_ok=True)

        all_files = []
        for jd in julian_days:
            all_files += glob.glob(os.path.join(Observation_dir, jd, '*.nc'))
        print('------------ Creating frames ----------------------')
        for dt in target_dts:
            timestamp_str = dt.strftime('%Y%m%d%H%M')
            full_png_path = os.path.join(out_frame_dir, f"{timestamp_str}.png")
            tracking_png_path = os.path.join(tracking_frame_dir, f"{timestamp_str}.png")
        
            closest_file = find_closest_file(dt, all_files)
            if closest_file:
                d_all = read_Tb(closest_file)
                dt_obj = datetime.strptime(timestamp_str, "%Y%m%d%H%M")            
            
                if (Full_domain and os.path.exists(full_png_path)):
                    print(f"Skipping {timestamp_str} — Full PNG already exists")
                elif Full_domain:
                    lat_min, lat_max = 1.82, 25.94 # Use ncdump -v LATX on wrfoutput file to get list of latitudes
                    lon_min, lon_max = -111.03, -85.97 # Use ncdump -v LONGX on wrfoutput file to get list of latitudes
                    plot_Tb(d_all, lat_min, lat_max, lon_min, lon_max, out_frame_dir, timestamp_str, Satellite, channel)
            
                if (Tracking_domain and os.path.exists(tracking_png_path)):
                    print(f"Skipping {timestamp_str} — Tracking PNG already exists")
                elif Tracking_domain and dt_obj in storm_center_df.index:
                    center_lon = storm_center_df.loc[dt_obj]["Center_Lon"]
                    center_lat = storm_center_df.loc[dt_obj]["Center_Lat"]
                    plot_tracking_Tb(d_all, tracking_frame_dir, timestamp_str, Satellite, channel, center_lat, center_lon)

                print(f"Saved frame: {timestamp_str}")
            else:
                print(f"No file found near {dt}")

        if Full_domain:
            print('------------ Creating Full Domain GIF ----------------------')
            Util_Vis.create_gif(out_frame_dir, gif_output_path, fps=2)
            zip_output_path = f"{out_frame_dir}/{Satellite}_{Storm}_{channel}_IR_frames.zip"
            Util_Vis.zip_frames(out_frame_dir, zip_output_path)

        if Tracking_domain:
            print('------------ Creating Tracking Domain GIF ----------------------')
            Util_Vis.create_gif(tracking_frame_dir, tracking_gif_output, fps=0.003)
            zip_output_path = f"{tracking_frame_dir}/{Satellite}_{Storm}_{channel}_IR_tracking_frames.zip"
            Util_Vis.zip_frames(tracking_frame_dir, zip_output_path)

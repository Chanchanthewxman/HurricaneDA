import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def read_ensemble_csvs(csv_dir, skip_rows, wind_threshold):
    csv_files = sorted(glob(os.path.join(csv_dir, "storm_center_*.csv")))
    data = {}
    for file in csv_files:
        member_id = os.path.splitext(os.path.basename(file))[0].split("_")[-1]
        df = pd.read_csv(file, parse_dates=["Datetime"])
        df = df.iloc[::skip_rows].reset_index(drop=True)
        df = df[df["Max_Wind"] >= wind_threshold].reset_index(drop=True)
        data[member_id] = df
    return data

def plot_mean_track(data_dict, track_plot_output, Storm, Exper_name):
    timestamps = data_dict[list(data_dict.keys())[0]]["Datetime"]

    all_lats = []
    all_lons = []

    for member_df in data_dict.values():
        all_lats.append(member_df["Center_Lat"].values)
        all_lons.append(member_df["Center_Lon"].values)

    all_lats = np.array(all_lats)
    all_lons = np.array(all_lons)

    mean_lats = np.mean(all_lats, axis=0)
    mean_lons = np.mean(all_lons, axis=0)

    # Mean start and end points
    start_lat = np.mean(all_lats[:, 0])
    start_lon = np.mean(all_lons[:, 0])
    end_lat = np.mean(all_lats[:, -1])
    end_lon = np.mean(all_lons[:, -1])

    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.gridlines(draw_labels=True)

    # Spaghetti tracks
    for i in range(all_lats.shape[0]):
        ax.plot(all_lons[i], all_lats[i], color="blue", linewidth=1, transform=ccrs.PlateCarree())

    # Mean track
    ax.plot(mean_lons, mean_lats, color="black", linewidth=4, label="Mean Track", transform=ccrs.PlateCarree())

    # Start and end markers
    ax.scatter(start_lon, start_lat, color='white', edgecolor='black', marker='*', s=250,
               transform=ccrs.PlateCarree(), zorder=5, label='Start')
    ax.scatter(end_lon, end_lat, color='white', edgecolor='black', marker='X', s=200,
               transform=ccrs.PlateCarree(), zorder=5, label='End')

    ax.set_title(f"{Storm} ({Exper_name}) Spaghetti Plot", fontsize=9, fontweight="bold")
    ax.legend()
    plt.savefig(track_plot_output, dpi=300, bbox_inches="tight")

def remove_outliers(data_array):
    """Remove outliers beyond 1.5 * IQR"""
    q1 = np.percentile(data_array, 25, axis=0)
    q3 = np.percentile(data_array, 75, axis=0)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    # Mask outliers
    mask = (data_array >= lower_bound) & (data_array <= upper_bound)
    masked = np.where(mask, data_array, np.nan)
    return masked

def plot_mean_intensity(data_dict, intensity_plot_output, Storm, Exper_name):
    timestamps = data_dict[list(data_dict.keys())[0]]["Datetime"]

    all_winds = []
    all_press = []

    for member_df in data_dict.values():
        all_winds.append(member_df["Max_Wind"].values)
        all_press.append(member_df["Min_Press"].values)

    all_winds = np.array(all_winds)
    all_press = np.array(all_press)

    # Remove outliers
    wind_clean = remove_outliers(all_winds)
    press_clean = remove_outliers(all_press)

    wind_mean = np.nanmean(wind_clean, axis=0)
    wind_q25 = np.nanpercentile(wind_clean, 25, axis=0)
    wind_q75 = np.nanpercentile(wind_clean, 75, axis=0)

    press_mean = np.nanmean(press_clean, axis=0)
    press_q25 = np.nanpercentile(press_clean, 25, axis=0)
    press_q75 = np.nanpercentile(press_clean, 75, axis=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Wind subplot
    ax1.plot(timestamps, wind_mean, color='blue', label='Mean Max Wind')
    ax1.fill_between(timestamps, wind_q25, wind_q75, color='blue', alpha=0.3, label='IQR')
    ax1.set_ylabel("Wind Speed (m/s)")
    ax1.grid(True)
    ax1.set_title(f"{Storm} ({Exper_name}) Mean Maximum Wind Speed", fontsize=9, fontweight="bold")
    ax1.legend()
    ax1.tick_params(axis='x', labelrotation=45)

    # --- Add dashed lines for Saffir-Simpson thresholds ---
    thresholds = [18, 33, 43, 50, 58, 70]
    ymin, ymax = ax1.get_ylim()
    for thresh in thresholds:
        if ymin <= thresh <= ymax:
            ax1.axhline(y=thresh, linestyle='--', color='gray', linewidth=0.8)
            ax1.text(timestamps[0], thresh + 0.5, f'{thresh} m/s',
                     fontsize=9, color='gray', va='bottom', ha='left')

    # Pressure subplot
    ax2.plot(timestamps, press_mean, color='red', label='Mean Min Pressure')
    ax2.fill_between(timestamps, press_q25, press_q75, color='red', alpha=0.3, label='IQR')
    ax2.set_ylabel("Min Pressure (hPa)")
    ax2.grid(True)
    ax2.set_title(f"{Storm} ({Exper_name}) Mean Minimum Pressure", fontsize=9, fontweight="bold")
    ax2.legend()
    ax2.tick_params(axis='x', labelrotation=45)

    plt.tight_layout()
    fig.autofmt_xdate()
    plt.savefig(intensity_plot_output, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    # ---------- Configuration -------------------------
    Storm = 'Otis'
    Exper_name = 'CONV_GTS'
    Storm_track_csv = '/expanse/lustre/projects/pen116/cpruett/Post_forecast_proc/Otis_CONV_GTS_2023-10-16_06:00:00_to_2023-10-21_06:00:00'
    wind_threshold = 0 # m/s
    skip_rows = 1 # 1 equals no skipping
    vis_dir = f'./{Storm}_visuals'
    exp_dir = f'{vis_dir}/{Exper_name}_visuals'
    TI_dir = f'{exp_dir}/Track_Intensity_visuals'
    plot_type_dir = f'{TI_dir}/Mean'
    track_plot_output = f'{plot_type_dir}/Mean_{Storm}_{Exper_name}_track.png'
    intensity_plot_output = f'{plot_type_dir}/Mean_{Storm}_{Exper_name}_intensity.png'
    # -------------------------------------------------------    
    
    os.makedirs(vis_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)
    os.makedirs(TI_dir, exist_ok=True)
    os.makedirs(plot_type_dir, exist_ok=True)

    data = read_ensemble_csvs(Storm_track_csv, skip_rows, wind_threshold)

    print("Generating mean track plot...")
    plot_mean_track(data, track_plot_output, Storm, Exper_name)

    print("Generating intensity time series plot...")
    plot_mean_intensity(data, intensity_plot_output, Storm, Exper_name)

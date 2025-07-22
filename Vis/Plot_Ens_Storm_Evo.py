import sys
import os
from glob import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
import Util_Vis as uv

def plot_storm_track_ens(csv_file, member_id, track_plot_output, wind_threshold, skip_rows, Storm, Exper_name):
    """Plots the storm track with colors based on wind speed."""
    df = pd.read_csv(csv_file)

    # Apply row skipping
    df = df.iloc[::skip_rows, :]

    # Apply wind speed threshold filter
    df = df[df["Max_Wind"] >= wind_threshold]

    if df.empty:
        print(f"No data points exceed the wind threshold of {wind_threshold} m/s.")
        return
    
    # Define the domain
    lat_min = df["Center_Lat"].min() - 3
    lat_max = df["Center_Lat"].max() + 3
    lon_min = df["Center_Lon"].min() - 3
    lon_max = df["Center_Lon"].max() + 3

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.COASTLINE, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # Gridlines
    lon_ticks = list(range(math.ceil(lon_min)-1, math.ceil(lon_max)+1, 3))
    lat_ticks = list(range(math.ceil(lat_min)-1, math.ceil(lat_max)+1, 3))

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

    # Get first and last points
    first_row = df.iloc[0]
    last_row = df.iloc[-1]

    # Plot storm track as a line
    ax.plot(df["Center_Lon"], df["Center_Lat"], linestyle='-', linewidth=1.5, color='black', alpha=0.7)

    # Plot storm track
    for _, row in df.iterrows():
        plt.scatter(row["Center_Lon"], row["Center_Lat"], 
                    color=uv.saffir_simpson_category(row["Max_Wind"]),
                    edgecolors='black', s=40)

    # Mark first and last points with distinct symbols
    plt.scatter(first_row["Center_Lon"], first_row["Center_Lat"], 
                color='white', edgecolors='black', marker='*', s=100, label="Start")
    plt.scatter(last_row["Center_Lon"], last_row["Center_Lat"], 
                color='white', edgecolors='black', marker='X', s=100, label="End")

    # Map of categories and their labels/colors
    category_labels = {
        "blue": "Tropical Depression (<18 m/s)",
        "green": "Tropical Storm (18-32 m/s)",
        "yellow": "Category 1 (33-42 m/s)",
        "orange": "Category 2 (43-49 m/s)",
        "red": "Category 3 (50-57 m/s)",
        "purple": "Category 4 (58-69 m/s)",
        "black": "Category 5 (≥70 m/s)"
    }

    # Define the order of colors from strongest to weakest
    color_order = ["black", "purple", "red", "orange", "yellow", "green", "blue"]

    # Determine which categories are used in the data
    used_colors = set(df["Max_Wind"].apply(uv.saffir_simpson_category))

    # Sort used colors by defined severity order
    sorted_used_colors = [color for color in color_order if color in used_colors]
    
    # Create legend only for used categories
    legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markersize=8,
                                 markerfacecolor=color, label=category_labels[color])
                      for color in sorted_used_colors] 
    
    # Add start and end markers to legend
    legend_patches.append(plt.Line2D([0], [0], marker='*', color='w', markersize=10, 
                                     markerfacecolor='white', markeredgecolor='black', label="Start"))
    legend_patches.append(plt.Line2D([0], [0], marker='X', color='w', markersize=10, 
                                     markerfacecolor='white', markeredgecolor='black', label="End"))

    plt.legend(handles=legend_patches, loc="upper left", fontsize=8, title="Saffir-Simpson Scale")

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title(f"Ensemble {member_id}: {Storm} ({Exper_name}) Track", fontsize=9, fontweight="bold")
    plt.grid()
    plt.savefig(track_plot_output, dpi=300, bbox_inches="tight")
    plt.close()

def plot_storm_intensity_ens(csv_file, member_id, intensity_plot_output, wind_threshold, skip_rows, Storm, Exper_name):
    """Plots max wind speed and min pressure over time."""
    df = pd.read_csv(csv_file)
    
    # Apply row skipping
    df = df.iloc[::skip_rows, :]

    # Apply wind speed threshold filter
    df = df[df["Max_Wind"] >= wind_threshold]

    df["Datetime"] = pd.to_datetime(df["Datetime"])
    
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax2 = ax1.twinx()
    
    ax1.plot(df["Datetime"], df["Max_Wind"], 'b-', label="Max Wind Speed")
    ax2.plot(df["Datetime"], df["Min_Press"], 'r-', label="Min Pressure")
    
    ax1.set_xlabel("Datetime", fontsize=9)
    ax1.set_ylabel("Wind Speed (m/s)", color='b', fontsize=9)
    ax2.set_ylabel("Surface Pressure (hPa)", color='r', fontsize=9)

    ax1.tick_params(axis='y', labelcolor='b', labelsize=9)
    ax2.tick_params(axis='y', labelcolor='r', labelsize=9)

    # --- Add dashed lines for Saffir-Simpson thresholds ---
    thresholds = [18, 33, 43, 50, 58, 70]
    ymin, ymax = ax1.get_ylim()
    for thresh in thresholds:
        if ymin <= thresh <= ymax:
            ax1.axhline(y=thresh, linestyle='--', color='gray', linewidth=0.8)
            ax1.text(df["Datetime"].iloc[0], thresh + 0.5, f'{thresh} m/s',
                     fontsize=9, color='gray', va='bottom', ha='left')

    plt.title(f"Ensemble {member_id}: {Storm} ({Exper_name}) Intensity Evolution", fontsize=9, fontweight="bold")
    fig.autofmt_xdate()
    for label in ax1.get_xticklabels():
        label.set_fontsize(9)
    plt.savefig(intensity_plot_output, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    # ---------- Configuration -------------------------
    Storm = 'Otis'
    Exper_name = 'CONV_GTS'
    csv_dir = '/expanse/lustre/projects/pen116/cpruett/Post_forecast_proc/Otis_CONV_GTS_2023-10-16_06:00:00_to_2023-10-21_06:00:00'
    wind_threshold = 0 # m/s
    skip_rows = 1 # 1 equals no skipping
    vis_dir = f'./{Storm}_visuals'
    exp_dir = f'{vis_dir}/{Exper_name}_visuals'
    TI_dir = f'{exp_dir}/Track_Intensity_visuals'
    plot_type_dir = f'{TI_dir}/Ensemble'
    # -------------------------------------------------------

    os.makedirs(vis_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)
    os.makedirs(TI_dir, exist_ok=True)
    os.makedirs(plot_type_dir, exist_ok=True)

    csv_files = sorted(glob(os.path.join(csv_dir, "storm_center_*.csv")))
    if not csv_files:
        print("No ensemble CSV files found.")
        sys.exit(1)

    for csv_file in csv_files:
        member_id = os.path.splitext(os.path.basename(csv_file))[0].split("_")[-1]

        track_plot_output = f'{plot_type_dir}/{member_id}_track.png'
        intensity_plot_output = f'{plot_type_dir}/{member_id}_intensity.png'

        print(f"Plotting ensemble member {member_id}...")

        plot_storm_track_ens(csv_file, member_id, track_plot_output, wind_threshold, skip_rows, Storm, Exper_name)
        plot_storm_intensity_ens(csv_file, member_id, intensity_plot_output, wind_threshold, skip_rows, Storm, Exper_name)
    
    # Zip all plots
    zip_output_path = f'{plot_type_dir}/{Storm}_{Exper_name}_Ensemble_Tracks_Intensities.zip'
    uv.zip_frames(plot_type_dir, zip_output_path)  

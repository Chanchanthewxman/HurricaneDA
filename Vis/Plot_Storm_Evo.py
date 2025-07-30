import os
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

def plot_storm_track(csv_file, track_plot_output, wind_threshold, skip_rows, Storm, Exper_name):
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
    
    start_date = first_row["Datetime"]
    end_date = last_row["Datetime"]

    # Add start and end markers to legend
    legend_patches.append(plt.Line2D([0], [0], marker='*', color='w', markersize=10, 
                                     markerfacecolor='white', markeredgecolor='black', label=f"Start {start_date}"))
    legend_patches.append(plt.Line2D([0], [0], marker='X', color='w', markersize=10, 
                                     markerfacecolor='white', markeredgecolor='black', label=f"End {end_date}"))

    plt.legend(handles=legend_patches, loc="upper right", fontsize=9, title="Saffir-Simpson Scale")

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title(f"{Storm} ({Exper_name}) Track", fontsize=9, fontweight="bold")
    plt.grid()
    plt.savefig(track_plot_output, dpi=300, bbox_inches="tight")

def plot_storm_track_exp(csv_file, track_plot_output, wind_threshold, skip_rows, Storm, Exper_name, exp_lat_min, exp_lat_max, exp_lon_min, exp_lon_max):
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
    track_lat_min = df["Center_Lat"].min() - 3
    track_lat_max = df["Center_Lat"].max() + 3
    track_lon_min = df["Center_Lon"].min() - 3
    track_lon_max = df["Center_Lon"].max() + 3

    exp_lat_min_adj = exp_lat_min - 3
    exp_lat_max_adj = exp_lat_max + 3
    exp_lon_min_adj = exp_lon_min - 3
    exp_lon_max_adj = exp_lon_max + 3

    lat_min = min(track_lat_min, exp_lat_min_adj)
    lat_max = max(track_lat_max, exp_lat_max_adj)
    lon_min = min(track_lon_min, exp_lon_min_adj)
    lon_max = max(track_lon_max, exp_lon_max_adj)

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

    start_date = first_row["Datetime"]
    end_date = last_row["Datetime"]

    # Add start and end markers to legend
    legend_patches.append(plt.Line2D([0], [0], marker='*', color='w', markersize=10,
                                     markerfacecolor='white', markeredgecolor='black', label=f"Start {start_date}"))
    legend_patches.append(plt.Line2D([0], [0], marker='X', color='w', markersize=10,
                                     markerfacecolor='white', markeredgecolor='black', label=f"End {end_date}"))

    # Plot experiment domain
    box_lats = [exp_lat_min, exp_lat_max, exp_lat_max, exp_lat_min, exp_lat_min]
    box_lons = [exp_lon_min, exp_lon_min, exp_lon_max, exp_lon_max, exp_lon_min]
    ax.plot(box_lons, box_lats, linestyle='--', color='red', linewidth=1.2, transform=ccrs.PlateCarree(), label='Experiment Domain')

    plt.legend(handles=legend_patches, loc="upper right", fontsize=9, title="Saffir-Simpson Scale")

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title(f"{Storm} ({Exper_name}) Track", fontsize=9, fontweight="bold")
    plt.grid()
    plt.savefig(track_plot_output, dpi=300, bbox_inches="tight")


def plot_storm_intensity(csv_file, intensity_plot_output,  wind_threshold, skip_rows, Storm, Exper_name):
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

    plt.title(f"{Storm} ({Exper_name}) Intensity Evolution", fontsize=9, fontweight="bold")
    fig.autofmt_xdate()
    for label in ax1.get_xticklabels():
        label.set_fontsize(9)
    plt.savefig(intensity_plot_output, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    # ---------- Configuration -------------------------
    Storm = 'Beryl'
    Exper_name = 'HURDAT'
    Storm_track_csv = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/HURDAT/BERYL_2024_best_track.csv'
    wind_threshold = 0 # m/s
    skip_rows = 1 # 1 equals no skipping
    AddExpDom = True  # Set to True to include experiment domain box
    exp_lat_min = -5
    exp_lat_max = 22
    exp_lon_min = -70
    exp_lon_max = -13
    vis_dir = f'./{Storm}_visuals'
    exp_dir = f'{vis_dir}/{Exper_name}_visuals'
    TI_dir = f'{exp_dir}/Track_Intensity_visuals'
    plot_type_dir = f'{TI_dir}/Archive'
    track_plot_output = f'{plot_type_dir}/{Storm}_{Exper_name}_track.png' 
    intensity_plot_output = f'{plot_type_dir}/{Storm}_{Exper_name}_intensity.png'
    # -------------------------------------------------------

    os.makedirs(vis_dir, exist_ok=True)
    os.makedirs(exp_dir, exist_ok=True)
    os.makedirs(TI_dir, exist_ok=True)
    os.makedirs(plot_type_dir, exist_ok=True)

    if AddExpDom:
        plot_storm_track_exp(Storm_track_csv, track_plot_output, wind_threshold, skip_rows, Storm, Exper_name, exp_lat_min, exp_lat_max, exp_lon_min, exp_lon_max)
        plot_storm_intensity(Storm_track_csv, intensity_plot_output, wind_threshold, skip_rows, Storm, Exper_name)
    else:
        plot_storm_track(Storm_track_csv, track_plot_output, wind_threshold, skip_rows, Storm, Exper_name)
        plot_storm_intensity(Storm_track_csv, intensity_plot_output, wind_threshold, skip_rows, Storm, Exper_name)

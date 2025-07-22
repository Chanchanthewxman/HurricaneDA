#!/usr/bin/env python3

# This module consists of functions for visualizing post-EnKF products.
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import imageio.v2 as imageio  # Use v2 for compatibility
import zipfile
import glob
import os
import subprocess
import pandas as pd
# --------------------------------------------------------------
# ---------------------- COLOR MAP -----------------------------
# --------------------------------------------------------------

# colormap for MW Tbs
def newJet(max_T, min_T, min_Jet):

    max_T = 300
    min_T = 80
    min_Jet = 150
    jetLength = max_T - min_Jet + 1
    notJetLength = (min_Jet - 1) - (min_T + 1) + 1 

    jet = cm.get_cmap('jet', max_T-min_T)
    Myjet = jet(np.linspace(0,1,jetLength))

    jet_red = Myjet[:,0]
    jet_green = Myjet[:,1]
    jet_blue = Myjet[:,2]

    jetAdd_red = np.linspace(1.0, jet_red[0], notJetLength)
    jetAdd_green = np.linspace(1.0, jet_green[0], notJetLength)
    jetAdd_blue  = np.linspace(1.0, jet_blue[0], notJetLength) 

    cm_red = np.concatenate([np.insert(jetAdd_red,0,0),np.append(jet_red,0)])
    cm_green = np.concatenate([np.insert(jetAdd_green,0,0),np.append(jet_green,0)])
    cm_blue = np.concatenate([np.insert(jetAdd_blue,0,0),np.append(jet_blue,0)])

    newJet_value = np.column_stack([cm_red,cm_green,cm_blue, np.ones(len(cm_red))])
    newJet = ListedColormap(newJet_value)
    
    return newJet

# Plot Tb differences
def newRWB(max_T, min_T, min_RWB):

    max_T = 10
    min_T = -10
    min_RWB = -8
    RWBLength = 18 # max_T - min_RWB + 1
    notRWBLength = 2 #(min_RWB - 1) - (min_T + 1) + 1

    RWB = cm.get_cmap('RdBu_r', max_T-min_T)
    MyRWB = RWB(np.linspace(0,1,RWBLength))

    RWB_red = MyRWB[:,0]
    RWB_green = MyRWB[:,1]
    RWB_blue = MyRWB[:,2]

    #RWBAdd_red = np.linspace(1.0, RWB_red[0], notRWBLength)
    #RWBAdd_green = np.linspace(1.0, RWB_green[0], notRWBLength)
    #RWBAdd_blue  = np.linspace(1.0, RWB_blue[0], notRWBLength)

    cm_red = np.append( np.insert(RWB_red,0,0), 0)
    cm_green = np.append( np.insert(RWB_green,0,0), 0)
    cm_blue = np.append( np.insert(RWB_blue,0,0), 0)

    newRWB_value = np.column_stack([cm_red,cm_green,cm_blue])
    #newRWB_value = np.column_stack([cm_red,cm_green,cm_blue, np.ones(len(cm_red))])
    newRWB = ListedColormap(newRWB_value)

    return newRWB

# colormap for IR Tbs
def IRcmap(c_int):
    
    cmap = np.zeros( (int(140/c_int),3) ) 

    cmap[0:int(5/c_int), 0] = np.append( np.arange(0.5, 1, c_int/(5-c_int)*0.5), 1.0)
    cmap[0:int(5/c_int), 1] = 0
    cmap[0:int(5/c_int), 2] = np.append( np.arange(0.5, 1, c_int/(5-c_int)*0.5), 1.0)
    

    cmap[int(5/c_int)-1:int(10/c_int):int(c_int/c_int),0] = 1
    cmap[int(5/c_int)-1:int(10/c_int):int(c_int/c_int),1] = np.append( np.arange(0, 0.5, c_int/5*0.5), 0.5)
    cmap[int(5/c_int)-1:int(10/c_int):int(c_int/c_int),2] = 1

    cmap[int(10/c_int)-1:int(20/c_int):int(c_int/c_int),0] = np.append( np.arange(0.8, 0.0, -c_int/10*0.8), 0.0)
    cmap[int(10/c_int)-1:int(20/c_int):int(c_int/c_int),1] = np.append( np.arange(0.8, 0.0, -c_int/10*0.8), 0.0)
    cmap[int(10/c_int)-1:int(20/c_int):int(c_int/c_int),2] = np.append( np.arange(0.8, 0.0, -c_int/10*0.8), 0.0)

    cmap[int(20/c_int)-1:int(30/c_int):int(c_int/c_int),0] = np.append( np.arange(0, 1.0, c_int/10), 1)
    cmap[int(20/c_int)-1:int(30/c_int):int(c_int/c_int),1] = 0
    cmap[int(20/c_int)-1:int(30/c_int):int(c_int/c_int),2] = 0

    cmap[int(30/c_int)-1:int(40/c_int):int(c_int/c_int),0] = 1
    cmap[int(30/c_int)-1:int(40/c_int):int(c_int/c_int),1] = np.append( np.arange(0, 1.0, c_int/10), 1)
    cmap[int(30/c_int)-1:int(40/c_int):int(c_int/c_int),2] = 0

    cmap[int(40/c_int)-1:int(50/c_int):int(c_int/c_int),0] = np.append( np.arange(1.0, 0, -c_int/10), 0) 
    cmap[int(40/c_int)-1:int(50/c_int):int(c_int/c_int),1] = 1
    cmap[int(40/c_int)-1:int(50/c_int):int(c_int/c_int),2] = 0

    cmap[int(50/c_int)-1:int(60/c_int):int(c_int/c_int),0] = 0
    cmap[int(50/c_int)-1:int(60/c_int):int(c_int/c_int),1] = np.append( np.arange(1.0, 0, -c_int/10), 0)  
    cmap[int(50/c_int)-1:int(60/c_int):int(c_int/c_int),2] = np.append( np.arange(0, 1.0, c_int/10), 1)

    cmap[int(60/c_int)-1:int(70/c_int):int(c_int/c_int),0] = 0
    cmap[int(60/c_int)-1:int(70/c_int):int(c_int/c_int),1] = np.append( np.arange(0, 1.0, c_int/10), 1)
    cmap[int(60/c_int)-1:int(70/c_int):int(c_int/c_int),2] = 1
    
    cmap[int(70/c_int)-1:int(140/c_int):int(c_int/c_int),0] = np.append( np.arange(1.0, 0, -c_int/70), 0)
    cmap[int(70/c_int)-1:int(140/c_int):int(c_int/c_int),1] = np.append( np.arange(1.0, 0, -c_int/70), 0)
    cmap[int(70/c_int)-1:int(140/c_int):int(c_int/c_int),2] = np.append( np.arange(1.0, 0, -c_int/70), 0)

    IRcmap = ListedColormap( cmap )
    
    return IRcmap

def create_gif_ffmpeg(frame_dir, gif_name, fps=2):
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-pattern_type", "glob",
        "-i", os.path.join(frame_dir, "*.png"),
        "-vf", "scale=800:-1",  # optional scaling
        gif_name
    ]
    subprocess.run(cmd)

def create_gifs_by_indices(frame_dir, gif_basename, indices, fps=2):
    """
    Create multiple GIFs from frames in `frame_dir`, based on `indices` defining segment breaks.
    
    Parameters:
    - frame_dir (str): Path to directory containing PNG frames.
    - gif_basename (str): Base name for the output GIFs (e.g., 'ERA5_Otis_700hPa_part').
    - indices (list of int): List of frame indices to use as boundaries for splitting the GIFs.
    - fps (int): Frames per second for the GIF.
    """
    frame_files = sorted(glob.glob(os.path.join(frame_dir, '*.png')))
    total_frames = len(frame_files)

    # Ensure indices are sorted and include 0 at the beginning and total_frames at the end
    all_indices = [0] + indices + [total_frames]

    for i in range(len(all_indices) - 1):
        start = all_indices[i]
        end = all_indices[i + 1]
        segment_frames = frame_files[start:end]
        if not segment_frames:
            continue
        images = [imageio.imread(f) for f in segment_frames]
        gif_name = f"{gif_basename}_part{i+1}.gif"
        imageio.mimsave(gif_name, images, duration=1/fps)
        start_print = start + 1
        print(f"GIF created: {gif_name} (frames {start_print} through {end})")

def create_gif(frame_dir, gif_name, fps):
    images = []
    frame_files = sorted(glob.glob(os.path.join(frame_dir, '*.png')))
    for filename in frame_files:
        images.append(imageio.imread(filename))
    imageio.mimsave(gif_name, images, duration=1/fps)
    print('GIF created:', gif_name)

def zip_frames(frame_dir, zip_name):
    """Zip all PNG frames in a directory."""
    with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for frame in sorted(os.listdir(frame_dir)):
            if frame.endswith('.png'):
                frame_path = os.path.join(frame_dir, frame)
                zipf.write(frame_path, arcname=frame)  # arcname avoids full path
    print(f"Zipped frames into: {zip_name}")

def get_storm_center(csv_path):
    df = pd.read_csv(csv_path, parse_dates=["Datetime"])
    df.set_index("Datetime", inplace=True)
    return df

def saffir_simpson_category(wind_speed):
    """Returns a color based on the Saffir-Simpson scale (m/s)."""
    if wind_speed < 18:
        return "blue"  # Tropical Depression
    elif wind_speed < 33:
        return "green"  # Tropical Storm
    elif wind_speed < 43:
        return "yellow"  # Category 1
    elif wind_speed < 50:
        return "orange"  # Category 2
    elif wind_speed < 58:
        return "red"  # Category 3
    elif wind_speed < 70:
        return "purple"  # Category 4
    else:
        return "black"  # Category 5

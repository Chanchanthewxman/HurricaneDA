import os
import imageio
from PIL import Image

def combine_multiple_frame_sets(frame_dirs, output_gif, layout='horizontal', duration=300):
    """
    Combine frames from multiple directories into one GIF.

    Parameters:
        frame_dirs (list): List of directories containing .png frames (2 or 3).
        output_gif (str): Output path for the combined GIF.
        layout (str): 'horizontal' or 'vertical' stacking.
        duration (int): Duration per frame in milliseconds.
    """
    if not (2 <= len(frame_dirs) <= 3):
        raise ValueError("Only 2 or 3 frame directories are supported.")

    # Get sorted common filenames across all sets
    common_files = None
    for directory in frame_dirs:
        png_files = sorted([f for f in os.listdir(directory) if f.endswith('.png')])
        if common_files is None:
            common_files = set(png_files)
        else:
            common_files &= set(png_files)
    matched_times = sorted(common_files)

    frames = []
    for fname in matched_times:
        imgs = [Image.open(os.path.join(d, fname)) for d in frame_dirs]

        if layout == 'horizontal':
            total_width = sum(img.width for img in imgs)
            max_height = max(img.height for img in imgs)
            new_img = Image.new('RGB', (total_width, max_height))
            x_offset = 0
            for img in imgs:
                new_img.paste(img, (x_offset, 0))
                x_offset += img.width
        elif layout == 'vertical':
            max_width = max(img.width for img in imgs)
            total_height = sum(img.height for img in imgs)
            new_img = Image.new('RGB', (max_width, total_height))
            y_offset = 0
            for img in imgs:
                new_img.paste(img, (0, y_offset))
                y_offset += img.height
        else:
            raise ValueError("Layout must be 'horizontal' or 'vertical'")

        frames.append(new_img)

    frames[0].save(output_gif, save_all=True, append_images=frames[1:], duration=duration, loop=0)
    print(f"Combined GIF saved to: {output_gif}")
    
# Example usage
combine_multiple_frame_sets(
    frame_dirs=[
        '/expanse/lustre/projects/pen116/cpruett/Vis/ERA5_Otis_700hPa_frames/tracking_frames',
        '/expanse/lustre/projects/pen116/cpruett/Vis/GOES-18_Channel-7_frames/tracking_frames',
        '/expanse/lustre/projects/pen116/cpruett/Vis/Otis_CONV_GTS_frames/tracking_frames'
    ],
    output_gif='Comparion_Otis_ERA5-GOES-CONV_GTS_tracking.gif',
    layout='horizontal',  # or 'vertical'
    duration=300  # ms per frame
)


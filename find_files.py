#!/usr/bin/env python3

# Example script to find and read satellite data files
#!/usr/bin/env python3
# I have a list of files from glob. They are all FCI netCDF files.
# I don't want to use satpy to read them, I want to read them directly with netCDF4.
# Create a python script that creates an image from the data in these files.
# The image should be a composite of the data from all files.
# The image should be saved as a PNG file. The resolution should be 1km (11136x11136 pixels).
# Use the band vir_04 (visible channel 4) for the composite.
import hdf5plugin


from datetime import datetime
import warnings
import numpy as np
import netCDF4 as nc
from PIL import Image
import os
from glob import glob

def read_fci_radiance(file_path):
    """
    Read VIS_04 effective radiance from FCI netCDF file.
    
    Args:
        file_path: Path to the netCDF file
        
    Returns:
        numpy array of radiance values
    """
    try:
        dataset = nc.Dataset(file_path, 'r')
        radiance = dataset['/data/vis_04/measured/effective_radiance'][:]
        dataset.close()
        return radiance
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def create_composite(file_list, output_size=(11136, 11136)):
    """
    Create a composite image from multiple FCI files.
    
    Args:
        file_list: List of netCDF file paths
        output_size: Tuple of (width, height) for output image
        
    Returns:
        numpy array containing the composite image
    """
    print(f"Creating composite from {len(file_list)} files...")
    
    # Initialize the composite array
    composite = np.zeros(output_size, dtype=np.float32)
    count = np.zeros(output_size, dtype=np.int32)
    
    for i, file_path in enumerate(file_list):
        print(f"Processing file {i+1}/{len(file_list)}: {os.path.basename(file_path)}")
        
        radiance = read_fci_radiance(file_path)
        
        if radiance is None:
            continue
            
        # Handle different array dimensions
        if radiance.ndim == 2:
            # If data is already 2D, use directly
            h, w = radiance.shape
            
            # Ensure the data fits within our composite size
            h_end = min(h, output_size[0])
            w_end = min(w, output_size[1])
            
            # Add to composite (averaging overlapping regions)
            composite[:h_end, :w_end] += radiance[:h_end, :w_end].astype(np.float32)
            count[:h_end, :w_end] += 1
        else:
            print(f"Unexpected data shape: {radiance.shape}")
            continue
    
    # Average overlapping regions
    mask = count > 0
    composite[mask] /= count[mask]
    
    return composite

def normalize_to_8bit(data):
    """
    Normalize radiance data to 8-bit range (0-255) for image display.
    
    Args:
        data: Input array with radiance values
        
    Returns:
        8-bit normalized array
    """
    # Remove invalid values
    valid_mask = np.isfinite(data) & (data > 0)
    
    if not np.any(valid_mask):
        print("Warning: No valid data found!")
        return np.zeros_like(data, dtype=np.uint8)
    
    # Calculate percentiles for contrast enhancement
    p2, p98 = np.percentile(data[valid_mask], [2, 98])
    
    # Clip and normalize
    data_clipped = np.clip(data, p2, p98)
    data_normalized = ((data_clipped - p2) / (p98 - p2) * 255)
    
    return data_normalized.astype(np.uint8)

def main():
    # Get list of FCI netCDF files
    # Modify this pattern to match your files

    file_pattern = "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_00[0-9][0-9].nc"  # Adjust pattern as needed
    file_list = sorted(glob(file_pattern))

    
    if not file_list:
        print(f"No files found matching pattern: {file_pattern}")
        print("Please update the file_pattern variable in the script.")
        return
    
    print(f"Found {len(file_list)} files")
    
    # Create composite
    composite = create_composite(file_list, output_size=(11136, 11136))
    
    # Normalize to 8-bit
    print("Normalizing data to 8-bit...")
    image_data = normalize_to_8bit(composite)
    
    # Create and save image
    print("Creating image...")
    img = Image.fromarray(image_data, mode='L')
    
    output_filename = "fci_vis04_composite.png"
    print(f"Saving image to {output_filename}...")
    img.save(output_filename)
    
    print(f"Done! Image saved as {output_filename}")
    print(f"Image size: {img.size}")

if __name__ == "__main__":
    main()

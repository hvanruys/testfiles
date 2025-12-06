import numpy as np
import xarray as xr
from PIL import Image
import glob
import hdf5plugin
import matplotlib.pyplot as plt
import time

# Start timer
start_time = time.time()

# Configuration
file_pattern = "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_00[0-9][0-9].nc"
files = sorted(glob.glob(file_pattern))
files.pop()  # Remove last file if needed
output_file = 'fci_composite_xarray.png'
resolution = 11136  # 1km resolution

# Initialize output array
composite_R = np.zeros((resolution, resolution), dtype=np.float32)
composite_G = np.zeros((resolution, resolution), dtype=np.float32)
composite_B = np.zeros((resolution, resolution), dtype=np.float32)

total_subsolar_lat = np.zeros((resolution, resolution), dtype=np.float32)
total_subsolar_lon = np.zeros((resolution, resolution), dtype=np.float32)

# Process each file
for file_path in sorted(files):
    print(f"Processing: {file_path}")

    try:
        # Open netCDF file with xarray - open groups explicitly
        
        index_offset = xr.open_dataset(file_path, engine='netcdf4')['index_offset'].values
        print(f"index_offset: {index_offset}")
        
        # Read radiance data from respective groups
        radiance_B = xr.open_dataset(file_path, engine='netcdf4', group='data/vis_04/measured')['effective_radiance'].values
        radiance_G = xr.open_dataset(file_path, engine='netcdf4', group='data/vis_05/measured')['effective_radiance'].values
        radiance_R = xr.open_dataset(file_path, engine='netcdf4', group='data/vis_06/measured')['effective_radiance'].values
        
        # Get index_map as raw integers
        index_map = xr.open_dataset(file_path, engine='netcdf4', group='data/vis_04/measured')['index_map']
        #index_map = index_map_var.values.astype(np.uint16)
        
        my_fill_value = index_map.attrs.get('_FillValue', 65535)
        
        # Read position information
        start_row = int(xr.open_dataset(file_path, engine='netcdf4', group='data/vis_06/measured')['start_position_row'].values)
        start_col = int(xr.open_dataset(file_path, engine='netcdf4', group='data/vis_06/measured')['start_position_column'].values)
        end_row = int(xr.open_dataset(file_path, engine='netcdf4', group='data/vis_06/measured')['end_position_row'].values)
        end_col = int(xr.open_dataset(file_path, engine='netcdf4', group='data/vis_06/measured')['end_position_column'].values)
        
        # Get the geometric parameter vectors
        subsolar_lat = xr.open_dataset(file_path, engine='netcdf4', group='state/celestial')['subsolar_latitude'].values
        subsolar_lon = xr.open_dataset(file_path, engine='netcdf4', group='state/celestial')['subsolar_longitude'].values
        subsat_lat = xr.open_dataset(file_path, engine='netcdf4', group='state/platform')['subsatellite_latitude'].values
        subsat_lon = xr.open_dataset(file_path, engine='netcdf4', group='state/platform')['subsatellite_longitude'].values

        # Place segment in composite (convert to 0-based indexing)
        composite_R[start_row-1:end_row, start_col-1:end_col] = radiance_R
        composite_G[start_row-1:end_row, start_col-1:end_col] = radiance_G
        composite_B[start_row-1:end_row, start_col-1:end_col] = radiance_B
        
        # Adjust index_map by subtracting index_offset
        index_map_adjusted = np.where(
            index_map == 65535,
            65535,
            index_map - index_offset
        )
        
        # Use ravel() to flatten the 2D indices and then reshape back
        indices_flat = index_map_adjusted.ravel()
        # Clamp indices to valid range before indexing
        valid_indices = np.clip(indices_flat, 0, len(subsolar_lat) - 1)
        
        # Create masked arrays to handle fill values
        subsolar_lat_values = np.where(
            index_map_adjusted == 65535,
            0,
            subsolar_lat[valid_indices].reshape(index_map_adjusted.shape)
        )
        total_subsolar_lat[start_row-1:end_row, start_col-1:end_col] = subsolar_lat_values
        
        subsolar_lon_values = np.where(
            index_map_adjusted == 65535,
            0,
            subsolar_lon[valid_indices].reshape(index_map_adjusted.shape)
        )
        total_subsolar_lon[start_row-1:end_row, start_col-1:end_col] = subsolar_lon_values
        
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        continue

# Normalize to 0-255 range for image output
print(f"composite_R raw range min : {composite_R.min()}")
print(f"composite_R raw range max : {composite_R.max()}")
print(f"composite_R unique values: {len(np.unique(composite_R))}")

valid_data_R = composite_R[composite_R != 65535.0]
valid_data_G = composite_G[composite_G != 65535.0]
valid_data_B = composite_B[composite_B != 65535.0]

print(f"Valid data R len: {len(valid_data_R)}")
print(f"Composite_R len: {len(composite_R.flatten())}")

if len(valid_data_R) > 0:
    vmin = np.percentile(valid_data_R, 0)
    vmax = np.percentile(valid_data_R, 100)
    
    composite_normalized_R = np.where(
        composite_R == 65535.0,
        0,
        np.clip((composite_R - vmin) / (vmax - vmin) * 255, 0, 255)
    )
else:
    composite_normalized_R = composite_R

if len(valid_data_G) > 0:
    vmin = np.percentile(valid_data_G, 0)
    vmax = np.percentile(valid_data_G, 100)
    
    composite_normalized_G = np.where(
        composite_G == 65535.0,
        0,
        np.clip((composite_G - vmin) / (vmax - vmin) * 255, 0, 255)
    )
else:
    composite_normalized_G = composite_G

if len(valid_data_B) > 0:
    vmin = np.percentile(valid_data_B, 1)
    vmax = np.percentile(valid_data_B, 100)
    
    composite_normalized_B = np.where(
        composite_B == 65535.0,
        0,
        np.clip((composite_B - vmin) / (vmax - vmin) * 255, 0, 255)
    )
else:
    composite_normalized_B = composite_B

# Create alpha channel
alpha_channel = np.where(composite_R == 65535.0, 0, 255).astype(np.uint8)

# Convert to uint8
composite_uint8_R = composite_normalized_R.astype(np.uint8)
composite_uint8_G = composite_normalized_G.astype(np.uint8)
composite_uint8_B = composite_normalized_B.astype(np.uint8)

# Stack RGBA channels
rgb_array = np.stack([composite_uint8_R, composite_uint8_G, composite_uint8_B, alpha_channel], axis=-1)

# Save as PNG
img = Image.fromarray(rgb_array, mode='RGBA')
imgout = img.transpose(method=Image.Transpose.FLIP_TOP_BOTTOM)
imgout.save(output_file)

print(f"Composite image saved to: {output_file}")
print(f"Image shape: {composite_uint8_R.shape}")
print(f"Value range: {composite_uint8_R.min()} - {composite_uint8_R.max()}")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nScript execution time: {elapsed_time:.2f} seconds")
print(f"Script execution time: {elapsed_time/60:.2f} minutes")
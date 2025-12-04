import numpy as np
from netCDF4 import Dataset
from PIL import Image
import glob
import hdf5plugin
import matplotlib.pyplot as plt
import time

# Start timer
start_time = time.time()
# Configuration
file_pattern = "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_00[0-9][0-9].nc"  # Adjust pattern as needed
files = sorted(glob.glob(file_pattern))
files.pop() # Remove last file if needed
output_file = 'fci_composite.png'
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
        # Open netCDF file
        nc = Dataset(file_path, 'r')
        index_offset = nc.variables['index_offset'][:]
        print(f"index_offset: {index_offset}")
          
        # Read radiance data
        #radiance_B = nc['/data/vis_04/measured/effective_radiance'][:]
        radiance_B = nc.groups['data'].groups['vis_04'].groups['measured'].variables['effective_radiance'][:]
        radiance_G = nc['/data/vis_05/measured/effective_radiance'][:]
        radiance_R = nc['/data/vis_06/measured/effective_radiance'][:]

        # Get index_map as raw integers
        index_map_var = nc.groups['data'].groups['vis_04'].groups['measured'].variables['index_map']
        index_map_var.set_auto_scale(False)
        index_map = index_map_var[:].astype(np.uint16)

        #print(f"index_map_var shape: {index_map_var.shape}")
        #print(f"index_map_var dtype: {index_map_var.dtype}")
                

        my_fill_value = index_map_var.getncattr('_FillValue')


        # Read position information
        start_row = int(nc['/data/vis_06/measured/start_position_row'][:])
        start_col = int(nc['/data/vis_06/measured/start_position_column'][:])
        end_row = int(nc['/data/vis_06/measured/end_position_row'][:])
        end_col = int(nc['/data/vis_06/measured/end_position_column'][:])
        
        # Get the geometric parameter vectors
        subsolar_lat = nc.groups['state'].groups['celestial'].variables['subsolar_latitude'][:]
        subsolar_lon = nc.groups['state'].groups['celestial'].variables['subsolar_longitude'][:]
        subsat_lat = nc.groups['state'].groups['platform'].variables['subsatellite_latitude'][:]
        subsat_lon = nc.groups['state'].groups['platform'].variables['subsatellite_longitude'][:]

        # Replace fill values (65535) with 0
        #radiance_R = np.where(radiance_R == 65535, 0, radiance_R)
        #radiance_G = np.where(radiance_G == 65535, 0, radiance_G)
        #radiance_b = np.where(radiance_B == 65535, 0, radiance_B)
        
        # Place segment in composite (convert to 0-based indexing)
        composite_R[start_row-1:end_row, start_col-1:end_col] = radiance_R
        composite_G[start_row-1:end_row, start_col-1:end_col] = radiance_G
        composite_B[start_row-1:end_row, start_col-1:end_col] = radiance_B

        index_map_adjusted = np.where(
            index_map == 65535,
            65535,
            index_map - index_offset
        )

        # Use ravel() to flatten the 2D indices and then reshape back
        #indices_flat = index_map_adjusted.ravel()
        #subsolar_lat_values = subsolar_lat[indices_flat].reshape(index_map_adjusted.shape)
        #total_subsolar_lat[start_row-1:end_row, start_col-1:end_col] = subsolar_lat_values

        #subsolar_lon_values = subsolar_lon[indices_flat].reshape(index_map_adjusted.shape)
        #total_subsolar_lon[start_row-1:end_row, start_col-1:end_col] = subsolar_lon_values

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

               
        #print(f"subsolar_lat length: {len(subsolar_lat)} subsolar_lon length: {len(subsolar_lon)}")
        #print(f"subsat_lat length: {len(subsat_lat)} subsat_lon length: {len(subsat_lon)}")
        #print(f"index_map_adjusted min: {index_map_adjusted.min()} max: {index_map_adjusted.max()}")
        
        nc.close()
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        continue

# Normalize to 0-255 range for image output
# Exclude zeros (fill values) from normalization
#valid_data_R = composite_R[composite_R > 0]
#valid_data_G = composite_G[composite_G > 0]
#valid_data_B = composite_B[composite_B > 0]
# Before normalization, check the data range
print(f"composite_R raw range min : {composite_R.min()}")
print(f"composite_R raw range max : {composite_R.max()}")
print(f"composite_R unique values: {len(np.unique(composite_R))}")

'''
valid_data_R = valid_data_R[~np.isnan(valid_data_R)]
valid_data_G = valid_data_G[~np.isnan(valid_data_G)]
valid_data_B = valid_data_B[~np.isnan(valid_data_B)]
'''
valid_data_R = composite_R[composite_R != 65535.0]
valid_data_G = composite_G[composite_G != 65535.0]
valid_data_B = composite_B[composite_B != 65535.0]

print(f"Valid data R len: {len(valid_data_R)}")
print(f"Composit_R len: {len(composite_R.flatten())}")

if len(valid_data_R) > 0:
    vmin = np.percentile(valid_data_R, 0)  # Use 1st percentile for better contrast
    vmax = np.percentile(valid_data_R, 100)  # Use 99th percentile
    
    composite_normalized_R = np.where(
        composite_R == 65535.0,
        0,
        np.clip((composite_R - vmin) / (vmax - vmin) * 255, 0, 255)
    )
else:
    composite_normalized_R = composite_R

if len(valid_data_G) > 0:
    vmin = np.percentile(valid_data_G, 0)  # Use 1st percentile for better contrast
    vmax = np.percentile(valid_data_G, 100)  # Use 99th percentile
    
    composite_normalized_G = np.where(
        composite_G == 65535.0,
        0,
        np.clip((composite_G - vmin) / (vmax - vmin) * 255, 0, 255)
    )
else:
    composite_normalized_G = composite_G

if len(valid_data_B) > 0:
    vmin = np.percentile(valid_data_B, 1)  # Use 1st percentile for better contrast
    vmax = np.percentile(valid_data_B, 100)  # Use 99th percentile
    
    composite_normalized_B = np.where(
        composite_B == 65535.0,
        0,
        np.clip((composite_B - vmin) / (vmax - vmin) * 255, 0, 255)
    )
else:
    composite_normalized_B = composite_B

# Create alpha channel (0 for fill values 65535, 255 for valid data)
alpha_channel = np.where(composite_R == 65535.0, 0, 255).astype(np.uint8)


# Convert to uint8
composite_uint8_R = composite_normalized_R.astype(np.uint8)
composite_uint8_G = composite_normalized_G.astype(np.uint8)
composite_uint8_B = composite_normalized_B.astype(np.uint8)
#print(f"After normalization unique values: {len(np.unique(composite_uint8_R[~np.isnan(composite_normalized_R)]))}")

# Create matplotlib figure
#plt.figure(figsize=(12, 10))
#plt.imshow(composite_uint8_R, cmap='gray')
#plt.colorbar(label='Radiance (8-bit)')
#plt.title('FCI VIS_04 Composite')
#plt.xlabel('Column')
#plt.ylabel('Row')
#plt.tight_layout()
#plt.savefig('fci_composite_plot.png', dpi=100)
#plt.show()
'''
plt.figure(figsize=(10, 6))
plt.hist(composite_uint8_R.flatten(), bins=256, range=(0, 256), edgecolor='black')
plt.xlabel('Pixel Value (0-255)')
plt.ylabel('Frequency')
plt.title('Histogram of FCI VIS_04 Composite (Red Channel)')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fci_composite_histogram.png', dpi=100)
plt.show()
'''

'''
rgb_array = np.stack([composite_uint8_R, composite_uint8_G, composite_uint8_B, alpha_channel], axis=-1)
# Save as PNG
img = Image.fromarray(rgb_array, mode='RGBA')
imgout = img.transpose(method=Image.Transpose.FLIP_TOP_BOTTOM)
imgout.save(output_file)

print(f"Composite image saved to: {output_file}")
print(f"Image shape: {composite_uint8_R.shape}")
print(f"Value range: {composite_uint8_R.min()} - {composite_uint8_R.max()}")
'''
#print(f"Region total_subsolar_lon [5000:5050, 5000:5050]:\n{total_subsolar_lon[5000:5010, 5000:5010]}")
#print(f"Region total_subsolar_lat [5000:5050, 5000:5050]:\n{total_subsolar_lat[5000:5010, 5000:5010]}")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nScript execution time: {elapsed_time:.2f} seconds")
print(f"Script execution time: {elapsed_time/60:.2f} minutes")


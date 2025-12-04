import netCDF4 as nc
import numpy as np
import hdf5plugin

# Open the file
#thefile = "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-FDHSI-FD--CHK-BODY--DIS-NC4E_C_EUMT_20251129113703_IDPFI_OPE_20251129113415_20251129113455_N_JLS_O_0070_0020.nc"
thefile = "W_XX2.nc"
fci_file = nc.Dataset(thefile, "r")

# Get index_map as raw integers
index_map_var = fci_file.groups['data'].groups['vis_06'].groups['measured'].variables['index_map']
index_map_var.set_auto_scale(False)
index_map = index_map_var[:].astype(np.int32)
#index_map_att = fci_file.groups['data'].groups['vis_06'].groups['measured'].variables['index_map']

# List all attributes
print("Attributes of index_map_var:")
for attr in index_map_var.ncattrs():
    print(f"  {attr}: {index_map_var.getncattr(attr)}")

# Get fill value
fill_value = index_map_var.getncattr("_FillValue")
fill_value1 = index_map_var.get_fill_value()
index_offset = fci_file.variables['index_offset'][:]
# Create mask for valid pixels
#valid_mask = (index_map != fill_value)

# Use only valid indices
#valid_idx = index_map[valid_mask]


# Get the geometric parameter vectors
subsolar_lat = fci_file.groups['state'].groups['celestial'].variables['subsolar_latitude'][:]
subsolar_lon = fci_file.groups['state'].groups['celestial'].variables['subsolar_longitude'][:]
subsat_lat = fci_file.groups['state'].groups['platform'].variables['subsatellite_latitude'][:]
subsat_lon = fci_file.groups['state'].groups['platform'].variables['subsatellite_longitude'][:]

# For pixel at row 100, column 200
row, col = 0, 350
idx = index_map[row, col]
print(f"idx = {idx}")

# Check dimensions
print(f"index_map shape: {index_map.shape}")
print(f"index_map min: {np.min(index_map[index_map != fill_value])}")
print(f"index_map max: {np.max(index_map[index_map != fill_value])}")
print(f"unique valid indices:  {np.max(index_map[index_map != fill_value]) - np.min(index_map[index_map != fill_value]) + 1}")
print(f"index offset: {index_offset}")
print(f"subsolar_lat length: {len(subsolar_lat)}")
print(f"subsolar_lon length: {len(subsolar_lon)}")
print(f"Fill value: {fill_value}")

#idxn = idx - np.min(index_map[index_map != fill_value])
idxn  = idx - index_offset

# Get the geometric parameters for this specific pixel
pixel_subsolar_lat = subsolar_lat[idxn]
pixel_subsolar_lon = subsolar_lon[idxn]
pixel_subsat_lat = subsat_lat[idxn]
pixel_subsat_lon = subsat_lon[idxn]
print(f"Pixel at row {row}, col {col} has:")
print(f"  Sub-solar Latitude: {pixel_subsolar_lat}")
print(f"  Sub-solar Longitude: {pixel_subsolar_lon}")
print(f"  Sub-satellite Latitude: {pixel_subsat_lat}")
print(f"  Sub-satellite Longitude: {pixel_subsat_lon}")

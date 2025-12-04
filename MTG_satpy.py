#!/usr/bin/env python3
import hdf5plugin
from satpy.scene import Scene
from satpy import find_files_and_readers
import glob
import time

# Start timer
start_time = time.time()

# define path to FCI test data folder
path_to_data = '/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_00[0-9][0-9].nc'
files = glob.glob('/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_00[0-9][0-9].nc')
# find files and assign the FCI reader
#files = find_files_and_readers(filter_parameters=path_to_data, reader='fci_l1c_nc')

# create an FCI scene from the selected files
scn = Scene(filenames=files, reader='fci_l1c_nc')

# print available dataset names for this scene (e.g. 'vis_04', 'vis_05','ir_38',...)
#print(scn.available_dataset_names())

# print available composite names for this scene (e.g. 'natural_color', 'airmass', 'convection',...)
print(scn.available_composite_names())

# load the datasets/composites of interest
#scn.load(['ndvi_hybrid_green_fully_sunzencorrected','true_color'], upper_right_corner='NE')
scn.load(['true_color_raw'], upper_right_corner='NE')
# note: the data inside the FCI files is stored upside down. The upper_right_corner='NE' argument
# flips it automatically in upright position.

# you can access the values of a dataset as a Numpy array with
#vis_04_values = scn['vis_04'].values

# resample the scene to a specified area (e.g. "eurol1" for Europe in 1km resolution)
#scn_resampled = scn.resample("eurol", resampler='nearest', radius_of_influence=5000)

# save the resampled dataset/composite to disk
scn.save_dataset("true_color_raw", filename='./true_color_raw.png')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nScript execution time: {elapsed_time:.2f} seconds")
print(f"Script execution time: {elapsed_time/60:.2f} minutes")

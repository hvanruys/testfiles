import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

print('Numpy version: ', np.__version__)
print('Xarray version: ', xr.__version__)

# Load a sample NetCDF file (replace 'sample.nc' with your file)
#ds = xr.open_dataset('W_XX2.nc')
#print(ds['index_offset'])
#for att in ds.attrs:
#   print(f"Attributes: {att} = {ds.attrs[att]}")
    #print(ds[var])
    #print(f"Attributes: {ds[var].attrs}")

#    ds = xr.open_dataset('xarray-data-master/air_temperature.nc')
#    da=ds['air']
#    with xr.set_options(display_style="text"):
#        print(ds.air.attrs)
#print(ds.coords['index'])
#print(ds.attrs)
#print(ds)

'''
x = np.linspace(-np.pi, np.pi, 19)
f = np.sin(x)
da_f = xr.DataArray(f, dims=['x'], coords={'x': x}) #, name='sine_wave', attrs={'description': 'Sine wave values'})
print(da_f[:10])
da_f.sel(x=slice(0, np.pi)).plot()
#da_f.plot(marker='o')
print(da_f.sel(x=0))
da_g = da_f**2 + 1
(da_f * da_g).sel(x=slice(-1, 1)).plot(marker='o')
da_f.plot()
plt.show()
'''

ds = xr.open_dataset('tutorial-data/sst/NOAA_NCDC_ERSST_v3b_SST-1960.nc')
xr.set_options(display_style="text")
print("opening NOAA dataset")
sst = ds['sst']
print(sst)
print(sst.sel(lon=180, lat=0, time='1960-01-16', method='nearest'))
#sst.sel(time='1960-06-15').plot(vmin=-2, vmax=30)
#sst.sel(lon=180).transpose().plot()
#sst.sel(lon=180, lat=40).plot()

plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
ax.coastlines()

sst[0].plot(transform=ccrs.PlateCarree(), vmin=-2, vmax=30,
            cbar_kwargs={'shrink': 0.4})

plt.show()
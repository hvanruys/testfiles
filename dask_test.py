import dask.array as da
x = da.random.random((100000, 100000), chunks=(1000, 1000))
result = x.mean().compute()
print(result)
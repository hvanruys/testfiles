#!/usr/bin/env python3
import sys
import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    print("Missing dependency: pip install netCDF4")
    sys.exit(1)

def summarize(path):
    ds = Dataset(path, "r")
    print("File:", path)
    # Global attributes
    print("\nGlobal attributes:")
    for k in ds.ncattrs():
        print(f"  {k} = {getattr(ds, k)}")

    # Dimensions
    print("\nDimensions:")
    for name, dim in ds.dimensions.items():
        print(f"  {name} = {len(dim)} (unlimited={dim.isunlimited()})")

    # Variables
    print("\nVariables:")
    for name, var in ds.variables.items():
        print(f"  {name} dtype={var.dtype} shape={var.shape}")
        attrs = {a: getattr(var, a) for a in var.ncattrs()}
        if attrs:
            print(f"    attrs: {attrs}")
        # small data preview
        try:
            if var.size == 0:
                print("    (empty)")
            elif var.size <= 20:
                print("    data:", var[:].tolist())
            else:
                # show a single-element sample (first index for each dim)
                idx = tuple(0 for _ in var.shape) if var.shape else (0,)
                sample = var[idx] if var.shape else var[:]
                print("    sample[0,...]:", np.array(sample).tolist())
        except Exception:
            pass

    # Optional: plot first 2D variable if matplotlib available
    try:
        import matplotlib.pyplot as plt
        two_d = next((v for v in ds.variables.values() if getattr(v, "ndim", 0) == 2), None)
        if two_d is not None:
            arr = two_d[:]
            plt.imshow(arr, origin="lower")
            plt.title(f"{two_d.name} (shape={arr.shape})")
            plt.colorbar()
            print("\nDisplaying plot for first 2D variable:", two_d.name)
            plt.show()
    except Exception:
        pass

    ds.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        thefile = "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-FDHSI-FD--CHK-BODY--DIS-NC4E_C_EUMT_20251129113703_IDPFI_OPE_20251129113415_20251129113455_N_JLS_O_0070_0020.nc"
    else:
    	thefile = sys.argv[1]
    	
    summarize(thefile)

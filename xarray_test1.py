import numpy as np
from matplotlib import pyplot as plt

im = np.copy(plt.imread('xarray-data-master/guinea_bissau.JPG'))
plt.imshow(im)
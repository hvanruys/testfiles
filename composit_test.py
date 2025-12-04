import numpy as np

composite_R = np.zeros((10, 10), dtype=np.float32)
rad_1 = np.random.rand(5, 10) * 1000
rad_2 = np.full((5, 10), 35)
composite_R[:5, :] = rad_1
composite_R[5:, :] = rad_2
print(composite_R)
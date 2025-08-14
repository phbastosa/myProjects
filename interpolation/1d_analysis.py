import numpy as np
import matplotlib.pyplot as plt

xi =-10.0
xf = 10.0

original_Nx = 1001
rescaled_Nx = 51

original_x = np.linspace(xi, xf, original_Nx)
rescaled_x = np.linspace(xi, xf, rescaled_Nx)

original_F = np.fromfile("original1d.bin", dtype = np.float32, count = original_Nx)
rescaled_F = np.fromfile("rescaled1d.bin", dtype = np.float32, count = rescaled_Nx)

cubic_F = np.fromfile("cubic1d.bin", dtype = np.float32, count = original_Nx)
linear_F = np.fromfile("linear1d.bin", dtype = np.float32, count = original_Nx)

plt.figure(figsize = (15,7))

plt.plot(original_x, original_F, "-k", label = "original function")
plt.plot(rescaled_x, rescaled_F, "or", label = "rescaled function")
plt.plot(original_x, cubic_F, "-g", label = "cubic interpolated function")
plt.plot(original_x, linear_F, "-b", label = "linear interpolated function")

plt.legend(loc = "lower right")
plt.xlim([xi-0.5, xf+0.5])
plt.tight_layout()
plt.show()

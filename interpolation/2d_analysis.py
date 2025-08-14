import numpy as np
import matplotlib.pyplot as plt

xi =-10.0
xf = 10.0
yi =-10.0
yf = 10.0

original_Nx = 1001
original_Ny = 1001

rescaled_Nx = 51
rescaled_Ny = 51

original_N = original_Nx*original_Ny
rescaled_N = rescaled_Nx*rescaled_Ny

xp, yp = np.meshgrid(np.linspace(xi, xf, rescaled_Nx), np.linspace(yi, yf, rescaled_Ny))

original_F = np.fromfile("original2d.bin", dtype = np.float32, count = original_N).reshape([original_Nx, original_Ny], order = "F")
rescaled_F = np.fromfile("rescaled2d.bin", dtype = np.float32, count = rescaled_N).reshape([rescaled_Nx, rescaled_Ny], order = "F")

cubic_F =  np.fromfile("cubic2d.bin", dtype = np.float32, count = original_N).reshape([original_Nx, original_Ny], order = "F")
linear_F =  np.fromfile("linear2d.bin", dtype = np.float32, count = original_N).reshape([original_Nx, original_Ny], order = "F")

scale = 0.01*np.max(np.abs(original_F))

diff_cubic = np.zeros_like(original_F)
diff_linear = np.zeros_like(original_F)

diff_cubic[25:-25,25:-25] = original_F[25:-25,25:-25] - cubic_F[25:-25,25:-25]
diff_linear[25:-25,25:-25] = original_F[25:-25,25:-25] - linear_F[25:-25,25:-25]

plt.figure(figsize = (15,7))

plt.subplot(231)
plt.imshow(original_F, aspect = "auto", cmap = "jet", extent = [xi,xf,yf,yi])
plt.colorbar()

plt.scatter(xp, yp, color = "black", s = 1)
plt.title("Original function")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(232)
plt.imshow(linear_F, aspect = "auto", cmap = "jet", extent = [xi,xf,yf,yi])
plt.colorbar()
plt.title("Linear interpolated function")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(233)
plt.imshow(cubic_F, aspect = "auto", cmap = "jet", extent = [xi,xf,yf,yi])
plt.colorbar()
plt.title("Cubic interpolated function")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(234)
plt.imshow(rescaled_F, aspect = "auto", cmap = "jet", extent = [xi,xf,yf,yi])
plt.colorbar()
plt.title("Rescaled function")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(235)
plt.imshow(diff_linear, aspect = "auto", cmap = "jet", vmax = scale, vmin = -scale, extent = [-10,10,10,-10])
plt.colorbar()
plt.title("Linear interpolation error")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(236)
plt.imshow(diff_cubic, aspect = "auto", cmap = "jet", vmax = scale, vmin = -scale, extent = [-10,10,10,-10])
plt.colorbar()
plt.title("Cubic interpolation error")
plt.xlabel("x")
plt.ylabel("y")

plt.tight_layout()
plt.show()

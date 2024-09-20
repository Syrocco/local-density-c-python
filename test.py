import cell_list
from functions2 import Data
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

def findMax(averageDensity, bins = 1000, plot = False):
    phiM = np.mean(averageDensity)
    y, x = np.histogram(average_density, bins = bins, density = True)
    x = 0.5*(x[1:] + x[:-1])
    Y = savgol_filter(y, 25, 3)

    xb = x[x < phiM]
    yb = Y[x < phiM]
    
    xa = x[x >= phiM]
    ya = Y[x >= phiM]
    
    # Find argmax of y below and above phiM
    liquid = xb[np.argmax(yb)]
    solid = xa[np.argmax(ya)]
    if plot:
        plt.plot(x, y)
        plt.plot(x, Y)
        plt.vlines([liquid, solid], 0, np.max(Y))
    return liquid, solid
    
    
data = Data("N_20000dtnoise_1.000res_1.000gamma_0.010T_1.000phi_0.600000rat_0.500vo_4.000ao_50.000delta_0.022Lx_323.604Ly_323.604q_0.000v_1.dump")

Lx_val = data.Lx
Ly_val = data.Ly
dump = data.dump
dump.jump_to_frame(dump.nframes - 1)

x_vals = dump.get_atompropf("x")
y_vals = dump.get_atompropf("y")
r_vals = dump.get_atompropf("radius")

cell_size_val = 2
rad = 4


cell_list_instance = cell_list.CellList(x_vals, y_vals, r_vals, Lx_val, Ly_val, cell_size_val)

average_density = cell_list_instance.get_rad_density_array(1001000, rad)
findMax(average_density, plot = True)
plt.xlabel(r"$\phi$")
plt.ylabel(r"$proba$")
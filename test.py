import cell_list
import random
from functions2 import Data
import matplotlib.pyplot as plt
from diptest import dip
import numpy as np


data = Data("/hdd/Documents/mips/data/LAMMPS/third tests/N_15000phi_0.6vo_600T_0m_0.457142857142857.dumpM")

Lx_val = data.Lx
Ly_val = data.Ly
dump = data.dump
dump.jump_to_frame(dump.nframes - 1)
x_vals = dump.get_atompropf("x")
y_vals = dump.get_atompropf("y")

cell_size_val = 2.5
rad = 2


cell_list_instance = cell_list.CellList(x_vals, y_vals, Lx_val, Ly_val, cell_size_val)

average_density = cell_list_instance.get_density_array(10000000, rad, 0.56123102415)
plt.hist(average_density, bins = 1000, density = True)
print(dip(average_density))
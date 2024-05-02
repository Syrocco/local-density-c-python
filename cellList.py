import numpy as np
import matplotlib.pyplot as plt
from functions2 import Data


class CellList:
    def toCellX(self, x):
        return min(int(x / self.cell_size), self.num_cells_x - 1)
    
    
    
    def toCellY(self, y):
        return min(int(y / self.cell_size), self.num_cells_y - 1)
    
    
    def pbcX(self, dx):
        if (dx >= self.Lx/2):
            return dx - self.Lx
        elif (dx <= -self.Lx/2):
            return dx + self.Lx
        return dx
    
    
    
    def pbcY(self, dy):
        if (dy >= self.Ly/2):
            return dy - self.Ly
        elif (dy <= -self.Ly/2):
            return dy + self.Ly
        return dy
    
    
    def pbcCellX(self, a):
        if a < 0:
            return a + self.num_cells_x
        elif a >= self.num_cells_x:
            return a - self.num_cells_x
        return a
    
    
    def pbcCellY(self, a):
        if a < 0:
            return a + self.num_cells_y
        elif a >= self.num_cells_y:
            return a - self.num_cells_y
        return a
    
    def __init__(self, x, y, Lx, Ly, cell_size):
        self.x = x
        self.y = y
        self.Lx = Lx
        self.Ly = Ly
        self.cell_size = cell_size
        self.cells = {}
        self.num_cells_x = int(self.Lx / self.cell_size)
        self.num_cells_y = int(self.Ly / self.cell_size)
        self.particleToCell = []
    
        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                self.cells[(i, j)] = []
        
       
        for idx, (px, py) in enumerate(zip(self.x, self.y)):
            cell_x = self.toCellX(px)
            cell_y = self.toCellY(py)
            self.cells[(cell_x, cell_y)].append(idx)
            self.particleToCell.append([cell_x, cell_y])
    
    def get_particles_in_cell(self, cell_x, cell_y):
        return self.cells.get((cell_x, cell_y), [])

    def get_density(self, x, y, rad, radius = 0.56123102415):
        mini = (rad - radius)**2
        maxi = (rad + radius)**2
        def f(X):
            if X < mini:
                return 1
            elif X > maxi:
                return 0
            else:
                return (1-0)/(mini - maxi)*(X - maxi)
                
        numCell = int(rad/self.cell_size) + 2
        print(numCell)
        cell_x_i, cell_y_i = self.toCellX(x), self.toCellY(y)
        c = 0
        XX = []
        YY = []
        l = []
        for i in range(-numCell, numCell):
            for j in range(-numCell, numCell):
                for k in self.cells[self.pbcCellX(cell_x_i + i), self.pbcCellY(cell_y_i + j)]:
                    if k == 4445:
                        print(i, j, self.pbcCellX(cell_x_i + i), self.pbcCellY(cell_y_i + j))
                    value = f(self.pbcX(x - self.x[k])**2 + self.pbcY(y - self.y[k])**2)
                    if value > 0:
                        l.append(k)
                        c += value

        return c/(np.pi*rad**2)*np.pi*radius**2, XX, YY, np.array(sorted(l))
        
data = Data("/hdd/Documents/mips/data/LAMMPS/second tests/N_5000phi_0.6vo_40T_0m_0.264444444444444.dumpM")

Lx = data.Lx
Ly = data.Ly
data.dump.jump_to_frame(data.dump.nframes - 1)
x = data.dump.get_atompropf("x")
y = data.dump.get_atompropf("y")


cellList = CellList(x, y, Lx, Ly, 10)
A = []
np.random.seed(18)
for i in range(1):
    print(i)
    X = np.random.random()*Lx
    Y = np.random.random()*Ly
    a, b, c, l = cellList.get_density(X, Y, 15)
    A.append(a)
  
plt.hist(A, bins = 50)
if 0:
    plt.figure(figsize = (10, 10*Ly/Lx))
    plt.scatter(x, y, s = 8)
    plt.scatter(b, c, s = 3)
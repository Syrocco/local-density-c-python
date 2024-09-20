from parserL import Dump
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib
from matplotlib.patches import Circle
import os
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 18})
from scipy.signal import correlate
from scipy.spatial.distance import pdist
import re


class Data():
    def __init__(self, loc, dump = True, intruderLenght = 0.1, jump = 1):
        self.loc = loc
        if dump:
            try:
                self.dump = Dump(loc)
            except:
                pass
        x = loc.split("/")[-1]
        if (x[-1] == 'p') or (x[-1] == "4"):
            self.lammps = False
            self.extract_variables()
            self.L = self.Lx
            try:
                self.sizeratio = round(self.dump.get_atompropf("radius")[-1], 4)
                self.Nb = np.sum(self.dump.get_atompropi("type"))
                self.Ns = self.N - self.Nb
                self.compo = self.Nb/self.N
            except:
                pass
            self.phi = float(x.split("phi_")[1][:8])
            #self.phi = round((np.pi*(self.sizeratio**2*self.Ns + 1*self.Nb))/(self.Ly*self.Lx), 9)
            try:
                path = loc.split(".dump")[0] + '.thermo'
                with open(path, 'r') as file:
                    first_line = file.readline().strip()
                    column_names = first_line.split()
                if column_names[0] == 't':


                    try:
                        data = np.loadtxt(path, skiprows = 1)
                    except:
                        data = np.loadtxt(path, skip_footer = 1, skiprows = 1)

                    if data.ndim == 1:
                        data = data.reshape(1, -1)

                    for i, col_name in enumerate(column_names):
                        setattr(self, col_name, data[:, i])
            except:
                pass

            else:
                try:
                    try:
                        temp = np.loadtxt(loc.split(".dump")[0] + '.thermo')
                    except:
                        temp = np.genfromtxt(loc.split(".dump")[0] + '.thermo', skip_footer = 1)
                    if np.shape(temp)[1] == 4:
                        self.t = temp[:, 0]
                        self.ncol = temp[:, 1]
                        self.E = temp[:, 2]
                        self.p = temp[:, 3]
                    elif np.shape(temp)[1] == 3:
                        self.t = temp[:, 0]
                        self.E = temp[:, 1]
                        self.p = temp[:, 3]

                    #print(temp)
                    else:
                        self.t = temp[:, 0]
                        self.ss = temp[:, 1]
                        self.bb = temp[:, 2]
                        self.sb = temp[:, 3]
                        self.freq = temp[:, 4]
                        self.E = temp[:, 5]
                        self.p = temp[:, 6]
                        """self.Ebb = temp[:, 7]
                    self.Esb = temp[:, 8]
                    self.vr = temp[:, 9]"""
                    try:
                        self.px = temp[:, 7]
                    except:
                        pass
                except (IndexError, FileNotFoundError, ValueError):
                    try:

                        temp = np.loadtxt(loc.split(".dump")[0] + '.thermo')
                        self.t = temp[0]
                        self.ss = temp[1]
                        self.bb = temp[2]
                        self.sb = temp[3]
                        self.freq = temp[4]
                        self.E = temp[5]
                        self.p = temp[6]
                        """self.Ebb = temp[:, 7]
                        self.Esb = temp[:, 8]
                        self.vr = temp[:, 9]"""

                    except (IndexError, FileNotFoundError, ValueError):
                        try:
                            temp = np.loadtxt(loc.split(".dump")[0] + '.thermo')
                            self.t = temp[:, 0]
                            self.p = temp[:, 1]
                            self.E = temp[:, 2]
                        except (IndexError, FileNotFoundError, ValueError):
                            pass
        elif x[-1] == 'L':
            self.lammps = True
            self.phi = float(x.split("phi_")[1].split('freq')[0])
            self.freq =  float(x.split("freq_")[1].split('T')[0])
            self.T = float(x.split("T_")[1].split('h_')[0])
            self.h = float(x.split("h_")[1].split('.d')[0])
            try:
                self.L = self.dump.get_boxx()[0]
                self.N = self.dump.get_natoms()
            except:
                pass
            """
            self.Nb = np.sum(self.dump.get_atompropi("type") == 1)
            self.Ns = self.N - self.Nb
            self.compo = self.Nb/self.N"""
            #typ = self.dump.get_atompropi("type")
            #rad = self.dump.get_atompropf("radius")
            #rb = rad[np.where(typ == 1)[0][0]]
            """
            try:
                rs = rad[np.where(typ == 2)[0][0]]
                self.sizeratio = round(rs/rb, 4)
            except Exception:
                self.sizeratio = 1"""
     
            try:

                temp = loc.split(".dump")[0] + '.thermo'           
                f = open(temp, 'r')
                lines = f.readlines()

               	#idx = lines.index('Step Atoms c_EkXY \n')
               	for i, line in enumerate(lines):

                    if "Atoms" in line:
                        idx = i
                        break

                #idx = lines.index('   Step        Atoms        c_EkXY    \n')
                #idx = lines.index('Step Atoms c_EkXY \n')
                ticks = []
                E = []
                D = []

                for i in range(idx + 1, len(lines)):
                    a = lines[i].split(" ")

                    try:
                        try:
                            a.remove("\n")
                        except:
                            pass
                        while '' in a:
                            a.remove('')

                        ticks.append(float(a[0]))
                        string = a[2]
                        string2 = a[-1]
                        if "\n" in string:
                            E.append(float(string[:-1]))
                        else:
                            E.append(float(string))
                        if "\n" in string2:
                            D.append(float(string2[:-1]))
                        else:
                            D.append(float(string2))

                    except:
                        break
                f.close()
                #print(E)
                self.ticks = np.array(ticks)
                self.E = np.array(E)
                self.D = np.array(D)
            except:
                pass
        elif x[-1] == 'M':
            self.extract_variables()
            self.Lx = self.dump.get_boxx()[0]
            self.Ly = self.dump.get_boxy()[0]
        try:
            data = np.loadtxt(loc.split(".dump")[0] + '.intruder')
            time = data[:, 0]
            time -= time[0]
            vx = data[:, 3]
            vy = data[:, 4]
            vx -= np.mean(vx)
            vy -= np.mean(vy)
            corr = autoV(vx, vy)
            self.time = np.copy(time[:int(len(time)*intruderLenght)][::jump])
            self.corr = np.copy(corr[:int(len(time)*intruderLenght)][::jump])
        except:
            try:
                data = np.genfromtxt(loc.split(".dump")[0] + '.intruder', skip_footer = 1)
                time = data[:, 0]
                time -= time[0]
                vx = data[:, 3]
                vy = data[:, 4]
                vx -= np.mean(vx)
                vy -= np.mean(vy)
                corr = autoV(vx, vy)
                self.time = np.copy(time[:int(len(time)*intruderLenght)][::jump])
                self.corr = np.copy(corr[:int(len(time)*intruderLenght)][::jump])
            except:
                pass
    def extract_variables(self):
        # Define the pattern to match variable names and values
        pattern = r'(\w+?)_(\d*\.?\d+)'

        # Find all matches
        matches = re.findall(pattern, self.loc)

        # Create a dictionary to store the variable names and values
        variables = {}

        # Iterate over matches and populate the dictionary
        for match in matches:
            variable_name = match[0]
            variable_value = float(match[1]) if '.' in match[1] else int(match[1])
            variables[variable_name] = variable_value

        # Assign values to variables with the same names
        for name, value in variables.items():
            setattr(self, name, value)



def autoV(vx, vy):
    return 0.5*(autocorr(vx)/np.var(vx) + autocorr(vy)/np.var(vy))

def autocorr(x):
    n = x.size
    result = correlate(x, x, mode='full')
    return result[result.size//2:]/np.arange(n, 0, -1)

def image(data):
    color = ["royalblue", "goldenrod"]
    dump = data.dump
    dump.jump_to_frame(dump.nframes - 2)
    x = dump.get_atompropf('x')
    y = dump.get_atompropf('y')
    r = dump.get_atompropf('radius')
    c = dump.get_atompropi("type")
    fig, ax = plt.subplots(figsize=(12,12))
    plt.xlim(0, data.L)
    plt.ylim(0, data.L)
    for i in range(len(x)):
        ax.add_patch(Circle((x[i], y[i]), r[i], zorder=10, color = color[c[i]]))
    
           
def dataArray(loc, dump = True, more = False, intruderLenght = 0.1, jump = 1):

    files = glob.glob(loc + "*.dump*")
    allDump = []
    for i in range(len(files)):
        print(str(i) + "/" + str(len(files))) 
        try:
            allDump.append(Data(files[i], dump, intruderLenght, jump))

        except IndexError:
            print(i)

    phi = []
    sizeratio = []
    delta = []
    amp = []
    height = []
    for i in range(len(allDump)):
        try:
            if not (allDump[i].phi in phi):
                phi.append(allDump[i].phi)
        except AttributeError:
            pass
        try:
            if not (allDump[i].sizeratio in sizeratio):
                sizeratio.append(allDump[i].sizeratio)
        except AttributeError:
            pass
        try:
            if not (allDump[i].T in amp):
                amp.append(allDump[i].T)
        except AttributeError:
            pass
        try:
            if not (allDump[i].h in height):
                height.append(allDump[i].h)
        except AttributeError:
            pass
        try:
            if not (allDump[i].delta in delta):
                delta.append(allDump[i].delta)
        except AttributeError:
            pass

    phi = np.sort(np.array(phi))
    sizeratio = np.sort(np.array(sizeratio))
    #delta = np.sort(np.array(delta))
    #gamma =  np.array([a.gamma for a in allDump])
    #res = np.array([a.res for a in allDump])  

    amp = np.sort(np.array(amp))
    height = np.sort(np.array(height))
    
    if more:
        return allDump, amp, height     
    return allDump


def energy(data, total = True):
    dump = data.dump
    n = dump.nframes - 2
    if data.lammps:
        big = np.abs(dump.get_atompropi("type") - 2).astype(bool)
        m = big*6.54498e-05 + np.logical_not(big)*6.37063e-06         
    else:
        m = dump.get_atompropf('m')
    if total:
        E = np.zeros(n - 1)
        t = np.zeros(n - 1)
        for i in range(0, n - 1):
            #print(i,"/", n - 1)
            dump.jump_to_frame(i)
            vx = dump.get_atompropf('vx')
            vy = dump.get_atompropf('vy')
            E[i] = 0.5*np.sum((m*(vx**2+vy**2)))/data.N
            t[i] = dump.get_timestep()
        return t, E
    else:
        if data.lammps:
            big = np.abs(dump.get_atompropi("type") - 2).astype(bool)
        else: 
            big = dump.get_atompropi("type").astype(bool)
        small = np.logical_not(big)
        Eb = np.zeros(n - 1)
        Es = np.zeros(n - 1)
        t = np.zeros(n - 1)
        for i in range(n - 1):
            #print(i,"/", n - 1)
            dump.jump_to_frame(i)
            vx = dump.get_atompropf('vx')
            vy = dump.get_atompropf('vy')
            Eb[i] = 0.5*np.sum(big*m*(vx**2+vy**2))/data.Nb
            Es[i] = 0.5*np.sum(small*m*(vx**2+vy**2))/data.Ns
            t[i] = dump.get_timestep()
        return t, Eb, Es
    

def pbc(a, L):
    if a < -L/2:
        return a + L
    if a > L/2:
        return a - L
    return a

def q4(data, fra = 0, num = 4):
    dump = data.dump
    L = data.L
    n = dump.nframes
    if fra > 0:
        a = [fra]
    else:
        a = range(n)
        q4 = np.zeros(n)
        t = np.zeros(n)
    for frame in a:
        dump.jump_to_frame(frame)
        if data.lammps:
            big = np.abs(dump.get_atompropi("type") - 2).astype(bool)
        else:
            big = dump.get_atompropi("type").astype(bool)
        x = dump.get_atompropf('x')[big]
        y = dump.get_atompropf('y')[big]
        radius = dump.get_atompropf('radius')[big][0]
        M = int(L/(3*radius))
        grid = -np.ones((M, M, 10)).astype(int)
        counter = np.zeros((M, M)).astype(int)
        for i in range(len(x)):
            X = int(M*(x[i]%L)/L)
            Y = int(M*(y[i]%L)/L)
            grid[X, Y, counter[X, Y]] = i
            counter[X, Y] += 1
            
        Q4 = np.zeros(len(x)).astype(complex)
        for X in range(M):
            for Y in range(M):
                for i in range(counter[X, Y]):
                    s = 0
                    xref = x[grid[X, Y, i]]
                    yref = y[grid[X, Y, i]]
                    for j in [-1, 0, 1]:
                        for k in [-1, 0, 1]:
                            Xbox = (X + j)%M
                            Ybox = (Y + k)%M
                            xtest = x[grid[Xbox, Ybox, :counter[Xbox, Ybox]]]
                            ytest = y[grid[Xbox, Ybox, :counter[Xbox, Ybox]]]
                            for z in range(len(xtest)):
                                dx = pbc(xref - xtest[z], L)
                                dy = pbc(yref - ytest[z], L)
    
                                if (dx != 0) and (np.sqrt(dx*dx + dy*dy) < 2.1*radius):
    
                                    Q4[grid[X, Y, i]] += np.exp(1j*num*np.arctan2(dy, dx))
                                    s += 1
                    #print(s)
                    if (s != 0):
                        Q4[grid[X, Y, i]]/=s
                    #print("--------", Q4[grid[X, Y, i]])
        if fra > 0:
            return np.mean(np.abs(Q4))
        q4[frame] = np.mean(np.abs(Q4))
        t[frame] = dump.get_timestep()
    return t, q4

def q4LxLy(data, fra = 0, num = 4):
    dump = data.dump
    Lx = data.Lx
    Ly = data.Ly
    n = dump.nframes
    if fra > 0:
        a = [fra]
    else:
        a = range(n)
        q4 = np.zeros(n)
        t = np.zeros(n)
    for frame in a:
        dump.jump_to_frame(frame)
        if data.lammps:
            big = np.abs(dump.get_atompropi("type") - 2).astype(bool)
        else:
            big = dump.get_atompropi("type").astype(bool)
        x = dump.get_atompropf('x')[big]
        y = dump.get_atompropf('y')[big]
        radius = dump.get_atompropf('radius')[big][0]
        Mx = int(Lx/(3*radius))
        My = int(Ly/(3*radius))
        grid = -np.ones((Mx, My, 10)).astype(int)
        counter = np.zeros((Mx, My)).astype(int)
        for i in range(len(x)):
            X = int(Mx*(x[i]%Lx)/Lx)
            Y = int(My*(y[i]%Ly)/Ly)
            grid[X, Y, counter[X, Y]] = i
            counter[X, Y] += 1
            
        Q4 = np.zeros(len(x)).astype(complex)
        for X in range(Mx):
            for Y in range(My):
                for i in range(counter[X, Y]):
                    s = 0
                    xref = x[grid[X, Y, i]]
                    yref = y[grid[X, Y, i]]
                    for j in [-1, 0, 1]:
                        for k in [-1, 0, 1]:
                            Xbox = (X + j)%Mx
                            Ybox = (Y + k)%My
                            xtest = x[grid[Xbox, Ybox, :counter[Xbox, Ybox]]]
                            ytest = y[grid[Xbox, Ybox, :counter[Xbox, Ybox]]]
                            for z in range(len(xtest)):
                                dx = pbc(xref - xtest[z], Lx)
                                dy = pbc(yref - ytest[z], Ly)
    
                                if (dx != 0) and (np.sqrt(dx*dx + dy*dy) < 2.1*radius):
    
                                    Q4[grid[X, Y, i]] += np.exp(1j*num*np.arctan2(dy, dx))
                                    s += 1
                    #print(s)
                    if (s != 0):
                        Q4[grid[X, Y, i]]/=s
                    #print("--------", Q4[grid[X, Y, i]])
        if fra > 0:
            return np.mean(np.abs(Q4))
        q4[frame] = np.mean(np.abs(Q4))
        t[frame] = dump.get_timestep()
    return t, q4


def Sscalar(data, frame = None):
    dump = data.dump
    N = 30
    K = 3000
    q = np.linspace(0.00001, K, N)
    if frame == None:
        dump.jump_to_frame(dump.nframes - 1)
    else:
        dump.jump_to_frame(frame)
    x = dump.get_atompropf('x')
    y = dump.get_atompropf('y')
    XY = np.abs(pdist(np.vstack((x, y)).T))
    Sres = np.zeros(N)
    for i, qiter in enumerate(q):
        print(i)
        Sres[i] = 1 + np.sum(np.sin(qiter*XY)/(qiter*XY))/data.N
    return q*dump.get_atompropf("radius")[0]/(2*np.pi), Sres






def S(data, frame = None, info = True):
    def structureFactor(x, y, kx, ky):
        dot = x*kx + y*ky
        return ((np.sum(np.cos(dot), axis = 0))**2 + (np.sum(np.sin(dot), axis = 0))**2)/len(x)

    dump = data.dump
    N = 300
    halfN = int(N/2)
    K = 15
    k = np.linspace(-K, K, N)
    if frame == None:
        dump.jump_to_frame(dump.nframes - 1)
    else:
        dump.jump_to_frame(frame)
    if data.lammps:
        big = np.abs(dump.get_atompropi("type") - 2).astype(bool)
    else:
        big = dump.get_atompropi("type").astype(bool)
    x = dump.get_atompropf('x')[big]
    y = dump.get_atompropf('y')[big]
    Lpost = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            Lpost[i, j] = structureFactor(x, y, k[i], k[j])
    Lpost[Lpost < 3] = 0
    Lpost = np.log(Lpost**2 + 1)
    Lpost = (Lpost - 1)/2.54
    p = 10
    Lpost[halfN - p:halfN + p, halfN - p:halfN + p] = Lpost[halfN + p + 1, halfN + p]
    #Lpost[np.where(Lpost>40)] = 0
    if info:
        x = dump.get_atompropf('x')
        y = dump.get_atompropf('y')

        r = dump.get_atompropf('radius')
        c = dump.get_atompropi("type")
        color = ["royalblue", "goldenrod"]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18,9))
        ax1.imshow(Lpost, extent=[-K, K, -K, K], cmap = plt.get_cmap('binary'))
        ax1.title.set_text(rf' $L =$ {data.L}, $\phi =$ {data.phi}, $N_b/N$ = {data.Nb/data.N} and $r_s/r_b = $ {data.sizeratio} ')

        for i in range(len(x)):
            ax2.add_patch(Circle((x[i], y[i]), r[i], zorder=10, color = color[c[i]]))
        plt.xlim(0, data.L)
        plt.ylim(0, data.L)
    else:
        return Lpost

def g2(data, frame, ave = 1):

    
    dump = data.dump
    L =  data.L

    N = data.N
    M = int(N/10)
    h = np.zeros(M)
    for j in range(ave):
        print(j)
        k = 0
        dump.jump_to_frame(frame - j)
        x = dump.get_atompropf('x')[dump.get_atompropi("type").astype(bool)]
        y = dump.get_atompropf('y')[dump.get_atompropi("type").astype(bool)]
        rij = np.zeros(int(N*(N-1)/2))
        for i in range(N):
            for j in range(i):
                dx = pbc(x[i] - x[j], L)
                dy = pbc(y[i] - y[j], L)
                rij[k] = np.sqrt(dx*dx + dy*dy)
                k += 1 
        htemp, r = np.histogram(rij, M, range = (0, L/10))
        h += htemp
    h = h/ave
    r = 0.5*(r[1:] + r[:-1])
    g = L*L*h/(np.pi*N*(N-1)*r*(r[1]-r[0]))
    return r, g
    

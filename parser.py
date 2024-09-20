"""
Python dump parser
"""

import numpy as np


class Dump:
    def __init__(self,path):
        """
        Scans through file to count frames and save their position.
        Uses first frame to get the list of available items.
        """
        self.path=path
        self.fo=open(path,"r")
        self.framepos=[]
        self.items=[]
        self.nframes=0
        self.curframe=0
        #Count frames and save their position
        line=self.fo.readline()
        while line!="":
            if line[:len("ITEM: TIMESTEP")]=="ITEM: TIMESTEP":
                self.nframes+=1
                self.framepos.append(self.fo.tell()-len(line))
            line=self.fo.readline()
        #Use first frame to determine items
        self.fo.seek(self.framepos[0],0)
        self.items.append("TIMESTEP")
        self.fo.readline()  #Skip first TIMESTEP item
        line=self.fo.readline()
        while line[:len("ITEM: TIMESTEP")]!="ITEM: TIMESTEP" and line!="":
            if line[:len("ITEM:")]=="ITEM:":
                self.items.append(line.split()[1])
            line=self.fo.readline()
        self.fo.seek(self.framepos[0],0)
    def __del__(self):
        try:
            self.fo.close()
        except:
            pass
    def next_frame(self):
        if self.curframe+1>=self.nframes:
            print("[Dump.next_frame]: No frame %d in a dump containing %d frames."%(self.curframe+1,self.nframes))
            return
        self.curframe+=1
        self.fo.seek(self.framepos[self.curframe],0)
    def jump_to_frame(self,frame):
        if frame>=self.nframes:
            print("[Dump.jump_to_frame]: No frame %d in a dump containing %d frames."%(frame,self.nframes))
            return
        self.curframe=frame
        self.fo.seek(self.framepos[self.curframe],0)
    def reach(self,item):
        """
        Sets file cursor at the beginning of the ITEM line in current frame.
        """
        if item not in self.items:
            print("[Dump.reach]: No item %s in this dump."%item)
            return
        self.fo.seek(self.framepos[self.curframe])
        line=self.fo.readline()
        while line[:len("ITEM: %s"%item)]!="ITEM: %s"%item:
            line=self.fo.readline()
        self.fo.seek(self.fo.tell()-len(line),0)
    def readcols(self,n,span,skiprows=0):
        """
        Returns the span values forming the nth column below the cursor, as a list of strings. Skiprows allows to offset the column with respect to current cursor position.
        """
        vals=np.zeros(span)
        for i in range(skiprows):
            self.fo.readline()
        for i in range(span):
            vals[i]=self.fo.readline().split()[n]
        return vals
    def readcoli(self,n,span,skiprows=0):
        """
        Returns the span values forming the nth column below the cursor, as a list of ints. Skiprows allows to offset the column with respect to current cursor position.
        """
        return self.readcols(n,span,skiprows).astype(int)
    def readcolf(self,n,span,skiprows=0):
        """
        Returns the span values forming the nth column below the cursor, as a list of floats. Skiprows allows to offset the column with respect to current cursor position.
        """
        return self.readcols(n,span,skiprows).astype(float)
    def get_timestep(self):
        self.reach("TIMESTEP")
        return self.readcolf(0,1,1)[0]
    def get_natoms(self):
        self.reach("NUMBER")
        return self.readcoli(0,1,1)[0]
    def get_boxx(self):
        self.reach("BOX")
        header=self.fo.readline().split()
        return self.readcolf(1,1,0)[0],header[3]
    def get_boxy(self):
        self.reach("BOX")
        header=self.fo.readline().split()
        return self.readcolf(1,1,1)[0],header[4]
    def get_atomheader(self):
        """
        Returns the header of available atom properties as a list of strings.
        """
        self.reach("ATOMS")
        return self.fo.readline().split()[2:]
    def get_atompropi(self,prop):
        """
        Returns an array of natoms values of atom property prop for the current frame.
        Returned values are ints. Use get_atompropf for float properties.
        """
        natoms=self.get_natoms()
        header=self.get_atomheader()
        if prop not in header:
            print("[Dump.get_atomprop]: No property %s in this dump."%prop)
            return
        index=header.index(prop)
        return self.readcoli(index,natoms,0)
    def get_atompropf(self,prop):
        """
        Returns an array of natoms values of atom property prop for the current frame.
        Returned values are float. Use get_atompropi for int properties.
        """
        natoms=self.get_natoms()
        header=self.get_atomheader()
        if prop not in header:
            print("[Dump.get_atomprop]: No property %s in this dump."%prop)
            return
        index=header.index(prop)
        return self.readcolf(index,natoms,0)


import numpy as np
import gena.readwrite as rw
import h5py
import gena.geometric as gmt
import gena.projection as proj
import gena.astro as ast

#func
def mid_ary(ary):
    return (ary[:-1]+ary[1:])/2

def names(self):
    return [attr for attr in dir(self) if not attr.startswith('__') and not callable(getattr(self, attr))]

class prob():
    def __init__(self,name):
        #prob name, A or B
        self.name = name
        #numer of events
        self.N = []
        self.E = []
        #timestamp UTC
        self.t = []
        #prob localtion in J2000
        self.s = []
        #q matrix 
        self.rot = []
        #moon location   
        self.moon = []
        #event vector :p(oint) s(tart) - p(oint) e(nd)
        self.vec = []
        #satellite z axis pointing at GSE Y = 0 plane
        self.pz = []
        #ena class
        self.ena = self.ena()
        #fov class
        self.fov = self.fov()

    def generate(self,csv):
        N,ti,s,Rot,vec,E = rw.preload(csv)
        ts = [t.timestamp() for t in ti]
        moon = ast.moon(ti)
        self.N = N
        self.E = E
        self.t = ts
        self.s = s
        self.rot = Rot
        self.moon = moon
        self.vec = vec
        self.pz = proj.projects(ts,s,np.array([0,0,-1]),Rot,moon)
        #self means prob_inst arg that passed to calcuate
        self.fov.calculate(self)
        self.ena.calculate(self)

    def get(self):
        t = self.t 
        s = self.s 
        r = self.rot
        m = self.moon
        return [t,s,r,m]


    #def dump(self, fn, group_name=None):
    #    # If no group name is provided, use the class name
    #    if group_name is None:
    #        group_name = self.__class__.__name__
    #    # Open the HDF5 file in append mode
    #    with h5py.File(fn, 'a') as f:
    #        # If the group already exists, delete it to avoid duplicates
    #        if group_name in f:
    #            del f[group_name]
    #        group = f.create_group(group_name)
    #        # Iterate over the attributes of the instance
    #        for attr_name in dir(self):
    #            # Skip special and callable attributes
    #            if attr_name.startswith('__') or callable(getattr(self, attr_name)):
    #                continue
    #            attr_value = getattr(self, attr_name)
    #            # Check if the attribute has a dump method (e.g., another class instance)
    #            if hasattr(attr_value, 'dump'):
    #                # Recursively call the dump method
    #                print(fn, f'{group_name}/{attr_name}')
    #                #attr_value.dump(fn, f'{group_name}/{attr_name}')
    #            else:
    #                # Save the attribute value as a dataset in the group
    #                #group.create_dataset(attr_name, data=attr_value)
    #                print('exceptions')
    def dump(self,fn):
        rw.dumph5(fn,self.name,['N','t','s','rot','moon','vec','pz'],\
                [self.N, self.t,self.s,self.rot,self.moon,self.vec,self.pz])       
        self.ena.dump(fn)
        self.fov.dump(fn)

    def load(self,fn):
        [self.N, self.t , self.s, self.rot, self.moon, self.vec, self.pz ] = \
            rw.loadh5(fn,self.name,['N','t','s','rot','moon','vec','pz'])       
        self.ena.load(fn)
        self.fov.load(fn)


    class ena:
        def __init__(self,t=[],xz=[],pt=[],E=[],r=[]):
            self.name='ena'
            self.t = t
            self.E = E
            self.xz = xz
            self.r = r
            self.pt = pt
            self.E_index = np.argsort(E)

        def calculate(self,prob_inst):
            t,s,r,m = prob_inst.get()
            self.t = prob_inst.t
            self.E = prob_inst.E
            vec = prob_inst.vec
            pname = prob_inst.name
            x,z = proj.projection(t,s,vec,r,m,pname)
            self.xz = [x,z]
            distance,pl,ps = gmt.vec_2_rphi(vec)
            self.pt = [pl,ps]
            self.r = distance
            self.E_index = np.argsort(self.E)

        def sub_E(self,El,Eu):
            E = self.E
            E_index = self.E_index
            sorted_E = E[E_index]
            sub_idx = (sorted_E >= El) & (sorted_E <= Eu)
            sub_ori_idx = E_index[sub_idx]
            sub_ena = prob.ena(self.t[sub_ori_idx],self.xz[:,sub_ori_idx],self.pt[:,sub_ori_idx],self.E[sub_ori_idx],self.r[sub_ori_idx])
            return sub_ena

        def sub_T(self,tl,tu):
            t = self.t
            sub_idx = (t >= tl) & (t <= tu)
            sub_ena = prob.ena(self.t[sub_idx],self.xz[:,sub_idx],self.pt[:,sub_idx],self.E[sub_idx],self.r[sub_idx])
            return sub_ena

        #slice ENA sequence according to their energy and time
        def slice(self,E_intervals,T_intervals):
            subsets=[]
            for tl, tu in zip(T_intervals[:-1],T_intervals[1:]):
                for el, eu in zip(E_intervals[:-1],E_intervals[1:]):
                    subsets.append(self.sub_E(el,eu).sub_T(tl,tu))
            return subsets

        def hist_xz(self,extend=[-30,30,-20,20],dx=2.5):
            xl,xr,zd,zu = extend
            grid_x = np.arange(xl,xr,dx)
            grid_z = np.arange(zd,zu,dx)
            x = self.xz[0]
            z = self.xz[1]
            H_ena = np.histogram2d(x,z,[grid_x,grid_z])
            Z,X = np.meshgrid(mid_ary(grid_z),mid_ary(grid_x))
            return X,Z,H_ena[0]

        def hist_pt(self,extend=[-30,30,-20,20],dx=2.5):
            xl,xr,zd,zu = extend
            grid_x = np.arange(xl,xr,dx)
            grid_z = np.arange(zd,zu,dx)
            x = self.pt[0]
            z = self.pt[1]
            H_ena = np.histogram2d(x,z,[grid_x,grid_z])
            Z,X = np.meshgrid(mid_ary(grid_z),mid_ary(grid_x))
            return X,Z,H_ena[0]


        def load(self,fn):
            [self.t,self.r,self.xz,self.pt,self.E,self.E_index] = \
                    rw.loadh5(fn,'ena',['t','r','xz','pt','e','e_idx'])

        def dump(self,fn):
            rw.dumph5(fn,'ena',['t','r','xz','pt','e','e_idx'],\
                    [self.t,self.r,self.xz,self.pt,self.E,self.E_index])
        #def dump(self, fn):
        #    # Call the dump method of the parent class
        #    self.prob_inst.dump(fn)


        #FOV class
        #           x
        #           |
        #4 ur       |
        ######################### 1 ur
        #           #           # 
        #           #           # 
        #           #           # 
        #           #           # 
        #     A     z      B    # ------y
        #           #           # 
        #           #           # 
        #           #           # 
        #3 dl       #           # 
        ######################### 2 dr
    class fov:
        def __init__(self,t=[],p=[],v=[]):
            self.name='fov'
            #init from file
            self.t = t
            self.p = p
            self.v = v

        #load from file
        def load(self,fn):
            self.p,self.v = rw.loadh5(fn,'fov',['p','v'])

        def dump(self,fn):
            rw.dumph5(fn,'fov',['p','v'],[self.p,self.v])

        #def dump(self, fn):
        #    # Call the dump method of the parent class
        #    self.prob_inst.dump(fn)

        def calculate(self,prob_inst,short_fov=5,long_fov=22.5):
            t,s,r,m = prob_inst.get()
            pname = prob_inst.name
            self.t = prob_inst.t
            p,v = proj.get_vertex(t,s,r,m,short_fov,long_fov,pname)
            #self.pz = P[0]
            self.p = p
            self.v = v

        def draw(self,ax,idx,c='k'):
            x=[]
            y=[]
            for p in self.v:
                x.append(p[0][idx])
                y.append(p[1][idx])
                ax.scatter(p[0][idx],p[1][idx],marker='+',color='k')
            #append the first one to get circ closed
            x.append(self.v[0][0][idx])
            y.append(self.v[0][1][idx])
            ax.plot(x,y,color=c)
            ax.scatter(self.p[0][idx],self.p[1][idx],marker='+',color='k')
            return 0



        #fov see objects
        #from j2000 to satellite, inv Rot

if __name__ == "main":
    print('unit test here')

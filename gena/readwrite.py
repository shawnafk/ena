import numpy as np
import pandas as pd
import gena.timefunc as tf
import gena.geometric as gmt
def load_csv(fn,column='',Nrange=''):
    #0: bj time  string
    #1: satllite time (s)
    #2: ms
    #3: angle_z_sun
    #4 - 6: q1 q2 q3 
    #7 - 9: mj2000 x y z (m) /1e3 to km
    #10 11: x1 y1 (mm) /1e6 to km
    #12 13: x2 y2 (mm)
    #14: distance
    #15: speed
    #16: species
    #17: energy
    out = []
    df = pd.read_csv(fn,encoding='gbk',header=0).iloc[:,1:]
    if Nrange == '':
        N = df.shape[0]
        Nrange=np.arange(N)
    if column == '':
        for i in range(18):
            out.append(np.array(df.iloc[Nrange,i]))
        return out
    else:
        for i in column:
            out.append(np.array(df.iloc[Nrange,i]))
        return out

def preload(fn):
    #load needed
    SatelliteTime,ms,q1,q2,q3,sx,sy,sz,x1,y1,x2,y2,E = load_csv(fn,[1,2,4,5,6,7,8,9,10,11,12,13,17])
    SatelliteTime = SatelliteTime + ms/1000
    
    N = SatelliteTime.shape[0]
    
    DateTime = [tf.SatelliteTime_2_Datetime(st) for st in SatelliteTime]
    
    Rot = gmt.q2rot(q1,q2,q3,N)
    
    x1,y1,x2,y2 = np.array([x1,y1,x2,y2]) * 1e-6
    
    z1,z2 = np.zeros(N), np.ones(N)*50 * 1e-6
    #event location at frame centered at MPC 0, first layer
    
    loc_s = np.array([x1,y1,z1]).transpose()
    
    loc_e = np.array([x2,y2,z2]).transpose()
    
    vec = loc_e - loc_s
    
    #s: loc_satellite_mj2000, km
    s = np.array([sx,sy,sz]).transpose()/1e3
    
    return  N,DateTime,s,Rot,vec,E



def dumph5(name,grp,varname,var):
    with h5py.File(name, 'a') as f:
        if grp in f:
            group = f[grp]
        else:
            #group like dict
            group = f.create_group(grp)
        #dateset is array
        for vn,v in zip(varname,var):
            if vn in f[grp]:
                del f[grp][vn]
                group.create_dataset(vn, data=v)
            else:
                group.create_dataset(vn, data=v)

#only load one eveal return
def loadh5(name,grp,varnames):
    with h5py.File(name, 'r') as f:
        dat_s = []
        for varname in varnames:
            dat = f[grp][varname]
            dim = len(dat.shape)
            if dim == 0:
                dat_s.append(dat[()])
            else:
                dat_s.append(dat[...])
        return dat_s


import h5py
def output(name,dict):
    for grp_name, grp_value in dict.items():
        for dataset_name, dataset_value in grp_value.items():
            dumph5(name,grp_name,dataset_name,dataset_value)
    return 0

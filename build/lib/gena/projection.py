import numpy as np
import gena.timefunc as tf
import gena.geometric as gmt
import gena.coord as cod
from gena.const import *
# ------------------------------ projection --------------------
#from statllite location at moon j2000 and staellite vec
#to Earth xz 
def projects(ts,s,v,rot,moon):
    s_ej2000 = cod.mj2000_2_ej2000(s,moon)
    s_gse = cod.ej2000_2_gse(ts,s_ej2000)
    #N * 3 * 3 @ 3 * 1 * 1  = N * 3 * 1 . => N * 3
    v_j2000 = (rot @ v[...,np.newaxis])[...,0]
    v_gse = cod.ej2000_2_gse(ts,v_j2000)
    t = (0 - s_gse[:,1])/v_gse[:,1]
    x_proj = s_gse[:,0] + t * v_gse[:,0]
    z_proj = s_gse[:,2] + t * v_gse[:,2]
    return x_proj/RE,z_proj/RE

def projection(ts,s,v,rot,moon,prob=''):
    if prob=='A':
        rp = gmt.rot_x(-5)
        newv = (rp@v[...,np.newaxis])[:,:,0]
    elif prob=='B':
        rp = gmt.rot_x(5)  
        newv = (rp@v[...,np.newaxis])[:,:,0]
    return projects(ts,s,newv,rot,moon)   


#method 2, use angle directly to get new vec
#angle in satellite fov
#-z is 0,0 
def projection2(ts,s,v,rot,moon,prob=''):
    if prob=='A':
        newv = gmt.rot_vec(v,+5,'short')
    elif prob=='B':
        newv = gmt.rot_vec(v,-5,'short')
    return projects(ts,s,newv,rot,moon)   

import magnetosphere as mag
def verify_angle(r,s,moon_bary,t):
    sun_ej2000 = mag.myget_body(t,body='sun')
    ts = [_.timestamp() for _ in t]
    sun_gse =  cod.ej2000_2_gse(ts,sun_ej2000)
    sun_vec_gse = sun_gse - cod.mj2000_2_ej2000(s,moon_bary)
    reversed_zaxis_j2000 = r @ np.array([0,0,-1])
    reversed_zaxis_gse = cod.ej2000_2_gse(ts,reversed_zaxis_j2000)
    ac = gmt.get_ang(reversed_zaxis_gse,sun_vec_gse)
    #verify angle between -z and sun vec
    return ac

def get_boundary_vec(short_fov,long_fov,vec_axis):
    R_u = gmt.rot_x(short_fov)
    R_d = gmt.rot_x(-short_fov)
    R_l = gmt.rot_y(long_fov)
    R_r = gmt.rot_y(-long_fov)
    vec_1 = R_r @ (R_u @ vec_axis)
    vec_2 = R_r @ (R_d @ vec_axis)
    vec_3 = R_l @ (R_d @ vec_axis)
    vec_4 = R_l @ (R_u @ vec_axis)
    return vec_1,vec_2,vec_3,vec_4

def get_vertex(t,s,rot,moon,short_fov,long_fov,prob):
        #satellite pointing
        v_axis_s = np.array([0,0,-1])
        if prob == 'A':
            rp = gmt.rot_x(-5)
            v = (rp@v_axis_s)
        elif prob== 'B':
            rp = gmt.rot_x(+5)  
            v = (rp@v_axis_s)
        v1,v2,v3,v4 = get_boundary_vec(short_fov,long_fov,v)
        #va1,va2,va3,va4 = get_boundary_vec(short_fov,long_fov,va)
        #vb1,vb2,vb3,vb4 = get_boundary_vec(short_fov,long_fov,vb)

        p = projects(t,s,v,rot,moon)
        #pb = projects(t,s,va,rot,moon)
        vtx = []
        for v in [v1,v2,v3,v4]:
            vtx.append(projects(t,s,v,rot,moon))
        #vtxb = []
        #for v in [vb1,vb2,vb3,vb4]:
        #    vtxb.append(projects(t,s,v,rot,moon))
        return p,vtx

def get_boundary_vec2(v,short_fov,long_fov):
        vec_1 = gmt.rot_vec(gmt.rot_vec(v,short_fov,'short'),long_fov,'long')
        vec_2 = gmt.rot_vec(gmt.rot_vec(v,short_fov,'short'),-long_fov,'long')
        vec_3 = gmt.rot_vec(gmt.rot_vec(v,-short_fov,'short'),-long_fov,'long')
        vec_4 = gmt.rot_vec(gmt.rot_vec(v,-short_fov,'short'),long_fov,'long')
        return vec_1,vec_2,vec_3,vec_4

def get_vertex2(t,s,rot,moon,prob='',short_fov=5,long_fov=22.5):
    #satellite pointing
    v_axis_s = np.array([0,0,-1])
    ps,pl = 0,0
    if prob=='A':
        ps = 5
        v = gmt.rphi_2_vec(ps,pl,1)
    elif prob=='B':
        ps = -5
        v = gmt.rphi_2_vec(ps,pl,1)
    v1,v2,v3,v4 = get_boundary_vec2(v,short_fov,long_fov)
    #ps = projects(t,s,v_axis_s,rot,moon)
    p = projects(t,s,v,rot,moon)
    vtx = []
    for v in [v1,v2,v3,v4]:
        vtx.append(projects(t,s,v,rot,moon))
    return p,vtx

if __name__ == "__main__":
    print("unit test here")

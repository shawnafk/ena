import gena.readwrite as rw
import gena.magnetosphere as obj
import numpy as np
from gena import GENA
import gena.timefunc as tf
from matplotlib import pyplot as plt
import time

prob = GENA.prob("B")
prob.load('test.h5')

ena = prob.ena
fov = prob.fov

fo = './out/'
#1. verify pointing and fov region
f,ax = plt.subplots()
obj.Earth().draw(ax)
ps = prob.pz
pc = fov.p
ax.plot(ps[0,:],ps[1,:])
ax.plot(pc[0,:],pc[1,:])
fov.draw(ax,0)
fov.draw(ax,-1)
ax.set_aspect('equal')
plt.savefig(fo+"z_pointing_fov"+'.pdf',bbox_inches='tight')

#verify particle counts flux v.s. time
#count particle in given intervals


ts = ena.t
bins = np.arange(ts[0], ts[-1], 1*60)
hist, bin_edges = np.histogram(ts, bins=bins)
t = [tf.Timestamp_2_Datetime(_) for _ in bin_edges]
#set time grid as hours
f,ax=plt.subplots()
ax.plot(t[1:],hist)
plt.savefig(fo+"ena_flux"+'.pdf',bbox_inches='tight')

#verify all ena events
f,ax=plt.subplots()
x,z,h = ena.hist_xz()
ax.pcolormesh(x,z,h,vmin=0,cmap='jet')

ts = time.time()
UT = tf.Str_2_Datetime('2024-10-11 00:00:00').timestamp()
SWDP = 10
DST = -300
BETA = 5
IMF_BY = 4
IMF_BZ = -10
M_MS=5


earth_obj = obj.Earth()
bs_obj = obj.Bowshock(IMF_BZ,BETA,M_MS,SWDP)
mp_obj = obj.Magnetopause(IMF_BZ,SWDP)
field_obj =  obj.T96_fieldlines(UT,SWDP, DST, IMF_BY, IMF_BZ)
earth_obj.draw(ax)
bs_obj.draw(ax)
mp_obj.draw(ax)
field_obj.draw(ax)
fov.draw(ax,0)
fov.draw(ax,-1)
te = time.time()
time_difference = (te - ts )*1000
print("time cost:", time_difference, "ms")
ax.set_aspect('equal')
ax.set_xlim(-30,30)
ax.set_ylim(-14,14)
plt.savefig(fo+"projection_all"+'.pdf',bbox_inches='tight')

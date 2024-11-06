import numpy as np
import gena.timefunc as tf
import gena.readwrite as rw
import gena.magnetosphere as obj
from gena import GENA

prob = GENA.prob("B")
prob.load('test.h5')

fn = 'test.h5'
#function lized
ena = prob.ena
fov = prob.fov

#given utc start time
t0 = tf.Str_2_Datetime('2024-10-10 23:00:00')
#time step
delta_t = 900
#steps
n = 21
#time_idx = tf.Get_Idx(t0,dt,n,ts)
target_time = t0.timestamp() + np.arange(0,n)*delta_t

E_min = 0
E_max = np.inf
ena_slices = ena.slice([E_min,E_max],target_time)

UT = tf.Str_2_Datetime('2024-10-11 00:00:00').timestamp()
SWDP = 10
DST = -300
BETA = 5
IMF_BY = 4
IMF_BZ = -10
M_MS=5
earth_obj = obj.Earth(color='w')
bs_obj = obj.Bowshock(IMF_BZ,BETA,M_MS,SWDP,color='w')
mp_obj = obj.Magnetopause(IMF_BZ,SWDP,color='w')
field_obj =  obj.T96_fieldlines(UT,SWDP, DST, IMF_BY, IMF_BZ,color='w')


from matplotlib import pyplot as plt
f,ax = plt.subplots(4,5)
ax = ax.ravel()
for i in range(n-1):
    x,z,h = ena_slices[i].hist_xz()
    im = ax[i].pcolormesh(x,z,h,vmin=0,vmax=20,cmap='jet')
    ax[i].set_xlim(-30,30)
    ax[i].set_ylim(-15,15)
    ax[i].set_aspect('equal')
    date = tf.Datetime_2_Str(tf.Timestamp_2_Datetime(target_time[i]))
    ax[i].set_title(date[5:-5])
 
    fov.draw(ax[i],0)
    earth_obj.draw(ax[i])
    bs_obj.draw(ax[i])
    mp_obj.draw(ax[i])
    field_obj.draw(ax[i])
cb = f.colorbar(im,ax=ax,location='right')
f.set_size_inches(15,15/4*3)
plt.savefig("time_ord"+'.pdf',bbox_inches='tight')

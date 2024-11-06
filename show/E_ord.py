import numpy as np
import gena.timefunc as tf
import gena.readwrite as rw
from gena import GENA
import gena.magnetosphere as obj


prob = GENA.prob("B")
prob.load('test.h5')

ena = prob.ena
fov = prob.fov

T_min = ena.t[0]
T_max = ena.t[-1]
target_E = [0,1,5,10,50,100,np.inf]
ena_slices = ena.slice(target_E,[T_min,T_max])


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


from matplotlib import pyplot as plt
f,ax = plt.subplots(2,3)
ax = ax.ravel()
for i in range(6):
    x,z,h = ena_slices[i].hist_xz()
    im = ax[i].pcolormesh(x,z,h,vmin=0,cmap='jet')
    ax[i].set_xlim(-30,30)
    ax[i].set_ylim(-15,15)
    ax[i].set_aspect('equal')
    ax[i].set_title(str(target_E[i])+'-'+str(target_E[i+1]))
    
    fov.draw(ax[i],0)
    earth_obj.draw(ax[i])
    bs_obj.draw(ax[i])
    mp_obj.draw(ax[i])
    field_obj.draw(ax[i])
    cb = f.colorbar(im,ax=ax[i],location='right')
f.set_size_inches(10,10/4*3)
plt.savefig("E_ord"+'.pdf',bbox_inches='tight')

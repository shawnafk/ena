#astroobj
import numpy as np
import gena.timefunc as tf
from jplephem.spk import SPK
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
kernel = SPK.open(current_dir+'/Ephemeris/de432s.bsp')
from astropy.coordinates import  get_body

#given timestamp get moon localtion in J2000
def moon(t):
    t_iso = [tf.Datetime_2_ISOT(_) for _ in t]
    t_jd =  [_.jd for _ in t_iso]
    moon_bary = kernel[3,301].compute(t_jd).transpose()
    return moon_bary

def myget_body(t,body='moon'):
    #m to km
    isot = [tf.Datetime_2_ISOT(_) for _ in t]
    AU = 149597870700/1e3
    N = isot.shape[0]
    body_ej2000_xyz = np.zeros((N,3))
    #t = Time(iso_format_string, scale='utc')
    #t = Time(iso_format_string, format='isot', scale='utc')
    for i in range(N):
        body_ej2000 = get_body(body,isot[i]).represent_as('cartesian')
        body_ej2000_xyz[i,0] = np.double(body_ej2000.x)*AU
        body_ej2000_xyz[i,1] = np.double(body_ej2000.y)*AU
        body_ej2000_xyz[i,2] = np.double(body_ej2000.z)*AU
    #scj2000 = sc.transform_to(GCRS(obstime=Time("J2000"))).represent_as('cartesian')
    return body_ej2000_xyz

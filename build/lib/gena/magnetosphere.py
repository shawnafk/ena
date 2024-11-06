import numpy as np
import gena.timefunc as tf
import gena.coord as cod
from geopack import geopack

class Magnetosphere:
    def __init__(self, color='k', lw=2):
        self.lines = []
        self.color = color
        self.lw = lw

    def draw(self, ax):
        for l in self.lines:
            ax.plot(l[0], l[1], color=self.color, lw=self.lw)

class Earth(Magnetosphere):
    def __init__(self, color='k', lw=1):
        super().__init__(color, lw)
        t = np.linspace(0, 2 * np.pi, 128)
        x = np.cos(t)
        y = np.sin(t)
        self.lines = [(x, y)]


class Bowshock(Magnetosphere):
    def __init__(self,B_z,beta,M_ms,D_p,color='k', lw=1):
    #def __init__(self,B_z,beta,M_ms,D_p,rlim,color,lw):
        super().__init__(color,lw)
        self.B_z  = B_z
        self.beta = beta
        self.M_ms = M_ms
        self.D_p  = D_p
        x,r = self.BS()
        #rbs = np.sqrt(bs[0,:]**2 + bs[1,:]**2)
        #idx_bs = rbs < rlim
        #self.lines = [bs[0,idx_bs],bs[1,idx_bs]]
        self.lines = [(x,r)]
    
    def BS_args(self):
        B_z  = self.B_z 
        beta = self.beta
        M_ms = self.M_ms
        D_p  = self.D_p 
 
        a_1=11.1266
        a_2=0.0010
        a_3=-0.0005
        a_4=2.5966
        a_5=0.8182
        a_6=-0.0170
        a_7=-0.0122
        a_8=1.3007
        a_9=-0.0049
        a_10=-0.0328
        a_11=6.047
        a_12=1.029
        a_13=0.0231
        a_14=-0.002
        if B_z >0:
            r_0=a_1 * (1+a_2*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1 / a_11)
            alpha = a_5 * (1+a_13 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
        else:
            r_0=a_1 * (1+a_3*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1/a_11)
            alpha=a_5*(1+a_6 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
        epsi = a_12
        return r_0, alpha, epsi

    def BS(self):
        r0,alp,epsi = self.BS_args()
        theta = np.linspace(-np.pi*0.9,np.pi*0.9,128)
        r = r0 * ( (1+epsi)/(1+epsi * np.cos(theta) ) ) ** alp
        x = r * np.cos(theta)
        R = r * np.sin(theta)
        return x,R
    
class Magnetopause(Magnetosphere):
    def __init__(self,B_z,D_p,color='k', lw=1):
        super().__init__(color,lw)
        self.B_z  = B_z
        self.D_p  = D_p
        x,r = self.MP()
        self.lines = [(x,r)]
    
    def MP_args(self):
        B_z  = self.B_z 
        D_p  = self.D_p 
        a_1 = 11.646
        a_2 = 0.216
        a_3 = 0.122
        a_4 = 6.215
        a_5 = 0.578
        a_6 = -0.009
        a_7 = 0.012
        if B_z > 0:
            r0  = a_1 * D_p**(-1/a_4)
        elif B_z < 0 and B_z > -8:
            r0 = (a_1 + a_2* B_z) * D_p**(-1/a_4)
        else:
            r0 =( a_1 + 8*a_3 - 8 * a_2 + a_3 * B_z ) * D_p**(-1/a_4)
        alpha = (a_5 + a_6 * B_z) * (1+ a_7 * D_p)
        return r0, alpha

    def MP(self):
        r0,alp = self.MP_args()
        theta = np.linspace(-np.pi*0.9,np.pi*0.9,128)
        r = r0 * ( 2/(1+np.cos(theta) ) ) ** alp
        x = r * np.cos(theta)
        R = r * np.sin(theta)
        return x,R

class T96_fieldlines(Magnetosphere):
    def __init__(self,ut,swdp,dst,by,bz,color='k', lw=1):
        super().__init__(color,lw)
        self.ut  = ut
        self.swdp = swdp
        self.dst = dst
        self.by  = by
        self.bz  = bz
        #rbs = np.sqrt(bs[0,:]**2 + bs[1,:]**2)
        #idx_bs = rbs < rlim
        #self.lines = [bs[0,idx_bs],bs[1,idx_bs]]
        self.lines = self.trace_lines()
    #solar wind pressure pdyn (nanopascals)
    #dst (nanotesla)
    #byimf (nanotesla)
    #bzimf (nanotesla)
    def trace_lines(self):
        ut = self.ut
        pm=[self.swdp,self.dst,self.by,self.bz]
        lines = []
        #x,z = earth_surface(20)
        geopack.recalc(ut)
        #mp
        for dir in [-1,1]:
            #magnetofront
            for _x in np.arange(2,5):
                res =  geopack.trace(_x,0,0,dir,rlim=20,parmod=pm,exname='t96')
                lines.append([res[3],res[5]])
            #magnetotail
            for _z in np.arange(2,4):
                #x=-5
                res =  geopack.trace(-10,0,_z,dir,rlim=25,parmod=pm,exname='t96')
                lines.append([res[3],res[5]])
            for _z in np.arange(-4,-2):
                #x=-5
                res =  geopack.trace(-10,0,_z,dir,rlim=25,parmod=pm,exname='t96')
                lines.append([res[3],res[5]])
            #upper cusp
            for _x in np.arange(-1,3):
                res =  geopack.trace(_x,0,2,dir,rlim=10,parmod=pm,exname='t96')
                lines.append([res[3],res[5]])
            #lower cusp
            for _x in np.arange(-1,3):
                res =  geopack.trace(_x,0,-2,dir,rlim=10,parmod=pm,exname='t96')
                lines.append([res[3],res[5]])
            #gse_lines=[]
            #for l in lines:
            #    N = len(l[0])
            #    xyz = np.array([l[0],np.zeros(N),l[1]])
            #    gse_xyz = cod.gsm_2_gse(np.ones(N)*ut,xyz.transpose())
            #    gse_lines.append([gse_xyz[:,0],gse_xyz[:,2]])
            #return gse_lines
            return lines

if __name__ == "__main__":
    print("unit test here")
    import timefunc as tf
    UT = tf.Str_2_Datetime('2024-10-11 00:00:00').timestamp()
    SWDP = 10
    DST = -300
    BETA = 5
    IMF_BY = 4
    IMF_BZ = -10
    M_MS=5
    earth_obj = Earth()
    bs_obj = Bowshock(IMF_BZ,BETA,M_MS,SWDP)
    mp_obj = Magnetopause(IMF_BZ,SWDP)
    field_obj =  T96_fieldlines(UT,SWDP, DST, IMF_BY, IMF_BZ)

    from matplotlib import pyplot as plt
    f,ax=plt.subplots()
    earth_obj.draw(ax)
    bs_obj.draw(ax)
    mp_obj.draw(ax)
    field_obj.draw(ax)




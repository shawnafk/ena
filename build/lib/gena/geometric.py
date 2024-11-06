import numpy as np
def q2rot(q1,q2,q3,N):
    mat = np.zeros((N,3,3))
    q0 = np.sqrt(1-q1**2-q2**2-q3**2)
    mat[:,0,0] = 1 - 2*q2**2 - 2* q3**2
    mat[:,0,1] = 2*q1*q2 - 2* q0 * q3
    mat[:,0,2] = 2 * q1 * q3 + 2 * q0 * q2
         
    mat[:,1,0] = 2*q1*q2 + 2* q0 * q3
    mat[:,1,1] = 1 - 2*q1**2 - 2* q3**2
    mat[:,1,2] = 2 * q2 * q3 - 2 * q0 * q1
         
    mat[:,2,0] = 2*q1*q3 - 2* q0 * q2
    mat[:,2,1] = 2*q2*q3 + 2* q0 * q1
    mat[:,2,2] = 1 - 2*q1**2 - 2*q2**2
    return mat

#only one
def rot_x(deg):
    deg = deg/180*np.pi
    r11 = 1
    r12 = 0
    r13 = 0
    
    r21 = 0
    r22 = np.cos(deg)
    r23 = -np.sin(deg)
    
    r31 = 0
    r32 = np.sin(deg)
    r33 = np.cos(deg)
    
    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat

def rot_y(deg):
    deg = deg/180*np.pi
    r11 = np.cos(deg)
    r12 = 0
    r13 = np.sin(deg)
    
    r21 = 0
    r22 = 1
    r23 = 0
    
    r31 = -np.sin(deg)
    r32 = 0
    r33 = np.cos(deg)
    
    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat

def rot_z(deg):
    deg = deg/180*np.pi
    r11 = np.cos(deg)
    r12 = -np.sin(deg)
    r13 = 0
    
    r21 = np.sin(deg)
    r22 = np.cos(deg)
    r23 = 0
    
    r31 = 0
    r32 = 0
    r33 = 1
    
    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat
#--------------------------------------------------------------- angle manipulation -------------------------
def vec_2_rphi(vec):
    #return the angle along long axis: x , and short axis: y  
    phi_l = np.arctan(vec[:,0]/vec[:,2])/np.pi*180
    phi_s = np.arctan(vec[:,1]/vec[:,2])/np.pi*180
    r = np.linalg.norm(vec,axis=1)
    return r,phi_l,phi_s

#unit r = 1
def rphi_2_vec(ps,pl,N=1):
    z_a = np.ones(N)
    x_a = np.tan(pl/180*np.pi)
    y_a = np.tan(ps/180*np.pi)
    r_a = np.sqrt(x_a**2 + y_a**2 + z_a**2)
    v = np.zeros((N,3))
    v[:,0] = x_a/r_a
    v[:,1] = y_a/r_a
    v[:,2] = z_a/r_a
    return v

#rotation in the long (x), wide direction (along short axis)
#V has shape N * 3
def rot_vec(vec,phi,ax):
    N=vec.shape[0]
    _,pl,ps = vec_2_rphi(vec)
    if ax == "long": 
        pl = pl + phi
    elif ax == "short":
        ps = ps + phi
    return rphi_2_vec(ps,pl,N)

def get_ang(a,b):
    adb = np.sum(np.multiply(a,b),axis=1)
    na = np.linalg.norm(a,axis=1)
    nb = np.linalg.norm(b,axis=1)
    ang = np.arccos(adb/na/nb)/np.pi*180
    return  ang



if __name__ == "__main__":
    print("unit test here")
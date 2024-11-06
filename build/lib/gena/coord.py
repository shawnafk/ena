from pyspedas import cotrans
import numpy as np
from gena.const import *
#
#vec true
#mcp + offset 
#offset is 0
def mcp_2_sc(in_mcp,offset):
    return in_mcp + offset

#sc -> mj2000  
def sc_2_mj2000(in_sc,rot,sc_center_mj2000):
    return rot @ in_sc + sc_center_mj2000 

#mj2000  -> ej2000
def mj2000_2_ej2000(in_mj2000,loc_m_ej2000):
    return in_mj2000 + loc_m_ej2000

#sequential
#ej2000  -> gse
def ej2000_2_gse(tin,ej2000): 
    gse = cotrans(time_in=tin,data_in=ej2000,coord_in = 'j2000',coord_out='gse')
    #ej2000 to gei
    #Rgei = cotrans_lib.subj20002gei(time_in],[ej2000])[0]
    #gei to gse
    #Rgse = cotrans_lib.subgei2gse([time_in],[Rgei])[0]
    return gse
def gsm_2_gse(tin,gsm): 
    gse = cotrans(time_in=tin,data_in=gsm,coord_in = 'gsm',coord_out='gse')
    return gse
if __name__ == "__main__":
    print("unit test here")

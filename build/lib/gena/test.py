from gena import GENA
#path = './data/output/Detector_B/time_sorted/2A_matched/'
path = './'
date = '20241011'
surfix = '_Detector_B.csv'
fn = path+date+surfix

prob_B = GENA.prob("B")
prob_B.generate(fn)
prob_B.dump('test.h5')
prob_B.load('test.h5')

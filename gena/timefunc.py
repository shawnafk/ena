from datetime import datetime, timedelta
import pytz
import numpy as np

#dateime object has timezone information
#datetime.strptime(time_str, FORMAT)
#datetime.fromtimestamp(timestamp)
#dtobj.timestamp()
#dtobj.strftime('%Y-%m-%d %H:%M:%S')

format = '%Y-%m-%d %H:%M:%S%z'

#create Datetime object from str, timestamp to utc 1970-01-01 00:00:00
def Str_2_Datetime(str,tz='+0000'):
    return datetime.strptime(str+tz, format)

def Timestamp_2_Datetime(ts,tz='UTC'):
    return datetime.fromtimestamp(ts,pytz.timezone(tz))

def Datetime_2_Str(dt):
    return dt.strftime(format)

#get timestamp of the event from satellite time and gives a UTC datetime and string
def SatelliteTime_2_Datetime(ST):
    SatelliteTime0 = datetime.strptime('2021-01-01 00:00:00+0800', format)
    TS = SatelliteTime0.timestamp() + ST
    DT = datetime.fromtimestamp(TS,pytz.timezone('UTC'))
    return DT

from astropy.time import Time
#convert timestamp to iso time
def Datetime_2_ISOT(dt):
    isot = Time(dt, scale='utc')
    return isot


#time function (useless)
def Find_Closest_Idx(timestamp_list, t0):
    left, right = 0, len(timestamp_list) - 1
    closest_index = -1
    min_diff = float('inf') 

    while left <= right:
        mid = (left + right) // 2
        current_diff = abs(timestamp_list[mid] - t0)

        # 如果当前的差值小于之前的最小差值，则更新最接近的索引和最小差值
        if current_diff < min_diff:
            min_diff = current_diff
            closest_index = mid

        # 如果当前时间戳小于t0，则在右半部分继续查找
        if timestamp_list[mid] < t0:
            left = mid + 1
        # 如果当前时间戳大于t0，则在左半部分继续查找
        else:
            right = mid - 1

    return closest_index

#given time t0  in string in YYYY-MM-DD hh:mm:ss
#time step dt in (s)
#number of intervals n
#extract n index of the left values of those specific times: t0+dt*n, n=0 to n-1
#thus np can split to n+1 (useless)
def Get_Idx(dt,deltat,n,timestamp_list):
    ts = dt.timestamp()
    target_time = ts + np.arange(0,n)*deltat
    idx = []
    for i in range(n):
        idx.append(Find_Closest_Idx(timestamp_list,target_time[i]))
    return idx

if __name__ == "__main__":

    print("unit test here")

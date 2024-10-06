import math
from math import radians, cos, sin, asin, sqrt
import numpy as np
import pandas as pd
import os
import sys
def geodistance(lng1,lat1,lng2,lat2):
        lng1, lat1, lng2, lat2 = map(radians, [lng1, lat1, lng2, lat2])
        dlon=lng2-lng1
        dlat=lat2-lat1
        a=sin(dlat/2)**2+cos(lat1)*cos(lat2)*sin(dlon/2)**2
        dis=2*asin(sqrt(a))*6371
        return dis

sta_list = pd.read_table('station_sort.lst',sep = ' ', header = None)
stn = sta_list.iloc[:,0]
stlon = sta_list.iloc[:,1]
stlat = sta_list.iloc[:,2]
stlon = np.array(stlon)
stlat = np.array(stlat)
bc = pd.read_table('station_sort.lst',sep = ' ', header = None)
bcst = bc.iloc[:,0]
bclo = bc.iloc[:,1]
bcla = bc.iloc[:,2]
bclo = np.array(bclo)
bcla = np.array(bcla)
period = [3]

wavel_ref = 2*3
for per in period:
	refer_vel  = 3
	wavel = per*refer_vel
	dir0 = 'db_period_'+str(per)
	if not os.path.isdir(dir0):
		os.mkdir(dir0)
	for i in range(len(bclo)):
		loc_s_stn = bcst[i]
		loc_s_lon = bclo[i]
		loc_s_lat = bcla[i]	
		dir1 = str(dir0)+'/'+loc_s_stn
		if not os.path.isdir(dir1):
			os.mkdir(dir1)
		f1 = open(dir1+'/source_list','w+')
		for x in range(0,len(sta_list)):
			if round(geodistance(loc_s_lon,loc_s_lat,stlon[x],stlat[x]),3) <= 1.0*wavel_ref:
				distance = round(geodistance(loc_s_lon,loc_s_lat,stlon[x],stlat[x]),3)
				f1.writelines(stn[x]+' '+str(stlon[x])+' '+str(stlat[x])+' '+str(distance)+'\n')
		f1.close()
		for j in range(len(bclo)):
			loc_r_stn = bcst[j]
			loc_r_lon = bclo[j]
			loc_r_lat = bcla[j]
			if round(geodistance(loc_s_lon,loc_s_lat,loc_r_lon,loc_r_lat),3) < 3.0*wavel:
				continue
			dir2 = str(dir1)+'/'+loc_r_stn
			if not os.path.isdir(dir2):
				os.mkdir(dir2)
			f2 = open(dir2+'/receiver_list','w+')
			for k in range(len(stn)):
				distance = round(geodistance(stlon[k],stlat[k],loc_r_lon,loc_r_lat),3)
				if distance <= 1.0*wavel_ref:
					f2.writelines(stn[k]+' '+str(stlon[k])+' '+str(stlat[k])+' '+str(distance)+'\n')
			f2.close()

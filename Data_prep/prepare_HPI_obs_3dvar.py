# Script to generate tcvitals record to file that is written in 3Dvar format &
# which is readable to WRFDA

# Credit to Scott Sieron
# Author: Zhu (Judy) Yao. Aug 19, 2024

import numpy as np
import datetime
import subprocess

#### parameters ####
storm_name = 'Otis'
tint = 60
obs_error = 300.0
OUTPUT_DIR = '/expanse/lustre/scratch/cpruett/temp_project/Otis/HPI'
TCvital_dir = '/expanse/lustre/scratch/cpruett/temp_project/Otis/TCVitals'
####################

### functions ###
def run_script(cmd):
    print(cmd)
    subprocess.getoutput(cmd)
    return
###################

## make 3dvar formatted HPI files from TCvitals
#time = datetime.datetime(int(line[19:23]),int(line[23:25]),int(line[25:27]),int(line[28:30]),int(line[30:32]),0)
date, mslp, lon, lat = [], [], [], []
itime = datetime.datetime(2023,10,16,00,0,0)
ftime = datetime.datetime(2023,10,22,12,0,0)

time = itime
while time <= ftime:
    filename = TCvital_dir+'/'+time.strftime('%Y%m%d%H%M')+'.'+storm_name+'-tcvitals.dat'
    infile = open(filename)
    line = infile.read()
    lon_slp, lat_slp = float(line[38:42])/10., float(line[33:36])/10.
    if line[42] == 'W':
       lon_slp = lon_slp * -1.
    minslp = float(line[52:56])
    print('original TCvitals: '+time.strftime('%Y%m%d%H%M')+', ('+str(lon_slp)+'W '+str(lat_slp)+'N) '+str(minslp)+'hPa')
    date.append(time)
    mslp.append(minslp)
    lat.append(lat_slp)
    lon.append(lon_slp)
    time = time + datetime.timedelta(minutes = 360)

time = max(itime,np.amin(date))
while time <= min(ftime,np.amax(date)):
    rec = max(1, np.searchsorted(date,time))
    dt1 = float((time-date[rec-1]).seconds)
    dt2 = float((date[rec]-time).seconds)
    lat_slp = ( lat[rec-1]*dt2  + lat[rec]*dt1)/(dt1+dt2)
    lon_slp = ( lon[rec-1]*dt2  + lon[rec]*dt1)/(dt1+dt2)
    minslp  = (mslp[rec-1]*dt2 + mslp[rec]*dt1)/(dt1+dt2)*100.

    # output
    # template
    file_output = OUTPUT_DIR+'/HPI_obs_gts_'+time.strftime('%Y-%m-%d_%H:%M')+':00.3DVAR'
    print('writing... '+file_output)
    print('interpolated TCvitals: ('+str(lon_slp)+'W '+str(lat_slp)+'N) '+str(minslp)+'hPa')

    f = open(file_output,'w')

    f.write('TOTAL =      1, MISS. =-888888.,'+'\n')
    f.write('SYNOP =      0, METAR =      1, SHIP  =      0, BUOY  =      0, BOGUS =      0, TEMP  =      0,'+'\n')
    f.write('AMDAR =      0, AIREP =      0, TAMDAR=      0, PILOT =      0, SATEM =      0, SATOB =      0,'+'\n')
    f.write('GPSPW =      0, GPSZD =      0, GPSRF =      0, GPSEP =      0, SSMT1 =      0, SSMT2 =      0,'+'\n')
    f.write('TOVS  =      0, QSCAT =      0, PROFL =      0, AIRSR =      0, OTHER =      0,'+'\n')
    f.write('PHIC  =  17.50, XLONC = -73.00, TRUE1 =   0.00, TRUE2 =  30.00, XIM11 =   1.00, XJM11 =   1.00,'+'\n')
    f.write('base_temp= 290.00, base_lapse=  50.00, PTOP  =  1000., base_pres=100000., base_tropo_pres= 20000., base_strat_temp=   215.,'+'\n')
    f.write('IXC   =    550, JXC   =    300, IPROJ =      1, IDD   =      1, MAXNES=      1,'+'\n')
    f.write('NESTIX=    550,'+'\n')
    f.write('NESTJX=    300,'+'\n')
    f.write('NUMC  =      1,'+'\n')
    f.write('DIS   =   27.0,'+'\n')
    f.write('NESTI =      1,'+'\n')
    f.write('NESTJ =      1,'+'\n')
    f.write('INFO  = PLATFORM, DATE, NAME, LEVELS, LATITUDE, LONGITUDE, ELEVATION, ID.'+'\n')
    f.write('SRFC  = SLP, PW (DATA,QC,ERROR).'+'\n')
    f.write('EACH  = PRES, SPEED, DIR, HEIGHT, TEMP, DEW PT, HUMID (DATA,QC,ERROR)*LEVELS.'+'\n')
    f.write('INFO_FMT = (A12,1X,A19,1X,A40,1X,I6,3(F12.3,11X),6X,A40)'+'\n')
    f.write('SRFC_FMT = (F12.3,I4,F7.2,F12.3,I4,F7.3)'+'\n')
    f.write('EACH_FMT = (3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2))'+'\n')
    f.write('#------------------------------------------------------------------------------#'+'\n')

# file information
    dataset = 'FM-15        '+time.strftime('%Y-%m-%d_%H:%M:%S')+' MADIS                                         1'+"{0:12.3f}".format(float(lat_slp))+"{0:23.3f}".format(float(lon_slp))+"{0:23.3f}".format(0.0)+"{0:>21}".format('CWZQ')+'\n'
    f.write(dataset)
    dataset = "{0:12.3f}".format(float(minslp))+"{0:4d}".format(0)+"{0:7.2f}".format(obs_error)+' -888888.000 -88  0.200'+'\n'
    f.write(dataset)
    dataset = ' -888888.000  -5 100.00 -888888.000  -5 100.00 -888888.000  -5 100.00            -888888.000  -5 100.00 -888888.000  -5 100.00 -888888.000 -11   2.00            -888888.000 -11  10.00'+'\n'
    f.write(dataset)

    time = time + datetime.timedelta(minutes = tint) 
run_script('chmod 755 '+OUTPUT_DIR+'/'+time.strftime('%Y')+'/obs_gts*')



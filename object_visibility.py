#!/usr/bin/env python
# This script the calculates the visibility of stars from a particular position on Earth. It also calulates the average topocentric position of the moon for that particular date given in conf file.
# As output it gives a plot of IST vs. alt. of sources. The distances of the sources from the moon are also displayed in terminal.
from numpy import *
import matplotlib.pyplot as plt 
import jdcal
import time
import sys
from datetime import date

# Calculates LST for a given date, time and Longitute of the place.
def lst(yr, mn, day, ut, lon):
    jd1 = sum(jdcal.gcal2jd(yr, mn, day))
    d_2000 = jd1 + ut/24. - 2451545.0
    LST = 100.46 + 0.985647*d_2000 + lon + 15.0*ut
    if LST < 0:
        LST = LST + 360.
    resd, LST = divmod(LST, 360.)
    return LST/15.

# Calculation of ALT-AZM from RA-Dec
def alt_azm(rat, dect, almst, lat):
    ha = almst - rat
    sinal = sin(radians(dect))*sin(radians(lat)) + cos(radians(dect))*cos(radians(lat))*cos(radians(ha*15.))
    alt = degrees(arcsin(sinal))
    sinaz = -cos(radians(dect))*sin(radians(ha*15.))/cos(radians(alt))
    az = 180.0 - degrees(arcsin(sinaz))
    return alt, az

# Calculation of Lunar position
def moon_pos(yr, mn, day, ut, lat, lon):
    D= (date(yr, mn, day) - date(yr, 1, 1)).days
    jd1 = sum(jdcal.gcal2jd(yr, 1, 0))
    D2000 = jd1 - 2451545.0
    temp, gmst = divmod((280.461+360.98564737*D2000)/15.0, 24)
    gmst=gmst+0.0657098243*D+1.00273791*ut
    almst= gmst + lon/15.
    if almst > 24:
        almst = almst -24.

    t = (sum(jdcal.gcal2jd(yr,mn,day)) + ut/24. -2451545.0)/36525.
 
    theta=100.46+36000.77*t+lon+15.*ut
    al=218.32000 + 481267.88300*t + 6.29000*sin(radians(134.90000 + 477198.85000*t))
    al=al-1.27*sin(radians(259.2-413335.38*t))+0.66*sin(radians(235.7+890534.23*t))
    al=al+0.21*sin(radians(269.9+954397.7*t))-0.19*sin(radians(357.5+35999.05*t))-0.11*sin(radians(186.6+966404.05*t))

    ab=5.13*sin(radians(93.3+483202.03*t))+0.28*sin(radians(228.2+960400.87*t))-0.28*sin(radians(318.3+6003.18*t))
    ab=ab-0.17*sin(radians(217.6-407332.2*t))

    apar=0.9508+0.0518*cos(radians(134.9+477198.85*t))+0.0095*cos(radians(259.2-413335.38*t))
    apar=apar+0.0078*cos(radians(235.7+890534.23*t))+0.0028*cos(radians(269.9+954397.70*t))

    rr=1./sin(radians(apar))

    dl=cos(radians(ab))*cos(radians(al))
    dm=0.9175*cos(radians(ab))*sin(radians(al))
    dm=dm-0.3978*sin(radians(ab))
    dn=0.3978*cos(radians(ab))*sin(radians(al))+0.9175*sin(radians(ab))

    rag=abs(degrees(arctan(dm/dl)))

    if (dm > 0) & (dl > 0):
        rag = rag
    if (dm < 0) & (dl > 0):
        rag = 360. - rag
    if (dm > 0) & (dl < 0):
        rag = 180. - rag
    if (dm < 0) & (dl < 0):
        rag = 180. + rag

# rag, decg: GEOCENTRIC LUNAR COORDINATES
    rag = rag/15.                                                       
    decg = degrees(arcsin(dn))

    x1=rr*dl
    y1=rr*dm                                                      
    z1=rr*dn
    xt=x1-cos(radians(lat))*cos(radians(theta))
    yt=y1-cos(radians(lat))*sin(radians(theta))
    zt=z1-sin(radians(lat))
    rrt=sqrt(xt**2 + yt**2 + zt**2)
    rat=abs(degrees(arctan(yt/xt)))
    
    if (xt > 0) & (yt > 0):
        rat = rat
    if (xt > 0) & (yt < 0):
        rat = 360. - rat
    if (xt < 0) & (yt > 0):
        rat = 180. - rat
    if (xt < 0) & (yt < 0):
        rat = 180. + rat

    rat = rat/15.
    dect = degrees(arcsin(zt/rrt))
    apart= degrees(arcsin(1./rrt))
    sdm = 3600.*0.2725*apart

    altm, azm = alt_azm(rat, dect, almst, lat)
    return rat, dect, altm, azm

###Calculation for Sun rise and set
def sun_rise_set(yr, mn, day, lat, lon, alt):
    local_noon = 12.0 + ((82.5 - lon)/15.0)
    D= (date(yr, mn, day) - date(yr, 1, 1)).days
    t1 = sin(radians(360.0*(D+284.0)/365.0))
    temp = -1*tan(radians(lat))*tan(radians(23.44*t1))

    H = abs((1./15.)*degrees(arccos(temp)))
    srise = local_noon - H  - (alt/1500.)/60.
    sset = local_noon + H + (alt/1500.)/60.
    return srise, sset

#####------------- Definitions are over ------------------


##------------------------------------------------------##
##------------- Main Program Begins here ---------------##

##----------- Reading Configuration file ----------------
if len(sys.argv)<2 :
    print('\nUsage : {0} input.conf'.format(sys.argv[0]))
    print('input.conf contains the user inputs.\n')
    sys.exit(1)

try : 
    configfile=open(sys.argv[1],'r')
except IOError :
    print ("Error: Cannot open the file "+sys.argv[1]+". Setup the config file correctly, before running the script.")
    sys.exit(1)

for con in configfile:
    con=con.rstrip()
    if len(con.split()) >= 2 :
        if con.split()[0] == "DATE=" :
            yr=int(con.split()[1])
            mn=int(con.split()[2])
            day=int(con.split()[3])
        elif con.split()[0] == "OBSERVATORY=" :
            lat=float(con.split()[1])
            lon=float(con.split()[2])
            altitude=float(con.split()[3])
        elif con.split()[0] == "UTDIFF=" :
            utd=float(con.split()[1])
        elif con.split()[0] == "STARCOORDFILE=" :
            fi=con.split()[1]
        elif con.split()[0] == "OUTFILE=" :
            fo=con.split()[1]
            
configfile.close()
##----------- END of Reading Configuration file ----------------

##------------ Calculating Sunrise/Sunset time -----------------
srise, sset = sun_rise_set(yr, mn, day, lat, lon, altitude)

srise_ut = srise - utd
sset_ut = sset - utd
if srise_ut < 3 :
    srise_ut = 24 + srise_ut

mor_twi_ut = srise_ut - 1.0 + utd
eve_twi_ut = sset_ut + 1.0 + utd

plt.plot([mor_twi_ut,mor_twi_ut],[-5,95], 'r--')
plt.plot([eve_twi_ut,eve_twi_ut],[-5,95], 'r--')

plt.figtext(0.21,0.91, 'Eve Twil.', color='r', size='small' , ha='center')
plt.figtext(0.82,0.91, 'Mor Twil.', color='r', size='small', ha='center')
plt.figtext(0.12,0.92, 'Sun set', color='b', size='small' , ha='center')
plt.figtext(0.91,0.92, 'Sun Rise', color='b',  size='small', ha='center')

ut = sset_ut
ut0 = ut
##----------- END of Calculating Sunrise/Sunset time -----------


##-------------- Plotting Moon's positions ------------------

# Sampling time at which lunar positions will be calculated
samp = 5.      # in min 
stp = samp/60. # in Hr
nsamp = int(abs(sset_ut - srise_ut)/stp)

ut_p = []
al_m = []
al_s = []
ut = ut0
for i in range(nsamp):
    ut = ut + stp
    rat, dect, altm, azm = moon_pos(yr, mn, day,ut, lat, lon)
    ut_p.append(ut + utd)
    al_m.append(altm)

plt.plot(ut_p,al_m, 'b--')

##------------ END of Plotting Moon's positions -------------


##----------- Plotting Star's positions --------------------
moon_position_time = 24.0    # IST
moon_position_time = moon_position_time - utd 
print '\nDistance from the moon: '
starfile=open(fi,'r')
j=1
for con in starfile:
    con=con.rstrip()
    srah=float(con.split()[0])
    sram=float(con.split()[1])
    sras=float(con.split()[2])
    sded=float(con.split()[3])
    sdem=float(con.split()[4])
    sdes=float(con.split()[5])
    
    sra = srah + sram/60. + sras/3600.
    sdec = sded + sdem/60. + sdes/3600.
   
    mra, mdec, altm, azm = moon_pos(yr, mn, day, moon_position_time, lat, lon)
#    print mra, mdec, sra, sdec
    dra = 15.0*(sra - mra)
    ddec = sdec - mdec

    dist_moon = degrees(arccos((sin(radians(sdec))*sin(radians(mdec)))+ (cos(radians(sdec))*cos(radians(mdec))*cos(radians(dra)))))
    print ('Star%d:  %5.1f'%(j, dist_moon))

    al_s = []
    ut_p = []
    ut = ut0
    ut_s = ut0
    for i in range(nsamp):
        ut = ut + stp
        lst1 = lst(yr, mn, day, ut, lon)
        sal, saz = alt_azm(sra, sdec, lst1, lat)
        ut_p.append(ut + utd)
        al_s.append(sal)
        if ((ut - ut_s) > 1.5) & (sal > 0.0):
            ut_s = ut
            plt.text(ut_s + utd, sal, str(j), fontsize=12)
    plt.plot(ut_p, al_s, 'k')
    j = j+1
##----------- END of Plotting Star's positions --------------
if lat >= 0:
    lat_mod = str(lat)+'N'
else:
    lat_mod = str(lat)+'S'
if lon >= 0:
    lon_mod = str(lon)+'E'
else:
    lon_mod = str(lon)+'W'


plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48],\
          [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
plt.yticks([10,20,30,40,50,60,70,80,90])  
plt.xlabel('UT + '+ str(utd)+' hrs')
plt.ylabel('Altitude')
plt.title('Object visibility on '+str(yr)+'/'+str(mn)+'/'+str(day)+'\n', color='b')
plt.figtext(0.51,0.91, 'Observatory coords: '+lon_mod+'  '+lat_mod, color='k', size='medium', ha='center')#,  weight='roman')
plt.xlim( (sset_ut+utd,srise_ut+utd) ) 
plt.ylim( (0,90) ) 
plt.grid()
plt.savefig(fo)
plt.show()

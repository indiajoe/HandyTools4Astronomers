#!/usr/bin/env python
# This is a crude Rise and Set time calculator for HCT-Hanle and IGO-Girawali Observatory.
# IMP: This is approximate calculation. No effect due to atmosphere refraction included
#------------------------------------------------------------indiajoe@gmail.com
# Default Elevation limit is assumed to be 30 degree.. (Change at LINE below)
Limiting_Altitude=25             # Set Default visibility elevation limit HERE <-----#######

from math import *
import time
import datetime
import sys
import readline # for interative editing of input line

def is_Number(numb) :
    try:
        float(numb)
        return True
    except ValueError:
        return False

if len(sys.argv) < 2 :
    print "--------------------------------------------------------"
    print " Usage  : "+sys.argv[0]+" Observatory  [Altitude_Limit]" 
    print " Observatories allowed are "
    print " HCT "
    print " IGO"
    print " GMRT"
    print " Example:   "+sys.argv[0]+" HCT 30"
    print " Altitude_Limit is optional (in degrees). :- Default Value =" +str(Limiting_Altitude)
    print "----------------------------------------------------------------"
    exit(1)

if (sys.argv[1] == "HCT") : 
    Longitude=78.9641
    Latitude=32.77944
elif (sys.argv[1] == "IGO") : 
    Longitude=73.666
    Latitude=19.0833
elif (sys.argv[1] == "GMRT") : 
    Longitude=74.04974
    Latitude=19.09651
else : 
    print "Not able to recognise Observatory ",sys.argv[1],". Enter only HCT,IGO or GMRT"
    exit(1)

if ( len(sys.argv) == 3 and is_Number(sys.argv[2]) ) :
    Limiting_Altitude=float(sys.argv[2])
#Longitude=74.05
IST_Longitude=82.58
Correction_Long=(IST_Longitude-Longitude)*24/360.0

Elevation_Limit= radians(Limiting_Altitude)    

#Year = input('Year of observation ')
Year=2010
Month = input('MONTH of Observation (Enter in number:- Eg 1 for January ) : ')
Day = input('DAY of Observation : ')

diff = datetime.date(2010, Month, Day) - datetime.date(2010, 9, 23)

if (diff.days > 0) : No_days=diff.days
else: No_days=365 + diff.days

#At Midnight the following RA start will be at the meridian
RA_at_Mid=24.0*No_days/365.0 -Correction_Long
print "RA at Midnight will be ", RA_at_Mid
#print "Press Cntrl-C to exit after done with all sources"
while 1 :
    try:
        rainput=raw_input('RA  (Eg: hr min sec) : ').strip(' ')  #strips away leading and trailing whitespaces.
        decinput=raw_input('Dec (Eg: Deg min sec) : ').strip(' ')
    except (KeyboardInterrupt, EOFError):
        print("\n Exiting.. \n ")
        exit(0)

    RAdelimiter=' '
    DECdelimiter=' '
    if rainput.find(':') > 0 : RAdelimiter=':'
    if decinput.find(':') > 0 : DECdelimiter=':'
    try:
        h,m,s =map(float, rainput.split(RAdelimiter))
        d,dm,ds =map(float, decinput.split(DECdelimiter))
    except (ValueError):
        print("Error: Wrong input given. \n---------------------------------------------")
        continue
    
    RA=h + m/60.00 + s/3600.00
    if  (decinput.split(DECdelimiter)[0][0] == "-") :  #Negative declination
        dm=-dm
        ds=-ds
    DEC=d + dm/60.00 + ds/3600.00
    RA_diff=RA-RA_at_Mid
    Hours=floor(RA_diff)
    Minut=(RA_diff-Hours)*60
    Transit_min= floor(Minut)
    if (RA_diff > 0):
        Transit_hr=Hours
    else:
        Transit_hr=24+Hours
#    print "Transit of the star will occur at ",int(Transit_hr),":",int(Transit_min) ,"hours"
    # Rise and Set time calculation
    try: 
        H_rad= acos(sin(Elevation_Limit)/(cos(radians(Latitude))*cos(radians(DEC))) - tan(radians(Latitude))*tan(radians(DEC)))
    except (ValueError):
        print("---------------------------------------------")
        print("MathError: WARNING: The object probably NEVER rises nor sets")
        print("Object's Transit time = "+str(int(Transit_hr))+":"+str(int(Transit_min)))
        print("---------------------------------------------")
        continue
    H=H_rad*24.0/(2*3.14)  #in hours
    Rise_Hours= floor(RA_diff-H)
    Rise_Minut=(RA_diff-H-Rise_Hours)*60
    Rise_min= floor(Rise_Minut)
    if (RA_diff - H > 0):
        Rise_hr=Rise_Hours
    else:
        Rise_hr=24+Rise_Hours
#    print "Rise of the star will occur at ",int(Rise_hr),":",int(Rise_min) ,"hours"
    Set_Hours= floor(RA_diff+H)
    Set_Minut=(RA_diff+H-Set_Hours)*60
    Set_min= floor(Set_Minut)
    if (RA_diff + H > 0):
        Set_hr=Set_Hours
    else:
        Set_hr=24+Set_Hours
#    print "Setting of the star will occur at ",int(Set_hr),":",int(Set_min) ,"hours"

    print "*Rise* --------- *Transit* ---------- *Set*"
    print int(Rise_hr),":",int(Rise_min),"----------",int(Transit_hr),":",int(Transit_min) ,"-----------",int(Set_hr),":",int(Set_min)
    print "---------------------------------------------"
    


#!/usr/bin/env python
# This is a simple JD date calculator 
# The input is two columns of YearMonthDate and TimeinSeconds OR Time in HR:MIN:SEC 
# Example: cat 2columnTextfile | ./JD_calc_pipeInput.py 
# Where the 2columnTextfile is in format 20120831 5436  <OR> 20120831 16:25:34
#...................................................indiajoe@gmail.com
import datetime
import sys
while 1 :
    try :
        Input=raw_input()
    except EOFError :
        break
    YearMonthDate,Time=Input.split()

    Year=int(YearMonthDate[0:4])
    Month=int(YearMonthDate[4:6])
    Day=int(YearMonthDate[6:8])
    JD19900101=2447893
    diff = datetime.date(Year, Month, Day) - datetime.date(1990, 1, 1)
    try :
        timefrac=eval(Time)/(24*60*60.0)
    except SyntaxError:
        t=Time.split(':')
        timefrac=(int(t[0])*60*60.0+int(t[1])*60+int(t[2]))/(24*60*60.0)
    print(str(JD19900101+diff.days+timefrac))

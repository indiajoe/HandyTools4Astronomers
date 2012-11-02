#!/usr/bin/env python
# This is a simple JD date calculator
# Useage: ./JD_calc.py 2010 08 31
#..........................indiajoe@gmail.com
import datetime
import sys
Year=int(sys.argv[1])
Month=int(sys.argv[2])
Day=int(sys.argv[3])
JD19900101=2447893
diff = datetime.date(Year, Month, Day) - datetime.date(1990, 1, 1)
print JD19900101+diff.days

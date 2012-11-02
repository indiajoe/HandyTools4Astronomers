#!/usr/bin/env python
# This is a simple JD date to Ordinary Calander
import datetime
import sys
JDinput=int(sys.argv[1])
#Month=int(sys.argv[2])
#Day=int(sys.argv[3])
JD19900101=2447893
JDdiff=JDinput-JD19900101
New_date = datetime.timedelta(JDdiff) + datetime.date(1990, 1, 1)
print New_date

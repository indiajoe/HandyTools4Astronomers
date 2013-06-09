#!/usr/bin/env python
# Based on example in  http://aladin.u-strasbg.fr/java/FAQ.htx#ToC112
#-------------
import subprocess

imagesize=6 #in arminutes
OBJECTLIST=('M1', 'M104', 'NGC2024')
AladinLocation='/usr/bin/Aladin.jar'

p = subprocess.Popen(['java -jar '+AladinLocation],
                     shell=True,
                     stdin=subprocess.PIPE
                       )
p.stdin.write('grid on\n')
for obj in OBJECTLIST :
#    p.stdin.write('reset; get aladin,VizieR(GSC1.2),simbad '+obj+';\n')
    for filt in ('J','H','K'):  # For 2mass all sky survey images
        p.stdin.write('reset; get allsky("2MASS '+filt+'",jpeg) '+obj+' ;\n')
        p.stdin.write('zoom '+str(imagesize)+'arcmin;save '+obj+filt+'.jpg\n')

p.stdin.write('quit\n')
p.wait()

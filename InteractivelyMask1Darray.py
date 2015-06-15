#/usr/bin/env python
""" This script is to interactively mask points in a 1 D numpy array.
Usefull for masking bad pixels in 1D spectrum """

import matplotlib.pyplot as plt
import numpy as np

def InteractivelyMask(InputArray,Mask=None):
    """ Interactively mask points in input Array and return the mask array.
    An optional Initial Mask can also be provided which will be overwitten."""

    if Mask is None:
        Mask = np.zeros(len(InputArray),dtype=bool) 

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('m to mask the nearest point, u to unmask a masked point')
    ax.plot(InputArray,linestyle='--', drawstyle='steps-mid',marker='.',alpha=0.5)
    maskedLine, = ax.plot(np.ma.array(InputArray,mask=Mask),'.',color='g')
    PlotedLines = [maskedLine]
    # Define the function to run while key is pressed
    def on_key(event):
        if event.key == 'm' :
            Nearestpoint = int(round(event.xdata))
            print('Masking index = {0}'.format(Nearestpoint))
            Mask[Nearestpoint] = True
        elif event.key == 'u' :
            Nearestpoint = int(round(event.xdata))
            print('Unmasking index = {0}'.format(Nearestpoint))
            Mask[Nearestpoint] = False

        ax.lines.remove(PlotedLines[0])  #Removing previous line plots
        PlotedLines[0], = ax.plot(np.ma.array(InputArray,mask=Mask),'.',color='g')
        ax.figure.canvas.draw()

    cid = fig.canvas.mpl_connect('key_press_event', on_key)
    plt.show()
    return Mask


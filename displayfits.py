#!/usr/bin/env python
# This example shows how to display a fits image in matplotlib
#-------------------------------------------indiajoe
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

def display(image,scale='histeq',title=''):
    """ Display 2D fits images or ndarrays using imshow of matplotlib """
    #Input Sanity check
    if type(image)==str :  #User has given filename
        try :
            img=pyfits.convenience.getdata(image)
        except IOError as err:
            print("Cannot open file to display:{0} \n{1}".format(image,err))
            return
    elif type(image)==np.ndarray :
        img=image
    else :
        print("Unknown Input format to display")
        return
    # Now we proceed to displaying the image
    if scale=='histeq' : #We have to do Histogram Equalisation of image data. See http://en.wikipedia.org/wiki/Histogram_equalization
        img4disp=np.empty(img.size,dtype=img.dtype)  #Creating a 1d array of same size
        img4disp[np.argsort(img,axis=None)]=np.arange(img.size,dtype=img.dtype)  #Creating the 1d hist equalised image
        img4disp=img4disp.reshape(img.shape)  # Reshaping to previous 2d image size

    #More scale options can be added here later.    
    else:    
        img4disp=img
    
    plt.imshow(img4disp,origin='lower',cmap=plt.cm.gray)  #Display the image in gray scale in ds9 style
    plt.title(title)           # Set title if any
    plt.show()   # Show the image window
    return

    

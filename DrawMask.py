#!/usr/bin/env python
""" This tool is to apply a polygon mask by interatively drawing any polynomial on image 

    Adapted from http://matplotlib.org/1.4.3/examples/event_handling/poly_editor.html
"""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
from matplotlib.patches import Polygon
from matplotlib.pyplot import cm

class MaskDrawer(object):
    """An interactive polygon mask drawer on an image.
    Parameters
    ----------
    ax     : matplotlib plot axis
    
    Inpimg : 2d numpy array
          Input image to overlay for drawing mask
    Mask : Boolean numpy array same size of inpimg
           A Mask which will be used as initial mask and updated upon confirmation
    max_ds : float
           Max pixel distance to count as a vertex hit.
    PolyAtStart : List of vertices
           A list of square vertices to draw the initial polygon
    Key-bindings
    ------------
    't' : toggle vertex markers on and off. When vertex markers are on,
          you can move them, delete them
    'd' : delete the vertex under point
    'i' : insert a vertex at point. You must be within max_ds of the
          line connecting two existing vertices
    'n' : Invert the region selected by polynomial for masking
    'c' : Confirm the polygon and update the mask
    """ 
    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, ax, Inpimg, Mask, max_ds=10,PolyAtStart = [(50,50),(100,50),(100,100),(50,100)]):
        self.showverts = True
        self.max_ds = max_ds
        self.Mask = Mask
        self.img = Inpimg
        self.maskinvert = False
        # imshow the image
        self.imgplot = ax.imshow(np.ma.filled(np.ma.array(self.img,mask=self.Mask,fill_value=np.nan)), cmap=cm.gray)
         
        self.poly = Polygon(PolyAtStart, animated=True,
                            fc='y', ec='none', alpha=0.5)
 
        ax.add_patch(self.poly)
        ax.set_clip_on(False)
        ax.set_title("Click and drag a point to move it; "
                     "'i' to insert; 'd' to delete.\n"
                     "'n' to invert the region for masking, 'c' to confirm & apply the mask.")
        self.ax = ax
         
        x, y = zip(*self.poly.xy)
        self.line = Line2D(x, y, color='none', marker='o', mfc='r',
                               alpha=0.7, animated=True)
#        self._update_line()
        self.ax.add_line(self.line)
         
        self.poly.add_callback(self.poly_changed)
        self._ind = None # the active vert
         
        canvas = self.poly.figure.canvas
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback) 
        self.canvas = canvas

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state


    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]

        if d[ind]>=self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        if event.key=='t':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts: self._ind = None
        elif event.key=='d':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.poly.xy = [tup for i,tup in enumerate(self.poly.xy) if i!=ind]
                self.line.set_data(zip(*self.poly.xy))
        elif event.key=='i':
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y # display coords
            for i in range(len(xys)-1):
                s0 = xys[i]
                s1 = xys[i+1]
                d = dist_point_to_segment(p, s0, s1)
                if d<=self.epsilon:
                    self.poly.xy = np.array(
                        list(self.poly.xy[:i]) +
                        [(event.xdata, event.ydata)] +
                        list(self.poly.xy[i:]))
                    self.line.set_data(zip(*self.poly.xy))
                    break
                    
        elif event.key=='n':
            """ Flips the region inside out of the polygon to be masked """
            print('Inverting the mask region')
            self.maskinvert = not self.maskinvert


        elif event.key=='c':
            """ Confirm the drawn polynomial shape and add update the mask """
            self.UpdateMask()
            #Update the imshowed image with new mask
            self.imgplot.set_data(np.ma.filled(np.ma.array(self.img,mask=self.Mask,fill_value=np.nan)))
            self.imgplot.figure.canvas.draw()

        self.canvas.draw()

    def UpdateMask(self):
        """ Updates the maks array with points insied the polygon """
        print('Updating the original Mask..')
        Path = self.poly.get_path()
        h, w = self.Mask.shape
        y, x = np.mgrid[:h,:w]
        XYpoints = np.transpose((x.ravel(), y.ravel()))
        NewMask = Path.contains_points(XYpoints)
        if self.maskinvert :
            NewMask = ~NewMask
        # Combine the two mask by taing an elemet wise or
        self.Mask = self.Mask | NewMask.reshape((h,w))

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata

        self.poly.xy[self._ind] = x,y
        self.line.set_data(zip(*self.poly.xy))

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)


def main():
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import sys
    imgfname = sys.argv[1]
    outimgname = sys.argv[2]

    ax = plt.subplot(111)

    with fits.open(imgfname) as hdulist:
        Image = hdulist[0].data
        Mask = np.zeros(Image.shape, dtype=bool)  # An array of all false
        Img4Disp = HistEquilize(Image)
        MC = MaskDrawer(ax,Img4Disp,Mask)
        plt.show()
        MaskedImage = np.ma.filled(np.ma.array(Image,mask=MC.Mask,fill_value=np.nan))
        hdulist[0].data = MaskedImage
        hdulist.writeto(outimgname)
        np.save(imgfname[:-5]+'_Mask.npy',MC.Mask)


def HistEquilize(img):
    """Returns a histogram equlised array for display"""
    img4disp=np.empty(img.size,dtype=img.dtype)  #Creating a 1d array of same size
    img4disp[np.argsort(img,axis=None)]=np.arange(img.size,dtype=img.dtype)  #Creating the 1d hist equalised image
    img4disp=img4disp.reshape(img.shape)  # Reshaping to previous 2d image size
    return img4disp

if __name__ == '__main__':
    main()

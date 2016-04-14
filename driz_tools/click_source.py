from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
import glob
from matplotlib.colors import PowerNorm, LogNorm
import numpy as np
from photutils.morphology import centroid_com, centroid_2dg
from scipy import ndimage

def parse_args():
    """Parse command line arguments.

    Parameters:
        Nothing

    Returns:
        arguments: argparse.Namespace object
            An object containing all of the added arguments.

    Outputs:
        Nothing
    """

    input_help = 'Images to register.  Default is all visit level drizzle images (f*dr?.fits)'
    ncoords_help = 'Number of sources to click on?  Default 3.'

    inp = ''

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default=inp)
    parser.add_argument('-n', type=int, help=ncoords_help, action='store',
        required=False, default=3)
    arguments = parser.parse_args()
    return arguments

def onpress(event):
    global ix, iy
    if event.key == 'a':
        ix, iy = event.xdata, event.ydata
        print 'x = %d, y = %d'%(
            ix, iy)

        global coords
        coords.append((ix, iy))

        if len(coords) == ncoords:
            fig.canvas.mpl_disconnect(cid)
        return coords

def mark(im):
    data = fits.getdata(im,1)
    print im
    global fig
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(111)
    disp = ax1.imshow(data, cmap='afmhot', origin='lower',\
    interpolation='none', norm=PowerNorm(0.5), vmin=0.,vmax=1.)
    plt.colorbar(disp)
    print 'Hover over source location and press the  \'a\' key'
    global cid
    cid = fig.canvas.mpl_connect('key_press_event', onpress)
    plt.show()

def calc_centroids(section):
    y2, x2 = centroid_2dg(section)
    print 'VALUE ', section[int(y2+.5),int(x2+.5)]
    # print y2,x2
    return x2+.5, y2+.5

def show_centroids(coords, im):
    data = fits.getdata(im,1)
    f = open(im.replace('.fits','_sci1_xy_catalog.coo'),'w')
    for c in coords:
        x, y = c[1], c[0]
        section = data[x-10:x+10,y-10:y+10]
        cx,cy = calc_centroids(section)
        actual_x, actual_y = x+cx-10, y+cy-10
        print actual_x, actual_y
        f.write('{} {}\n'.format(actual_y,actual_x))
        # xout,yout = ndimage.measurements.center_of_mass(section)
        # print section[int(xout),int(yout)]
        # fig = plt.figure(figsize=(10,10))
        # ax1 = fig.add_subplot(111)
        # cutout = ax1.imshow(section,cmap='afmhot', origin='lower',\
        # interpolation='none', norm=PowerNorm(.4))
        # plt.colorbar(cutout)
        # ax1.plot([cy],[cx],'kx',scalex=False, scaley=False)
        # plt.show()
    f.close()

if __name__ == '__main__':
    options = parse_args()
    assert (len(glob.glob(options.i)) == 1), 'Need exactly 1 image to manual source finding'
    coords = []
    global ncoords
    ncoords = options.n
    im = glob.glob(options.i)[0]
    mark(im)
    # coords = [(1616.7190236317433, 3128.1874459665414)]
    print coords
    show_centroids(coords, im)

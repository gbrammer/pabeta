from astropy.io import fits
from drizzlepac import pixtosky
import glob, os, shutil
from multiprocessing import Pool
import argparse
import numpy as np
import math


def parse_args():
    """Parse command line arguements.

    Parameters:
        Nothing

    Returns:
        arguments: argparse.Namespace object
            An object containing all of the added arguments.

    Outputs:
        Nothing
    """

    input_help = 'Images to calculate bounds of: *dr?.fits *fl?.fits etc'
    plot_help = 'Plot bounding boxes of silhouettes?'

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default='f*dr?.fits')
    parser.add_argument('-p', help=plot_help, action='store_true',
        required=False)
    arguments = parser.parse_args()

    return arguments

def get_corners(im):
    hdr = fits.getheader(im,1)
    xmax = hdr['NAXIS1']-1
    ymax = hdr['NAXIS2']-1
    # print xmax,ymax
    data = fits.getdata(im)
    # print data.shape
    # print data[ymax,xmax]
    r,d = pixtosky.xy2rd(im+'[sci,1]',hms=False,x=[0,0,xmax,xmax],y=[0,ymax,ymax,0])
    if 'flt.fits' in im or 'flc.fits' in im:
        try:
            r2,d2 = pixtosky.xy2rd(im+'[sci,2]',hms=False,x=[0,0,xmax,xmax],y=[0,ymax,ymax,0])
            r[1] = r2[1]
            d[1] = d2[1]
            r[2] = r2[2]
            d[2] = d2[2]
        except:
            print 'FLT/FLC is not multichip'
    return r,d

def bounds(ra,dec):
    delta_ra = (max(ra)-min(ra))*math.cos(math.radians(min(np.absolute(dec))))
    delta_dec = max(dec)-min(dec)
    delta_ra_arcsec = delta_ra*3600.
    delta_dec_arcsec = delta_dec*3600.
    print 'OUTPUT IMAGE IS {}\" x {}\"'.format(delta_ra_arcsec, delta_dec_arcsec)
    print 'IR: {} x {} pix'.format(delta_ra_arcsec/0.1, delta_dec_arcsec/0.1)
    print 'UVIS+WFC: {} x {} pix'.format(delta_ra_arcsec/0.05, delta_dec_arcsec/0.05)
    return

def make_silhouette_plot(ra,dec,res, ims):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.path import Path
    import matplotlib.patches as patches
    colors = cm.rainbow(np.linspace(0, 1, len(res)))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n_images = len(ims)
    alpha = 1./float(n_images+1.)+.1
    for i in range(len(res)):
        verts = np.array(res[i]).T
        codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor=colors[i], lw=1, alpha=alpha)
        ax.add_patch(patch)
    ax.set_xlim(max(ra),min(ra))
    ax.set_ylim(min(dec),max(dec))
    plt.show()

if __name__ == '__main__':
    options = parse_args()
    ims = glob.glob(options.i)
    p = Pool(8)
    res = p.map(get_corners,ims)
    print ' ___________ '
    coords = np.concatenate([np.array(cor_poss).T for cor_poss in res])
    ra = coords[:,0]
    dec = coords[:,1]
    bounds(ra,dec)
    if options.p:
        make_silhouette_plot(ra,de,res,ims)

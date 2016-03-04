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

    input_help = 'Images to calculate bounds of. Default *dr?.fits'
    plot_help = 'Plot bounding boxes of silhouettes? Default False'
    write_help = 'Write out file with ra/dec size/midpt? Default False'

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help=input_help, action='store',
        required=False, default='f*dr?.fits', type=str)
    parser.add_argument('-p', help=plot_help, action='store_true',
        required=False)
    parser.add_argument('-w', help=write_help, action='store_true',
        required=False)
    arguments = parser.parse_args()
    return arguments

def get_corners(im):
    """Calculate ra/dec positions of image corners

    Parameters:
        im : string
            The name of the image on which to calculate corner positions

    Returns:
        (r,d) : tuple of ndarrays
            first element is array of corners' RAs, second is array of corners' Decs

    Outputs:
        Nothing
    """
    hdr = fits.getheader(im,1)
    xmax = hdr['NAXIS1']-1
    ymax = hdr['NAXIS2']-1
    data = fits.getdata(im)
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

def bounds(ra,dec,write=False):
    """Calculate RA/Dec bounding box properties from multiple RA/Dec points

    Parameters:
        ra : ndarray
            array of right ascension positions
        dec : ndarray
            array of declination positions
        write : bool
            switch controlling writing of file with bounding box parameters

    Returns:
        Nothing

    Outputs:
        dimensions.txt
            output file containing RA/Dec dimensions in arcsec, as well as
            RA/dec midpts in decimal degrees
    """
    delta_ra = (max(ra)-min(ra))*math.cos(math.radians(min(np.absolute(dec))))
    delta_dec = max(dec)-min(dec)
    delta_ra_arcsec = delta_ra*3600.
    delta_dec_arcsec = delta_dec*3600.
    ra_midpt = (max(ra)+min(ra))/2.
    dec_midpt = (max(dec)+min(dec))/2.
    if write:
        f = open('dimensions.txt','w')
        f.write('# final drizzle image sizing information\n')
        f.write('# first row is ra and dec size in arcsec (for outnx/outny param)\n')
        f.write('# second row is ra and dec midpts in deg (for final_ra/dec)\n')
        f.write('{} {}\n'.format(delta_ra_arcsec,delta_dec_arcsec))
        f.write('{} {}\n'.format(ra_midpt,dec_midpt))
        f.close()
    print 'OUTPUT IMAGE IS {}\" x {}\"'.format(delta_ra_arcsec, delta_dec_arcsec)
    print 'IR: {} x {} pix'.format(delta_ra_arcsec/0.1, delta_dec_arcsec/0.1)
    print 'UVIS+WFC: {} x {} pix'.format(delta_ra_arcsec/0.05, delta_dec_arcsec/0.05)
    return

def make_silhouette_plot(ra,dec,res):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.path import Path
    import matplotlib.patches as patches
    colors = cm.rainbow(np.linspace(0, 1, len(res)))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n_images = len(res)
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
    bounds(ra,dec,options.w)
    if options.p:
        make_silhouette_plot(ra,dec,res,ims)

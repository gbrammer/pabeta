import argparse
import glob
import math
import os
import shutil

from astropy.io import fits
from astropy.wcs import WCS
from drizzlepac import pixtosky
from multiprocessing import Pool
import numpy as np

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

    input_help = 'Images to calculate bounds of. Default f*dr?.fits'
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

def get_footprints(im):
    fps = []
    hdu = fits.open(im)
    flt_flag = 'flt.fits' in im or 'flc.fits' in im
    for ext in hdu:
        if 'SCI' in ext.name:
            hdr = ext.header
            wcs = WCS(hdr, hdu)
            fp = wcs.calc_footprint(hdr, undistort=flt_flag)
            fps.append(fp)
    return fps


def bounds(fp_list,write=False):
    """Calculate RA/Dec bounding box properties from multiple RA/Dec points

    Parameters:
        fp_list : list of list of ndarrays
            each ndarray is a 4x2 array corresponding to the sky positions
            of the 4 corners of a science extension.
            [[im1 ext1 corners, im1 ext2 corners], [im2 ext1 corners]]
        write : bool
            switch controlling writing of file with bounding box parameters

    Returns:
        Nothing

    Outputs:
        dimensions.txt
            output file containing RA/Dec dimensions in arcsec, as well as
            RA/dec midpts in decimal degrees
    """
    # flatten list of extensions into numpy array of all corner positions
    merged = []
    for im in fp_list:
        for ext in im:
            merged.append(ext)
    merged = np.vstack(merged)
    print merged.shape
    ra = merged[:,0]
    dec = merged[:,1]
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
    return ra, dec

def make_silhouette_plot(fp_list, ra, dec):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.path import Path
    import matplotlib.patches as patches

    colors = cm.rainbow(np.linspace(0, 1, len(fp_list)))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n_images = len(fp_list)
    alpha = 1./float(n_images+1.)+.1

    for i in range(len(fp_list)):
        im = fp_list[i]
        for ext in im:
            plot_patch(ax,ext,colors[i],alpha)
    ax.set_xlim(min(ra),max(ra))
    ax.set_ylim(min(dec),max(dec))
    plt.show()

def plot_patch(ax ,corners, color, alpha):
    from matplotlib.path import Path
    import matplotlib.patches as patches
    codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO]
    path = Path(corners, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=1, alpha=alpha)
    ax.add_patch(patch)


if __name__ == '__main__':
    options = parse_args()
    ims = glob.glob(options.i)
    p = Pool(8)
    fp_list = p.map(get_footprints, ims)
    ra, dec = bounds(fp_list,options.w)
    if options.p:
        make_silhouette_plot(fp_list, ra, dec)

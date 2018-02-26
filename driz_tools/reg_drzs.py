from drizzlepac import tweakreg
import glob, os
from astropy.io import fits
from multiprocessing import Pool
import numpy as np
from stsci.tools import teal
import argparse

import astropy.units as u
from astropy.table import Table
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia


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

    ref_help = 'Reference image to align other drizzled images to.  Default is first image'
    input_help = 'Images to register.  Default is all visit level drizzle images (f*dr?.fits)'
    refcat_help = 'Reference catalog to use.  Defaults is \'\' (Makes refcat from input images).'
    teal_help = 'Show teal interface for TweakReg?  Default False'
    gaia_help = 'Get source catalog from Gaia?  Default False'

    im = ''
    inp = 'f*dr?.fits'
    refcat = ''

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, help=ref_help, action='store',
        required=False, default=im)
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default=inp)
    parser.add_argument('-c', type=str, help=refcat_help, action='store',
        required=False, default=refcat)
    parser.add_argument('-t', help=teal_help, action='store_true',
        required=False)
    parser.add_argument('-g', help=gaia_help, action='store_true',
        required=False)
    arguments = parser.parse_args()
    return arguments

def get_gaia_cat(ims, cat_name='gaia'):
    """Get the Gaia catalog for the area of input images"""
    from calc_bounds import bounds, get_footprints

    print 'Calculating coordinate ranges for Gaia query:'
    footprint_list = map(get_footprints, ims)
    ras, decs = bounds(footprint_list)

    ra_midpt = (np.amax(ras)+np.amin(ras))/2.
    dec_midpt = (np.amax(decs)+np.amin(decs))/2.
    ra_width = (np.amax(ras)-np.amin(ras))
    dec_height = (np.amax(decs)-np.amin(decs))

    print '\nPerforming Gaia query:'
    coord = SkyCoord(ra=ra_midpt, dec=dec_midpt, unit=(u.degree, u.degree), frame='icrs')
    width = Quantity(ra_width, u.deg)
    height = Quantity(dec_height, u.deg)
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    print 'Sources returned: {}'.format(len(r))

    assert len(r) > 0, 'No sources found in Gaia query\n'

    cat_file_name = '{}.cat'.format(cat_name)
    print 'Writing Gaia source catalog: {}\n'.format(cat_file_name)
    Table([r['ra'], r['dec']]).write(cat_file_name, format='ascii.fast_commented_header')
    return cat_file_name

if __name__ == '__main__':
    ir_ims, vis_ims = [], []
    options = parse_args()

    # Gets the desired input images to align
    drzs = sorted(glob.glob(options.i))
    for drz in drzs: print drz

    assert len(drzs) > 0, 'No input images!\n'

    # Choose reference image for WCS tranforms
    if options.r != '':
        if os.path.exists(options.r):
            print 'Refence image is {}'.format(options.r)
            ref = options.r
        else:
            raise IOError('Reference input {} does not exist'.format(options.r))
    else:
        ref = drzs[0]

    # Generate lists of IR and optical images
    for f in drzs:
        hdr = fits.getheader(f)
        if hdr['DETECTOR'] == 'IR':
            ir_ims.append(f)
        elif fits.getval(f, 'DETECTOR') in ['UVIS', 'WFC', 'SBC', 'PC', 'HRC', 'FUV-MAMA']:
            vis_ims.append(f)

    # Set source finding parameters in reference images
    # only needed if NOT using external catalog
    if fits.getval(ref, 'DETECTOR') == 'IR':
        thresh = 15.
        cw = 2.5
    if fits.getval(ref, 'DETECTOR') == 'HRC':
        thresh = 4.
        cw = 4.
    elif fits.getval(ref, 'DETECTOR') in ['UVIS', 'WFC', 'SBC', 'PC']:
        thresh = 5.
        cw = 3.5

    cat = options.c

    # Determine WCS name
    if '_hsc.radec' in options.c:
        wcsname = 'HSC'
    elif options.c == 'gaia.cat':
        wcsname = 'GAIA'
    else:
        wcsname = 'TWEAK'

    if options.g:
        print '\nUsing Gaia source catalog as reference'
        cat = get_gaia_cat(drzs)
        wcsname = 'GAIA'

    # show teal interface if desired
    if options.t: teal.teal('tweakreg')


    # Set source detection Parameters
    # Set the conv_width values for optical/ir images
    vis_cw = 3.5
    ir_cw = 2.5

    # Set the thresholds for source detection for optical/ir images
    vis_thresh = 5.
    ir_thresh = 3.

    if len(vis_ims)>0:
        print 'Aligning optical images'
        tweakreg.TweakReg(sorted(vis_ims), updatehdr=True, expand_refcat=False,
        enforce_user_order=True,refimage=ref,
        imagefindcfg={'threshold':vis_thresh.,'conv_width':vis_cw},
        refimagefindcfg={'threshold':thresh,'conv_width':cw},
        refcat=cat,shiftfile=True,outshifts='vis_shifts.txt',
        wcsname=wcsname, interactive=False,see2dplot=False,
        fitgeometry='general',reusename=True,separation=0.0)

    if len(ir_ims)>0:
        print 'Aligning IR images'
        tweakreg.TweakReg(sorted(ir_ims), updatehdr=True, expand_refcat=False,
        enforce_user_order=True,refimage=ref,
        imagefindcfg={'threshold':ir_thresh,'conv_width':ir_cw},
        refimagefindcfg={'threshold':thresh,'conv_width':cw},
        refcat=cat,shiftfile=True,outshifts='ir_shifts.txt',
        wcsname=wcsname, interactive=False, see2dplot=False,
        fitgeometry='general',reusename=True,separation=0.0)

from drizzlepac import astrodrizzle, tweakreg, tweakback
from stwcs import updatewcs
import glob, os
from astropy.io import fits
from multiprocessing import Pool
from stsci.tools import teal
import argparse


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
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    ir_ims, vis_ims = [], []
    options = parse_args()
    drzs = glob.glob(options.i)
    for drz in drzs: print drz
    if len(drzs) == 0:
        print 'No input images'
        raise
    if options.r != '':
        if options.r in drzs:
            print 'Refence image is {}'.format(options.r)
            # drzs.remove(options.r)
        ref = options.r
    else:
        ref = drzs[0]
        # drzs.remove(drzs[0])
    for f in drzs:
        hdr = fits.getheader(f)
        if hdr['DETECTOR'] == 'IR':
            ir_ims.append(f)
        elif hdr['DETECTOR'] == 'UVIS' or hdr['DETECTOR'] == 'WFC':
            vis_ims.append(f)
    if fits.getval(ref, 'DETECTOR') == 'IR':
        thresh = 5.
        cw = 2.5
    elif fits.getval(ref, 'DETECTOR') == 'UVIS' or fits.getval(ref, 'DETECTOR') == 'WFC':
        thresh = 20.
        cw = 3.5

    # Determine WCS name
    if '_hsc.radec' in options.c:
        wcsname = 'HSC'
    else:
        wcsname = 'TWEAK'

    if options.t: teal.teal('tweakreg')

    if len(vis_ims)>0:
        tweakreg.TweakReg(vis_ims, updatehdr=True, expand_refcat=False,enforce_user_order=False,refimage=ref,
        imagefindcfg={'threshold':5.,'conv_width':3.5}, refimagefindcfg={'threshold':thresh,'conv_width':cw},
        refcat=options.c,shiftfile=True,outshifts='vis_shifts.txt', wcsname=wcsname)
    if len(ir_ims)>0:
        tweakreg.TweakReg(ir_ims, updatehdr=True, expand_refcat=False,enforce_user_order=True,refimage=ref,
        imagefindcfg={'threshold':5.,'conv_width':2.5}, refimagefindcfg={'threshold':thresh,'conv_width':cw},
        refcat=options.c,shiftfile=True,outshifts='ir_shifts.txt', wcsname=wcsname)

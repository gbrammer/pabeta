from drizzlepac import astrodrizzle, tweakreg, tweakback
from stwcs import updatewcs
import glob, os
from astropy.io import fits
from multiprocessing import Pool
from stsci.tools import teal
import argparse


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

    ref_help = 'Reference image to align other drizzled images to.  Default is first image'
    input_help = 'Images to register.  Default is all visit level drizzle images (f*dr?.fits)'

    im = ''
    inp = 'f*dr?.fits'

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, help=ref_help, action='store',
        required=False, default=im)
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default=inp)
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    teal.teal('tweakreg')
    ir_ims, vis_ims = [], []
    options = parse_args()
    drzs = glob.glob(options.i)
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
        thresh = 10.
        cw = 2.5
    elif fits.getval(ref, 'DETECTOR') == 'UVIS' or fits.getval(ref, 'DETECTOR') == 'WFC':
        thresh = 100.
        cw = 3.5
    if len(vis_ims)>0:
        tweakreg.TweakReg(vis_ims, updatehdr=True, expand_refcat=True,enforce_user_order=False,refimage=ref,
        imagefindcfg={'threshold':35.,'conv_width':3.5}, refimagefindcfg={'threshold':thresh,'conv_width':cw})
    if len(ir_ims)>0:
        tweakreg.TweakReg(ir_ims, updatehdr=True, expand_refcat=True,enforce_user_order=False,refimage=ref,
        imagefindcfg={'threshold':25.,'conv_width':2.5}, refimagefindcfg={'threshold':thresh,'conv_width':cw})

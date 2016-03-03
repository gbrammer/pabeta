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

    filts_help = 'Reference image to align other drizzled images to.  Default is first image'

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help=filts_help, nargs='+',
        required=False, default=None)
    # parser.add_argument('-i', type=str, help=input_help, action='store',
    #     required=False, default=inp)
    arguments = parser.parse_args()
    return arguments

if __name__ == '__main__':
    options = parse_args()
    if options.f != None:
        filts = [filt.upper() for filt in filts]
    filts=['F606W']
    fs = glob.glob('*fl?.fits')
    teal.teal('astrodrizzle')
    for filt in filts:
        exps = []
        for f in fs:
            expfil = fits.getheader(f)['FILTER']
            if expfil == filt:
                print f
                exps.append(f)
        hdr = fits.getheader(exps[0])
        inst = hdr['INSTRUME']
        targ = os.getcwd().split('/')[-1]
        det = hdr['DETECTOR']
        out = '_'.join([filt, targ, inst, det])
        astrodrizzle.AstroDrizzle(exps,output=out, mdriztab=False, num_cores=12, in_memory=False)

    # teal.teal('astrodrizzle')
    #
    # drzs = glob.glob('f[2-8]*dr?.fits')
    # exps = []
    # for im in drzs:
    #     exps += tweakback.extract_input_filenames(im)
    #     astrodrizzle.AstroDrizzle(exps,output='master', mdriztab=False, num_cores=12, in_memory=False, final_pixfrac=0.85, final_scale=0.03)

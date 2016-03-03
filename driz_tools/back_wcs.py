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

    input_help = 'Images to tweak back to flts. Default is f*dr?.fits'
    inp = 'f*dr?.fits'

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default=inp)
    arguments = parser.parse_args()
    return arguments

def tback(drz):
    flts = tweakback.extract_input_filenames(drz)
    print 'Tweaking back exposures for {} with input ims:'.format(drz)
    for f in flts:
        print f
    try:
        tweakback.tweakback(drz)
    except:
        return

if __name__ =='__main__':
    options = parse_args()
    drzs=glob.glob(options.i)
    p = Pool(32)
    p.map(tback, drzs)

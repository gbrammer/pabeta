import argparse
import glob
import numpy as np
import os

from astropy.io import fits
from multiprocessing import Pool
from PIL import Image

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

    input_help = 'Images to make WHT jpgs for. Default is F*dr?.fits'
    inp = 'F*dr?.fits'

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default=inp)
    arguments = parser.parse_args()
    return arguments

def bitrange(data):
    """Scale data from [min, max] to [0, 255] for 8 bit images"""
    lo = np.amin(data)
    hi = np.amax(data)
    data = (((data-lo)/(hi-lo))*254.999).astype('uint8')
    return data

def make_wht_jpg(im):
    """Open WHT FITS Image, scale to 8 bit and save ouput array as jpg"""
    wht_im = fits.getdata(im, 'WHT')
    wht_im = bitrange(wht_im)
    img = Image.fromarray(wht_im)
    im_root = '_'.join(im.split('_')[:-1])
    if img.mode != 'RGB':
        img = img.convert('RGB')
    img.save('wht_{}_.jpg'.format(im_root))
    print 'made wht image jpg for {}'.format(im_root)


if __name__ =='__main__':
    options = parse_args()
    drzs=glob.glob(options.i)
    p = Pool(32)
    p.map(make_wht_jpg, drzs)

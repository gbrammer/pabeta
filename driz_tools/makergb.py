from astropy.io import fits
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob
from PIL import Image
from multiprocessing import Pool
from ginga.util import zscale


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

    red_help = 'Red Image name'
    green_help = 'Green Image name'
    blue_help = 'Blue Image name'
    scale_help = 'Scale to use, either \"log\" (default), \"sqrt"\" or \"asinh\"'
    lo_help = 'Lower percentage threshold for clipping?  Default 0.1%%'
    hi_help = 'Upper percentage threshold for clipping?  Default 99.0%%'
    pseudog_help = 'Make pseudo green from red and blue ims?  Default False'
    deactivate_mp_help = 'Deactivate multiprocessing (recommended for large images)? Default False'
    zscale_help= 'Compute scale boundaries from zscasle algorithm (z1,100.*z2)?  Default False'

    parser = argparse.ArgumentParser()
    parser.add_argument('r', type=str, help=red_help, action='store')
    parser.add_argument('g', type=str, help=green_help, action='store')
    parser.add_argument('b', type=str, help=blue_help, action='store')

    parser.add_argument('-s', type=str, help=scale_help, action='store',
        required=False, default='log')
    parser.add_argument('-l', type=float, help=lo_help, action='store',
        required=False, default=.1)
    parser.add_argument('-u', type=float, help=hi_help, action='store',
        required=False, default=99.0)
    parser.add_argument('-p', help=pseudog_help, action='store_true',
        required=False)
    parser.add_argument('-d', help=deactivate_mp_help, action='store_true',
        required=False)
    parser.add_argument('-z', help=zscale_help, action='store_true',
        required=False)
    arguments = parser.parse_args()
    return arguments



def scl(data,lo=None,hi=None):
    if lo: print lo
    else: print 'No params provided, percent scaling'

def pctscale(data,lo=0.1,hi=99.5):
    '''Compute thresholds at lo/hi percent, and clip data below/above'''
    lv, hv = np.nanpercentile(data, [lo,hi])
    print 'Getting nanmask'
    nanmask = np.isnan(data)

    print 'Low val at {}%% level is {}'.format(lo,lv)
    print 'High val at {}%% level is {}'.format(hi,hv)
    data[nanmask] = lv
    data[data<lv] = lv
    data[data>hv] = hv
    return data

def do_zscale(data):
    z1, z2 = zscale.zscale(data)
    print 'Getting nanmask'
    nanmask = np.isnan(data)

    print 'Low val at z1 is {}'.format(z1)
    print 'High val at 100*z2 level is {}'.format(100.*z2)

    data[nanmask] = z1
    data = data.clip(z1, 100.*z2)
    return data

def bitrange(data):
    '''Remap data from [min,max] to [0,255] for 8 bit writing'''
    lo = np.amin(data)
    hi = np.amax(data)
    data = (((data-lo)/(hi-lo))*254.999).astype('uint8')
    return data

def asinh(data):
    '''Remap data from [min,max] to [0,1.0] and apply asinh scale'''
    lo = np.amin(data)
    hi = np.amax(data)
    data = ((data-lo)/(hi-lo))
    data = np.arcsinh(10.0*data)/3.0
    return data

def logscl(data):
    '''Remap data from [min,max] to [0,1.0] and apply log scale'''
    a = 1000.
    lo = np.amin(data)
    hi = np.amax(data)
    data = ((data-lo)/(hi-lo))
    data = np.log(a*data+1.)/np.log(a)
    return data

def sqrt(data):
    '''Remap data from [min,max] to [0,1.0] and apply sqrt scale'''
    lo = np.amin(data)
    hi = np.amax(data)
    data = ((data-lo)/(hi-lo))
    data = np.sqrt(data)
    return data

def wrap(im):
    if options.z:
        rscl = scale_method(do_zscale(im))
    else:
        rscl = scale_method(pctscale(im,lo,hi))
    r8bit = bitrange(rscl)
    # print im, 'DONE'
    print '\n'
    return r8bit


if __name__ == '__main__':
    # print drzs

    options = parse_args()

    rim, gim, bim = options.r, options.g, options.b

    lo = options.l
    hi = options.u
    ims = [fits.getdata(rim), fits.getdata(gim), fits.getdata(bim)]
    if not ims[0].shape == ims[1].shape == ims[2].shape:
        print 'DIMENSIONS OF IMAGES ARE NOT EQUAL'
        print 'This code does not yet support images of different sizes or pixel scales'
        raise

    if options.s == 'log':
        scale_method = logscl
    elif options.s == 'asinh':
        scale_method = asinh
    elif options.s == 'sqrt':
        scale_method = sqrt
    else:
        print 'Invalid scale method given, use -h for help.'
        raise

    if options.d:
        if options.p:  # skip making the green if using pseudogreen
            scls = [wrap(ims[0]),None,wrap(ims[2])]
        else:
            scls = map(wrap, ims)
    else:
        scls = Pool(3).map(wrap, ims)

    filt_names = [os.path.split(x)[-1].split('_')[0] for x in [rim, gim, bim]]
    if options.p:
        filt_names[1] = 'PSEUDOGREEN'
    filt_string = '_'.join(filt_names)

    if options.p:
        print 'Making pseudogreen from average red and blue\n'
        scls[1] = np.nanmean([scls[0],scls[2]],axis=0).astype('uint8')

    print 'Making image object'
    img = Image.fromarray(np.flipud(np.dstack(scls)))
    print 'Saving image'

    targ = os.path.split(os.getcwd())[-1]
    if options.z:
        jpg_name = '{}_{}_{}_zscale.jpeg'.format(targ,filt_string,scale_method.func_name)
    else:
        jpg_name = '{}_{}_{}_{}_{}.jpeg'.format(targ,filt_string,scale_method.func_name,lo,hi)
    img.save(jpg_name)
    print 'Image saved: {}'.format(jpg_name)

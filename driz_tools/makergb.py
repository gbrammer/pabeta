from astropy.io import fits
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob
from PIL import Image
from multiprocessing import Pool


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
    lo_help = 'Lower percentage threshold for clipping?  Deafult 0.1%%'
    hi_help = 'Upper percentage threshold for clipping?  Deafult 99.0%%'

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
    # print 'fraction of pixels as nan:', np.mean(nanmask.astype('float'))
#     plt.imshow(nanmask.astype('float'))
#     plt.show()
    if lv<0:
        lv = 0.
    print 'Low val at {}%% level is {}'.format(lo,lv)
    print 'High val at {}%% level is {}'.format(hi,hv)
    data[nanmask] = lv
    data[data<lv] = lv
    data[data>hv] = hv
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
    a = 50.
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
    rscl = scale_method(pctscale(im,lo,hi))
    r8bit = bitrange(rscl)
    # print im, 'DONE'
    return r8bit


if __name__ == '__main__':
    drzs = sorted(glob.glob('F*_dr?.fits'))
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
    scls = Pool(3).map(wrap, ims)
    # scls = map(wrap, ims)
    img = Image.fromarray(np.fliplr(np.dstack(scls)))
    print 'Making image'
    img.save('final_{}_{}_{}.jpeg'.format(scale_method.func_name,lo,hi))
    print 'Image made'

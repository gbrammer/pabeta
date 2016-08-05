from drizzlepac import astrodrizzle, tweakreg, tweakback
from stwcs import updatewcs
import glob, os
from astropy.io import fits
from multiprocessing import Pool
from stsci.tools import teal
import argparse
import math
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

    filts_help = 'Filters to make final products. Default is all filters'
    hipeec_help = 'USE Hi-PEEC scales: UVIS/WFC 0.04 IR 0.12 SBC 0.025. Default False'

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help=filts_help, nargs='+',
        required=False, default=None)
    parser.add_argument('-hipeec', help=hipeec_help, action='store_true',
        required=False)
    # parser.add_argument('-i', type=str, help=input_help, action='store',
    #     required=False, default=inp)
    arguments = parser.parse_args()
    return arguments

def wfc3_exp_filt(im):
    """Return filter used in WFC3 image"""
    filt = fits.getval(im,'FILTER')
    return (im,filt)

def acs_exp_filt(im):
    """Return real (non CLEAR*) filter used in ACS image"""
    hdr = fits.getheader(im)
    filt1 = hdr['FILTER1']
    filt2 = hdr['FILTER2']
    filt = filt1
    if filt1 == 'CLEAR1L' or filt1 == 'CLEAR1S':
        filt = filt2
    return (im,filt)

def make_filt_dict(im_filts):
    """Converts tuple list from (im,filter) to filter:[im1,im2,...]"""
    im_to_filt = {}
    for ln in im_filts:
        im_to_filt[ln[0]] = ln[1]
    filt_to_im_list = {}
    for filt in set(im_to_filt.values()):
        tmp = []
        for ln in im_filts:
            if ln[1] == filt:
                tmp.append(ln[0])
        filt_to_im_list[filt] = tmp
    return filt_to_im_list

def parse_filters(filts=None):
    """Structures exposures into list of lists of exposures for each filter

    Makes list of lists of exposures. Each individual list is the list of
    exposures for a single ACS or WFC3 filter.  The argument filts dictates
    which filters to return the exposures for.  This does not overcome filter
    name degeneracy, so the list of lists will contain a list of WFC3 F814W
    images as well as a list of ACS F814W images (though the exposures
    for each are inseparate lists).  If filts is None will make list of lists
    for all filters of data in current directory

    Parameters:
        filts : list of strings
            list of filter names to get the exposures for.  If None (default)
            all filters in current directory are used.

    Returns:
        exps_by_filt : list of lists of strings
            Each individual list is a list of image filenames for exposures
            taken in a single filter

    Outputs: Nothing
    """
    p = Pool(32)
    acs_ims = glob.glob('j*fl?.fits')
    wfc3_ims = glob.glob('i*fl?.fits')
    exps_by_filt = []
    if len(wfc3_ims) > 0:
        wfc3_im_filts = p.map(wfc3_exp_filt,wfc3_ims)
        wfc3_dict = make_filt_dict(wfc3_im_filts)
        if filts != None:
            for filt in wfc3_dict.keys():
                if filt not in filts: del wfc3_dict[filt]
        exps_by_filt += wfc3_dict.values()
    if len(acs_ims) > 0:
        acs_im_filts = p.map(acs_exp_filt,acs_ims)
        acs_dict = make_filt_dict(acs_im_filts)
        if filts != None:
            for filt in acs_dict.keys():
                if filt not in filts: del acs_dict[filt]
        exps_by_filt += acs_dict.values()

    return exps_by_filt

def final_drizzle(exps):
    hdr = fits.getheader(exps[0])
    inst = hdr['INSTRUME']
    if inst == 'ACS':
        filt1 = hdr['FILTER1']
        filt2 = hdr['FILTER2']
        filt = filt1
        if filt1 == 'CLEAR1L' or filt1 == 'CLEAR1S':
            filt = filt2
    else: filt = hdr['FILTER']
    targ = os.getcwd().split('/')[-1]
    det = hdr['DETECTOR']
    out = '_'.join([filt, targ, inst, det])

    if os.path.exists(out+'_drz.fits') or os.path.exists(out+'_drc.fits'): return

    print out, len(exps)

    combine_nhigh = 1
    med_alg = 'iminmed'
    if len(exps) > 4:
        med_alg = 'imedian'
        if len(exps) > 7:
            combine_nhigh = 3


    if det == 'IR':
        scl = 0.1
    elif det == 'WFC' or det == 'UVIS':
        scl = 0.05
    elif det == 'SBC':
        scl = 0.025

    if options.hipeec:
        print 'USING HIPEEC PIXEL SCALES'
        if det == 'IR':
            scl = 0.12
        elif det == 'WFC' or det == 'UVIS':
            scl = 0.04
        elif det == 'SBC':
            scl = 0.025

    if os.path.exists('dimensions.txt'):
        dims = np.loadtxt('dimensions.txt')
        outnx = int(math.ceil(dims[0,0])/scl)
        outny = int(math.ceil(dims[0,1])/scl)
        ra = dims[1,0]
        dec = dims[1,1]
        rot=0.
    else:
        outnx, outny, ra, dec, rot = None, None, None, None, None
    if det == 'IR' or det == 'SBC' or len(exps)==1:
        astrodrizzle.AstroDrizzle(exps,output=out, mdriztab=False, num_cores=1,
                                in_memory=False,final_wcs=True,final_rot=rot,
                                final_outnx=outnx,final_outny=outny, final_ra=ra,
                                final_dec=dec,final_scale=scl,median=False,
                                blot=False,driz_cr=False,runfile='ADRIZ_{}'.format(out),
                                clean=True,build=True)
    else:
        astrodrizzle.AstroDrizzle(exps,output=out, mdriztab=False, num_cores=1,
                                in_memory=False,final_wcs=True,final_rot=rot,
                                final_outnx=outnx,final_outny=outny, final_ra=ra,
                                final_dec=dec,final_scale=scl,combine_type=med_alg,
                                combine_nhigh=combine_nhigh,runfile='ADRIZ_{}'.format(out),
                                clean=True,build=True)

    input_wcs = fits.getval(exps[0],'wcsname',1)
    if input_wcs == 'HSC':
        prod_name = glob.glob('{}_dr?.fits'.format(out))[0]
        fits.setval(prod_name,keyword='wcsname',value='HSC',extname='sci')


if __name__ == '__main__':
    options = parse_args()
    filts = options.f
    if options.f != None:
        filts = [filt.upper() for filt in options.f]
    exps_by_filt = parse_filters(filts)
    teal.teal('astrodrizzle')
    p = Pool(24)
    p.map(final_drizzle,exps_by_filt)
    # map(final_drizzle,exps_by_filt)

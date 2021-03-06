from drizzlepac import astrodrizzle
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
    hipeec_help = 'USE Hi-PEEC scales: UVIS/WFC 0.04 IR 0.12 SBC 0.03. Default False'
    teal_help = 'Show teal interface for AstroDrizzle?  Default False'
    scl_help = 'Use default camera pixel scales?  Default False.'
    ncore_help = 'Number of cores to use with parallel processing? Default 32'

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help=filts_help, nargs='+',
        required=False, default=None)
    parser.add_argument('-hipeec', help=hipeec_help, action='store_true',
        required=False)
    parser.add_argument('-t', help=teal_help, action='store_true',
        required=False)
    parser.add_argument('-d', help=scl_help, action='store_true',
        required=False)
    parser.add_argument('-n', help=ncore_help, action='store',
        required=False, default=32, type=int)
    # parser.add_argument('-i', type=str, help=input_help, action='store',
    #     required=False, default=inp)
    arguments = parser.parse_args()
    return arguments

def wfc3_exp_filt(im):
    """Return filter used in WFC3 image"""
    filt = fits.getval(im,'FILTER')
    return (im,filt)

def wfpc2_exp_filt(im):
    """Return filter used in WFPC2 image"""
    filt = fits.getval(im,'FILTNAM1')
    return (im,filt)

def wfpc2_hdr_corr(im):
    """Add information to WFPC2 Data Global Header"""
    hdu = fits.open(im,mode='update')
    det = hdu[1].header['DETECTOR']
    if det == 1:
        hdu[0].header['DETECTOR'] = 'PC'
    elif det in [2, 3, 4]:
        hdu[0].header['DETECTOR'] = 'WFC'
    hdu.close()

def acs_exp_filt(im):
    """Return real (non CLEAR*) filter used in ACS image"""
    hdr = fits.getheader(im)
    filt1 = hdr['FILTER1']
    filt2 = hdr['FILTER2']
    det = hdr['DETECTOR']
    filt = filt1
    if filt1 == 'CLEAR1L' or filt1 == 'CLEAR1S':
        filt = filt2
    filt+=det
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
    wfpc2_ims = glob.glob('u*c0m.fits')
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
    if len(wfpc2_ims) > 0:
        wfpc2_im_filts = p.map(wfpc2_exp_filt,wfpc2_ims)
        p.map(wfpc2_hdr_corr,wfpc2_ims)
        wfpc2_dict = make_filt_dict(wfpc2_im_filts)
        if filts != None:
            for filt in wfpc2_dict.keys():
                if filt not in filts: del wfpc2_dict[filt]
        exps_by_filt += wfpc2_dict.values()

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
    elif inst == 'WFPC2':
        filt = hdr['FILTNAM1']
    elif inst == 'STIS':
        filt = hdr['OPT_ELEM']
    else: filt = hdr['FILTER']
    targ = os.getcwd().split('/')[-1]
    det = hdr['DETECTOR']
    out = '_'.join([filt, targ, inst, det])

    if os.path.exists(out+'_drz.fits') or os.path.exists(out+'_drc.fits'):
        print 'file exists: ', out+'_drz.fits'
        return

    print out, len(exps)

    combine_nhigh = 1
    med_alg = 'iminmed'
    if len(exps) > 150:
        med_alg = 'imedian'
        if len(exps) > 7:
            combine_nhigh = 3


    if det == 'IR':
        scl = 0.1
    elif det == 'WFC' or det == 'UVIS' or det == 'PC':
        scl = 0.05
    elif det == 'SBC' or det == 'HRC' or det == 'FUV-MAMA':
        scl = 0.025

    if options.hipeec:
        print 'USING HIPEEC PIXEL SCALES'
        if det == 'IR':
            scl = 0.12
        elif det == 'WFC' or det == 'UVIS' or det == 'PC':
            scl = 0.04
        elif det == 'SBC' or det == 'HRC' or det == 'FUV-MAMA':
            scl = 0.03

    if options.d:
        if det == 'IR':
            scl = 0.13
        elif det == 'WFC' or det == 'PC':
            scl = 0.05
        elif det == 'UVIS':
            scl = 0.04
        elif det == 'SBC' or det == 'HRC' or det == 'FUV-MAMA':
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

    mem, nc = True, 1
    if options.n == 0:
        nc = None
        mem = False
    if det == 'IR' or det == 'SBC' or det == 'FUV-MAMA' or len(exps)==1:
        astrodrizzle.AstroDrizzle(exps,output=out, mdriztab=False, num_cores=nc,
                                in_memory=mem,final_wcs=True,final_rot=rot,
                                final_outnx=outnx,final_outny=outny, final_ra=ra,
                                final_dec=dec,final_scale=scl,median=False,
                                blot=False,driz_cr=False,runfile='ADRIZ_{}'.format(out),
                                clean=True,build=True, context=False)
    elif det == 'PC':
        astrodrizzle.AstroDrizzle(exps,output=out, mdriztab=False, num_cores=nc,
                                in_memory=mem,final_wcs=True,final_rot=rot,
                                final_outnx=outnx,final_outny=outny, final_ra=ra,
                                final_dec=dec,final_scale=scl,combine_type=med_alg,
                                combine_nhigh=combine_nhigh,runfile='ADRIZ_{}'.format(out),
                                clean=True,build=True, driz_cr_snr='5.5 3.5',
                                driz_cr_scale='2.0 1.5', context=False)
    else:
        astrodrizzle.AstroDrizzle(exps,output=out, mdriztab=False, num_cores=nc,
                                in_memory=mem,final_wcs=True,final_rot=rot,
                                final_outnx=outnx,final_outny=outny, final_ra=ra,
                                final_dec=dec,final_scale=scl,combine_type=med_alg,
                                combine_nhigh=combine_nhigh,runfile='ADRIZ_{}'.format(out),
                                clean=True,build=True, context=False)

    input_wcs = fits.getval(exps[0],'wcsname',1)
    if input_wcs == 'HSC':
        prod_name = glob.glob('{}_dr?.fits'.format(out))[0]
        fits.setval(prod_name,keyword='wcsname',value='HSC',extname='sci')


if __name__ == '__main__':
    options = parse_args()
    filts = options.f
    if options.f != None:
        filts = [filt.upper() for filt in filts]
    exps_by_filt = parse_filters(filts)
    if options.t:
        teal.teal('astrodrizzle')

    if options.n != 0:
        p = Pool(options.n)
        p.map(final_drizzle,exps_by_filt)
    else:
        map(final_drizzle,exps_by_filt)

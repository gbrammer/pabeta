from astropy.io import fits
import glob, os, shutil
from multiprocessing import Pool


def up_keywords(im):
    with fits.open(im, mode='update') as image:
        hdr = image[0].header
        inst = hdr['INSTRUME']
        dete = hdr['DETECTOR']
        if inst != 'ACS' or dete != 'WFC':
            print '{} IS NOT AN ACS/WFC IMAGE'.format(im)
            return
        print 'Updating image for {} {} {}'.format(im,inst,dete)
        path = '/astro/pabeta/distortionfiles/'
        npols = {'F435W': 'wfc_update3_CLEAR1L_F435W_npl.fits',
             'F475W': 'wfc_update3_F475W_CLEAR2L_npl.fits',
             'F502N': 'wfc_update3_F502N_CLEAR2L_npl.fits',
             'F550M': 'wfc_update3_F550M_CLEAR2L_npl.fits',
             'F555W': 'wfc_update3_F555W_CLEAR2L_npl.fits',
             'F606W': 'wfc_update3_F606W_CLEAR2L_npl.fits',
             'F625W': 'wfc_update3_F625W_CLEAR2L_npl.fits',
             'F658N': 'wfc_update3_F658N_CLEAR2L_npl.fits',
             'F660N': 'wfc_update3_CLEAR1L_F660N_npl.fits',
             'F775W': 'wfc_update3_F775W_CLEAR2L_npl.fits',
             'F814W': 'wfc_update3_CLEAR1L_F814W_npl.fits',
             'F850LP': 'wfc_update3_F850LP_CLEAR2L_npl.fits'}
        d2i = 'wfc_update3_d2im.fits'
        presm4 = 'wfc_update3_pre_idc.fits'
        postsm4 = 'wfc_update3_post_idc.fits'

        filt1 = hdr['FILTER1']
        filt2 = hdr['FILTER2']
        filt = filt1
        if filt1 == 'CLEAR1L' or filt1 == 'CLEAR1S':
            filt = filt2
        hdr['NPOLFILE'] = path + npols[filt]
        date = hdr['DATE-OBS']
        year = int(date.split('-')[0])
        if year >= 2009:
            hdr['IDCTAB'] = path + postsm4
        else:
            hdr['IDCTAB'] = path + presm4
        hdr['D2IMFILE'] = path + d2i
    return

if __name__ == '__main__':
    ims = glob.glob('j*fl?.fits')
    Pool(8).map(up_keywords,ims)

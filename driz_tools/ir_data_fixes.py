import glob, os, shutil
from astropy.io import fits
import astropy.io.fits as pyfits # for compatibility...
from multiprocessing import Pool
import argparse
import shutil
import numpy as np
import numpy.ma
import wfc3tools
from datetime import datetime

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

    copy_help = 'Copy files back to original target directory? Default False'

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', help=copy_help, action='store_true',
        required=False)
    arguments = parser.parse_args()
    return arguments

def fix_ramp(raw):
    if fits.getval(raw,'filter') != 'F110W' or os.path.exists(raw.replace('raw.fits','flb.fits')) or fits.getval(raw,'subarray') == 'T':
        return
    else:
        print 'Flattening Ramp for: {}'.format(raw)
        make_IMA_FLT(raw)
        # shutil.move(raw.replace('raw.fits','ima.tra'),'ex/')
        # shutil.move(raw.replace('raw.fits','.tra'),'ex/')
        flt = raw.replace('raw.fits','flt.fits')
        shutil.move(flt,flt.replace('flt.fits','flb.fits'))
        print 'FLATTENING COMPLETE'
        return


def update_flt(raw):
    # Copy Pipeline flt to ir_aux_files directory
    # Add Gabe's extra bad pixels to FLT DQ
    # Flag persistence pixels (>0.6*ERR) in FLT DQ if persist.fits exists
    # Replace FLT data with flattened ramp from make_IMA_FLT
    # Add header keywords indicating when this script was run
    # Write fixes into original FLT to preserve WCS transformations
    flt = raw.replace('raw.fits', 'flt.fits')
    if not os.path.exists(flt):
        try:
            orig_flt = glob.glob('../targets/*/{}'.format(flt))[0]
            shutil.copy(orig_flt, '.')
        except:
            print raw, fits.getval(raw,'TARGNAME'), 'NOT FOUND, REMOVING'
            os.remove(raw)
            return
        # print 'BEGINNING FLT UPDATE'
        # print '________________________________________________________'

    else:
        print 'FLT {} already in directory, skipping'.format(flt)
        return flt
    hdu = fits.open(flt, mode='update')
    orig_dq = hdu['DQ'].data
    new_flags = fits.getdata('/astro/pabeta/wfc3/data/badpix_spars200_Nov9.fits')
    if orig_dq.shape == (1014,1014):
        new_dq = np.bitwise_or(orig_dq,new_flags)
    else:
        hdr1 = hdu[1].header
        ltv1, ltv2 = abs(hdr1['LTV1']), abs(hdr1['LTV2'])
        naxis1, naxis2 = hdr1['NAXIS1'], hdr1['NAXIS2']
        flags_sub = new_flags[ltv2:ltv2+naxis2,ltv1:ltv1+naxis1]
        new_dq = np.bitwise_or(orig_dq,flags_sub)
    today = datetime.today()
    date = '{}-{}-{}'.format(today.year,today.month,today.day)
    hdu[0].header['BPIXFLAG'] = date

    proposid = hdu[0].header['PROPOSID']
    root = hdu[0].header['ROOTNAME']
    vis = root[4:6]
    pers_path_potential = '/grp/hst/wfc3a/GO_Links/{}/Visit{}/Persist/{}'.format(proposid,vis,root.lower()+'_persist.fits')
    pers_path = glob.glob(pers_path_potential)
    if len(pers_path)==1:
        pers_path = pers_path[0]
        err = hdu['ERR'].data
        pers_flags = np.zeros(err.shape, dtype=np.uint16)
        pers_data = fits.getdata(pers_path)
        pers_flags[pers_data>0.6*err] = 1024
        new_dq = np.bitwise_or(new_dq,pers_flags)
        hdu[0].header['PERSFLAG'] = date
        print 'FLAGGING PERSISTENCE'
    else: hdu[0].header['PERSFLAG'] = 'NONE'
    print 'ADDING BAD PIXELS TO DQ ARRAY FOR {}'.format(flt)
    hdu['DQ'].data = new_dq
    flb = flt.replace('flt.fits','flb.fits')
    if os.path.exists(flb):
        flb_data = fits.getdata(flb,1)
        hdu['SCI'].data = flb_data
        print 'REPLACING PIPELINE FLT DATA WITH FLATTENED RAMP FOR {}'.format(flt)
        hdu[0].header['FLATRAMP'] = date
    hdu.close()
    return flt


def check_flt(flt_path):
    flt = os.path.split(flt_path)[-1]
    det = fits.getval(flt_path,'DETECTOR')
    filt = fits.getval(flt_path,'FILTER')
    if det == 'IR' and filt == 'F110W':
        if not os.path.exists(flt.replace('flt.fits','raw.fits')):
            print 'NO RAW FOR {}'.format(flt)
            asn_id = fits.getval(flt_path, 'asn_id')
            if asn_id == 'NONE':
                asn_id = fits.getval(flt_path, 'rootname').upper()
            return asn_id
            # raise
    return None

def get_targ_dir(flt):
    try: orig_flt = glob.glob('../targets/*/{}'.format(flt))[0]
    except: print 'FLT IS:', flt
    targ_dir = os.path.split(orig_flt)[0]
    return targ_dir

def backup_flts(targ_dir):
    # Copies 'unfixed' flts into backup folder within target directory
    print targ_dir
    backup_dir = '{}/backup_flts'.format(targ_dir)
    if not os.path.exists(backup_dir):
        os.mkdir(backup_dir)
    for flt_path in glob.glob('{}/i*flt.fits'.format(targ_dir)):
        flt = os.path.split(flt_path)[-1]
        if os.path.exists(flt):
            if os.path.exists('{}/{}'.format(backup_dir,flt)):
                print 'BACKUP EXISTS OR NON-ORIGINAL FOR {}'.format(flt_path)
                continue
            else:
                shutil.move(flt_path,backup_dir)
                shutil.copy(flt,targ_dir)
        else: continue
    return

def copy_fixed_flts(flts):
    p = Pool(32)
    print flts
    targ_dirs = set(p.map(get_targ_dir,flts))
    p.map(backup_flts, targ_dirs)

if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/astro/pabeta/wfc3/')
    from reprocess_wfc3 import make_IMA_FLT, fetch_calibs
    raws = glob.glob('*raw.fits')
    p = Pool(32)
    missing_ids = p.map(check_flt, glob.glob('/astro/pabeta/targets/*/i*flt.fits'))
    p.map(fix_ramp, raws)
    fixed_flts = set(p.map(update_flt, raws))
    if None in fixed_flts:
        fixed_flts.remove(None)
    copy_fixed_flts(fixed_flts)

    if missing_ids.count(None) == len(missing_ids):
        print 'ALL RAWS HERE'
    else:
        missing_ids = set(missing_ids)
        missing_ids.remove(None)
        print 'MISSING RAWS:'
        for mid in missing_ids: print '{}, '.format(mid)

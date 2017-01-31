from drizzlepac import tweakreg, tweakback, updatehdr
from astropy.io import fits
from astropy import wcs
from multiprocessing import Pool
import argparse
import glob, os, shutil
from back_wcs import tback
from stsci.tools import teal

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

    input_help = 'Image to register.  Shifts calculated will apply to all mosaics.  Default is F814W*drc.fits'
    align_help = 'Calculate alignments using TweakReg?  Default False.'
    teal_help = 'Show teal interface for TweakReg?  Default False'

    inp = 'F814W*drc.fits'


    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help=input_help, action='store',
        required=False, default=inp)
    parser.add_argument('-a', help=align_help, action='store_true',
        required=False)
    parser.add_argument('-t', help=teal_help, action='store_true',
            required=False)
    arguments = parser.parse_args()
    return arguments


def calc_shift(im, radec):
    hdr = fits.getheader(im)
    det = hdr['DETECTOR']
    if det == 'WFC':
        cw = 3.5
    elif det == 'UVIS':
        cw = 3.5*(0.04/0.05)
    elif det == 'IR':
        cw = 2.5*(0.13/.1)
    tweakreg.TweakReg(im, updatehdr=False, expand_refcat=False,enforce_user_order=True,
    imagefindcfg={'threshold':5.,'conv_width':cw},refcat=radec,shiftfile=True,
    outshifts='single_shift.txt',outwcs='hsc_shifts_wcs.fits', refimage=None)

def shift_wrap():
    good_shift = raw_input('Is shift good? y/[n]  ')
    if good_shift == 'y':
        untweaked, shiftfile = write_full_shiftfile()
        updatehdr.update_from_shiftfile(shiftfile,wcsname='HSC')
        # Pool(40).map(shift_wcs,all_ims) # Doesn't work because of bug in updatewcs_with_shift :(
        Pool(32).map(tback, untweaked)
    else:
        print 'Shift solution rejected/invalid character'
        return

# Workaround for updatewcs_with_shift errors, write shiftfile and use update_from shiftfile
def write_full_shiftfile(input_shiftfile='single_shift.txt',output='hsc_shifts.txt'):
    # Writes shiftfile for all top level mosaics, using transforms from HSC TweakReg
    f = open(input_shiftfile)
    shift_lines = [x.strip() for x in f.readlines()]
    f.close()

    last_line = shift_lines[-1]
    filename = last_line.split()[0]
    mosaics = glob.glob('F*dr?.fits')

    g = open(output, 'w')
    for ln in shift_lines[:-1]: g.write('{}\n'.format(ln))

    untweaked = []
    for im in glob.glob('F*dr?.fits'):
        current_wcs = wcs.WCS(fits.getheader(im,1)).wcs.name
        if current_wcs == 'HSC':
            print '{} already aligned to HSC, skipping'.format(im)
            continue
        entry = last_line.replace(filename, im)
        g.write('{}\n'.format(entry))
        untweaked.append(im)
    g.close()
    return untweaked, output

def shift_wcs(im):
    try:
        im_wcs = wcs.WCS(fits.getheader(im,1))
    except:
        raise AssertionError('COULDN\'T FIND WCS FOR {}'.format(im))
    if im_wcs.wcs.name == 'HSC':
        print '{} already aligned to HSC, skipping'.format(im)
        return
    f = [ln.strip('\n') for ln in open('hsc_shifts.txt').readlines()][-1]
    xsh, ysh, rot, scl, xrms, yrms = f.split()[1:]
    ref = 'hsc_shifts_wcs.fits'
    print im
    try:
        updatehdr.updatewcs_with_shift(im,ref,xsh=xsh,ysh=ysh,rot=rot,scale=scl,wcsname='HSC',force=False,sciext='SCI')
        tweakback.tweakback(im,force=False,origwcs='DRZWCS',newname='HSC',wcsname='HSC')
    except:
        print 'COULDN\'T UPDATE {}, ERRORING'.format(im)
        raise

if __name__ == '__main__':
    options = parse_args()
    ims = glob.glob(options.i)
    if len(ims) != 1:
        print 'ERROR MULTIPLE OR NO IMAGES SELECTED FOR CALCULATING SHIFTS'
        raise ValueError('Too many input arguments')
    targ = os.path.split(os.getcwd())[-1]
    cat_path = '/astro/pabeta/HSC_Catalogs'
    all_cats = glob.glob('{}/*_hsc.radec'.format(cat_path))
    cats = []

    if not os.path.exists(targ+'_hsc.radec'):
        for derp in all_cats:
            if targ in derp.upper():
                print derp
                cats.append(derp)
        if len(cats) > 1:
            print 'TOO MANY CATALOGS MATCHED'
            raise
        else:
            cat = cats[0]
            shutil.copy(cat,'.')
            cat = os.path.split(cat)[-1]
            print cat

    else: cat = targ+'_hsc.radec'

    if options.t: teal.teal('tweakreg')

    if os.path.exists(cat):
        if options.a: calc_shift(ims[0],cat)
        shift_wrap()
    else:
        print 'NO HSC ENTRY IN CATALOG DIRECTORY FOR {}'.format(targ)

from drizzlepac import astrodrizzle, tweakreg, tweakback
from stwcs import updatewcs
import glob, os, shutil
from astropy.io import fits
from multiprocessing import Pool
from stsci.tools import teal

def upwcs(f):
  print f
  updatewcs.updatewcs(f)

def driz(exps):
    hdr = fits.getheader(exps[0])
    if hdr['INSTRUME'] == 'WFC3':
        filt = hdr['FILTER']
    elif hdr['INSTRUME'] == 'ACS':
        filt = hdr['FILTER1']
        if filt == 'CLEAR1L' or filt == 'CLEAR1S':
            filt = hdr['FILTER2']
    targ = hdr['TARGNAME']
    prop = hdr['PROPOSID']
    vis = exps[0][4:6]
    inst = hdr['INSTRUME']
    dete = hdr['DETECTOR']
    pvis = exps[0][1:6]
    print targ, filt, len(exps)
    fname = '_'.join([filt,targ,inst,dete,pvis]).lower()

    med_alg = 'minmed'
    if len(exps) > 4:
        med_alg = 'median'

    if os.path.exists(fname+'_drz.fits') or os.path.exists(fname+'_drc.fits'): return
    if fname == 'f814w_ngc3034_acs_wfc_8mx19':
        # astrodrizzle.AstroDrizzle(exps,output=fname,num_cores=1, in_memory=False, clean=True, build=True, combine_type='median', runfile=fname)
        return
    if hdr['DETECTOR'] == 'WFC':
        astrodrizzle.AstroDrizzle(exps,output=fname,num_cores=1, in_memory=False, clean=True, build=True, combine_type=med_alg, runfile=fname)
    else:
        astrodrizzle.AstroDrizzle(exps,output=fname,num_cores=1, mdriztab=True, in_memory=False, clean=True, build=True, combine_type=med_alg, runfile=fname)

def parse_wfc3_visit(visit):
    # print 'Processing orbit {} of {}'.format(i+1,len(wfc3_orbs))
    # i+=1
    tmp = glob.glob(visit+'*fl?.fits')
    hdr = fits.getheader(tmp[0])
    print 'Target appears to be {} for visit {} with {} exposures'.format(hdr['TARGNAME'], visit, len(tmp))
    filts = set([fits.getheader(f)['FILTER'] for f in tmp])
    exps_by_filt = []
    for filt in filts:
        exps = []
        for f in tmp:
            if fits.getheader(f)['FILTER'] == filt:
                exps.append(f)
        exps_by_filt.append(exps)
    return exps_by_filt

def parse_acs_visit(visit):
    tmp = glob.glob(visit+'*fl?.fits')
    hdr = fits.getheader(tmp[0])
    print 'Target appears to be {} for visit {} with {} exposures'.format(hdr['TARGNAME'], visit, len(tmp))
    acs_filts = []
    for f in tmp:
        fhdr = fits.getheader(f)
        acs_filts.append(fhdr['FILTER1'])
        acs_filts.append(fhdr['FILTER2'])
    filts = set([filt for filt in acs_filts if filt[0]=='F'])

    exps_by_filt = []
    for filt in filts:
        exps = []
        for f in tmp:
            if fits.getheader(f)['FILTER1'] == filt or fits.getheader(f)['FILTER2'] == filt:
                exps.append(f)
        exps_by_filt.append(exps)
    return exps_by_filt

if __name__ == '__main__':
    for flc in glob.glob('*flc.fits'):
        if os.path.exists(flc.replace('flc.fits','flt.fits')):
            if not os.path.exists('extra_flts'):
                os.mkdir('extra_flts')
            shutil.move(flc.replace('flc.fits','flt.fits'),'extra_flts')
    ims = sorted(glob.glob('*fl?.fits'))

    wfc3ims = sorted(glob.glob('i*fl?.fits'))
    wfc3_orbs = list(set([i[:6] for i in wfc3ims]))

    acsims = sorted(glob.glob('j*fl?.fits'))
    acs_orbs = list(set([i[:6] for i in acsims]))

    visbyfilt = []

    p = Pool(8)
    visbyfilt += p.map(parse_wfc3_visit,wfc3_orbs)
    visbyfilt += p.map(parse_acs_visit,acs_orbs)

    filter_visits = []
    for block in visbyfilt:
        for ln in block:
            filter_visits.append(ln)

    print '______________________________'
    print '______________________________'
    print 'Updating WCS'
    print '______________________________'
    print '______________________________'
    p.map(upwcs, ims)

    print '______________________________'
    print '______________________________'
    print 'Drizzling'
    print '______________________________'
    print '______________________________'
    print len(filter_visits)
    teal.teal('astrodrizzle')
    p.map(driz, filter_visits)

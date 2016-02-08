import unicorn
import glob
import os
import numpy as np
import threedhst
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

def prep():
    
    import glob
    import os
    
    import threedhst
    import threedhst.prep_flt_astrodrizzle as init
    from threedhst import catIO
    
    import unicorn
    import research.hawkiff
    
    ### Make ACS associations    
    unicorn.candels.make_asn_files(uniquename=True, translate={'-ROT':''})
    threedhst.options['FLT_PERSISTENCE_PATH'] = '../Persistence/'
    
    #### IR
    files=glob.glob('*-F1*asn.fits')
    radec=None
    for asn_file in files:
        init.prep_direct_grism_pair(direct_asn=asn_file, grism_asn=None, radec=radec, ACS=False, align_threshold=5, raw_path='../RAW/', order=0)
    
    ###
    for asn_file in files:
        drizzlepac.astrodrizzle.AstroDrizzle(asn_file, static=False, skysub=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_refimage='IRAS23436+5257-87-304-F814W_drc_sci.fits', final_pixfrac=1, context=False, resetbits=0, final_bits=576, preserve=False)
    
    ### Match WCS of NB observations
    init.copy_adriz_headerlets(glob.glob('*-F1*[MW]*asn.fits')[0], glob.glob('*-F1*[N]*asn.fits')[0], order=[1,0])
    
    ### Flag negative pixels    
    for asn_file in files:
        asn = threedhst.utils.ASNFile(asn_file)
        for exp in asn.exposures:
            im = pyfits.open('%s_flt.fits' %(exp), mode='update')
            neg = im['SCI'].data < -3*im['ERR'].data
            print '%s #negative = %d' %(exp, neg.sum())
            im['DQ'].data[neg] |= 4096
            im.flush()
        #
        drizzlepac.astrodrizzle.AstroDrizzle(asn_file, skysub=True, static=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_combine=True, clean=True, final_wcs=True, final_scale=0.1, final_rot=0, final_pixfrac=1, context=False, resetbits=0, final_bits=64, preserve=False)
        sci = pyfits.open(asn_file.replace('_asn', '_drz_sci'), mode='update')
        wht = pyfits.open(asn_file.replace('_asn', '_drz_wht'))
        sci[0].data[wht[0].data == 0] = 0
        sci.flush()
        
    #### ACS / UVIS
    files=glob.glob('*-F6*asn.fits')
    files=glob.glob('*-F[48]*asn.fits')
    files.sort()
    for asn_file in files[::-1]:
        asn = threedhst.utils.ASNFile(asn_file)
        for i in range(len(asn.exposures)):
            asn.exposures[i] = asn.exposures[i].split('_flc')[0]
    
        print asn_file, asn.exposures
        asn.write(asn_file)
        
        if ('814' in asn_file) | ('625' in asn_file):
            radec=None
        else:
            root = glob.glob('*-F[86]*W*sci.fits')[0].split('_sci')[0]
            if not os.path.exists('%s_sci.cat' %(root)):
                os.system('cp /user/brammer/software/share/gauss_5.0_9x9.conv .')
                research.hawkiff.make_catalog(root=root, sci_ext='sci', threshold=3, use_rms=False, subtract_background=True)
                cat = catIO.Table('%s_sci.cat' %(root))
                cat['X_WORLD', 'Y_WORLD'].write('ref_opt.radec', format='ascii.commented_header')
            
            radec='ref_opt.radec'
        
        radec = None
        init.prep_direct_grism_pair(direct_asn=asn_file, grism_asn=None, radec=radec, ACS=True, align_threshold=5, raw_path='../RAW/')
    #
    for asn_file in files:
        radec = None
        init.prep_direct_grism_pair(direct_asn=asn_file, grism_asn=None, radec=radec, ACS=True, align_threshold=5, raw_path='../RAW/')
    
    files=glob.glob('*sci.fits')
    ims = {}
    for file in files:
        im = pyfits.open(file)
        filter = file.split('-')[-1][:5]
        pix = np.sqrt(im[0].header['CD1_1']**2+im[0].header['CD1_2']**2)*3600
        scl = im[0].header['PHOTFLAM']/1.e-19*(0.1/pix)**2
        im[0].data *= scl
        ims[filter] = im
        
    ### Try Multidrizzle CR with combined image
    import drizzlepac
    
    asn_files=glob.glob('*-F[4-8]*asn.fits')
    exposures = []
    for file in asn_files:
        asn = threedhst.utils.ASNFile(file)
        exposures.extend(asn.exposures)
    
    target = '-'.join(file.split('-')[:-1])
    asn.exposures = exposures
    asn.product = target
    asn_file = '%s_asn.fits' %(asn.product)
    asn.write(asn_file, clobber=True)
    
    ### modify:  make copy of FLCs, make new mosaic and then copy FLCs back
    ##xxx 
    
    drizzlepac.astrodrizzle.AstroDrizzle(asn_file, skysub=False, clean=True, final_wcs=True, final_scale=0.05, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, preserve=False)
    
    #
    bits = 576
    
    drizzlepac.astrodrizzle.AstroDrizzle(asn_file, static=False, skysub=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_rot=0, final_scale=0.1, final_pixfrac=1, context=False, resetbits=0, final_bits=bits, preserve=False)
    
    #### rerun final mosaics
    refimage = '%s_drz_sci.fits' %(target)
    
    for file in asn_files:
        drizzlepac.astrodrizzle.AstroDrizzle(file, static=False, skysub=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_refimage=refimage, final_pixfrac=1, context=False, resetbits=0, final_bits=bits, preserve=False)
    
    for file in asn_files:
        print file
        sci = pyfits.open(glob.glob(file.replace('_asn', '_dr?_sci'))[0], mode='update')
        wht = pyfits.open(sci.filename().replace('_sci','_wht'))
        sci[0].data[wht[0].data == 0] = 0
        sci.flush()
        
    #### Try own driz CR
    import stwcs
    from drizzlepac import astrodrizzle, quickDeriv
    import scipy.ndimage as nd
    
    ref_file = 'ESO550-IG02-13-083-F435W_drc_sci.fits'
    ref_file = 'ESO550-IG02-13-083-F814W_drc_sci.fits'
    ref_file = 'ESO550-IG02-13-083_drc_sci.fits'
    ref_ext = 0
    ref = pyfits.open(ref_file)
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=ref_ext)
    
    asn_file = 'ESO550-IG02-13-083-F814W_asn.fits'
    asn = threedhst.utils.ASNFile(asn_file)
    for exp in asn.exposures:
        flt = pyfits.open('%s_flc.fits' %(exp)) #, mode='update')
        for ext in [1,2]:
            flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=('sci',ext))
            blotted_ref = astrodrizzle.ablot.do_blot(ref[ref_ext].data, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
            #
            blotDeriv = quickDeriv.qderiv(blotted_ref)
            
            scl = flt['sci', ext].header['PHOTFLAM']/ref[0].header['PHOTFLAM']/flt[0].header['EXPTIME']
            noise = flt['err', ext].data
            
            driz_scale = 1.2
            driz_snr = 3.5

            driz_scale = 1.0
            driz_snr = 6
            driz_f = 2
            
            # |data_image - blotted_image| > scale x deriv + SNR x noise
            abs = np.abs(flt['sci', ext].data*scl - blotted_ref*driz_f)
            xderiv = driz_scale*blotDeriv
            xsnr = driz_snr*scl*noise
            
            mask = abs > (xderiv + xsnr)
            new_dq = mask & ((flt['dq', ext].data & 4096) == 0)
            
    ### Try "cosmics" LA Cosmic
    im = pyfits.open('ESO550-IG02-13-083-F814W_drc_sci.fits')
    slx, sly = slice(1500, 2600), slice(1800, 2800)
    ext = 0
    subim = im[ext].data[sly, slx]
    h = im[ext].header
    subim *= im[0].header['EXPTIME']

    #### FLT
    im = pyfits.open('j9cv13pcq_flc.fits')
    ext = 0
    slx, sly = slice(1500, 2600), slice(0, 400)
    subim = im[4].data[sly,slx]
    h = im[0].header
        
    import cosmics
    c = cosmics.cosmicsimage(subim, pssl=20.0, gain=1, readnoise=h['READNSEA'], sigclip=4.0, sigfrac=0.3, objlim=3.8, satlevel=84700.0, verbose=True)
    c.run(maxiter=4)
    
    crflux = subim*c.mask
    plt.hist(np.log10(crflux[c.mask]), range=[1,6], bins=100, alpha=0.5)
    
    ####
    import research.pab.pab
    files=glob.glob('*-F*asn.fits')
    #files=glob.glob('IC*-F814*asn.fits')
    for file in files:
        asn = threedhst.utils.ASNFile(file)
        for exp in asn.exposures:
            flc = pyfits.open('%s_flc.fits' %(exp))
            research.pab.pab.run_lacosmic('%s_flc.fits' %(exp), split=2048, sigclip=5, pssl=flc[1].header['MDRIZSK0'])
        
        #
        drizzlepac.updatehdr.update_from_shiftfile(file.replace('asn.fits', 'shifts.txt'))
        drizzlepac.astrodrizzle.AstroDrizzle(file, skysub=True, clean=True, final_wcs=True, final_scale=0.05, final_pixfrac=0.8, context=False, resetbits=0, final_bits=576, preserve=False)

        drizzlepac.astrodrizzle.AstroDrizzle(file, output='sub_'+file.split('_asn')[0], skysub=True, static=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, clean=True, final_wcs=True, final_scale=0.08, final_ra=60.635233, final_dec=60.344597, final_outnx=800, final_outny=800, final_rot=0, final_pixfrac=0.8, context=False, resetbits=0, final_bits=576, preserve=False)
        
    ### Final mosaics
    f110w = glob.glob('*110W*asn.fits')[0]
    drizzlepac.astrodrizzle.AstroDrizzle(f110w, output='sub_'+f110w.split('_asn')[0], skysub=True, static=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, clean=True, final_wcs=True, final_scale=0.1, final_rot=0, final_pixfrac=1, context=False, resetbits=0, final_bits=576-512, preserve=False)
    
    files=glob.glob('*-F*asn.fits')
    for file in files:
        if 'F110W' in file:
            continue
        #
        drizzlepac.astrodrizzle.AstroDrizzle(file, output='sub_'+file.split('_asn')[0], skysub=True, static=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, clean=True, final_wcs=True, final_refimage='sub_'+f110w.replace('_asn', '_drz_sci'), final_pixfrac=1, context=False, resetbits=0, final_bits=576-512*('-F1' in file), preserve=False)
    
    for sci_im in glob.glob('sub*sci.fits'):
        print sci_im
        sci = pyfits.open(sci_im, mode='update')
        wht = pyfits.open(sci_im.replace('_sci','_wht'))
        sci[0].data[wht[0].data == 0] = 0
        sci.flush()
        
    #### Test alignment
    import align
    ra_list, de_list = np.loadtxt('f814w.radec', unpack=True)
    for geom in ['shift', 'shift', 'rxyscale', 'shift', 'shift']:
        align.get_align_to_subaru(sci='sub_ESO550-IG025-41-115-F673N_drc_sci.fits', wht='sub_ESO550-IG025-41-115-F673N_drc_wht.fits', verbose=False, fitgeometry=geom, align_data=(ra_list, de_list, ra_list*0), THRESH=1.2)
    
    #drizzlepac.tweakback.tweakback('ESO550-IG025-41-115-F673N_drc_sci.fits', force=True, verbose=True) #, origwcs='DRZWCS')

    ra_list, de_list = np.loadtxt('f110w.radec', unpack=True)
    for geom in ['shift', 'shift', 'rxyscale', 'shift', 'shift']:
        align.get_align_to_subaru(sci='sub_ESO550-IG025-41-115-F132N_drz_sci.fits', wht='sub_ESO550-IG025-41-115-F132N_drz_wht.fits', verbose=False, fitgeometry=geom, align_data=(ra_list, de_list, ra_list*0), THRESH=1.2)
    
    #drizzlepac.tweakback.tweakback('ESO550-IG025-41-115-F132N_drz_sci.fits', force=True, verbose=True, origwcs='DRZWCS')
    
    #align.get_align_to_subaru(sci='ESO550-IG025-41-115-F673N_drc_sci.fits', wht='ESO550-IG025-41-115-F673N_drc_wht.fits', verbose=True, fitgeometry='rxyscale', align_data=(ra_list, de_list, ra_list*0))

def line_ratios():
    """
    Testing for making the decrement maps
    """
    import glob
    files=glob.glob('*sci.fits')
    ims = {}
    for file in files:
        im = pyfits.open(file)
        if 'FILTER' in im[0].header.keys():
            filt = im[0].header['FILTER']
        else:
            filt = im[0].header['FILTER2']
        
        scl = im[0].header['PHOTFLAM']/1.e-19#*(5500./im[0].header['PHOTPLAM'])
        print '%s %s %.3f' %(file, filt, scl)
        ims[filt] = im[0].data*scl
            
def run_lacosmic(flc='j9cv13pcq_flc.fits', split=1024, pssl=0.0, sigclip=4.0, objlim=3.8, sigfrac=0.3):
    """
    
    drizzlepac.astrodrizzle.AstroDrizzle('ESO550-IG02-13-083-F814W_asn.fits', skysub=True, clean=True, final_wcs=True, final_scale=0.05, final_pixfrac=0.8, context=False, resetbits=0, final_bits=576, preserve=False)
    
    """
    
    import cosmics
    import threedhst
    
    im = pyfits.open(flc, mode='update')
    h = im[0].header
    
    nx, ny = 4096/split, 2048/split
    
    for ext in [1,2]:
        if flc.startswith('i'):
            dq = np.zeros((2051, 4096), dtype=int)
        else:
            dq = np.zeros((2048, 4096), dtype=int)
            
        for i in range(nx):
            for j in range(ny):
                threedhst.showMessage('ext: %d, (i,j)=%d,%d' %(ext, i,j))
                subim = im['sci',ext].data[j*split:(j+1)*split, i*split:(i+1)*split]
                c = cosmics.cosmicsimage(subim, pssl=pssl, gain=1, readnoise=h['READNSEA'], sigclip=objlim, sigfrac=sigfrac, objlim=objlim, satlevel=84700.0, verbose=True)
                c.run(maxiter=4)
                dq[j*split:(j+1)*split, i*split:(i+1)*split] += c.mask*1
        
        im['dq',ext].data |= dq*4096
        
        pyfits.writeto(im.filename().replace('.fits', '_lac_%d.fits' %(ext)), data=dq*4096, header=im['dq', ext].header, clobber=True)
        
    im.flush()
    
def min_drizzle(root='ESO550-IG02-13-083-F814W'):
    from drizzlepac.astrodrizzle import adrizzle
    
    asn = threedhst.utils.ASNFile('%s_asn.fits' %(root))
    
    ref = pyfits.open('%s_drc_sci.fits' %(root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)
    
    outsci = ref[0].data*0.
    outwht = ref[0].data*0.
    outcon = ref[0].data*0.
    min_sci = outsci+1e10
    
    for exp in asn.exposures:
        flt = pyfits.open('%s_flc.fits' %(exp)) #, mode='update')
        outsci*=0
        outwht*=0
        outcon*=0
        for ext in [1,2]:
            print 'Drizzle %s[sci,%d]' %(exp, ext)
            flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=('sci',ext))
            dq = flt['dq',ext].data*1
            for bit in [64,512]:
                dq -= dq & bit
            
            mask = (dq > 0) | (flt['sci',ext].data < -5*flt['err',ext].data)
            flt['sci', ext].data[mask] = 1.e30    
            adrizzle.do_driz(flt['sci',ext].data, flt_wcs, flt['err',ext].data, ref_wcs, outsci, outwht, outcon, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=1, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
        
        has_data = outsci != 0
        
        min_sci[has_data] = np.minimum(min_sci[has_data], outsci[has_data])
        ds9.view(min_sci)
    
    min_sci[min_sci > 1.e5] = 0
    ds9.view(min_sci)
    
def ngc5256():
    
    research.pab.align.coarse_align('NGC5256-39-138-F673N_drc_sci.fits', ref='NGC5256-14-170-F110W_drz_sci.fits')    
    research.pab.align.coarse_align('NGC5256-48-152-F814W_drc_sci.fits', ref='NGC5256-14-170-F110W_drz_sci.fits')    
    
    shifts = {'F110W':[0,0],
              'F132N':[0,0],
              'F130N':[0,0],
              'F673N':[-1.137,5.2527],
              'F435W':[8.51,1.35],
              'F814W':[8.51,1.35]}
    #
    files=glob.glob('*sci.fits')
    ims = {}
    for file in files:
        print file
        im = pyfits.open(file)
        filter = file.split('-')[-1][:5]
        pix = np.sqrt(im[0].header['CD1_1']**2+im[0].header['CD1_2']**2)*3600
        scl = im[0].header['PHOTFLAM']/1.e-19*(0.1/pix)**2
        im[0].data *= scl
        shifted = nd.shift(im[0].data, shifts[filter][::-1])
        ims[filter] = shifted
    
    pab = ims['F132N'] - ims['F130N']
    ha = ims['F673N'] - ims['F814W']
    
def mcg_02_01_051():
    import numpy as np
    import astropy.io.fits as pyfits
    
    import research.pab.align
    import stwcs
    
    # research.pab.align.coarse_align('MCG-02-01-051-46-243-F673N_drc_sci.fits', ref='MGC-02-01-051-01-240-F132N_drz_sci.fits')    
    # research.pab.align.coarse_align('ARP256NED01-02-060-F435W_drc_sci.fits', ref='MGC-02-01-051-01-240-F132N_drz_sci.fits')    
    # research.pab.align.coarse_align('ARP256NED01-02-060-F814W_drc_sci.fits', ref='MGC-02-01-051-01-240-F132N_drz_sci.fits')    
    
    pab = pyfits.open('MGC-02-01-051-01-240-F132N_drz_sci.fits')
    con = pyfits.open('MGC-02-01-051-01-240-F130N_drz_sci.fits')
    bb = pyfits.open('MGC-02-01-051-01-240-F110W_drz_sci.fits')
    pab_wht = pyfits.open('MGC-02-01-051-01-240-F132N_drz_wht.fits')
    con_wht = pyfits.open('MGC-02-01-051-01-240-F130N_drz_wht.fits')
    bb_wht = pyfits.open('MGC-02-01-051-01-240-F110W_drz_wht.fits')
    slx, sly = slice(896, 1413), slice(703, 1695)
    
    ######### Voronoi bins
    wcs = stwcs.wcsutil.HSTWCS(pab)
    slwcs = wcs.slice([sly, slx])
    pab_sub = (pab[0].data-con[0].data)[sly, slx]
    rms_sub = np.sqrt(1/pab_wht[0].data + 1/con_wht[0].data)[sly, slx]
    rms_sub[~np.isfinite(rms_sub)] = 1e10

    pab_sub = (pab[0].data)[sly, slx]
    rms_sub = np.sqrt(1/pab_wht[0].data)[sly, slx]
    
    rms_scl = 1
    # SN = 5
    # omask = np.abs(pab_sub/rms_sub) > SN
    # pab_sm = pab_sub*1
    # pab_sm[omask] = 0.
    # sigma, rnd_scl = 2, 1./np.sqrt(2)
    
    #sigma = 1
    # rnd_noise = np.random.normal(size=(100,100))
    # rnd_scl = np.std(nd.gaussian_filter(rnd_noise, sigma=sigma))
    
    # pab_sm = nd.gaussian_filter(pab_sm, sigma=sigma)
    # pab_sm[omask] = SN*100
    
    #rnd_scl = threedhst.utils.nmad((pab_sub/rms_sub)[~omask])
    rnd_scl = 1.
       
    #pyfits.writeto('pab_sci.fits', data=pab_sm, header=slwcs.to_header(), clobber=True)
    pyfits.writeto('pab_sci.fits', data=pab_sub, header=slwcs.to_header(), clobber=True)
    pyfits.writeto('pab_rms.fits', data=rms_sub*rms_scl, header=slwcs.to_header(), clobber=True)
    
    import research.pab
    research.pab.voronoi.go(cfgfile=None, fluximagefile='pab_sci.fits', fluxerrorfile='pab_rms.fits', distmpc=0.1, snthresh=SN, cvt_iters=10, max_z_kpc=100000.0, max_bin_size_pix=1000, output_prefix='pab_voronoi')

    #research.pab.voronoi.go(cfgfile='voronoi.cfg', output_prefix='pab_voronoi_cfg')
    
    v = pyfits.open('pab_voronoi.fits')
    mask = (np.isfinite(pab_sub) & (rms_sub < 10))
    omask = np.abs(pab_sub/rms_sub) > SN
    
    files=glob.glob('[AM]*sci.fits')
    for file in files:
        im = pyfits.open(file)
        im_wht = pyfits.open(file.replace('_sci','_wht'))
        wcs_im = stwcs.wcsutil.HSTWCS(im, ext=0)
        xy = wcs_im.all_world2pix([pab[0].header['CRVAL1']], [pab[0].header['CRVAL2']], 1)
        dx = pab[0].header['CRPIX1']-xy[0][0]
        dy = pab[0].header['CRPIX2']-xy[1][0]
        
        im_sub = nd.shift(im[0].data, (dy, dx))[sly, slx]
        wht_sub = nd.shift(im_wht[0].data, (dy, dx))[sly, slx]
        im_rms = 1/np.sqrt(wht_sub)
        im_rms[~np.isfinite(im_rms)] = 1e10
        
        binned_sci, binned_var, binned_area = research.pab.voronoi.rebin_output(im_sub, (im_rms)**2, v[0].data+1, mask=mask, ds9=None)
        binned_sci[omask] = im_sub[omask]
        binned_var[omask] = im_rms[omask]**2
        h = slwcs.to_header()
        filter=file.split('_')[0][-5:]
        h['FILTER'] = filter
        h['PHOTFLAM'] = im[0].header['PHOTFLAM']
        
        pyfits.writeto('sub_%s_sci.fits' %(filter), data=binned_sci, header=h, clobber=True)
        pyfits.writeto('sub_%s_var.fits' %(filter), data=binned_var, header=h, clobber=True)
        
        #ds9.view(binned_sci)
    
    files=glob.glob('sub*sci.fits')
    sub = {}
    for file in files:
        print file
        im = pyfits.open(file)
        sub[im[0].header['FILTER']] = im[0].data*im[0].header['PHOTFLAM']/4.e-19
        
    pab = sub['F132N'] - sub['F130N']
    ha = sub['F673N'] - (673.-435.)/(814.-435)*(sub['F814W']/sub['F435W'])
    ha = sub['F673N'] - sub['F814W']
    decrement = ha/pab/10.
    dec_mask = (pab > 0.05) & np.isfinite(decrement)   
    decrement[~dec_mask] = -10

def sdss_j110501():
    
    import research.pab.align
    
    research.pab.align.coarse_align(image='SDSS-J110501.98+594103.5-60-345-F775W_drc_sci.fits', ref='SDSS-J110501.98+594103.5-1A-074-F110W_drz_sci.fits')
    
    research.pab.align.coarse_align(image='SDSS-J110504.41+593957.3-65-066-F775W_drc_sci.fits', ref='SDSS-J110501.98+594103.5-1A-074-F110W_drz_sci.fits')
    
    for root in ['SDSS-J110501.98+594103.5-60-345', 'SDSS-J110504.41+593957.3-65-066']:
        im_ref = pyfits.open('%s-F775W_drc_sci.fits' %(root))
        asn_files=glob.glob('%s*asn.fits' %(root))
        for asn_file in asn_files:
            asn = threedhst.utils.ASNFile(asn_file)
            for exp in asn.exposures:
                print '%s -> %s' %(asn_file, exp)
                flc = pyfits.open('%s_flc.fits' %(exp), mode='update')
                
                mask = flc['DQ'].data == 0
                slx = slice(733, 948)
                sly = slice(775, 943)
                window = flc['SCI'].data[sly, slx]
                bg = np.median(window[mask[sly, slx]])
                flc['SCI'].data -= bg
                
                flc[1].header['CRVAL1'] += im_ref[0].header['DELTARA']
                flc[1].header['CRVAL2'] += im_ref[0].header['DELTADE']
                flc[1].header['DELTARA'] = im_ref[0].header['DELTARA']
                flc[1].header['DELTADE'] = im_ref[0].header['DELTADE']
                
                flc.flush()
                
    ### Drizzle IR images
    drizzlepac.astrodrizzle.AstroDrizzle('SDSS-J110501.98+594103.5-1A-074-F110W_asn.fits', output='SDSS-J110501-F110W', static=False, skysub=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_refimage=None, final_rot=0, final_scale=0.128, final_pixfrac=1, context=False, resetbits=0, final_bits=576, preserve=False)
    
    drizzlepac.astrodrizzle.AstroDrizzle('SDSS-J110501.98+594103.5-1A-074-F132N_asn.fits', output='SDSS-J110501-F132N', static=False, skysub=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_refimage='SDSS-J110501-F110W_drz_sci.fits', final_pixfrac=1, context=False, resetbits=0, final_bits=576, preserve=False)
    
    for filter in ['F438W', 'F673N', 'F775W']:
        exposures = []
        asn_files=glob.glob('SDSS*%s_asn.fits' %(filter))
        for asn_file in asn_files:
            asn = threedhst.utils.ASNFile(asn_file)
            exposures.extend(asn.exposures)
        
        flc_files = ['%s_flc.fits' %(exp) for exp in exposures]
        
        drizzlepac.astrodrizzle.AstroDrizzle(flc_files, output='SDSS-J110501-%s' %(filter), static=False, skysub=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_refimage='SDSS-J110501-F110W_drz_sci.fits', final_pixfrac=1, context=False, resetbits=0, final_bits=576, preserve=False)
        
    ######### Voronoi bins
    pab = pyfits.open('SDSS-J110501-F132N_drz_sci.fits')
    pab_wht = pyfits.open('SDSS-J110501-F132N_drz_wht.fits')
    con = pyfits.open('SDSS-J110501-F110W_drz_sci.fits')
    con_wht = pyfits.open('SDSS-J110501-F110W_drz_wht.fits')
    
    ### top
    slx, sly = slice(600, 953), slice(805, 1054)

    ### bottom
    slx, sly = slice(430, 814), slice(334, 508)
    
    slx, sly = slice(420, 953), slice(330, 1054)

    root='NGC6670-15-167'
    slx, sly = slice(605, 1369), slice(641, 1080)

    root='NGC7592-51-247'
    slx, sly = slice(702, 1222), slice(688, 1316)

    pab = pyfits.open('%s-F132N_drz_sci.fits' %(root))
    pab_wht = pyfits.open('%s-F132N_drz_wht.fits' %(root))
    con = pyfits.open('%s-F110W_drz_sci.fits' %(root))
    con_wht = pyfits.open('%s-F110W_drz_wht.fits' %(root))
    
    wcs = stwcs.wcsutil.HSTWCS(pab)
    slwcs = wcs.slice([sly, slx])
    # pab_sub = (pab[0].data-con[0].data)[sly, slx]
    # rms_sub = np.sqrt(1/pab_wht[0].data + 1/con_wht[0].data)[sly, slx]
    # rms_sub[~np.isfinite(rms_sub)] = 1e10

    pab_sub = (pab[0].data)[sly, slx]
    rms_sub = np.sqrt(1/pab_wht[0].data)[sly, slx]
    
    rms_scl = 1
    # SN = 5
    rnd_scl = 1.
       
    pyfits.writeto('pab_sci.fits', data=pab_sub, header=slwcs.to_header(), clobber=True)
    pyfits.writeto('pab_rms.fits', data=rms_sub*rms_scl, header=slwcs.to_header(), clobber=True)
    
    import research.pab
    research.pab.voronoi.go(cfgfile=None, fluximagefile='pab_sci.fits', fluxerrorfile='pab_rms.fits', distmpc=0.1, snthresh=SN, cvt_iters=10, max_z_kpc=100000.0, max_bin_size_pix=1000, output_prefix='pab_voronoi')

    #research.pab.voronoi.go(cfgfile='voronoi.cfg', output_prefix='pab_voronoi_cfg')
    
    v = pyfits.open('pab_voronoi.fits')
    mask = (np.isfinite(pab_sub) & (rms_sub < 10))
    omask = np.abs(pab_sub/rms_sub) > SN
    
    files=glob.glob('SDSS-J110501-F*sci.fits')
    files=glob.glob('NGC6670*-F*sci.fits')
    files=glob.glob('NGC7592*-F*sci.fits')
    
    for file in files:
        im = pyfits.open(file)
        im_wht = pyfits.open(file.replace('_sci','_wht'))
        wcs_im = stwcs.wcsutil.HSTWCS(im, ext=0)
        xy = wcs_im.all_world2pix([pab[0].header['CRVAL1']], [pab[0].header['CRVAL2']], 1)
        dx = pab[0].header['CRPIX1']-xy[0][0]
        dy = pab[0].header['CRPIX2']-xy[1][0]
        
        im_sub = nd.shift(im[0].data, (dy, dx))[sly, slx]
        wht_sub = nd.shift(im_wht[0].data, (dy, dx))[sly, slx]
        im_rms = 1/np.sqrt(wht_sub)
        im_rms[~np.isfinite(im_rms)] = 1e10
        
        binned_sci, binned_var, binned_area = research.pab.voronoi.rebin_output(im_sub, (im_rms)**2, v[0].data+1, mask=mask, ds9=None)
        binned_sci[omask] = im_sub[omask]
        binned_var[omask] = im_rms[omask]**2
        h = slwcs.to_header()
        filter=file.split('_')[0][-5:]
        h['FILTER'] = filter
        h['PHOTFLAM'] = im[0].header['PHOTFLAM']
        
        pyfits.writeto('sub_%s_sci.fits' %(filter), data=binned_sci, header=h, clobber=True)
        pyfits.writeto('sub_%s_var.fits' %(filter), data=binned_var, header=h, clobber=True)
        print filter
        
    #
    files=glob.glob('sub*sci.fits')
    sub = {}
    for file in files:
        print file
        im = pyfits.open(file)
        sub[im[0].header['FILTER']] = im[0].data*im[0].header['PHOTFLAM']/4.e-19
        
    pab = sub['F132N'] - sub['F110W']*0.85
    ha = sub['F673N'] - (775.-438.)/(775.-438)*(sub['F775W']/sub['F438W'])
    ha = sub['F673N'] - sub['F775W']
    decrement = ha/pab/10.
    dec_mask = (pab > 0.05) & np.isfinite(decrement)   
    decrement[~dec_mask] = -10
    
    
def m82():
    
    import glob
    import os
    
    import threedhst
    import threedhst.prep_flt_astrodrizzle as init
    from threedhst import catIO
    
    import unicorn
    import research.hawkiff
    
    import stsci.convolve
    import research.pab.align as align
    
    ### Make ACS associations    
    unicorn.candels.make_asn_files(uniquename=True, translate={'-ROT':''})
    threedhst.options['FLT_PERSISTENCE_PATH'] = '../Persistence/'
    
    #### IR
    files=glob.glob('*-F1*asn.fits')
    radec=None
    for asn_file in files:
        #init.prep_direct_grism_pair(direct_asn=asn_file, grism_asn=None, radec=radec, ACS=False, align_threshold=5, raw_path='../RAW/', order=0)
        threedhst.process_grism.fresh_flt_files(asn_file)
        asn = threedhst.utils.ASNFile(asn_file)
        for exp in asn.exposures:
            init.subtract_fixed_background(flt_file='%s_flt.fits' %(exp))
            
        #drizzlepac.astrodrizzle.AstroDrizzle(asn_file, static=False, skysub=False, driz_separate=True, driz_sep_wcs=False, median=True, blot=True, driz_cr=True, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_rot=0, final_scale=0.128, final_pixfrac=1, context=False, resetbits=0, final_bits=64, preserve=False)
    
    #
    ims = []
    hdrs = []
    
    for asn_file in files:
        asn = threedhst.utils.ASNFile(asn_file)
    
        ### cross correlation
        for exp in asn.exposures:
            flt = pyfits.open('%s_flt.fits' %(exp))
            slx, sly = slice(400,600), slice(400,600)
            #slx, sly = slice(450,550), slice(450,550)
            
            ims.append((flt[1].data*(flt[3].data == 0))[sly, slx])
            hdrs.append(flt[1].header)
    
    shifts = [[0,0]]
        
    for i in range(1,len(ims)):
            
        ccorr = stsci.convolve.correlate2d(ims[0], ims[i], fft=1)
        from astropy.modeling import models, fitting
        poly = models.Polynomial2D(degree=4)
        lsq = fitting.LevMarLSQFitter()
        yp, xp = np.indices(ccorr.shape)
        pfit = lsq(poly, xp, yp, ccorr)
        
        #ix = np.unravel_index(np.argmax(ccorr-pfit(xp, yp)), ccorr.shape)
        sh = ccorr.shape
        ix = [sh[0]/2, sh[1]/2]
        ccorr[ix[0], ix[1]] = np.nan
        
        
        ds9.view(ccorr-pfit(xp, yp))
        status = raw_input('Mark center of ccorr peak')
        xcf = np.cast[float](ds9.get('pan image').split())
        xc = np.cast[int](np.round(xcf))
        
        N = 10
        sly2, slx2 = slice(xc[1]-N-1,xc[1]+N-1), slice(xc[0]-N-1, xc[0]+N-1)
        sub = (ccorr-pfit(xp, yp))[sly2, slx2]
        
        flat = models.Polynomial2D(degree=1)
        gauss = models.Gaussian2D(100000, xcf[0], xcf[1], 1, 1, 0, bounds={'x_mean':[xcf[0]-5,xcf[0]+5], 'y_mean':[xcf[1]-5,xcf[1]+5]})
        ok = np.isfinite(sub)
        gfit = lsq(flat+gauss, xp[sly2, slx2][ok], yp[sly2, slx2][ok], sub[ok])
        shifts.append([ix[1]-gfit.x_mean_1, ix[0]-gfit.y_mean_1])
        print shifts[-1]
    

    wcs = pywcs.WCS(hdrs[0])
    crpix = [[hdrs[0]['CRPIX1'], hdrs[0]['CRPIX2']]]
    
    for i in range(1, len(ims)):
        wcs = pywcs.WCS(hdrs[i])
        xy = wcs.all_world2pix([hdrs[0]['CRVAL1']], [hdrs[0]['CRVAL2']], 1)
        crpix.append([xy[0][0], xy[1][0]])
        print crpix[-1]
    
    dx = np.array(shifts) - (np.array(crpix)-crpix[0])
    
    ### M82
    dx = np.array([[ 0.        ,  0.        ],
           [ 0.1692337 , -0.03445068],
           [ 0.11051269,  0.26955819],
           [ 0.04652894,  0.30199132]])
           
    i0 = 0
    for asn_file in files:
        asn = threedhst.utils.ASNFile(asn_file)
        fp = open(asn_file.replace('_asn.fits', '_shifts.txt'),'w')
        fp.write("""# frame: output
# refimage: %s_flt.fits[1]
# form: delta
# units: pixels
""" %(asn.exposures[0]))

        for exp in asn.exposures:
            fp.write('%s_flt.fits  %.4f  %.4f  0, 1.0 0.1 0.1\n' %(exp, dx[i0][0], dx[i0][1]))
            flt = pyfits.open('%s_flt.fits' %(exp), mode='update')
            flt[1].header['CRPIX1'] += dx[i0][0]
            flt[1].header['CRPIX2'] += dx[i0][1]
            flt.flush()
            
            i0+=1
        
        fp.close()
        
    for asn_file in files:
        drizzlepac.astrodrizzle.AstroDrizzle(asn_file, static=False, skysub=True, driz_separate=True, driz_sep_wcs=False, median=True, blot=True, driz_cr=True, driz_combine=True, final_wht_type='IVM', clean=True, final_wcs=True, final_rot=0, final_scale=0.1, final_pixfrac=1, context=False, resetbits=4096, final_bits=64, preserve=False)
    
    
    for asn_file in files:
        align.coarse_align(image=asn_file.replace('_asn.fits', '_drz_sci.fits'), ref='h_m82_v_s20_drz_sci.fits')
        
      
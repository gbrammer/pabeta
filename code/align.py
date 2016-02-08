#!/usr/bin/env python
"""
Align to Subaru and/or HST astrometric catalogs.

GBB

"""
import os
import shutil

import numpy as np
import pyfits
import pywcs
#import astropy.wcs as pywcs
#import astropy.io.fits as pyfits

import threedhst
from threedhst import catIO

import pyraf
from pyraf import iraf
    
def get_align_to_subaru(sci='M0416_Ks_c1_mp_avg.fits', wht='M0416_Ks_c1_mp_exp.fits', field='', clean=True, toler=3, verbose=False, fitgeometry='shift', shift_max=20, rms_max=1.1, rot_max=2, rot_only=True, THRESH=2, align_data=None):
    """
    Align HAWK-I images to the FF Subaru astrometric reference catalogs
    """
    
    #sci='M0416_Ks_c1_mp_avg.fits'; wht='M0416_Ks_c1_mp_exp.fits'
    
    ### Make object catalog
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'NONE'
    if wht is None:
        se.options['WEIGHT_TYPE']     = 'NONE'
    else:
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = wht
    
    se.options['FILTER']    = 'Y'
               
    se.options['DETECT_THRESH']    = '%d' %(THRESH)
    se.options['ANALYSIS_THRESH']  = '%d' %(THRESH)
    se.options['MAG_ZEROPOINT'] = '26.0'

    #### Run SExtractor on direct and alignment images
    ## direct image
    se.options['CATALOG_NAME']    = 'direct.cat'
    status = se.sextractImage(sci)
    threedhst.sex.sexcatRegions('direct.cat', 'direct.reg', format=2)
    
    directCat = threedhst.sex.mySexCat('direct.cat')
    
    #### Get the X/Y coords of the reference catalog    
    #head = pyfits.getheader(sci, 0)
    #wcs = pywcs.WCS(head)
    if 'M0416' in sci:
        ra_list, dec_list, mag = np.loadtxt(os.getenv('HAWKI')+'/FrontierFields/HST/hlsp_frontier_subaru_suprimecam_macs0416-astrom_R_v1_cat.txt', unpack=True)
        if ('c4' in sci):
            ra_list, dec_list, mag = np.loadtxt(os.getenv('HAWKI')+'/FrontierFields/HST/M0416/macs0416_f814w_radec.cat', unpack=True)
    #
    if 'M0717' in sci:
        ra_list, dec_list, mag = np.loadtxt('subaru.radec', unpack=True)

    if ('M1149' in sci) | (field == 'M1149'):
        ra_list, dec_list, mag = np.loadtxt('/Users/brammer/Research/VLT/HAWKI/MACS1149/hlsp_frontier_subaru_suprimecam_macs1149-astrom_R_v1_cat.txt', unpack=True)
            
    if 'A2744' in sci:
        ra_list, dec_list, mag = np.loadtxt(os.getenv('HAWKI')+'/FrontierFields/HST/hlsp_frontier_subaru_suprimecam_abell2744-astrom_i_v1_cat.txt', unpack=True)
        if ('c1' in sci) | ('c4' in sci):
            ra_list, dec_list, mag = np.loadtxt(os.getenv('HAWKI')+'/FrontierFields/HST/abell2744_f814w_radec.cat', unpack=True)
    
    if align_data is not None:
        ra_list, dec_list, mag = align_data
            
    im = pyfits.open(sci)
    print sci
    
    sh = im[0].shape
    head = im[0].header
    head['CUNIT1'] = 'deg'; head['CUNIT2'] = 'deg'
    wcs = pywcs.WCS(head)

    x_image, y_image = wcs.wcs_sky2pix(ra_list, dec_list, 1)
    
    try:
        x_image, y_image = wcs.wcs_sky2pix(ra_list, dec_list, 1)
    except:
        x_image, y_image = wcs.wcs_world2pix(ra_list, dec_list, 1)
    
    ok = (x_image > 0) & (y_image > 0) & (x_image < sh[1]) & (y_image < sh[1])

    x_image, y_image = x_image[ok], y_image[ok]
    
    fpr = open('align.reg','w')
    fpr.write('image\n')
    for i in range(ok.sum()): fpr.write('circle(%.6f, %.6f,0.3") # color=magenta\n' %(x_image[i], y_image[i]))
    fpr.close()
    
    # x_image, y_image = [], []
    # 
    # for ra, dec in zip(ra_list, dec_list):
    #     x, y = wcs.wcs_sky2pix([[ra, dec]], 1)[0]
    #     if (x > 0) & (y > 0) & (x < sh[1]) & (y < sh[1]):
    #         x_image.append(x)
    #         y_image.append(y)
    
    alignCat = catIO.EmptyCat()
    alignCat['X_IMAGE'] = np.array(x_image)
    alignCat['Y_IMAGE'] = np.array(y_image)
    
    xshift = 0
    yshift = 0
    rot = 0
    scale = 1.
    
    xrms = 2
    yrms = 2
    
    NITER = 5
    IT = 0
    while (IT < NITER):
        IT = IT+1
        
        #### Get x,y coordinates of detected objects
        ## direct image
        fp = open('direct.xy','w')
        for i in range(len(directCat.X_IMAGE)):
            fp.write('%s  %s\n' %(directCat.X_IMAGE[i],directCat.Y_IMAGE[i]))
        fp.close()

        ## alignment image
        fp = open('align.xy','w')
        for i in range(len(alignCat.X_IMAGE)):
            fp.write('%s  %s\n' %(np.float(alignCat.X_IMAGE[i])+xshift,
                       np.float(alignCat.Y_IMAGE[i])+yshift))
        fp.close()

        iraf.flpr()
        iraf.flpr()
        iraf.flpr()
        #### iraf.xyxymatch to find matches between the two catalogs
        pow = toler*1.
        try:
            os.remove('align.match')
        except:
            pass
        status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                       output="align.match",
                       tolerance=2**pow, separation=0, verbose=iraf.yes, Stdout=1)
        
        nmatch = 0
        while status1[-1].startswith('0') | (nmatch < 10) | (float(status1[-3].split()[1]) > 40):
            pow+=1
            os.remove('align.match')
            status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                           output="align.match",
                           tolerance=2**pow, separation=0, verbose=iraf.yes, Stdout=1)
            #
            nmatch = 0
            for line in open('align.match').xreadlines(  ): nmatch += 1
            
        if verbose:
            for line in status1:
                print line
        
                
        #### Compute shifts with iraf.geomap
        iraf.flpr()
        iraf.flpr()
        iraf.flpr()
        try:
            os.remove("align.map")
        except:
            pass
            
        status2 = iraf.geomap(input="align.match", database="align.map",
                    fitgeometry=fitgeometry, interactive=iraf.no, 
                    xmin=iraf.INDEF, xmax=iraf.INDEF, ymin=iraf.INDEF, ymax=iraf.INDEF,
                    maxiter = 10, reject = 2.0, Stdout=1)
        if verbose:
            for line in status2:
                print line
        
        #fp = open(root+'.iraf.log','a')
        #fp.writelines(status1)
        #fp.writelines(status2)
        #fp.close()
                
        #### Parse geomap.output 
        fp = open("align.map","r")
        for line in fp.readlines():
            spl = line.split()
            if spl[0].startswith('xshift'):
                xshift += float(spl[1])    
            if spl[0].startswith('yshift'):
                yshift += float(spl[1])    
            if spl[0].startswith('xrotation'):
                rot = float(spl[1])    
            if spl[0].startswith('xmag'):
                scale = float(spl[1])    
            if spl[0].startswith('xrms'):
                xrms = float(spl[1])    
            if spl[0].startswith('yrms'):
                yrms = float(spl[1])    
            
        fp.close()
        
        #os.system('wc align.match')
        print 'Shift iteration #%d, xshift=%f, yshift=%f, rot=%f, scl=%f (rms: %5.2f,%5.2f)' %(IT, xshift, yshift, rot, scale, xrms, yrms)
    
    os.system('cat align.match | grep -v "\#" | grep [0-9] | awk \'{print "circle(", $1, ",", $2, ",4) # color=green"}\' > d.reg')
    os.system('cat align.match | grep -v "\#" | grep [0-9] | awk \'{print "circle(", $3, ",", $4, ",4) # color=magenta"}\' > a.reg')
    
    shutil.copy('align.map', sci.replace('.fits', '.align.map'))
    shutil.copy('align.match', sci.replace('.fits', '.align.match'))
    
    #### Cleanup
    if clean:
        rmfiles = ['align.cat', 'align.map','align.match','align.reg','align.xy', 'direct.cat','direct.reg','direct.xy']
        
        for file in rmfiles:
            try:
                os.remove(file)
            except:
                pass
    
    fp = open(sci.replace('.fits', '.align.info'), 'w')
    fp.write('# image xshift yshift rot scale xrms yrms\n')
    fp.write('%s %.3f %.3f %.4f %.4f %.3f %.3f\n' %(sci, xshift, yshift, rot, scale, xrms, yrms))
    
    if (np.abs(xshift) > shift_max) | (np.abs(yshift) > shift_max) | (xrms > rms_max) | (yrms > rms_max):
        print 'Shifts out of allowed range.  Run again with increased shift_max to accept.'
        #return xshift, yshift, rot, scale, xrms, yrms
        ## Add a small shift that should come out easily with another 
        ## shift iteration
        xshift, yshift, rot, scale, xrms, yrms = 2,2,0,1.0,-99,-99
        
    for file in [sci, wht]:
        if ('r' in fitgeometry) & rot_only:
            xshift, yshift = 0, 0
            
        #apply_offsets(file, [[xshift, yshift, rot, scale]])
        from drizzlepac import updatehdr
        updatehdr.updatewcs_with_shift(file, sci, wcsname='DRZWCS',
                        rot=rot,scale=scale,
                        xsh=xshift, ysh=yshift,
                        fit=None,
                        xrms=xrms, yrms = yrms,
                        verbose=False, force=True, sciext=0)
        
    if '_dr' in sci:
        im = pyfits.open(sci)
        h = im[0].header
        for i in range(h['NDRIZIM']):
            flt_str = h['D%03dDATA' %(i+1)]
            if 'sci,2' in flt_str:
                continue
            #
            flt_im = flt_str.split('[')[0]
            ext = int(flt_str.split('[')[1][:-1].split(',')[1])
            updatehdr.updatewcs_with_shift(flt_im, sci, wcsname='GTWEAK', rot=rot, scale=scale, xsh=xshift, ysh=yshift,
                            fit=None, xrms=xrms, yrms = yrms, verbose=False, force=True, sciext='SCI')
                
        # im = pyfits.open(file, mode='update')
        # wcs = pywcs.WCS(im[0].header)
        # wcs.rotateCD(-rot)
        # wcs.wcs.cd /= scale
        # #
        # im[0].header['CRPIX1'] += xshift
        # im[0].header['CRPIX2'] += yshift
        # #
        # for i in [0,1]:
        #     for j in [0,1]:
        #         im[0].header['CD%d_%d' %(i+1, j+1)] = wcs.wcs.cd[i,j]
        # #        
        # im.flush()
    
    return xshift, yshift, rot, scale, xrms, yrms

def apply_offsets(file, params):
    import astropy.wcs as pywcs
    import astropy.io.fits as pyfits
    
    im = pyfits.open(file, mode='update')
    for p in params:
        xshift, yshift, rot, scale = p
        #print p
        wcs = pywcs.WCS(im[0].header)
        wcs.rotateCD(-rot)
        wcs.wcs.cd /= scale
        im[0].header['CRPIX1'] += xshift
        im[0].header['CRPIX2'] += yshift
        #
        for i in [0,1]:
            for j in [0,1]:
                im[0].header['CD%d_%d' %(i+1, j+1)] = wcs.wcs.cd[i,j]
    #        
    im.flush()
    
def go_coarse():
    
    align.coarse_align(image='M1149_test_sci.fits', ref='/Users/brammer/Research/HST/FrontierFields/Subaru/hlsp_clash_subaru_suprimecam_macs1149_rc_2010-v20120820_sw.fits')
    
    files = glob.glob('20150224*sci.fits')
    files = glob.glob('*c4*sci.fits')
    for file in files:
        params = []
        for iter in range(4):
            out = align.get_align_to_subaru(sci=file, wht=file.replace('sci', 'wht'), clean=True, toler=4, verbose=False, fitgeometry='rxyscale', shift_max=40, rot_max=2, rot_only=False, THRESH=3, field='M1149')
            params.append(out[:-2])
        
        np.savetxt(file+'.align', params, fmt='%7.4f')
            
def coarse_align(image='M0416_Ks/M0416_Ks_c1_mp_avg.fits', ref='../20131103_60.01/M0416_Ks/M0416_Ks_c1_mp_avg.fits'):
    """
    
    Do rough alignment of two images with DS9.  The script will display both images 
    and prompt you to center the frame on a common object (just click to pan the frames).  
    
    """
    import threedhst.dq

    ds9 = threedhst.dq.myDS9()

    im = pyfits.open(image)
    im_ref = pyfits.open(ref)
    
    ds9.frame(1)
    ds9.set('file %s' %(image))
    ds9.scale(-0.1,10)
    ds9.set('scale log')

    ds9.frame(2)
    ds9.set('file %s' %(ref))
    ds9.scale(-0.1,10)
    ds9.set('scale log')
    
    ds9.set('tile yes')
    ds9.set('mode pan')
    
    x = raw_input('Center both frames on an object and type <ENTER>.')

    ds9.frame(1)
    im_rd = np.cast[float](ds9.get('pan fk5').split())
    im_xy = np.cast[float](ds9.get('pan image').split())
    
    ds9.frame(2)
    ref_rd = np.cast[float](ds9.get('pan fk5').split())
    ref_xy = np.cast[float](ds9.get('pan image').split())
    
    delta = ref_rd-im_rd
    im = pyfits.open(image, mode='update')
    im[0].header['CRVAL1'] += delta[0]
    im[0].header['CRVAL2'] += delta[1]
    im[0].header['DELTARA'] = delta[0]
    im[0].header['DELTADE'] = delta[1]
    
    print 'Delta (pix)', im_xy-ref_xy
    
    print 'Delta: %.4f, %.4f "' %(delta[0]*3600, delta[1]*3600)
    
    im.flush()
    
def check_alignment(match='A2744_Ks_c1_mp_avg.align.match'):
    x_ref, y_ref, x_in, y_in, i_ref, i_in = np.loadtxt(match, unpack=True)
    
    dx = x_in - x_ref
    dy = y_in - y_ref
    
    pix_scale = 0.1
    pix_scale = 0.06
    scl = 0.5
    
    plt.scatter(dx, dy, alpha=0.1)
    plt.xlim(-scl/pix_scale, scl/pix_scale); plt.ylim(-scl/pix_scale, scl/pix_scale)
    plt.xlabel(r'$\Delta x$ [%.02f" pix]' %(pix_scale)); plt.ylabel(r'$\Delta y$ [%.02f" pix]' %(pix_scale)) 
    
    dr = np.sqrt(threedhst.utils.biweight(dx)**2+threedhst.utils.biweight(dy)**2)
    ok = (np.abs(dx) < 2*dr) & (np.abs(dy) < 2*dr)
    s = 200
    
    plt.quiver(x_in[ok], y_in[ok], dx[ok]*s, dy[ok]*s, scale=1, units='xy', alpha=0.5)
    plt.quiver(0.05*x_in.max(), 0.05*y_in.max(), s, 0, scale=1, units='xy', alpha=0.8, label='1 pix', color='red')
    #plt.legend()
    
    
if __name__ == '__main__':
    
    import sys
    
    cwd = os.getcwd()
    if 'A2744' in cwd:
        field = 'A2744'
    
    if 'M0416' in cwd:
        field = 'M0416'
    
    #print sys.argv
    
    chips = [1,2,3,4]
    if len(sys.argv) > 1:
        chips = []
        for c in sys.argv[1:]:
            chips.append(int(c))
    
    #print chips
    
    #### Run in three steps to refine the shifts: shift, rxyscale, shift
    #### The first gets you close, the second takes out any rotation, and the last
    #### is needed because I don't know the shift reference of "rxyscale"....
    for chip in chips:
        get_align_to_subaru(sci='%s_Ks_c%d_mp_avg.fits' %(field, chip), wht='%s_Ks_c%d_mp_exp.fits' %(field, chip), clean=False, toler=3, verbose=False, fitgeometry='shift')
        get_align_to_subaru(sci='%s_Ks_c%d_mp_avg.fits' %(field, chip), wht='%s_Ks_c%d_mp_exp.fits' %(field, chip), clean=False, toler=3, verbose=False, fitgeometry='rxyscale', rot_only=True)
        xshift, yshift, rot, scale, xrms, yrms = get_align_to_subaru(sci='%s_Ks_c%d_mp_avg.fits' %(field, chip), wht='%s_Ks_c%d_mp_exp.fits' %(field, chip), clean=False, toler=3, verbose=False, fitgeometry='shift')
        if xrms < -90:
            xshift, yshift, rot, scale, xrms, yrms = get_align_to_subaru(sci='%s_Ks_c%d_mp_avg.fits' %(field, chip), wht='%s_Ks_c%d_mp_exp.fits' %(field, chip), clean=False, toler=3, verbose=False, fitgeometry='shift')
           
    # get_align_to_subaru(sci='M0416_Ks_c2_mp_avg.fits', wht='M0416_Ks_c2_mp_exp.fits', clean=False, toler=3, verbose=False, fitgeometry='shift')
    # get_align_to_subaru(sci='M0416_Ks_c3_mp_avg.fits', wht='M0416_Ks_c3_mp_exp.fits', clean=False, toler=3, verbose=False, fitgeometry='shift')
    # get_align_to_subaru(sci='M0416_Ks_c4_mp_avg.fits', wht='M0416_Ks_c4_mp_exp.fits', clean=False, toler=3, verbose=False, fitgeometry='shift')
    
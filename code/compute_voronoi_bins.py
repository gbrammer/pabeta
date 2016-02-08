#!/usr/bin/env python
"""
Taken from https://github.com/AndrewRook/astro/tree/master/voronoi
"""
import sys
import numpy as np
import ConfigParser
import time
import os
import astropy.io.fits as pyfits
from sklearn.neighbors import BallTree
import scipy.ndimage
#import bottleneck as bn

def writefits(image,filename,header=None):
    hdu = pyfits.PrimaryHDU(data=image,header=header)
    if os.path.isfile(filename):
        os.remove(filename)
    hdu.writeto(filename)
    
def readcfg(cfgfile):
    config=ConfigParser.ConfigParser()
    config.read(cfgfile)
    cfgdict = {}
    
    cfgdict['fluximagefile'] = process_cfg_param(config,'data','fluximagefile',str,False)
    cfgdict['fluxerrorfile'] = process_cfg_param(config,'data','fluxerrorfile',str,False)
    cfgdict['distmpc'] = process_cfg_param(config,'data','distmpc',float,False)

    cfgdict['cvt_iters'] = process_cfg_param(config,'constraints','cvt_iters',int,False)
    cfgdict['snthresh'] = process_cfg_param(config,'constraints','snthresh',float,False)
    cfgdict['max_z_kpc'] = process_cfg_param(config,'constraints','max_z_kpc',float,False)
    cfgdict['max_bin_size_pix'] = process_cfg_param(config,'constraints','max_bin_size_pix',int,False)

    cfgdict['output_prefix'] = process_cfg_param(config,'output','output_prefix',str,False)

    return cfgdict
    
def rebin_output(sci, var, bin_image, mask=None, ds9=None):
    """
    Given the bin image compute the average sci and variance within the bins
    """
    if mask is None:
        mask = np.isfinite(sci)
    
    #
    binned_sci = bin_image*0.
    binned_var = bin_image*0.
    binned_area = bin_image*0.

    bins = np.unique(bin_image[bin_image > 0])
    for i, b in enumerate(bins):
        bin_mask = (bin_image == b) & mask
        nbin = bin_mask.sum()
        binned_sci[bin_mask] = sci[bin_mask].sum()/nbin
        binned_var[bin_mask] = var[bin_mask].sum()/nbin**2
        binned_area[bin_mask] = nbin
        print '%d: %d, %d' %(i, b, nbin.sum())
        if ds9 is not None:
            ds9.view(binned_sci)
            
    return binned_sci, binned_var, binned_area

#Processes an individual configuration parameter:
def process_cfg_param(config,section,parametername,vartype,canbenull=False):
    if canbenull:
        paramval = config.get(section,parametername)
        if paramval.lower() == 'none' and vartype != list:
            return None
    paramval = None
    if vartype == int:
        paramval = config.getint(section,parametername)
    elif vartype == float:
        paramval = config.getfloat(section,parametername)
    elif vartype == bool:
        paramval = config.getboolean(section,parametername)
    else:
        paramval = config.get(section,parametername)
        if vartype == list:
            paramval = paramval.split(",")
            for i in range(len(paramval)):
                tempparamval = paramval[i].lower()
                if tempparamval == 'true':
                    paramval[i] = True
                elif tempparamval == 'false':
                    paramval[i] = False
                if canbenull:
                    if tempparamval == 'none':
                        paramval[i] = None
    return paramval

def compute_rz(cfgdict,fluxhdu):
    #Get the center of the galaxy in pixel space:
    #galcentx,galcenty = fluxhdu.header['GALCENTX'],fluxhdu.header['GALCENTY']
    galcentx = fluxhdu.header['NAXIS1']/2
    galcenty = fluxhdu.header['NAXIS2']/2
    # pixscalex = np.sqrt(fluxhdu.header['CD1_1']**2+fluxhdu.header['CD1_2']**2)*3600. #Arcseconds
    # pixscaley = np.sqrt(fluxhdu.header['CD1_2']**2+fluxhdu.header['CD2_2']**2)*3600. #Arcseconds
    pixscalex = np.sqrt(fluxhdu.header['PC1_1']**2+fluxhdu.header['PC1_2']**2)*3600. #Arcseconds
    pixscaley = np.sqrt(fluxhdu.header['PC1_2']**2+fluxhdu.header['PC2_2']**2)*3600. #Arcseconds

    #Get the indices of each pixel as 2d arrays:
    i_x = np.arange(fluxhdu.data.shape[1])
    i_y = np.arange(fluxhdu.data.shape[0])
    i_x2d,i_y2d = np.meshgrid(i_x,i_y)
    #Convert to arcseconds around the galaxy center:
    asec_x2d = (i_x2d-galcentx)*pixscalex
    asec_y2d = (i_y2d-galcenty)*pixscaley
    #Convert to kpc in radius and height:
    kpc_conversion = cfgdict['distmpc']*1000./206265. #kpc per arcsecond
    kpc_r2d = asec_x2d*kpc_conversion
    kpc_z2d = asec_y2d*kpc_conversion

    return i_x2d,i_y2d,kpc_r2d,kpc_z2d
    
def bin_accrete(rkpc,zkpc,snthresh,flux,fluxerror,snapprox_factor=4,dist_factor=1.2,bin_start=0):
    sn = flux/fluxerror
    squarefluxerror = fluxerror**2
    #Make the ball tree:
    combrzarr = (np.vstack((rkpc,zkpc))).T
    start = time.time()
    tree = BallTree(combrzarr,leaf_size=3,metric='euclidean')
    print "Making BallTree: {0:.2f}".format(time.time()-start)
    #Compute the minimum distance between two points:
    start = time.time()
    dists,indices = tree.query(combrzarr,k=2)
    mindist = np.min(dists[:,1])
    print "Computing Minimum Distance: {0:.2f} s".format(time.time()-start)

    bin_id = np.zeros(len(rkpc),dtype=np.int)-1
    reached_sn_thresh = np.ones(len(rkpc),dtype=np.bool)
    binned_bool = np.zeros(len(rkpc),dtype=np.bool)
    good_sns = sn.copy()
    
    currpixel = np.argmax(good_sns)
    maxbins = len(rkpc)
    start = time.time()
    for i in range(maxbins):
        #Get enough nearest-neighbors to cover enough to reach the S/N threshold:
        approx_numpix = int(round(snapprox_factor*(snthresh/sn[currpixel])**2))
        if approx_numpix > len(rkpc):
            approx_numpix = len(rkpc)
        #If you have a pixel above the threshold, just add it and move on:
        if approx_numpix <= 0:
            bin_id[currpixel] = i+bin_start
            binned_bool[currpixel] = True
            good_sns[currpixel] = -100
            currpixel = np.argmax(good_sns)
            if good_sns[currpixel] < 0:
                break#you're done, all pixels are binned
        else:
            dists,indices = tree.query(combrzarr[currpixel,:],k=approx_numpix)
            dists = dists[0,:]
            indices = indices[0,:]
        
            #Screen out the ones that have been binned and/or are not close enough:
            not_binned = (binned_bool[indices] == False)
            not_too_far = np.ones(len(dists),dtype=np.bool)
            if len(dists) > 1:
                not_too_far[1:] = ((dists[1:]-dists[:-1]) < dist_factor*mindist)
            dists_good = dists[(not_binned & not_too_far)]
            indices_good = indices[(not_binned & not_too_far)]
            
            #Get the cumulative S/N:
            cum_sn = np.cumsum(flux[indices_good])/np.sqrt(np.cumsum(squarefluxerror[indices_good]))
            #Shift cum_sn by 1 to make sure to capture the point where the S/N crosses the threshold too:
            max_cum_sn = cum_sn[-1]
            cum_sn[1:] = cum_sn[:-1]
            cum_sn[0] = 0
            #Make a cutoff where cumulative S/N = threshold:
            pixels_to_bin = (cum_sn <= snthresh)
            if i%1000 == 0:
                print i,approx_numpix,np.sum(pixels_to_bin),np.sum(binned_bool),"{0:.2f}".format(time.time()-start)
            #Bin up all those pixels:
            bin_id[indices_good[pixels_to_bin]] = i+bin_start
            binned_bool[indices_good[pixels_to_bin]] = True
            good_sns[indices_good[pixels_to_bin]] = -100
            if max_cum_sn < snthresh:
                reached_sn_thresh[indices_good[pixels_to_bin]] = False
            #Choose new pixel:
            currpixel = np.argmax(good_sns)
            # if np.sum((pixels_to_bin) == False) > 0:
            #     #If there is at least one pixel that was chosen that is past the S/N threshold, use it to start the next bin:
            #     currpixel = indices_good[(pixels_to_bin==False)][0]
            # else:
            #     currpixel = np.argmax(good_sns)
            if good_sns[currpixel] < 0:
                break; #you're done, all pixels are binned
                
    bin_id[reached_sn_thresh == False] = -1
    return bin_id

def fix_bins(rkpc,zkpc,bin_ids):
    #Get the centroids of all the successful bins:
    good_bins = np.unique(bin_ids[bin_ids >= 0])
    bad_rkpc = rkpc[bin_ids < 0]
    bad_zkpc = zkpc[bin_ids < 0]
    comb_badrzarr = (np.vstack((bad_rkpc,bad_zkpc))).T

    rcentroid = scipy.ndimage.mean(rkpc,labels=bin_ids,index=good_bins)
    zcentroid = scipy.ndimage.mean(zkpc,labels=bin_ids,index=good_bins)

    #Ball-tree the centroids:
    combrzarr = (np.vstack((rcentroid,zcentroid))).T
    tree = BallTree(combrzarr,leaf_size=5,metric='euclidean')

    #Query all the bad pixels for their nearest centroid:
    dists,indices = tree.query(comb_badrzarr,k=1)
    indices = indices[:,0]
    bin_ids[bin_ids < 0] = good_bins[indices]
    return bin_ids

def cvt_bins(rkpc,zkpc,flux,fluxerror,bin_ids,niters = 4,bad_bin_id = None):

    combrzarr = (np.vstack((rkpc,zkpc))).T
    weights = (flux/fluxerror)**2
    weighted_rkpc = rkpc*weights
    weighted_zkpc = zkpc*weights
    old_rcentroids = None
    old_zcentroids = None
    unique_bins = np.unique(bin_ids)
    if bad_bin_id != None:
        unique_bins = unique_bins[unique_bins != bad_bin_id]
    max_distance = -1
    for i in range(niters):
        start = time.time()
        #Compute centroids:
        meanweights = scipy.ndimage.mean(weights,labels=bin_ids,index=unique_bins)
        numperbin = scipy.ndimage.sum(np.ones(len(weights)),labels=bin_ids,index=unique_bins)
        mean_weighted_rkpc = scipy.ndimage.mean(weighted_rkpc,labels=bin_ids,index=unique_bins)
        mean_weighted_zkpc = scipy.ndimage.mean(weighted_zkpc,labels=bin_ids,index=unique_bins)
        #print numperbin.min(),numperbin.max()
        #np.savetxt('test_weights.txt',zip(meanweights,numperbin))

        #print "Debug: ",meanweights.min(),meanweights.max()

        rcentroids = (mean_weighted_rkpc/meanweights)[numperbin > 0]
        zcentroids = (mean_weighted_zkpc/meanweights)[numperbin > 0]

        #Compute which bin each pixel belongs in:
        combrzarr_centroids = (np.vstack((rcentroids,zcentroids))).T
        tree = BallTree(combrzarr_centroids,leaf_size=10,metric='euclidean')
        dists,indices = tree.query(combrzarr,k=1)
        indices = indices[:,0]
        bin_ids = unique_bins[indices]
        print "CVT iteration {0:d}: Number of bins = {1:d}, Time = {2:.2f}s".format(i,len(rcentroids),time.time()-start)

        # if old_rcentroids != None:
        #     old_rcentroids = old_rcentroids[numperbin > 0]
        #     old_zcentroids = old_zcentroids[numperbin > 0]
        #     distances = np.sqrt((rcentroids-old_rcentroids)**2+(zcentroids-old_zcentroids)**2)
        #     max_distance = distances.max()
        #     print "CVT iteration {0:d}: Number of bins = {1:d}, Max Centroid Change = {2:.3e} kpc, Time = {3:.2f}s".format(i,len(rcentroids),max_distance,time.time()-start)
        # old_rcentroids = rcentroids
        # old_zcentroids = zcentroids

    return bin_ids,max_distance

#Get rid of bins with too little flux:
def remove_bad_bins(bin_ids,rkpc,zkpc,flux,fluxerror,snthresh):
    #Compute the total S/N of each bin:
    unique_bins = np.unique(bin_ids)
    bin_fluxes = scipy.ndimage.mean(flux,labels=bin_ids,index=unique_bins)
    num_per_bin = scipy.ndimage.sum(np.ones(bin_ids.shape),labels=bin_ids,index=unique_bins)
    bin_error = np.sqrt(scipy.ndimage.sum(fluxerror**2,labels=bin_ids,index=unique_bins))/num_per_bin
    bin_sn = bin_fluxes/bin_error

    #Re-assign too small bins to nearest good bin:
    bad_bins_bool = np.in1d(bin_ids,unique_bins[(bin_sn < snthresh)])
    bad_bin_id = unique_bins.min()-1
    bin_ids[bad_bins_bool] = bad_bin_id
    new_bin_ids,junk = cvt_bins(rkpc,zkpc,flux,fluxerror,bin_ids,niters = 1,bad_bin_id = bad_bin_id)

    return new_bin_ids

def remove_large_bins(bin_ids,max_bin_size_pix,bad_bin_id = -1):
    unique_bins = np.unique(bin_ids)
    num_per_bin = scipy.ndimage.sum(np.ones(bin_ids.shape),labels=bin_ids,index=unique_bins)
    bad_bins_bool = np.in1d(bin_ids,unique_bins[(num_per_bin > max_bin_size_pix)])
    bin_ids[bad_bins_bool] = bad_bin_id
    return bin_ids
    
def go(cfgfile=None, fluximagefile='sci.fits', fluxerrorfile='rms.fits', distmpc=1, snthresh=5, cvt_iters=20, max_z_kpc=1e5, max_bin_size_pix=500, output_prefix='voronoibins'):
    
    if cfgfile is not None:
        cfgdict = readcfg(cfgfile)
    else:
        cfgdict = {}

        cfgdict['fluximagefile'] = fluximagefile
        cfgdict['fluxerrorfile'] = fluxerrorfile
        cfgdict['distmpc'] = distmpc

        cfgdict['cvt_iters'] = cvt_iters
        cfgdict['snthresh'] = snthresh
        cfgdict['max_z_kpc'] = max_z_kpc
        cfgdict['max_bin_size_pix'] = max_bin_size_pix

        cfgdict['output_prefix'] = output_prefix
        
    #load in the files:
    fluxhdu = pyfits.open(cfgdict['fluximagefile'])[0]
    fluxerrorhdu = pyfits.open(cfgdict['fluxerrorfile'])[0]

    if fluxhdu.data.shape != fluxerrorhdu.data.shape:
        raise ValueError("Images must have the same shape!")

    
    #Make 2d arrays showing index values of each pixel as well as r and z values in kpc:
    i_x2d,i_y2d,kpc_r2d,kpc_z2d = compute_rz(cfgdict,fluxhdu)

    #Cut down arrays based on z and regions outside of the frame, simultaneously making them 1d:
    goodzbool = (np.abs(kpc_z2d) <= cfgdict['max_z_kpc']) & (fluxerrorhdu.data > 0)
    i_x_zcut = i_x2d[goodzbool]
    i_y_zcut = i_y2d[goodzbool]
    kpc_r_zcut = kpc_r2d[goodzbool]
    kpc_z_zcut = kpc_z2d[goodzbool]
    flux_zcut = fluxhdu.data[goodzbool].astype(np.float32)
    fluxerror_zcut = fluxerrorhdu.data[goodzbool].astype(np.float32)

    #Do a couple of basic tests suggested by Cappellari:
    if np.sum(flux_zcut)/np.sqrt(np.sum(fluxerror_zcut**2)) < cfgdict['snthresh']:
        raise ValueError("Not enough S/N in entire ROI!")
    if np.min(flux_zcut/fluxerror_zcut) >= cfgdict['snthresh']:
        raise ValueError("All pixels in ROI are above S/N threshold!")

    #Remove points that are already above threshold from bin centering computation:
    indivgoodpoints_bool = (flux_zcut/fluxerror_zcut >= cfgdict['snthresh']) 
    i_x_zcut_indivgoodpoints = i_x_zcut[indivgoodpoints_bool]
    i_y_zcut_indivgoodpoints = i_y_zcut[indivgoodpoints_bool]
    kpc_r_zcut_indivgoodpoints = kpc_r_zcut[indivgoodpoints_bool]
    kpc_z_zcut_indivgoodpoints = kpc_z_zcut[indivgoodpoints_bool]
    flux_zcut_indivgoodpoints = flux_zcut[indivgoodpoints_bool]
    fluxerror_zcut_indivgoodpoints = fluxerror_zcut[indivgoodpoints_bool]
    num_indivgoodpoints = np.sum(indivgoodpoints_bool)
    bin_id_indivgoodpoints = np.arange(len(i_x_zcut_indivgoodpoints))

    i_x_zcut_pointstobin = i_x_zcut[indivgoodpoints_bool == False]
    i_y_zcut_pointstobin = i_y_zcut[indivgoodpoints_bool == False]
    kpc_r_zcut_pointstobin = kpc_r_zcut[indivgoodpoints_bool == False]
    kpc_z_zcut_pointstobin = kpc_z_zcut[indivgoodpoints_bool == False]
    flux_zcut_pointstobin = flux_zcut[indivgoodpoints_bool == False]
    fluxerror_zcut_pointstobin = fluxerror_zcut[indivgoodpoints_bool == False]
    
    accreted_bin_id = bin_accrete(kpc_r_zcut_pointstobin,kpc_z_zcut_pointstobin,cfgdict['snthresh'],
                                  flux_zcut_pointstobin,fluxerror_zcut_pointstobin,bin_start = num_indivgoodpoints)
    #np.savetxt('compute_voronoi_bins_temporararysave.txt',zip(kpc_r_zcut_pointstobin,kpc_z_zcut_pointstobin,accreted_bin_id))
    #kpc_r_zcut_pointstobin,kpc_z_zcut_pointstobin,accreted_bin_id = np.loadtxt('compute_voronoi_bins_temporararysave.txt',unpack=True)
    fixed_bin_ids = fix_bins(kpc_r_zcut_pointstobin,kpc_z_zcut_pointstobin,accreted_bin_id)
    cvt_bin_ids,cvt_err = cvt_bins(kpc_r_zcut_pointstobin,kpc_z_zcut_pointstobin,flux_zcut_pointstobin,fluxerror_zcut_pointstobin,fixed_bin_ids,niters=cfgdict['cvt_iters'])
    small_bin_fixed_bin_ids = remove_bad_bins(cvt_bin_ids,kpc_r_zcut_pointstobin,kpc_z_zcut_pointstobin,flux_zcut_pointstobin,fluxerror_zcut_pointstobin,0.75*cfgdict['snthresh'])
    final_bin_ids = remove_large_bins(small_bin_fixed_bin_ids,cfgdict['max_bin_size_pix'],bad_bin_id = -1)
    # cvt_bin_image = np.zeros(i_x2d.shape)-1
    # cvt_bin_image[i_y_zcut_indivgoodpoints,i_x_zcut_indivgoodpoints] = bin_id_indivgoodpoints
    # cvt_bin_image[i_y_zcut_pointstobin,i_x_zcut_pointstobin] = cvt_bin_ids
    # output_bin_image = np.zeros((3,i_x2d.shape[0], i_x2d.shape[1]), dtype=int)-1
    # #output_bin_image[i_y_zcut_indivgoodpoints,i_x_zcut_indivgoodpoints] = bin_id_indivgoodpoints
    # output_bin_image[0,i_y_zcut_pointstobin,i_x_zcut_pointstobin] = final_bin_ids
    # output_bin_image[1,i_y_zcut_pointstobin,i_x_zcut_pointstobin] = flux_zcut_pointstobin
    # output_bin_image[2,i_y_zcut_pointstobin,i_x_zcut_pointstobin] = fluxerror_zcut_pointstobin

    output_bin_image = np.zeros((i_x2d.shape[0], i_x2d.shape[1]), dtype=int)-1
    #output_bin_image[i_y_zcut_indivgoodpoints,i_x_zcut_indivgoodpoints] = bin_id_indivgoodpoints
    output_bin_image[i_y_zcut_pointstobin,i_x_zcut_pointstobin] = final_bin_ids

    output_header = fluxhdu.header
    output_header.set('cvt_err',cvt_err,'Max centroid shift in final CVT (kpc)')
    output_header.set('cvt_iter',cfgdict['cvt_iters'],'Number of CVT Iterations')
    output_header.set('cvt_sn',cfgdict['snthresh'],'S/N threshold for binning')
    output_header.set('cvt_maxz',cfgdict['max_z_kpc'],'Max height to be binned (kpc)')
    # writefits(cvt_bin_image-output_bin_image,cfgdict['output_prefix']+"_cvtdiff.fits",output_header)
    # writefits(cvt_bin_image,cfgdict['output_prefix']+"_cvt.fits",output_header)
    writefits(output_bin_image,cfgdict['output_prefix']+".fits",output_header)

#This program supersedes (hopefully) mc_voronoi.py:
if __name__ =='__main__':
    if len(sys.argv) != 2:
        sys.exit("Syntax: [configuration file]")

    cfgfile = sys.argv[1]
    go(cfgfile)

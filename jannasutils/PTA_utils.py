#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 13:39:21 2018

@author: jgoldstein
"""
import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# import from this package
from .astro_utils import radec_to_thetaphi

def radec_location_to_ang(location):
    """
    Convert Ra and Dec as read from the pulsar data to theta, phi coordinates
    
    Parameters
    ----------
    location: list
        location should contain 4 strings, encoding Ra in hours, Ra in minutes,
        Dec in degrees, Dec in arcminutes
    
    Returns
    -------
    float
        theta (polar) coordinate in radians, between 0 and pi
    float
        phi (azimuthal) coordinate in radians, between 0 and 2 pi
    """
    try:
        rah = float(location[0])
        ram = float(location[1])
        decd = float(location[2])
        decm = float(location[3])
    except:
        raise ValueError('Could not convert location {} to four floats'.format(location))
        
    # minus sign is attached to the hours of Dec only,
    # but the minutes should also be negative if the hours are
    if decd < 0:
        decm = -decm

    return radec_to_thetaphi(rah, ram, decd, decm)

def sky_localisation(logLs, sourceloc=[np.pi/2.0, 0.0], frac=0.9):
    """
    Calculate 90%(or other) likelihood area of sky map
    
    Parameters
    ----------
    logLs: NumPy Array
        sky map of loglikelihoods
        can also be 2D array with samples over likelihood
    sourceloc: list or NumPy array
        default: [np.pi/2.0, 0.0]
        location of the injected source in theta, phi coordinates
    frac: float
        default: 0.9
        fraction of the likelihood to include in the localisation area
    
    Returns
    -------
    float
        total likelihood
    float
        fraction of the sky containing >= frac of the likelihood
    list
        list of pixels in likelihood area
    bool
        True if injected source location in the likelihood area
    """
    
    Ls = np.exp(logLs)
    # reduce dimensionality if need (if amplitude is sampled it is now
    # marginalised over)
    if Ls.ndim > 1:
        Ls = np.sum(Ls, axis=1)
    # normalise likelihood
    Ltot = np.sum(Ls)
    Lsnorm = Ls / Ltot
    
    npix = len(Ls) # size of the map
    # variables used to check if source in likelihood area
    nside = hp.npix2nside(npix)
    sourcepix = hp.ang2pix(nside, sourceloc[0], sourceloc[1])
    source_in_loc = False
    
    # count max likelihood pixels until a fraction of frac is obtained
    pixels = []
    Lcount = 0.0
    while Lcount < frac:
        maxLpix = np.argmax(Lsnorm)
        # check if we count the source pix at any point
        if maxLpix == sourcepix:
            source_in_loc = True
        pixels.append(maxLpix)
        Lcount += Lsnorm[maxLpix]
        # set the max likelihood pixel to zero once it is counted
        Lsnorm[maxLpix] = 0.0

    loc_fraction = len(pixels) / float(npix)
    if np.allclose(loc_fraction, 1.0 / len(Ls)):
        print("sky localisation likely unresolved due to pixel size")
    return Ltot, loc_fraction, pixels, source_in_loc

class LikelihoodData:
    """
    Likelihood data for 2D (theta, phi) or 3D (A, theta, phi)
        
    """    
    def __init__(self, datapath, fname, GWfreq=None):
        """
        Initialise the likelihood data by loading from a logL file as is the
        output from the nullstream calculations
        
        Parameters
        ----------
        datapath: string
            'path/to/directory/with/data'
        fname: string
            filename without starting 'logL_' and extension
            (we need the same name later to load the amplitude samples if any)
        GWfreq: float
            optional GWfreq of the signal in the data
        """
        # get the loglikelhood data
        loglpath = os.path.join(datapath, 'logL_' + fname + '.txt')
        self._logl = np.loadtxt(loglpath)
        self._ls = np.exp(self._logl) 
        
        self._npix = len(self._logl)
        self._nside = hp.npix2nside(self._npix)
        
        if self._logl.ndim > 1:
            # amplitudes are sampled, not marginalised over
            # so get the values of the amp samples as well
            self.amp_marginalised = False
            ampspath = os.path.join(datapath, 'amps_' + fname + '.txt')
            self.amps = np.loadtxt(ampspath)
            self._namps = len(self.amps)
        else:
            self.amp_marginalised = True
            
        # get the pulsars
        pulsarpath = os.path.join(datapath, 'pulsars_' + fname + '.txt')
        self._pulsars = np.loadtxt(pulsarpath)
        
        # save GWfreq if given
        self.GWfreq = GWfreq
        
        
    def get_RaDec_likelihood(self, ra_hours, ra_mins, dec_degs, dec_mins):
        """
        Get the amplitude likelihood distribution at a specific sky location.
        
        Parameters
        ----------
        ra_hours: int, float
            whole hours of the Right Ascension
        ra_mins: int, float
            whole minutes of the Right Ascension
        dec_degs: int, float
            whole degrees of the Declination
        dec_mins: int, float
            whole arcminutes of the Declination

        """
        theta, phi = radec_to_thetaphi(ra_hours, ra_mins, dec_degs, dec_mins)
        return self.get_skyloc_likelihood(theta, phi)
            
    def get_skyloc_likelihood(self, theta, phi):
        """
        Get the amplitude likelihood distribution at a specific sky locations.
        
        Parameters
        ----------
        theta: float
            theta coordinate in radians (0 - pi)
        phi: float
            phi coordinate in radians (0 - 2pi)
            
        Returns
        float or array-like: likelihood at (theta, phi)
            if array-like, likelihood at each amplitude sampled
        """
        # select log likelihoods at location p
        p = hp.ang2pix(self._nside, theta, phi)
        return self._ls[p]
    
    def get_skymap(self):           
        """
        Get the likelihood sky map (marginalised over amplitude).
        """
        if self.amp_marginalised:
            return np.exp(self._logl - np.max(self._logl))
        
        # marginalise over amplitude assuming a log flat prior and 
        # evenly distributed log(A) samples
        ls_sum = np.sum(self._ls, axis=1)
        logl_sum = np.log(ls_sum)
        m = np.exp(logl_sum - np.max(logl_sum))
        # set any values that were masked in the logl map again to UNSEEN
        m[np.all(self._logl == hp.UNSEEN, axis=1)] = hp.UNSEEN
        return m
    
    def plot_skymap(self, source=None, save=None, mask=False, 
                    transparent=False, title=None, vmin=None, vmax=None):
        m = self.get_skymap()
        if mask:
            loc, pixels, _ = self.calc_localisation()
            for p in range(len(m)):
                if p not in pixels:
                    m[p] = hp.UNSEEN
        hp.mollview(m, title=title, min=vmin, max=vmax)
        for i, p in enumerate(self._pulsars):
            hp.projplot(p[0], p[1], c='w', marker='*', markersize=8)
        if source is not None:
            hp.projplot(source[0], source[1], c='w', marker='o', markersize=6)
        #if title is not None:
        #    plt.title(str(title))
        if save is None:
            plt.show()
        else:
            plt.savefig(str(save), transparent=transparent)
    
    def plot_likelihood_v_amp(self, theta, phi, ax=None, true_A=None, log=True):
        """
        Plot likelihood vs amplitude sampels at sky location (theta, phi)
        
        Only for non-marginalised amplitude
        """
        if ax is None:
            fig = plt.figure()
            return_fig = True
            ax = fig.add_subplot(111)
        else:
            return_fig = False
            
        if log:
            ax.set_xscale('log')
        
        if self.amp_marginalised:
            print("Can not plot amplitude samples because data is " +
                  "marginalised over amplitude.")
            return None
        
        ls = self.get_skyloc_likelihood(theta, phi)
        ax.plot(self.amps, ls)
        if true_A is not None:
            ax.vlines(true_A, 0, 1)
        
        if return_fig:
            return fig, ax
        else:
            return ax
    
    def calc_localisation(self, frac=0.9, source=[np.pi/2.0, 0.0]):
        Ltot, locfrac, pixels, hit = sky_localisation(self._logl, 
                                                frac=frac, sourceloc=source)
        if hit:
            hittxt = ''
        else:
            hittxt = 'not '
        print('sky localisation A{0} = {1}, with source {2}in loc'.format(
                int(100*frac), locfrac, hittxt))
        return locfrac, pixels, hit
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 17:08:52 2018

@author: jgoldstein
"""
import numpy as np
import numpy.random as rd
from scipy.integrate import quad
from scipy.interpolate import interp1d

# imports from this package
from .constants import RAD_TOT, HOUR_TOT, MIN_TOT, DEG_TOT, ARCMIN_TOT
from .generic_utils import isIterable


def radec_to_thetaphi(ra_hours, ra_mins, dec_degs, dec_mins):
    """
    Convert Ra, Dec to theta, phi coordinates
    
    Parameters
    ----------
    ra_hours: scalar
        whole hours of the right ascension
    ra_mins: scalar
        minutes of the right ascension
    dec_degs: scalar
        whole degrees of the declination
    dec_mins: scalar
        arcminutes of the declination
        
    Returns
    -------
    float
        theta (polar) coordinate in radians, between 0 and pi
    float
        phi (azimuthal) coordinate in radians, between 0 and 2 pi
    """
    #convert Ra and Dec to radians
    ra_rad = RAD_TOT * (ra_hours / HOUR_TOT + ra_mins / MIN_TOT)
    dec_rad = RAD_TOT * (dec_degs / DEG_TOT + dec_mins / ARCMIN_TOT)
    
    # Convert Ra and Dec angles to theta and phi angles
    return radec_reduced_to_thetaphi(ra_rad, dec_rad, unit='rad')

def thetaphi_to_radec(theta, phi, return_radians=False): 
    """
    Convert theta, phi to Ra, Dec
    
    Parameters
    ----------
    theta: float
        theta (polar) coordinate in radians, between 0 and pi
    phi: float
        phi (polar) coordinate in radians, between 0 and 2 pi
    return_radians: bool
        default = False
        if True, return ra and dec in radians
    
    Returns
    -------
    float: whole hours of the right ascension
    float: minutes of the right ascencion
    float: whole degrees of the declination
    float: arcminutes of the declination
        Note, for negative declination both degrees and arcminutes are negative
    if return_radians, return float, float being ra and dec in radians respectively
    """
    dec_phase = (0.5 * np.pi - theta) / (RAD_TOT)
    ra_phase = (2 * np.pi - phi) / (RAD_TOT)
    
    if return_radians:
        return ra_phase * RAD_TOT, dec_phase * RAD_TOT
    
    ra_hours = int(ra_phase * HOUR_TOT)
    ra_mins = (ra_phase - (ra_hours / HOUR_TOT)) * MIN_TOT
    dec_degs = int(dec_phase * DEG_TOT)
    dec_mins = (dec_phase - (dec_degs / DEG_TOT)) * ARCMIN_TOT

    return ra_hours, ra_mins, dec_degs, dec_mins

def radec_reduced_to_thetaphi(ra, dec, unit='rad'):
    """
    Convert right ascencion and declination (single numbers) to theta, phi.
    
    Parameters
    ----------
    ra: float
        Right ascencion in radians or degrees
    dec: float
        Declination in radians or degrees
    unit: {'degrees', 'deg', 'radians', 'rad'}
        default = 'rad'
        If 'degrees' or 'deg', ra and dec are interpreted as degrees
        If 'radians' or 'rad', ra and dec are interpreted as radians
        
    Returns
    -------
    float
        theta (polar) coordinate in radians, between 0 and pi
    float
        phi (azimuthal) coordinate in radians, between 0 and 2 pi
    """
    if unit not in ['degrees', 'deg', 'radians', 'rad']:
        raise ValueError("Unknown unit {} has to be one of 'degrees',"
                         " 'deg', 'radians', 'rad'".format(unit))
    
    # if ra dec given in degrees, convert to radians
    if unit in ['degrees', 'deg']:
        ra = ra * (RAD_TOT/DEG_TOT)
        dec = dec * (RAD_TOT/DEG_TOT)
    
    # convert ra dec angels (in radians) to theta and phi angles
    theta = 0.5 * np.pi - dec
    phi = 2 * np.pi - ra
    return theta, phi
    
def ang_distance(A, B, small_distance_treshold=0.05):
    """
    Calculate angular distance in radians between A and B sky locations
    
    Parameters
    ----------
    A, B: list or NumPy Array
        theta, phi coordinates (in radians)
        with theta between 0 an pi and phi between 0 and 2 pi
    
    Returns
    -------
    float: angular distance in radians
    """
    raA, decA = thetaphi_to_radec(*A, return_radians=True)
    raB, decB = thetaphi_to_radec(*B, return_radians=True)

    dist = np.arccos(np.sin(decA) * np.sin(decB) +
                     np.cos(decA) * np.cos(decB) * np.cos(raA - raB))
    
    # if points are very close, the above formula is inaccurate so we use
    # the quadrature sum approximation (which only works for small distances)
    if dist < small_distance_treshold:
        dist = ((raA - raB)**2. + (decA - decB)**2.)**0.5
    
    return dist


def random_location(r=None):
    """
    Get a random sky location within radius r of middle point (pi/2, 0)
    
    Parameters
    ----------
    r: float or None
        default = None
        if None, use the whole sky to pick from
        if float, pick a position within a distance r from (pi/2, 0)
        since r = pi is the maximum distance, any bigger r uses the whole
        sky as well
    
    Returns
    -------
    float
        theta (polar) coordinate in radians, between 0 and pi
    float
        phi (azimuthal) coordinate in radians, between 0 and 2 pi
    """
    # since maximum distance is pi, any bigger isn't needed
    if not r or r > np.pi:
        r = np.pi
    
    # for r bigger than the horizon to pole distance, use a rejection
    # technique to see if a random point falls in the disk
    # because maths is difficult and it's not too expensive here
    if r > 0.5 * np.pi:
        while True:
            phi = rd.uniform(0.0, 2 * np.pi)
            u = rd.uniform(-1.0, 1.0)
            theta = np.arccos(u)
            y = 0.5 * np.pi - theta
            if phi > np.pi:
                x = 2.0 * np.pi - phi
            else:
                x = phi
            # reject point not in disk
            if np.sqrt(x**2.0 + y**2.0) < r:
               break
           
    # use picking points on a disk if the disk is small enough to fit
    # between the north and south pole (so for r < pi/2)
    else:   
        D = r**2.0
        a = rd.uniform(0.0, D)
        psi = rd.uniform(0.0, 2*np.pi)
        
        x = np.sqrt(a) * np.cos(psi)
        y = np.sqrt(a) * np.sin(psi)
    
        phi = x%(2 * np.pi)
        theta = 0.5 * np.pi - y
    
    return theta, phi

class Cosmology:
    def __init__(self, h0=0.73, Omatter=0.25, Olambda=0.75, 
                 max_z=1., num_interp=1000):
        """
        Set up LCDM cosmology
        
        h0: float (default 0.73)
            dimensionless Hubble parameter
        Omatter: float (default 0.25)
            Matter content of the universe
        Olambda: float (default 0.75)
            Dark energy content of the universe
        """
        self.h0 = h0
        self.Omatter = Omatter
        self.Olambda = Olambda
        self.max_z = max_z
    
        # create a set of redshifts and associated luminosity distances
        # only using redshift 0-1 but can go just a bit over 1
        redshifts = np.linspace(0., self.max_z, num_interp)
        lumDistances = [self._lum_distance_gridpopulator(r) for r in redshifts]
        # set up interpolators to use in methods to calculate luminosity 
        # distance from redshift or vice versa
        self.interpLumDistance = interp1d(redshifts, lumDistances, assume_sorted=True)
        self.interpRedshift = interp1d(lumDistances, redshifts, assume_sorted=True)
        self.max_lumD = float(self.lumDistance(self.max_z))
    
    def _lum_distance_gridpopulator(self, redshift):
        """
        Calculate luminosity distance from redshift (slow)
        """
        d = lambda z: 1. / (self.Omatter*(1 + z)**3.0 + self.Olambda)**0.5    
        res, err = quad(d, 0, redshift)
        como_distance =  (2.99792458 / self.h0) * 1000 * res
        return (1 + redshift) * como_distance
    
    def lumDistance(self, redshift):
        """
        Calculate luminosity distance from  redshift (interpolation)
    
        Parameters
        ----------
        redshift: float or array
        
        Returns
        -------
        float or array (same as redshift)
            distance in Mpc
        """
        try:
            lumD = self.interpLumDistance(redshift)
            try:
                return float(lumD)
            except TypeError:
                return lumD
        except: 
            print('[warning] Failed to interpolate the redshift {}, using maximum'.format(redshift))
            return self.max_lumD #6535.962093540928
    
    def redshift(self, lumDistance):
        """
        Calculate redshift from luminosity distance (interpolation)
    
        Parameters
        ----------
        lumDistance: float or array
            luminosity distance in Mpc

        Returns
        -------
        float or array (same as lumDistance)
            redshift
        """
        try:
            z = self.interpRedshift(lumDistance)
            try:
                return float(z)
            except TypeError:
                return z
        except:
            print('[warning] Failed to interpolate for luminosity distance {} Mpc, using maximum'.format(lumDistance))
            return self.max_z
        
def GW_amplitude(chirp_mass, redshift, frequency, cosmology):
    """
    (dimensionless) GW amplitdue from a BBH
    
    parameters
    ----------
    chirp_mass: float
        chirp_mass in solar masses
    redshift: float
        cosmological redshift of the source
    frequency: float
        frequency of the GW (not orbit) in s^-1
    """ 
    # calculate luminosity distance in natural units (seconds)
    lum_dist_Mpc = cosmology.lumDistance(redshift)
    lum_dist_sec = lum_dist_Mpc * 1.029e14
    # calculate redshifted chirp mass in natural units (seconds)
    cm_z_msun = (1 + redshift) * chirp_mass
    cm_z_sec = cm_z_msun * 4.925e-6
    
    return (4. / lum_dist_sec) * (cm_z_sec**(5./3) * 
                              (np.pi * frequency)**(2./3))
    
def decorate_chirp_mass(func):
    def func_wrapper(mBH, q):
        if isIterable(mBH):
            mBH_exp = np.expand_dims(mBH, axis=1)
            return func(mBH_exp, q).squeeze()
        else:
            return func(mBH, q)
    return func_wrapper
            
@decorate_chirp_mass
def chirp_mass(mBH, q):
    """
    Calculate chirp mass
    
    Parameters
    ----------
    mBH: float or array
        total mass of the binary black hole
    q: float or array
        mass ratio of the binary, between 0 and 1
    
    Returns
    -------
    float or array: chirp mass
        If mBH and q are both arrays of length N and M, respectively, 
        return an array of shape (N, M)
    """
    return (mBH * q**(3./5)) / ((1+q)**(6./5))


def decorate_logexp_MMbulge(func):
    def func_wrapper(logmBulge, alpha, beta):
        if isIterable(alpha):
            alpha_exp = np.expand_dims(alpha, axis=1)
            return func(logmBulge, alpha_exp, beta).squeeze()
        else:
            return func(logmBulge, alpha, beta)
    return func_wrapper

@decorate_logexp_MMbulge
def logexp_M_Mbulge_linear(logmBulge, alpha, beta):
    """
    Expectted log black hole mass given linear M-Mbulge relation 
    (such as Kormendy & Ho 2013 or Haring & Rix 2004)
    
    KorHo should have a scatter of 0.29 and HaRix of 0.3 around the
    expected log BH mass 
    
    Parameters
    ----------
    logmBulge: float
        10log of galaxy bulge mass in solar masses
    alpha, beta: float
        parameters of the M-Mbulge relation
        KorHo: alpha=8.69, beta=1.17
        HaRix: alpha=8.2, beta=1.12
    
    Returns
    -------
    float: expected log black hole mass in solar masses
    i.e. log10(mBh / Msun)
    """
    return alpha + beta * (logmBulge - 11)

def logexpKorHo(logmBulge, alpha=8.69, beta=1.17):
    return logexp_M_Mbulge_linear(logmBulge, alpha, beta)
logexpKorHo.__doc__ = logexp_M_Mbulge_linear.__doc__

def logexpHaRix(logmBulge, alpha=8.2, beta=1.12):
    return logexp_M_Mbulge_linear(logmBulge, alpha, beta)
logexpHaRix.__doc__ = logexp_M_Mbulge_linear.__doc__

def logexpShankar(logmBulge, a=7.574, b=1.946, c=-0.306, d=-0.011):
    """
    Expected log black hole mass given Shankar et al (2016) relation, 
    which non-linear (power 3 polynomial).
    
    Scatter on the expected value is 0.4 dex
    
    Paramters
    ---------
    logmBulge: float
        10log of galaxy bulge mass in solar masses
    a, b, c, d: float
        paramters of the M-Mbulge relation
        Sharkar: a=7.574, b=1.946, c=-0.306, d=-0.011
    """
    x = logmBulge - 11
    return a + b*x + c*x**2.0 + d*x**3.0


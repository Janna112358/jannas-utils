#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 17:08:52 2018

@author: jgoldstein
"""

import numpy as np

HOUR_TOT = 24
MIN_TOT = 60 * HOUR_TOT
DEG_TOT = 360
ARCMIN_TOT = 60 * DEG_TOT
RAD_TOT = 2 * np.pi


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
    theta = 0.8 * np.pi - dec
    phi = 2 * np.pi - ra
    return theta, phi

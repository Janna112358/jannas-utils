#!/usr/bin/env python3

import unittest
import jannasutils as ju
import numpy as np
from numpy import random as rd

        
class Test_radec_thetaphi_conversions(unittest.TestCase):
    
    def test_declination(self):
        theta1, phi1 = ju.radec_to_thetaphi(0, 0, 90, 0)
        self.assertTrue(np.allclose(theta1, 0))
        theta2, phi2 = ju.radec_to_thetaphi(0, 0, 0, 0)
        self.assertTrue(np.allclose(theta2, np.pi/2.0))
        theta3, phi3 = ju.radec_to_thetaphi(0, 0, -90, 0)
        self.assertTrue(np.allclose(theta3, np.pi))
    
    def test_backforth_conversion(self):
        coordinates = np.array([[0.01, 0], [0.01, 0.8], [0.3, 1.9], [0.5, 0], 
                              [0.5, 1.0], [0.99, 0], [0.99, 1.99]]) * np.pi
        for c in coordinates:
            self.assertTrue(np.allclose(c, 
                        ju.radec_to_thetaphi(*ju.thetaphi_to_radec(*c))))
            
        locations = np.array([[1, 10, 89, 50], [19, 9, -37, -44], 
                              [10, 4, 0, 0], [0, 0, -89, -55]])
        for l in locations:
            self.assertTrue(np.allclose(l, 
                        ju.thetaphi_to_radec(*ju.radec_to_thetaphi(*l))))

    def test_location_string(self):
        location = ['19', '09', '-37', '44']
        theta, phi = ju.radec_location_to_ang(location)
        self.assertTrue(theta >= 0 and theta < np.pi)
        self.assertTrue(phi >= 0 and phi < 2 * np.pi)
        self.assertTrue(np.allclose([19, 9, -37, -44], 
                                    ju.thetaphi_to_radec(theta, phi)))
        
class Test_angular_distance(unittest.TestCase):

    def test_straight_distance(self):
        # for fixed phi (right ascension), the distance must be delta theta
        self.assertTrue(np.allclose(1.8, 
                    ju.ang_distance([0.0, 0.1], [1.8, 0.1])))
        # on the equation theta=pi/2, the distance must be delta phi (mod 2 pi)
        self.assertTrue(np.allclose(3.0, 
                    ju.ang_distance([np.pi/2.0, 0.0], [np.pi/2.0, 3.0])))
        self.assertTrue(np.allclose(1.0, 
                    ju.ang_distance([np.pi/2.0, 0.0], [np.pi/2.0, 2*np.pi-1.0])))
        
    def test_smaller_than_pi(self):
        for i in range(10):
            randomA = [rd.random() * np.pi, rd.random() * 2*np.pi]
            randomB = [rd.random() * np.pi, rd.random() * 2*np.pi]
            self.assertTrue(ju.ang_distance(randomA, randomB) <= np.pi)
    
    def test_small_distance(self):
        # test some distances that are just bigger than the treshold for using
        # the small distance approximation (default 0.05), and some that are 
        # small (in which case the approximation IS used)
        deltas = [0.1, 0.06, 0.05, 0.01, 0.001]
        for delta in deltas:
            self.assertTrue(np.allclose(delta, 
                    ju.ang_distance([0.0, 1.2], [0.0, 1.2+delta])))
            self.assertTrue(np.allclose(delta, 
                    ju.ang_distance([np.pi/2.0, 3.5], [np.pi/2.0, 3.5+delta])))
        
class Test_random_location(unittest.TestCase):
    
    def test_output_range(self):
        for i in range(20):
            theta, phi = ju.random_location()
            self.assertTrue(theta >= 0 and theta < np.pi)
            self.assertTrue(phi >= 0 and phi < 2 * np.pi)
        
    # test this separately since a different technique is used to pick points
    # when r < 0.5 pi
    def test_output_range_smallr(self):
        for j in range(20):
            theta, phi = ju.random_location(r=0.3)
            self.assertTrue(theta >= 0 and theta < np.pi)
            self.assertTrue(phi >= 0 and phi < 2 * np.pi)
        
    # actual radius constraint is stronger, but this should do as a test
    def test_radius(self):
        radius = [0.01 * np.pi, 0.1 * np.pi, np.pi]
        for r in radius: 
            for k in range(20):
                theta, phi = ju.random_location(r=r)
                self.assertTrue(abs(theta - np.pi/2.0) <= r)
                self.assertTrue(min(phi, 2*np.pi - phi) <= r)

class Test_cosmology(unittest.TestCase):
    
    def test_backforth_conversion(self):
        COSMOLOGY = ju.Cosmology(max_z=1.0)
        
        redshifts = np.array([0, 0.1, 0.5, 1.0])
        self.assertTrue(np.allclose(redshifts, 
                        COSMOLOGY.redshift(COSMOLOGY.lumDistance(redshifts))))
        
        distances = np.array([0, 10, 1000, 5000, COSMOLOGY.max_lumD])
        self.assertTrue(np.allclose(distances, 
                        COSMOLOGY.lumDistance(COSMOLOGY.redshift(distances))))

    def test_output_format(self):
        COSMOLOGY = ju.Cosmology(max_z=2.0)
        
        self.assertTrue(type(COSMOLOGY.lumDistance(1.0)) == float)
        self.assertTrue(type(COSMOLOGY.redshift(6000)) == float)

class Test_GWamplitude(unittest.TestCase):
    
    def test_value(self):
        COSMOLOGY = ju.Cosmology(max_z=1.1)
        self.assertTrue(np.allclose(5.7742474481883122e-17, 
                        ju.GW_amplitude(10**9, 1.0, 1.e-9, COSMOLOGY)))

class Test_chirpmass(unittest.TestCase):
    
    def test_value(self):
        self.assertTrue(np.allclose(418742239.1639287, ju.chirp_mass(1e9, 0.6)))
    
    def test_output_shape(self):
        # for two floats, output should be float
        self.assertTrue(type(ju.chirp_mass(1e9, 0.6)) == float)
        # for two arrays of length N, M, output should be of shape (N, M)
        mBH = np.array([1e9, 1.1e9, 1.2e9])
        q = np.array([0.5, 0.6])
        self.assertTrue(ju.chirp_mass(mBH, 0.6).shape == (3,))
        self.assertTrue(ju.chirp_mass(1e9, q).shape == (2,))
        self.assertTrue(ju.chirp_mass(mBH, q).shape == (3, 2))

class Test_MMbule(unittest.TestCase):
    
    def test_output_type(self):
        self.assertTrue(type(ju.logexp_M_Mbulge_linear(11., 8., 1.1)) == float)
        for f in [ju.logexpHaRix, ju.logexpKorHo, ju.logexpShankar]:
            self.assertTrue(type(f(11.)) == float)
    

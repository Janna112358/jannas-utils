#!/usr/bin/env python3

import unittest
import jannasutils as ju
import numpy as np

        
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
        self.assertTrue(np.allclose([19, 9, -37, -44], ju.thetaphi_to_radec(theta, phi)))
        

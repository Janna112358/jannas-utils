#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 17:08:52 2018

@author: jgoldstein
"""

import six

def isIterable(x):
    return (not isinstance(x, six.string_types)) and hasattr(x, "__iter__")
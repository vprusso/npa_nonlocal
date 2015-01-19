# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        util.py
# Purpose:     Various utility functions for npa_nonlocal.
#
# Author:      Vincent Russo (vrusso@cs.uwaterloo.ca)
#
# Created:     1/13/2015
# Copyright:   (c) Vincent Russo 2015
# Licence:     GNU
#------------------------------------------------------------------------------
'''

import os
import itertools


def generate_bit_strings(n, basis):
    '''
    Generates all bit strings of length n.
    '''
    return ["".join(seq) for seq in itertools.product(basis, repeat=n)]
    
    
def list_2_str(_list):
    '''
    Converts a list of objects into a concatenation of strings.
    '''
    return ' '.join(map(str, _list))
    
def check_equal(iterator):
    '''
    Checks if elements in an iterable object are all equal to each other.
    '''
    return len(set(iterator)) <= 1
    
    
def clear():
    '''
    Clears the shell of the spyder application. Use either clear() or cls()
    '''
    os.system('cls')
    return None


def clear_all():
    '''
    Clears all the variables from the workspace of the spyder application,
    '''
    cls()
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue
        
        del globals()[var]

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

import itertools

def generate_bit_strings(n):
    '''
    Generates all bit strings of length n
    '''
    return ["".join(seq) for seq in itertools.product("01", repeat=n)]
    
    
def list_2_str(_list):
    '''
    Converts a list of objects into a concatenation of strings
    '''
    return ''.join(map(str, _list))
# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        chsh.py
# Purpose:     This example  
#
# Author:      Vincent Russo (vrusso@cs.uwaterloo.ca)
#
# References: [1] Navascues, M. and Pironio, S. and A. Acin. A convergent 
#                 hierarchy of semidefinite programs characterizing the set of 
#                 quantum correlations. New Journal of Physics, 2008, 073013. 
#
# Created:     1/11/2015
# Copyright:   (c) Vincent Russo 2015
# Licence:     GNU
#------------------------------------------------------------------------------
'''

# TODO:

ops = generate_measurement_operators(2,2)
seq = generate_sequence(ops, "1")
M = generate_moment_matrix(seq)

A00 = ops[1]; A10 = ops[2]
B00 = ops[3]; B10 = ops[4]
n = len(seq)

chsh_exp = A00*B00 + A10*B00 + A00*B10 - A10*B10

bell_mat = bell_operator_matrix(chsh_exp, M)
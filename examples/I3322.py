# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        I3322.py
# Purpose:     This example computes the maximal violation of the I3322
#              inequality via the NPA hierarchy method.  
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

import npa_io
import moment_matrix
import bell_violation


# I3322 inequality has 3 inputs / 2 outputs
num_inputs = 3
num_outputs = 2
npa_level = "1"

# Generate measurement operators of specified number of inputs / outputs.
meas_ops = moment_matrix.generate_measurement_operators(num_inputs, num_outputs)

## Generate sequence of operators for specified level of hierarchy.
seq = moment_matrix.generate_sequence(meas_ops, npa_level)

## Generate the corresponding moment matrix of sequence.
M = moment_matrix.generate_moment_matrix(seq)

# I3322 operators 
A00 = meas_ops[0]; A01 = meas_ops[1]; A10 = meas_ops[2]; A11 = meas_ops[3]; 
A20 = meas_ops[4]; A21 = meas_ops[5]; 

B00 = meas_ops[6]; B01 = meas_ops[7]; B10 = meas_ops[8]; B11 = meas_ops[9]; 
B20 = meas_ops[10]; B21 = meas_ops[11]; 

# I3322 inequality:
# A_1 B_1 + A_2 B_1 + A_3 B_1 + A_1 B_2 + A_2 B_2 - A_3 B_2  + A_1 B_3 
#- A_2 B_2 - A_1 - 2*B_1 - B_2
I3322_exp = A00*B00 + A00*B10 + A00*B20 + A10*B00 + A10*B10 - A10*B20 + \
A20*B00 - A20*B10 - A00 - 2*B00 - B10

# Returns a matrix that is weighted from the coefficients of the Bell 
# expression corresponding to the entries in the moment matrix.
bell_mat = bell_violation.bell_operator_matrix(I3322_exp, M)

print npa_io.convert_python_matrix_to_matlab(bell_mat)
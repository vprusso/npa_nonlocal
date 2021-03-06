# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        chsh.py
# Purpose:     This example computes the maximal violation of the CHSH 
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

# CHSH inequality has 2 inputs / 2 outputs
num_inputs = 2
num_outputs = 2
npa_level = "1"

# Generate measurement operators of specified number of inputs / outputs.
meas_ops = moment_matrix.generate_measurement_operators(num_inputs, num_outputs)

# Generate sequence of operators for specified level of hierarchy.
seq = moment_matrix.generate_sequence(meas_ops, npa_level)

# Generate the corresponding moment matrix of sequence.
M = moment_matrix.generate_moment_matrix(seq,False)

# CHSH operators { A_0^0, A_0^1, B_0^0, B_0^1 }
A00 = meas_ops[0]; A11 = meas_ops[3];
B00 = meas_ops[4]; B11 = meas_ops[7]

# CHSH inequality:
# A_1 B_1 + A_2 B_1 + A_1 B_2 - A_2 B_2 - A_1 - B_1
chsh_exp = A00*B00 + A11*B00 + A00*B11 - A11*B11 - A00 - B00

# Returns a matrix that is weighted from the coefficients of the Bell 
# expression corresponding to the entries in the moment matrix.
bell_mat = bell_violation.bell_operator_matrix(chsh_exp, M)

print npa_io.convert_python_matrix_to_matlab(bell_mat)
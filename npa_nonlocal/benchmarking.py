# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        benchmarking.py
# Purpose:     File for benchmarking the time of various functions in the
#              npa_nonlocal project. 
#
# Author:      Vincent Russo (vrusso@cs.uwaterloo.ca)
#
# Created:     1/19/2015
# Copyright:   (c) Vincent Russo 2015
# Licence:     GNU
#------------------------------------------------------------------------------
'''

import util
import moment_matrix

from timeit import timeit, Timer

num_inputs = 2
num_outputs = 2
npa_level = 1
parallel_reps = 1
M = moment_matrix.MomentMatrix(num_inputs, num_outputs, npa_level, parallel_reps)


''' Time trials for: moment_matrix.'''
# 0.0375638008118
meas_ops_input_2_output_2_level_1_reps_1 = Timer(\
    "M = moment_matrix.MomentMatrix(2,2,1,1)",\
    setup="import moment_matrix")

# 0.334450960159
meas_ops_input_2_output_2_level_1_reps_2 = Timer(\
    "M = moment_matrix.MomentMatrix(2,2,1,2)",\
    setup="import moment_matrix")
    
# 65.810256958    
meas_ops_input_2_output_2_level_1_A_reps_2 = Timer(\
    "M = moment_matrix.MomentMatrix(2,2,'1+A',2)",\
    setup="import moment_matrix")

# 98.2671711445
meas_ops_input_2_output_2_level_1_AB_reps_2 = Timer(\
    "M = moment_matrix.MomentMatrix(2,2,'1+AB',2)",\
    setup="import moment_matrix")
    
# 83.6279420853    
meas_ops_input_2_output_2_level_1_A_reps_2_noclass = Timer(\
    "M = moment_matrix.generate_moment_matrix(seq)",\
    setup="import moment_matrix;\
    ops = moment_matrix.generate_measurement_operators(2,2,False,2);\
    seq = moment_matrix.generate_sequence(ops, '1+A')")

# 92.7540969849
meas_ops_input_2_output_2_level_1_AB_reps_2_noclass = Timer(\
    "M = moment_matrix.generate_moment_matrix(seq)",\
    setup="import moment_matrix;\
    ops = moment_matrix.generate_measurement_operators(2,2,False,2);\
    seq = moment_matrix.generate_sequence(ops, '1+AB')")    

'''Display time trial results'''
#print meas_ops_input_2_output_2_level_1_reps_1.timeit(1)

#print meas_ops_input_2_output_2_level_1_reps_2.timeit(1)

print meas_ops_input_2_output_2_level_1_A_reps_2.timeit(1)
#print meas_ops_input_2_output_2_level_1_AB_reps_2.timeit(1)

###
print meas_ops_input_2_output_2_level_1_A_reps_2_noclass.timeit(1)
#print meas_ops_input_2_output_2_level_1_AB_reps_2_noclass.timeit(1)
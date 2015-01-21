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

''' Time trials for: moment_matrix.generate_measurement_operators'''
meas_ops_input_2_output_2_reps_1 = Timer(\
    "moment_matrix.generate_measurement_operators(inputs,outputs,False,reps)",\
    setup="import moment_matrix; inputs = 2; outputs = 2; reps = 1")

meas_ops_input_2_output_2_reps_2 = Timer(\
    "moment_matrix.generate_measurement_operators(inputs,outputs,False,reps)",\
    setup="import moment_matrix; inputs = 2; outputs = 2; reps = 2")
            
''' Time trials for: moment_matrix.'''

'''Display time trial results'''
print meas_ops_input_2_output_2_reps_1.timeit(10)
print meas_ops_input_2_output_2_reps_2.timeit(10) 

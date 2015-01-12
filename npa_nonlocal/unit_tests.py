# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        unit_tests.py
# Purpose:     Unit testing functions for npa nonlocal functions.
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

import unittest
from moment_matrix import *

###############################################################################
##  MOMENT_MATRIX.PY UNIT TESTS
###############################################################################

class TestMomentMatrixFunctions(unittest.TestCase):
    '''
    '''
    def setUp(self):
        
        # Refer to "Matrix Size" column in Table-1 under I_1 in [1]
        self.seq_len_input_2_output_2_level_1 = 5
        self.seq_len_input_3_output_3_level_1 = 9      
        self.seq_len_input_4_output_4_level_1 = 13
        self.seq_len_input_5_output_5_level_1 = 17
        self.seq_len_input_6_output_6_level_1 = 21
        self.seq_len_input_7_output_7_level_1 = 25
        self.seq_len_input_8_output_8_level_1 = 29
        
        # Refer to "Matrix Size" column in Table-1 under I_{1+AB} in [1]
        self.seq_len_input_2_output_2_level_1_AB = 9
        self.seq_len_input_3_output_3_level_1_AB = 25
        self.seq_len_input_4_output_4_level_1_AB = 49
        self.seq_len_input_5_output_5_level_1_AB = 81
        self.seq_len_input_6_output_6_level_1_AB = 121
        self.seq_len_input_7_output_7_level_1_AB = 169
        self.seq_len_input_8_output_8_level_1_AB = 225       

        # Refer to "Matrix Size" column in Table-2
        self.seq_len_input_3_output_2_level_1 = 7
        self.seq_len_input_3_output_2_level_1_AB = 16
        self.seq_len_input_3_output_2_level_1_A_AB = 22

        
    def test_generate_sequence(self):
        '''
        Tests for generate_sequence function in moment_matrix.py
        '''
        meas_ops_input_2_output_2 = generate_measurement_operators(2,2)
        meas_ops_input_3_output_3 = generate_measurement_operators(3,3)
        meas_ops_input_4_output_4 = generate_measurement_operators(4,4)
        meas_ops_input_5_output_5 = generate_measurement_operators(5,5)
        meas_ops_input_6_output_6 = generate_measurement_operators(6,6)
        meas_ops_input_7_output_7 = generate_measurement_operators(7,7)
        meas_ops_input_8_output_8 = generate_measurement_operators(8,8)
        
        seq_ops_input_2_output_2_level_1 = \
            generate_sequence(meas_ops_input_2_output_2, 1)

        seq_ops_input_3_output_3_level_1 = \
            generate_sequence(meas_ops_input_3_output_3, 1)

        seq_ops_input_4_output_4_level_1 = \
            generate_sequence(meas_ops_input_4_output_4, 1)
 
        seq_ops_input_5_output_5_level_1 = \
            generate_sequence(meas_ops_input_5_output_5, 1)
            
        seq_ops_input_6_output_6_level_1 = \
            generate_sequence(meas_ops_input_6_output_6, 1)
            
        seq_ops_input_7_output_7_level_1 = \
            generate_sequence(meas_ops_input_7_output_7, 1)

        seq_ops_input_8_output_8_level_1 = \
            generate_sequence(meas_ops_input_8_output_8, 1)
                     
        # Ensure the length of the sequence generated agrees with the results
        # of the length in reference [1]. 
        self.assertEqual(len(seq_ops_input_2_output_2_level_1), \
                         self.seq_len_input_2_output_2_level_1) 
                         
        self.assertEqual(len(seq_ops_input_3_output_3_level_1), \
                         self.seq_len_input_3_output_3_level_1)  

        self.assertEqual(len(seq_ops_input_4_output_4_level_1), \
                         self.seq_len_input_4_output_4_level_1)  

        self.assertEqual(len(seq_ops_input_5_output_5_level_1), \
                         self.seq_len_input_5_output_5_level_1)  

        self.assertEqual(len(seq_ops_input_6_output_6_level_1), \
                         self.seq_len_input_6_output_6_level_1)  

        self.assertEqual(len(seq_ops_input_7_output_7_level_1), \
                         self.seq_len_input_7_output_7_level_1)  

        self.assertEqual(len(seq_ops_input_8_output_8_level_1), \
                         self.seq_len_input_8_output_8_level_1)                          
        
################################################################################
## MAIN UNIT TEST DRIVER
################################################################################
if __name__ == '__main__':

    # run unit tests for util.py
    moment_matrix_suite = unittest.TestLoader().loadTestsFromTestCase(TestMomentMatrixFunctions)
    unittest.TextTestRunner(verbosity=2).run(moment_matrix_suite)

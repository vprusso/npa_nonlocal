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
from bell_violation import *


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

        # Refer to "Matrix Size" column in Table-2 in [1]
        self.seq_len_input_3_output_2_level_1 = 7
        self.seq_len_input_3_output_2_level_1_AB = 16
        self.seq_len_input_3_output_2_level_1_A_AB = 22
        
        # Generate measurement operators of specified input / output length
        self.meas_ops_input_2_output_2 = generate_measurement_operators(2,2)
        self.meas_ops_input_3_output_3 = generate_measurement_operators(3,3)
        self.meas_ops_input_4_output_4 = generate_measurement_operators(4,4)
        self.meas_ops_input_5_output_5 = generate_measurement_operators(5,5)
        self.meas_ops_input_6_output_6 = generate_measurement_operators(6,6)
        self.meas_ops_input_7_output_7 = generate_measurement_operators(7,7)
        self.meas_ops_input_8_output_8 = generate_measurement_operators(8,8)
        
        self.meas_ops_input_3_output_2 = generate_measurement_operators(3,2)
    
        # Check that despite level input, permutations of order do not matter.
        self.level_1 = 1
        self.level_1_AB = "1+AB"
        self.level_1_A_AB = "1+A+AB"
        self.level_1_A_B_AB = "1+A+B+AB"
        self.level_1_B_A_AB = "1+B+A+AB"
        self.level_1_AB_A_B = "1+AB_A_B"
        self.level_1_AB_B_A = "1+AB_B_A"
        
        # Generate sequence operators of specified input / output level 1:
        self.seq_ops_input_2_output_2_level_1 = \
            generate_sequence(self.meas_ops_input_2_output_2, self.level_1)
        self.seq_ops_input_3_output_3_level_1 = \
            generate_sequence(self.meas_ops_input_3_output_3, self.level_1)
        self.seq_ops_input_4_output_4_level_1 = \
            generate_sequence(self.meas_ops_input_4_output_4, self.level_1)            
        self.seq_ops_input_5_output_5_level_1 = \
            generate_sequence(self.meas_ops_input_5_output_5, self.level_1)        
        self.seq_ops_input_6_output_6_level_1 = \
            generate_sequence(self.meas_ops_input_6_output_6, self.level_1)            
        self.seq_ops_input_7_output_7_level_1 = \
            generate_sequence(self.meas_ops_input_7_output_7, self.level_1)            
        self.seq_ops_input_8_output_8_level_1 = \
            generate_sequence(self.meas_ops_input_8_output_8, self.level_1)

        # Generate sequence operators of specified input / output level 1+AB:
        self.seq_ops_input_2_output_2_level_1_AB = \
            generate_sequence(self.meas_ops_input_2_output_2, self.level_1_AB)
        self.seq_ops_input_3_output_3_level_1_AB = \
            generate_sequence(self.meas_ops_input_3_output_3, self.level_1_AB)
        self.seq_ops_input_4_output_4_level_1_AB = \
            generate_sequence(self.meas_ops_input_4_output_4, self.level_1_AB)            
        self.seq_ops_input_5_output_5_level_1_AB = \
            generate_sequence(self.meas_ops_input_5_output_5, self.level_1_AB)        
        self.seq_ops_input_6_output_6_level_1_AB = \
            generate_sequence(self.meas_ops_input_6_output_6, self.level_1_AB)            
        self.seq_ops_input_7_output_7_level_1_AB = \
            generate_sequence(self.meas_ops_input_7_output_7, self.level_1_AB)            
        self.seq_ops_input_8_output_8_level_1_AB = \
            generate_sequence(self.meas_ops_input_8_output_8, self.level_1_AB)
            
        # Generate sequence operators of 3 input / 2 output level 1, 1+AB, and
        # level 1+A+AB
        self.seq_ops_input_3_output_2_level_1 = \
            generate_sequence(self.meas_ops_input_3_output_2, self.level_1)
        self.seq_ops_input_3_output_2_level_1_AB = \
            generate_sequence(self.meas_ops_input_3_output_2, self.level_1_AB)
        self.seq_ops_input_3_output_2_level_1_A_AB = \
            generate_sequence(self.meas_ops_input_3_output_2, self.level_1_A_AB)
            
        # Generate moment matrices of various input / output / level
        self.moment_matrix_input_2_output_2_level_1 = \
            generate_moment_matrix(self.seq_ops_input_2_output_2_level_1)  
            
        # Moment matrix dimensions variables
        self.moment_matrix_dim_input_2_output_2_level_1 = \
            int(math.sqrt(len(self.moment_matrix_input_2_output_2_level_1))) 
        
        
    def test_generate_sequence(self):
        '''
        Tests for generate_sequence function in moment_matrix.py
        '''                             
        # Ensure the length of the sequence generated agrees with the results
        # of the length in Table-1 under I_1 in reference [1]. 
        self.assertEqual(len(self.seq_ops_input_2_output_2_level_1), \
                         self.seq_len_input_2_output_2_level_1)                        
        self.assertEqual(len(self.seq_ops_input_3_output_3_level_1), \
                         self.seq_len_input_3_output_3_level_1)  
        self.assertEqual(len(self.seq_ops_input_4_output_4_level_1), \
                         self.seq_len_input_4_output_4_level_1)  
        self.assertEqual(len(self.seq_ops_input_5_output_5_level_1), \
                         self.seq_len_input_5_output_5_level_1)  
        self.assertEqual(len(self.seq_ops_input_6_output_6_level_1), \
                         self.seq_len_input_6_output_6_level_1)  
        self.assertEqual(len(self.seq_ops_input_7_output_7_level_1), \
                         self.seq_len_input_7_output_7_level_1)  
        self.assertEqual(len(self.seq_ops_input_8_output_8_level_1), \
                         self.seq_len_input_8_output_8_level_1)    
                         
        # Ensure the length of the sequence generated agrees with the results
        # of the length in Table-1 under I_1+AB in reference [1].
        self.assertEqual(len(self.seq_ops_input_2_output_2_level_1_AB), \
                         self.seq_len_input_2_output_2_level_1_AB) 
        self.assertEqual(len(self.seq_ops_input_3_output_3_level_1_AB), \
                         self.seq_len_input_3_output_3_level_1_AB) 
        self.assertEqual(len(self.seq_ops_input_4_output_4_level_1_AB), \
                         self.seq_len_input_4_output_4_level_1_AB) 
        self.assertEqual(len(self.seq_ops_input_5_output_5_level_1_AB), \
                         self.seq_len_input_5_output_5_level_1_AB) 
        self.assertEqual(len(self.seq_ops_input_6_output_6_level_1_AB), \
                         self.seq_len_input_6_output_6_level_1_AB) 
        self.assertEqual(len(self.seq_ops_input_7_output_7_level_1_AB), \
                         self.seq_len_input_7_output_7_level_1_AB) 
        self.assertEqual(len(self.seq_ops_input_8_output_8_level_1_AB), \
                         self.seq_len_input_8_output_8_level_1_AB) 
                         
        # Ensure the length of the sequence generated agrees with the results
        # of the length in Table-2 in reference [1]
        self.assertEqual(len(self.seq_ops_input_3_output_2_level_1), \
                         self.seq_len_input_3_output_2_level_1)
        self.assertEqual(len(self.seq_ops_input_3_output_2_level_1_AB), \
                         self.seq_len_input_3_output_2_level_1_AB)
        
        self.assertEqual(len(self.seq_ops_input_3_output_2_level_1_A_AB), \
                         self.seq_len_input_3_output_2_level_1_A_AB)
                         
        def test_check_moment_matrix_entry_equiv(self):
            '''
            Tests for check_moment_matrix_entry_equiv function in 
            moment_matrix.py
            '''
            # For 2 inputs and 2 outputs at level 1, expected that 
            # M[i,j] = M[j,i] for all i,j
            for i in range(self.moment_matrix_dim_input_2_output_2_level_1):
                for j in range(self.moment_matrix_dim_input_2_output_2_level_1):
                    assertTrue(check_moment_matrix_entry_equiv(\
                        self.moment_matrix_input_2_output_2_level_1[i,j]),\
                        self.moment_matrix_dim_input_2_output_2_level_1[j,i] )

            # For 2 input and 2 outputs at level 1 expected that 
            # M[i,i] = M[1,i] and M[i,i] = M[i,1] for all i.
            for i in range(self.moment_matrix_dim_input_2_output_2_level_1):
                for j in range(self.moment_matrix_dim_input_2_output_2_level_1):
                    assertTrue(check_moment_matrix_entry_equiv(\
                        self.moment_matrix_input_2_output_2_level_1[i,i]),\
                        self.moment_matrix_dim_input_2_output_2_level_1[1,i] )
                    assertTrue(check_moment_matrix_entry_equiv(\
                        self.moment_matrix_input_2_output_2_level_1[i,i]),\
                        self.moment_matrix_dim_input_2_output_2_level_1[i,1] )

###############################################################################
##  BELL_VIOLATION.PY UNIT TESTS
###############################################################################

class TestBellViolationFunctions(unittest.TestCase):
    def setUp(self):
        pass
    

###############################################################################
##  NPA_IO.PY UNIT TESTS
###############################################################################

class TestNPAIOFunctions(unittest.TestCase):
    def setUp(self):
        pass
    
################################################################################
## MAIN UNIT TEST DRIVER
################################################################################
if __name__ == '__main__':

    # run unit tests for moment_matrix.py
    moment_matrix_suite = unittest.TestLoader().loadTestsFromTestCase(TestMomentMatrixFunctions)
    unittest.TextTestRunner(verbosity=2).run(moment_matrix_suite)

    # run unit tests for bell_violation.py
    bell_violation_suite = unittest.TestLoader().loadTestsFromTestCase(TestBellViolationFunctions)
    unittest.TextTestRunner(verbosity=2).run(bell_violation_suite)
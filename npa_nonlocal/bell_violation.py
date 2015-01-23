# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        bell_violation.py
# Purpose:     This file contains functions for generating  
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

import util
import npa_io
import moment_matrix

import math

from sympy.physics.quantum import Dagger, HermitianOperator, IdentityOperator
from sympy.core.numbers import Infinity, Integer, NegativeOne
from sympy.matrices import zeros

class BellViolation(object):
    
    def __init__(self):
        pass


    def parse_bell_exp(self):
        '''
        Function that parses through a bell expression string and converts it to
        a sympy object. 
        '''
        pass    

def bell_operator_matrix(bell_exp, M):
    '''
    Given a Bell expression (bell_expr) and a moment matrix, (M) 
    this function returns a matrix where the entries corresponding to the Bell
    expression are weighted in the positions in the moment matrix.
    '''
    
    n = int(math.sqrt(len(M)))

    bell_terms = []
    bell_coeffs = []
    bell_args = list(bell_exp.args)
    for i in range(len(bell_args)):
        
        # Split Bell expression up into factors
        factors = bell_args[i].as_ordered_factors()
    
        negate = False
        terms = []
        coeffs = []
    
        for j in range(len(factors)):
            # Factors will have a NegativeOne instance if term is negative
            if isinstance(factors[j], NegativeOne):
                negate = True        
            # Any scalars are represented as Integers
            elif isinstance(factors[j], Integer):
                if negate == True:
                    coeffs.append(-factors[j])
                else:
                    coeffs.append(factors[j])
            # The remaining terms are Hermitian Operators (meas. ops.)
            elif isinstance(factors[j], HermitianOperator):
                terms.append(factors[j])    
        
        # If a term doesn't have a scalar, we put a "1" to denote its scalar is
        # just a factor of 1.
        if len(coeffs) == 0:
            coeffs.append(1)
            if negate == True:
                coeffs[0] = -1
            
        terms = reduce(lambda x,y : x*y, terms)     
        bell_terms.append(terms)
        bell_coeffs.append(coeffs)
    
    # Go through matrix and properly weight entries with the coefficients
    # derived from the Bell expression.
    bell_mat = zeros(n,n)
    
    # find all entries in matrix that correspond to term
    for i in range(len(bell_terms)):
        equiv_entries = moment_matrix.find_all_equiv_moment_matrix_entries(bell_terms[i], M)
        # weight the term appropriately in new bell matrix
        for j in range(len(equiv_entries)):
            # XXX NOTE: Not sure why this line is here. Needs to be here to work
            # but why is it required to block out the diagonal entries??
            #if not check_equal(equiv_entries[j]):
            bell_mat[ equiv_entries[j] ] = bell_coeffs[i]

    return bell_mat




#ops = generate_measurement_operators(2,2,False,1)
#seq = generate_sequence(ops, "1")
#M = generate_moment_matrix(seq)
#
#A00 = ops[0]; A01 = ops[1]; A10 = ops[2]; A11 = ops[3];
#B00 = ops[4]; B01 = ops[5]; B10 = ops[6]; B11 = ops[7]
#
#n = len(seq)
#
##chsh_exp = A00*B00 + A11*B00 + A00*B11 - A11*B11 - A00 - B00
#chsh_exp = A00*B00 + A10*B00 + A00*B10 - A10*B10 - A00 - B00 
#
##bell_mat = bell_operator_matrix(chsh_exp, M)
#
#chsh_args = list(chsh_exp.args)
#
#chsh_terms = []
#chsh_coeffs = []
#for i in range(len(chsh_args)):
#    
#    # Split Bell expression up into factors
#    factors = chsh_args[i].as_ordered_factors()
#
#    negate = False
#    terms = []
#    coeffs = []
#
#    for j in range(len(factors)):
#        #
#        if isinstance(factors[j], NegativeOne):
#            negate = True        
#        #
#        elif isinstance(factors[j], Integer):
#            if negate == True:
#                coeffs.append(-factors[j])
#            else:
#                coeffs.append(factors[j])
#        #
#        elif isinstance(factors[j], HermitianOperator):
#            terms.append(factors[j])    
#    
#    #        
#    if len(coeffs) == 0:
#        coeffs.append(1)
#        if negate == True:
#            coeffs[0] = -1
#        
#    terms = reduce(lambda x,y : x*y, terms)     
#    chsh_terms.append(terms)
#    chsh_coeffs.append(coeffs)
#
#bell_mat = zeros(n,n)
## find all entries in matrix that correspond to term
#for i in range(len(chsh_terms)):
#    equiv_entries = find_all_equiv_moment_matrix_entries(chsh_terms[i], M)
#    # weight the term appropriately in new bell matrix
#    for j in range(len(equiv_entries)):
#        if not check_equal(equiv_entries[j]):
#            bell_mat[ equiv_entries[j] ] = chsh_coeffs[i]

# 2,2
# 6,6
#0 -1 0 0 0 -1 0 0 0
#-1 X 0 0 0 1 0 0 1
#0 0 0 0 0 0 0 0 0
#0 0 0 0 0 0 0 0 0        
#0 0 0 0 0 1 0 0 -1
#-1 1 0 0 1 X 0 0 0
#0 0 0 0 0 0 0 0 0 
#0 0 0 0 0 0 0 0 0
#0 1 0 0 -1 0 0 0 0

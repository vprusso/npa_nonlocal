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


from moment_matrix import *
from npa_io import *

from sympy import pprint

def bell_operator_matrix(bell_exp, moment_matrix):
    '''
    Given a Bell expression (bell_expr) and a moment matrix, (moment_matrix) 
    this function returns a matrix where the entries corresponding to the Bell
    expression are weighted in the positions in the moment matrix.
    '''
    bell_args = list(bell_exp.args)

    bell_terms = []
    bell_coeffs = []
    for i in range(len(bell_args)):
        
        # Split Bell expression up into factors
        factors = bell_args[i].as_ordered_factors()
    
        negate = False
        terms = []
        coeffs = []
    
        for j in range(len(factors)):
            #
            if isinstance(factors[j], NegativeOne):
                negate = True        
            #
            elif isinstance(factors[j], Integer):
                if negate == True:
                    coeffs.append(-factors[j])
                else:
                    coeffs.append(factors[j])
            #
            elif isinstance(factors[j], HermitianOperator):
                terms.append(factors[j])    
        
        #        
        if len(coeffs) == 0:
            coeffs.append(1)
            if negate == True:
                coeffs[0] = -1
            
        terms = reduce(lambda x,y : x*y, terms)     
        bell_terms.append(terms)
        bell_coeffs.append(coeffs)
    
    bell_mat = zeros(n,n)
    for i in range(n):
        for j in range(n):
            if M[i,j] in bell_terms and i != j:
                print M[i,j], i,j


#ops = generate_measurement_operators(3,2,False)
#seq = generate_sequence(ops, "1")
#
#A00 = ops[1]; A01 = ops[2]; A10 = ops[3]; 
#A11 = ops[4]; A20 = ops[5]; A21 = ops[6]
#
#B00 = ops[7]; B01 = ops[8]; B10 = ops[9];
#B11 = ops[10]; B20 = ops[11]; B21 = ops[12]

ops = generate_measurement_operators(2,2)
seq = generate_sequence(ops, "1")
M = generate_moment_matrix(seq)

A00 = ops[1]; A10 = ops[2]
B00 = ops[3]; B10 = ops[4]
n = len(seq)

chsh_exp = A00*B00 + A10*B00 + A00*B10 - A10*B10

bell_mat = bell_operator_matrix(chsh_exp, M)

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
#print chsh_terms
#print flatten(chsh_coeffs)
#
#for i in range(n):
#    for j in range(n):
#        if M[i,j] in chsh_terms and i != j:
#            print M[i,j], i,j
#pprint(M)
#content = generate_latex_matrix(M)
#write_file("TEST",".tex",content)
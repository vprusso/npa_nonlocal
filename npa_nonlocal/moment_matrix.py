# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        moment_matrix.py
# Purpose:     This file contains functions for creating the moment matrix 
#              introduced in [1].  
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

from sympy import *
from sympy.physics.quantum import Dagger, HermitianOperator, IdentityOperator
from sympy.core.numbers import Infinity, Integer, NegativeOne
from sympy.matrices import zeros

from npa_io import *
from util import *


def generate_moment_matrix(seq):
    '''
    Given a sequence of level l (denoted S^l), the n x n moment matrix 
    corresponding to S^l may be written in the form:
                M^l(u,v) = <psi| U^* V |psi>
    for any entry
    '''    
    n = len(seq)
    M = zeros(n,n)
    for i in range(n):
        for j in range(n):            
            simp_seq = simplify_moment_matrix_entry( Dagger(seq[i]) * seq[j] )                 
            M[i,j] = simp_seq
            #M[i,j] = Dagger(seq[i]) * seq[j]
    return M


def generate_measurement_operators(num_inputs, num_outputs, \
                                   short_meas=True, parallel_reps=1):
    '''
    Measurement operators for Alice and Bob.

        short_meas: One can reduce the number of entries in the measurement
        operators by noting that they sum to the identity. For larger
        computations, this form is ideal, however for other cases, one may
        wish to generate the most general list of measurements.

        parallel_reps: Number of repetitions carried out by party.         
    '''    
    meas_ops = []    

    # Assuming that measurement labels are {0,1}-valued.
    meas_labels = generate_bit_strings(parallel_reps)    
    num_meas = len(meas_labels)    
    
    # Longer form of measurements are generated.
    if short_meas == False:
        for i in range(num_meas):
            for j in range(num_meas):
                alice_label = "A^" + meas_labels[i] + "_" + meas_labels[j]
                bob_label = "B^" + meas_labels[i] + "_" + meas_labels[j]

                alice_meas_op = HermitianOperator(alice_label)
                alice_meas_op.is_commutative = False

                bob_meas_op = HermitianOperator(bob_label)
                bob_meas_op.is_commutative = False                
                
                meas_ops.append(HermitianOperator(alice_meas_op))
                meas_ops.append(HermitianOperator(bob_meas_op))

    # Shorter form of measurements are generated
    if short_meas == True:
        for i in range(num_inputs-1 + num_outputs-1):

            alice_label = "A^" + str(i) + "_" + "0"
            bob_label = "B^" + str(i) + "_" + "0"

            alice_meas_op = HermitianOperator(alice_label)
            alice_meas_op.is_commutative = False
            #alice_meas_op.is_real = True

            bob_meas_op = HermitianOperator(bob_label)
            bob_meas_op.is_commutative = False            
            #bob_meas_op.is_real = True
            
            meas_ops.append(HermitianOperator(alice_meas_op))
            meas_ops.append(HermitianOperator(bob_meas_op))
      
    return sorted(meas_ops, key=default_sort_key)


def generate_sequence(meas_ops, level):
    '''
    A sequence is generated from the list of meas_ops from function 
    generate_measurement_operators. The level can be either an integer, or an
    intermediate level. Appropriate intermediate levels are "l+X" where "l" is
    the level of the sequence, and "X" is the intermediate steps. For instance:
        "l+A", "l+B", "l+AB", "l+A+B", "l+AB+A", etc.
    are all appropriate intermediate levels. 
    '''    
    seq = meas_ops[:]
    inter_med_seq = False
    
    # If the level is not an integer, but instead an intermediate level value, 
    # the format is expected to be "l+AB", "l+A", or "l+B" where "l" is an 
    # integer value corresponding to the level, and "AB", "A", or "B" are the
    # intermediate levels between l and l+1.
    if isinstance(level, str):
        l_str = level.split('+')
        level = int(l_str[0]) 
        inter_med = l_str[1:]

        inter_med = map(lambda x: x.strip(), inter_med)
        inter_med_seq = True    
    
    for i in range(1,level):        
        
        n = len(seq)                  

        # Process all A_a^x
        for j in range(n/2):
            for k in range(n/2):        
                if j != k:
                    seq.append( seq[j] * seq[k] )
                    
        # Process all B_b^y
        for j in range(n/2, n):
            for k in range(n/2, n):
                if j != k:
                    seq.append( seq[j] * seq[k] )
        
        # Process all A_a^x B_b^y
        for j in range(n/2):
            for k in range(n/2):
                seq.append( seq[j] * seq[k+(n/2)] )

    # If the sequence is intermediate, process the last bit 
    if inter_med_seq == True:
        n = len(seq)
        
        if "A" in inter_med:
            for j in range(n/2):
                for k in range(n/2):        
                    if j != k:
                        seq.append( seq[j] * seq[k] )
        if "B" in inter_med:
            for j in range(n/2, n):
                for k in range(n/2, n):
                    if j != k:
                        seq.append( seq[j] * seq[k] )
        if "AB" in inter_med:
            for j in range(n/2):
                for k in range(n/2):
                    seq.append( seq[j] * seq[k+(n/2)] )

    # Add in Identity operator to the front of the sequence
    I = IdentityOperator()
    #I = HermitianOperator("I")
    seq[0:0] = [I]
    
    return seq
    

def find_all_equiv_moment_matrix_entries(entry, mat):
    '''
    Given an entry in the moment matrix, this function finds all other entries
    that are equal to the entry in question. The indices are returned as a
    list of tuples.
    '''
    n = int(math.sqrt(len(mat))) 
    
    equiv_indices = []
    for i in range(n):
        for j in range(n):
            if check_moment_matrix_entry_equiv(entry, mat[i,j]):
                equiv_indices.append( (i,j) )
    return equiv_indices
    
   
def check_moment_matrix_entry_equiv(entry_1, entry_2):
    '''
    Given two entries in the moment matrix, this function checks whether or not
    they are within the same equivalence class by performing various checks 
    based on the properties of the projective measurement operators. 
    '''
    
    # First ensure the length of the entry is the same.
    if len(str(entry_1)) == len(str(entry_2)):
        
        flip_entry_1 = entry_1.args[::-1]
        flip_entry_2 = entry_2.args[::-1]
        
        # If entries are identical strings:
        if str(entry_1) == str(entry_2):
            return True
        # If entries are mirrored strings:
        elif ( str(entry_1) == str(flip_entry_2) ) or \
             ( str(entry_2) == str(flip_entry_1) ):
            return True
            
    # Otherwise, the entries are not equal 
    else:
        return False


def generate_moment_matrix_equivalence_dict(mat):
    '''
    Given a moment matrix, this function returns a dictionary of all respective
    equivalent entries in the matrix. Note, depending on the size of the 
    moment matrix, this may take a while. 
    '''
    n = int(math.sqrt(len(mat))) 
    
    # Go through each entry in the moment matrix and compare against every
    # other entry. Store the result in a dictionary. 
    equiv_dict = {}
    
    for i in range(n):
        for j in range(n):
            equiv_dict[(i,j)] = find_all_equiv_moment_matrix_entries(mat[i,j], mat)   
    return equiv_dict


def simplify_moment_matrix_entry(entry):
    '''
    Since the measurement operators pair-wise commute, i.e. [A_a^x, B_b^y] = 0, 
    and since they are also projection operators, i.e. P^2 = P, we can possibly
    reduce the number of terms in certain entries in the moment matrix.
    '''
           
    if isinstance(entry, IdentityOperator):
        pass
    else:
        # Measurement operators satisfy [A_a^x, B_b^y] = 0 for all operators 
        # A_a^x and B_b^y
        args = list(entry.args)
        for k in range(len(args)-1):
            # Comparing A vs. B lexographically to determine commutation swap
            if str(args[k])[0] > str(args[k+1])[0]:
                tmp = args[k]
                args[k] = args[k+1]
                args[k+1] = tmp
        args = Mul(tuple(args))
                   
        entry = reduce(lambda x,y : x*y, args)   
                        
        # Measurement operators are projective, so enforce that P^2 = P for 
        # any collection of measurement operators in sequence.
        args = list(entry.args)
        new_args = []
        for k in range(len(args)):               
            if isinstance(args[k], Integer):
                pass
            # Remove identities since they are simply absorbed by the term.
            elif isinstance(args[k], IdentityOperator):
                pass 
            elif isinstance(args[k], Pow):
                new_args.append(args[k].base)
            else:
             new_args.append(args[k])                    
        new_args = Mul(tuple(new_args))
                  
        entry = reduce(lambda x,y : x*y, new_args)   
                                            
    return entry
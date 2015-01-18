# -*- coding: utf-8 -*-
'''
#------------------------------------------------------------------------------
# Name:        npa_io.py
# Purpose:     This file contains functions for various I/O processing. 
#              Specifically, functionality is provided for: 
#                   - outputtng large moment matrices to LaTex format
#                   - user prompt command line functions
#                   - reading / writing various text / tex files
#
# Author:      Vincent Russo (vrusso@cs.uwaterloo.ca)
#
# Created:     1/11/2015
# Copyright:   (c) Vincent Russo 2015
# Licence:     GNU
#------------------------------------------------------------------------------
'''

import math

from sympy import *
from sympy import pprint

from util import *


###############################################################################
#   Prompt functions
###############################################################################
def disp_prompt():
    '''
    Displays the main prompt / option list to the user. 
    '''
    # TODO:
    #   Enter numeber of inputs / outputs 
    #   Option to:
    #       Output moment matrix
    #       Run SDP
    pass


###############################################################################
#   File I/O functions
###############################################################################
def write_file(file_name, file_ext, content):
    '''
    Write file of specified extension in local directory. 
    '''
    if "." not in file_ext:
        file_ext = "." + file_ext
    print ("Writing file %s%s...") % (file_name, file_ext)
    with open(file_name+file_ext, 'w') as out_file:
        out_file.write(content)
    print ("Done.")


###############################################################################
#   LaTeX functions
###############################################################################
def generate_latex_matrix(mat, block_mat_format=False):
    '''
    Generate source for a .tex file to output very large matrices. The variable
    block_mat_format allows the user to specify if they wish to output the 
    matrix in a block format in LaTeX.    
    '''
    
    tex_src = """
    \\documentclass[10pt]{article}
    \\usepackage[landscape,left=1cm,right=1cm,top=1cm,bottom=1cm,a3paper]{geometry}
    \\usepackage{graphicx}
    \\usepackage{etoolbox}
    \\usepackage{calc}
    \\usepackage{mathpazo}
    \\usepackage{amsmath}

    \\def\I{\\mathbb{1}}

    \\newlength{\\myx}
    \\newlength{\\myy}

    \\newcommand{\\microspace}{\\mspace{0.5mu}}

    \\newcommand{\\ket}[1]{
        \\lvert\\microspace #1 \\microspace \\rangle}

    \\newcommand{\\bra}[1]{
        \\langle\\microspace #1 \\microspace \\rvert}   

    \\newcommand{\\resizetopage}[1]{%
    \\settowidth{\\myx}{#1}%
    \\settototalheight{\\myy}{#1}%
    \\ifdimcomp{\myx}{<}{\\myy}{%
        \\resizebox*{!}{\\textheight}{#1}}{%
        \\resizebox*{\\textwidth}{!}{#1}}}
        \\newcommand{\\pagematrix}[1]{
        \\resizetopage{%
            $
            #1%
            $
        }}

    \\begin{document}

    \\pagematrix{
    """

    dim = int(math.sqrt(len(mat))) 
    if block_mat_format == True:
        block_sz = str("c"*dim)    
        start_mat_tag = "\\left(\\begin{array}{" + block_sz + "|" + block_sz + "} \n" 
        end_mat_tag = "\n \\end{array} \\right) }"

    else:
        start_mat_tag = "\\left(\\begin{array}{" + "c"*dim + "}\n"
        end_mat_tag = "\n \\end{array} \\right) }"
    
    mat = MutableDenseMatrix(mat)
           
    # Write the actual matrix in LaTeX format
    output = start_mat_tag
    for i in range(dim):
        if i > 0:
            output += "\n"
        for j in range(dim):
            if block_mat_format == True and i == dim/2:
                output += "\\hline \n"                
            if j == dim - 1:
                delim = " \\\ "
            else:
                delim = " & "
            entry = (str(mat[i,j]).replace("*"," ")).replace("I", " \\I ")
            s = entry.split()
            for k in range(len(s)):
                if "^" in s[k]:
                    s[k] = s[k].replace("^","^{").replace("_","}_{")
                    s[k] = s[k] + "} "
            entry = list_2_str(s)
            
            #output += "\\bra{\\psi}" + entry + "\\ket{\\psi}" + delim
            output += entry + delim
    output += end_mat_tag

    tex_src += output + "\n"     
    tex_src += "\\end{document}"
    
    return tex_src
    
    
###############################################################################
#   MATLAB functions
###############################################################################

def convert_python_matrix_to_matlab(mat):
    '''
    Takes a python matrix and converts it one that can be used in MATLAB.
    '''
    pass
    
#    return matlab_mat 


def output_matlab_script(mat, bell_exp):
    '''
    Given a moment matrix and Bell expression, this function writes a MATLAB
    script that uses CVX to solve the SDP. The script 
    '''

    matlab_bell_exp = convert_python_matrix_to_matlab(bell_exp)    
    
    output = """
    cvx_begin sdp
    \t %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
    \t %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX 
    \t variable M(dim,dim) semidefinite symmetric
    \t maximize trace(B * M)      
    \t subject to 
    \t \t    % entry M(1,1) = <psi| I I |psi> = 1
    \t \t    M(1,1) == 1;
    """
    
    
#ops = generate_measurement_operators(2,2)
#print ops

#seq = generate_sequence(ops, "1+AB")
#print seq
#mat = generate_moment_matrix(seq)
#print generate_latex_matrix(mat)


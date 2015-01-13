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


###############################################################################
#   Prompt functions
###############################################################################
def disp_prompt():
    '''
    '''
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

    \\def\I{\\mathbb{1}}

    \\newlength{\\myx}
    \\newlength{\\myy}

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
            output += entry + delim
    output += end_mat_tag

    tex_src += output + "\n"     
    tex_src += "\\end{document}"
    
    return tex_src
    

#ops = generate_measurement_operators(2,2)
#print ops

#seq = generate_sequence(ops, "1+AB")
#print seq
#mat = generate_moment_matrix(seq)
#print generate_latex_matrix(mat)


"""
NPA_Nonlocal: A suite of functions useful for calculating and viewing the 
              moment matrices described in [1]. Functionality is also provided
              for calculating the maximal violation of Bell inequalities via
              the NPA hierarchy. 
              
References: [1] Navascues, M. and Pironio, S. and A. Acin. A convergent 
                hierarchy of semidefinite programs characterizing the set of 
                quantum correlations. New Journal of Physics, 2008, 073013.
"""

from distutils.core import setup

setup(
    name='npa_nonlocal',
    version='1.0',
    author='Vincent Russo',
    author_email='vincentrusso1@gmail.com',
    packages=['npa_nonlocal'],
    url='http://vprusso.github.io/',
    keywords=[
        'sdp',
        'semidefinite programming',
        'sdp hierarchy',
        'moment matrix',
        'NPA',
        'noncommuting variable',
        'sdpa'],
    license='LICENSE',
    description='A suite of functions useful for calculating the moment matrix\
                 and Bell violation matrices covered in the NPA paper.',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python'
    ],
    install_requires=[
        "sympy >= 0.7.2",
        "numpy"
    ],
)
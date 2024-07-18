# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 03:19:35 2024

@author: saxen
"""

from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    ext_modules=cythonize("CCGEngine.pyx"),
    include_dirs=[np.get_include()]
)

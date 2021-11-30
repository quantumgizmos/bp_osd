from .bposd import bposd_decoder

import os
from . import __file__
def get_include():
    path = os.path.dirname(__file__)
    return path

f=open(get_include()+"/VERSION")
__version__=f.read()
f.close()
<<<<<<< HEAD
from .bposd import bposd_decoder

import os
from . import __file__
def get_include():
    path = os.path.dirname(__file__)
    return path

f=open(get_include()+"/VERSION")
__version__=f.read()
f.close()
=======
from ldpc import bposd_decoder
import os
import bposd
def get_include():
    path = os.path.dirname(bposd.__file__)
    return path

with open(f"{get_include()}/VERSION","r") as f:
    __version__=f.read()


>>>>>>> dev2

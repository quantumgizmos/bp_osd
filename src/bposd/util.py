import os
import bposd
def get_include():
    path = os.path.dirname(bposd.__file__)
    return path
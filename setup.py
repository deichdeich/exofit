from distutils.core import setup, Extension
import numpy as np
 
module1 = Extension('rsky', sources = ['c_code/rsky.c'])#, extra_compile_args=["-fopenmp"], extra_link_args=["-lgomp"])
 
setup (name = 'rsky',
        version = '1.0',
        description = 'This is a demo package',
	include_dirs = [np.get_include()],
        ext_modules = [module1])

module2 = Extension('occultquad', sources = ['c_code/occultquad.c'])#, extra_compile_args=["-fopenmp"], extra_link_args=["-lgomp"])
 
setup (name = 'occultquad',
        version = '1.0',
        description = 'This is a demo package',
	include_dirs = [np.get_include()],
        ext_modules = [module2])

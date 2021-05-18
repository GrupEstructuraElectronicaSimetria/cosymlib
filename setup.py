from numpy.distutils.core import setup, Extension
from numpy.distutils.command.install import install as _install
from distutils.sysconfig import get_python_lib
from distutils.dir_util import copy_tree
import os, sys

def get_version_number():
    main_ns = {}
    for line in open('cosymlib/__init__.py', 'r').readlines():
        if not(line.find('__version__')):
            exec(line, main_ns)
            return main_ns['__version__']

shape = Extension('cosymlib.shape.shp',
                  # extra_compile_args=['-std=c99'],
                  # include_dirs=include_dirs_numpy,
                  sources=['fortran/shp.pyf', 'fortran/shp.f90'])

on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    ext_modules = []
else:
    ext_modules = [shape]


class PostInstallCommand(_install):
    def run(self):
        _install.run(self)
        from shutil import copyfile
        dir = os.path.dirname(__file__)
        files = os.listdir(dir + '/cosymlib/.libs')
        for file in files:
            filename = os.path.join(dir, 'cosymlib', '.libs', file)
            copyfile(filename, os.path.join(dir, 'cosymlib', 'shape', file))

        site_dir = get_python_lib()
        files = [f for f in os.listdir('./cosymlib/shape/') if os.path.isfile('./cosymlib/shape/' + f)]
        for file in files:
            filename = os.path.join(dir, 'cosymlib', 'shape', file)
            copyfile(filename, os.path.join(site_dir, 'cosymlib', 'shape', file))


setup(name='cosymlib',
      version=get_version_number(),
      description='Continuous measures of shape and symmetry',
      author='Efrem Bernuz',
      author_email='komuisan@gmail.com',
      packages=['cosymlib',
                'cosymlib.shape',
                'cosymlib.molecule',
                'cosymlib.molecule.geometry',
                'cosymlib.molecule.electronic_structure',
                'cosymlib.file_io',
                'cosymlib.symmetry',
                'cosymlib.symmetry.pointgroup',
                'cosymlib.tools',
                'cosymlib.simulation',
                'cosymlib.gui'],
      package_data={'': ['ideal_structures_center.yaml',
                         'periodic_table.yaml']},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib', 'symgroupy', 'wfnsympy', 'PyYAML', 'huckelpy'],
      scripts=['scripts/cosym',
               'scripts/shape',
               'scripts/gsym',
               'scripts/cchir',
               'scripts/esym',
               'scripts/shape_map',
               'scripts/shape_classic'],
      cmdclass={'install': PostInstallCommand} if sys.platform.startswith('win') else {},
      ext_modules=ext_modules)

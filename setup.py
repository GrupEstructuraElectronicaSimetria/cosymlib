from numpy.distutils.core import setup, Extension

def get_version_number():
    for l in open('cosymlib/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__


shape = Extension('cosymlib.shape.shp',
                  # extra_compile_args=['-std=c99'],
                  #include_dirs=include_dirs_numpy,
                  sources=['fortran/shp.pyf', 'fortran/shp.f90'])

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
               'scripts/shape_map',
               'scripts/shape_classic'],
      ext_modules=[shape])

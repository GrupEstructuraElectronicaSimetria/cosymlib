from numpy.distutils.core import setup, Extension

def get_version_number():
    for l in open('cosym/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__


shape = Extension('cosym.shape.shp',
                  # extra_compile_args=['-std=c99'],
                  #include_dirs=include_dirs_numpy,
                  sources=['fortran/shp.pyf', 'fortran/shp.f90'])

setup(name='cosym',
      version=get_version_number(),
      description='Continuous measures of shape and symmetry',
      author='Efrem Bernuz',
      author_email='komuisan@gmail.com',
      packages=['cosym',
                'cosym.shape',
                'cosym.molecule',
                'cosym.molecule.geometry',
                'cosym.molecule.electronic_structure',
                'cosym.file_io',
                'cosym.symmetry',
                'cosym.tools',
                'cosym.simulation',
                'cosym.gui'],
      package_data={'': ['ideal_structures_center.yaml']},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib', 'symgroupy', 'wfnsympy', 'PyYAML'],
      scripts=['scripts/cosym'],
      ext_modules=[shape])

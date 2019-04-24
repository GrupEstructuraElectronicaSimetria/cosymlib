from numpy.distutils.core import setup, Extension

def get_version_number():
    for l in open('symeess/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__


shape = Extension('symeess.shape.shp',
                  # extra_compile_args=['-std=c99'],
                  #include_dirs=include_dirs_numpy,
                  sources=['fortran/shp.pyf', 'fortran/shp.f90'])

setup(name='symeess',
      version=get_version_number(),
      description='Continuous measures of shape and symmetry',
      author='Efrem Bernuz',
      author_email='komuisan@gmail.com',
      packages=['symeess',
                'symeess.shape',
                'symeess.molecule',
                'symeess.molecule.geometry',
                'symeess.molecule.electronic_structure',
                'symeess.file_io',
                'symeess.symmetry',
                'symeess.tools',
                'symeess.gui'],
      package_data={'': ['ideal_structures_center.yaml']},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib', 'symgroupy', 'wfnsympy', 'PyYAML'],
      scripts=['scripts/symeess',
               # 'scripts/shape_script.py',
               # 'scripts/shape_old_script.py',
               # 'scripts/symgroup_script.py',
               # 'scripts/wfnsym_script.py'
               ],
      ext_modules=[shape])

from numpy.distutils.core import setup, Extension


shape = Extension('symeess.shape.shp',
                  # extra_compile_args=['-std=c99'],
                  #include_dirs=include_dirs_numpy,
                  sources=['fortran/shp.pyf', 'fortran/shp.f90'])

setup(name='symeess',
      version='0.6',
      description='Programa de simetria',
      author='Efrem Bernuz',
      author_email='komuisan@gmail.com',
      packages=['symeess',
                'symeess.shape',
                'symeess.molecule',
                'symeess.molecule.geometry',
                'symeess.molecule.electronic_structure',
                'symeess.file_io',
                'symeess.symmetry',
                'symeess.tools'],
      package_data={'': ['ideal_structures_center.yaml']},
      include_package_data=True,
      scripts=['scripts/shape_script.py',
               'scripts/shape_old_script.py',
               'scripts/symgroup_script.py',
               'scripts/wfnsym_script.py'],
      ext_modules=[shape])

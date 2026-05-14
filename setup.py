from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import sys, os
import subprocess
import pathlib
import shutil
import os


try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
except ModuleNotFoundError:
    _bdist_wheel = object


def get_version_number():
    main_ns = {}
    for line in open('cosymlib/__init__.py', 'r').readlines():
        if not(line.find('__version__')):
            exec(line, main_ns)
            return main_ns['__version__']


ext = Extension(
    "cosymlib.shape.shp",   # must be inside your package namespace
    sources=["src/shpmodule.c"],  # your actual sources
)


class MesonBuildExt(build_ext):
    def run(self):

        # make compilation dir if needed
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print('self.build_lib:', self.build_lib)

        # define module dir to place fortran extension
        workdir = os.path.dirname(os.path.abspath(__file__))
        # workdir = self.build_lib
        install_dir = pathlib.Path(workdir, 'cosymlib', 'shape')

        # build with meson
        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', install_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        if '--inplace' in sys.argv:
            subprocess.check_call(['meson', 'install', '-C', self.build_temp])


class InstallWithBuildExt(install):
    def run(self):

        self.build_temp = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'build/temp')

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # self.install_lib = self.install_lib.replace('purelib/', '')

        # define module dir to install fortran extension
        install_dir = pathlib.Path(self.install_lib, 'cosymlib', 'shape')
        install_dir = os.path.abspath(install_dir)

        # build with meson and install
        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', str(install_dir)])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        subprocess.check_call(['meson', 'install', '-C', self.build_temp])

        # install
        import distutils.command.install as orig
        orig.install.run(self)


class MesonBdistWheel(_bdist_wheel):
    def run(self):

        self.build_temp = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'build/temp')

        # create compile directory (overwrite if exists)
        if os.path.exists(self.dist_dir):
            shutil.rmtree(self.dist_dir)
        os.makedirs(self.dist_dir)

        # define project root dir
        workdir = os.path.dirname(os.path.abspath(__file__))

        # define distribution dir
        dist_dir = pathlib.Path(workdir, self.dist_dir)
        dist_dir = os.path.abspath(dist_dir)

        # build with meson
        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', dist_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])

        self.root_is_pure = False
        super().run()

    def has_ext_modules(self):   # <-- extra insurance
        return True

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False


on_rtd = os.environ.get('READTHEDOCS') == 'True'
ext_modules = []

setup(name='cosymlib',
      version=get_version_number(),
      description='Continuous measures of shape and symmetry',
      author='Efrem Bernuz & Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=find_packages(where="."),
      package_data={'': ['ideal_structures_center.yaml',
                         'periodic_table.yaml']},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib', 'symgroupy', 'wfnsympy', 'PyYAML', 'huckelpy', 'pointgroup'],
      scripts=['scripts/cosym',
               'scripts/shape',
               'scripts/gsym',
               'scripts/cchir',
               'scripts/esym',
               'scripts/mosym',
               'scripts/shape_map',
               'scripts/shape_classic'],
      ext_modules=[ext],
      cmdclass={'build_ext': MesonBuildExt,
                'install': InstallWithBuildExt,
                'bdist_wheel': MesonBdistWheel,
                },
      url='https://github.com/GrupEstructuraElectronicaSimetria/cosymlib',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )

.. highlight:: rst

Installation
============

:program:`cosymlib` is available in both GitHub and PyPI repositories (https://pypi.org/project/cosymlib/).
Installation via PyPI is simpler and it is recommended for most users.

from PyPI (recommended)
-----------------------

This installation requires :program:`pip`  ( https://pip.pypa.io/en/stable/installing/) to be installed
in your system. We strongly recommend the use of python environments (https://docs.python.org/3/library/venv.html).
For most users the basic installation instructions are:

1. Create a virtual environment at path <venv>::

    python3 -m venv <venv>

2. Activate environment ::

    # on MAC / Linux (bash shell)
    source <venv>/bin/activate

    # on windows (powershell)
    C:\> <venv>\Scripts\Activate.ps1

3. Install cosymlib ::

    pip install numpy
    pip install cosymlib

4. Deactivate environmen ::

    deactivate

.. note::
    To use :program:`cosymlib` it is necessary to activate the environment every time a new shell is open.
    All the scripts contained in :program:`cosymlib` will be accessible in this environment. On windows
    it is necessary to type *python* before the script name ::

        python <script_name> <script_options>

.. note::
   On windows it may be necessary to add user execution permissions to activate the environment.
    To do this open a poweshell as administrator and type ::

      Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser


from source code
----------------

Alternatively you can download the latest version of :program:`cosymlib` from github using :program:`git` (https://git-scm.com)
and install it manually through :file:`setup.py` file using :program:`setuptools` (https://setuptools.readthedocs.io/).

:program:`cosymlib` contains libraries written in Fortran that require a compiler to build them.
Before installing cosymlib make sure you have a working Fortran compiler installed in your system.
For UNIX based systems you can install GNU Fortran Compiler from package repositories by opening a terminal and
typing the following commands:

- **Linux**

  On YUM-based systems (Fedora/RedHat/CentOS) ::

    sudo yum install yum-utils

  On APT-based systems (Debian/Ubuntu) ::

    sudo apt-get build-dep

- **Mac**

 1. Install command-line tools: ::

     xcode-select --install

 2. Get Homebrew following the instructions at https://brew.sh, and install GCC formula by: ::

     brew install gcc

- **Windows**

 1. Install windows development environment :program:`Visual Studio` (https://developer.microsoft.com/en-us/windows/downloads/)

 2. Install C/Fortran compiler for windows. We have tested and recommend  :program:`mingw` (https://www.mingw-w64.org)


To install :program:`cosymlib` download the source code using :program:`git` in your computer by typing: ::

    git clone https://github.com/GrupEstructuraElectronicaSimetria/cosymlib.git

This creates a copy of the repository in your computer. You can keep it updated by synchronizing it
with GitHub repository by using the command: ::

    git pull

Once this is done, move to the repository root directory (where :file:`setup.py` is found) and type the
following command to install :program:`cosymlib` : ::

    python setup.py install --user

.. note::
    :file:`requirements.txt` file located at the repository root directory contains a list of all dependency
    python modules needed for :program:`cosymlib` to run. If any of them are missing in your system you will
    need to install them before running :program:`cosymlib`.

In both cases (PyPI & Github installations) the code will be installed as a :program:`python` module. To check that it is properly
installed you can run the :program:`python` interpreter and execute: ::

   import cosymlib

if the execution do not show any errors :program:`cosymlib` has been installed successfully.

.. note::
    For users with Apple M1, scipy library might not properly install when following the next instructions,
    to solve this, install manually: ::

     brew install openblas
     brew install lapack
     brew install python
     pip install cython pybind11 pythran numpy
     OPENBLAS=$(brew --prefix openblas) CFLAGS="-falign-functions=8 ${CFLAGS}" pip install --no-use-pep517 scipy==1.7.0

.. note::
    When using an IDE remember to select the python interpreter in the hombrew path, to find it: ::

     which python3
     >> /opt/homebrew/bin/python3



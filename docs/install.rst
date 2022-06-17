.. highlight:: rst

Installation
============

:program:`cosymlib` is available in both the GitHub and PyPI repositories (https://pypi.org/project/cosymlib/).
Installation via PyPI is simpler and it is recommended for most users. Follow the instructions below to
install :program:`cosymlib` in your computer.

Installing cosymlib from PyPI (recommended)
-------------------------------------------

This installation requires :program:`pip` (https://pip.pypa.io/en/stable/installing/) to be installed
in your system. We strongly recommend the use of python environments, for more details on this, refer to
https://docs.python.org/3/library/venv.html. For most users the basic installation should proceed as follows:

1. Create a virtual environment at path <venv> ::

    $ python3 -m venv <venv>

2. Activate this virtual environment ::

    # on MAC / Linux
    $ source <venv>/bin/activate

    # on windows (powershell) [see note below]
    C:\> <venv>\Scripts\Activate.ps1

3. Install :program:`cosymlib` ::

    $ pip install numpy
    $ pip install cosymlib

4. Deactivate the virtual environment ::

    $ deactivate


To use :program:`cosymlib` you will need to activate the virtual environment every time that you open a new shell.
On Linux/MAC all the scripts contained in :program:`cosymlib` will be accessible in this environment: ::

    $ source <venv>/bin/activate
    $ <script_name> <script_options>
    $ deactivate

On Windows, to execute the scripts you should type *python* followed by the full path of the script name: ::

    C:\> python <venv>\Scripts\<script_name> <script_options>

.. note::
    On Windows it may be necessary to add user execution permissions to activate the environment.
    To do this, open a poweshell as administrator and type::

      Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

    You should do this only once in order to gain execution permissions.

Installing cosymlib's source code
---------------------------------

Alternatively, you can download the latest version of :program:`cosymlib` from github using :program:`git` (https://git-scm.com)
and install it manually through the :file:`setup.py` file using :program:`setuptools` (https://setuptools.readthedocs.io/).

:program:`cosymlib` contains libraries written in Fortran that require a compiler to build them.
Before installing :program:`cosymlib` make sure you have a working Fortran compiler installed in your system.
For UNIX based systems you can install the GNU Fortran Compiler from package repositories by opening a terminal and
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

 1. Install the Windows development environment :program:`Visual Studio` (https://developer.microsoft.com/en-us/windows/downloads/)

 2. Install C/Fortran compiler for Windows. We have tested and recommend  :program:`mingw` (https://www.mingw-w64.org)


To install :program:`cosymlib`, download the source code using :program:`git` in your computer by typing: ::

    git clone https://github.com/GrupEstructuraElectronicaSimetria/cosymlib.git

This creates a copy of the repository in your computer. You can keep it updated by synchronizing it
with the GitHub repository by using the command: ::

    git pull

Once this is done, move to the repository root directory (where :file:`setup.py` is found) and type the
following command to install :program:`cosymlib` : ::

    python setup.py install --user

.. note::
    The :file:`requirements.txt` file located at the repository root directory contains a list of all dependency
    python modules needed for :program:`cosymlib` to run. If any of them are missing in your system you will
    need to install them before running :program:`cosymlib`.

In both cases (PyPI & Github installations) the code will be installed as a :program:`python` module. To check that it is properly
installed you can run the :program:`python` interpreter and execute: ::

   import cosymlib

If the execution does not show any errors, then :program:`cosymlib` has been installed successfully.

.. note::
    For users with Apple M1, the :program:`scipy` library might not properly install when following the
    instructions above. To solve this, install it manually: ::

     brew install openblas
     brew install lapack
     brew install python
     pip install cython pybind11 pythran numpy
     OPENBLAS=$(brew --prefix openblas) CFLAGS="-falign-functions=8 ${CFLAGS}" pip install --no-use-pep517 scipy==1.7.0

.. note::
    When using an IDE, remember to select the python interpreter in the hombrew path. To find it: ::

     which python3
     >> /opt/homebrew/bin/python3



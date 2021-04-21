.. highlight:: rst

Installation
============
:program:`cosymlib` contains libraries written in Fortran that require a compiler to build them.
Before installing cosymlib make sure you have a working Fortran compiler installed in your system.

Requirements
____________

For UNIX based systems you can install it from package repositories by opening a terminal and
typing the following commands:

- **Linux**

  On YUM-based systems (Feodora/RedHat/CentOS) ::

    sudo yum install yum-utils

  On APT-based systems (Debian/Ubuntu) ::

    sudo apt-get build-dep

- **Mac**

 1. Install command-line tools: ::

     xcode-select --install

 2. Get Homebrew following the instructions at https://brew.sh, and install GCC formula by: ::

     brew install gcc

Install
_______

:program:`cosymlib` is available in both GitHub and PyPI repositories (https://pypi.org/project/cosymlib/).
Installation via PyPI is simpler and it is recommended for most users.

from PyPI
---------

This installation requires :program:`pip`  ( https://pip.pypa.io/en/stable/installing/) to be installed
in your system. Once pip is properly installed you can install :program:`cosymlib` by typing: ::

    pip install cosymlib --user

if your system contains both :program:`python2` and :program:`python3` installed and you intend to install :program:`cosymlib`
for :program:`python3` use: ::

    pip3 install cosymlib --user

from GitHub
-----------

Alternatively you can download the latest version of :program:`cosymlib` from github using :program:`git` (https://git-scm.com)
and install it manually through :file:`setup.py` file using :program:`setuptools` (https://setuptools.readthedocs.io/).

First, download the code using :program:`git` in your computer by typing: ::

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

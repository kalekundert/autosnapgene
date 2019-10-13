************
Installation
************
AutoSnapGene can be installed from PyPI::

   $ pip install autosnapgene

You can check that the installation worked::

   $ autosnapgene
   Perform batch operations on SnapGene files...

Staden io_lib
=============
If you want to be able to use AutoSnapGene to :doc:`add sequencing traces to 
files <traces>`, you also need to install `io_lib from the Staden package 
<https://github.com/jkbonfield/io_lib>`_.  Specifically, the ``convert_traces`` 
command must be on your ``$PATH``.  This program is used to convert traces to 
the ZTR format used internally by SnapGene.  If you're using Linux, you can 
probably install Staden io_lib from your package manager:

Ubuntu::

   $ sudo apt-get install staden-io-lib-utils

Fedora::

   $ sudo dnf install staden-io_lib

Arch (via AUR)::

   $ https://aur.archlinux.org/staden-io_lib.git
   $ cd staden-io_lib
   $ less PKGBUILD
   $ makepkg -si

Otherwise, you'll need to build Staden io_lib from its source code.  Refer to 
the above link for the most up-to-date instructions.


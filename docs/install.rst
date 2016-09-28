Dependency installation
=======================

Mandatory:

- `samtools-v1.3.1 <https://github.com/samtools/samtools/releases>`_
- bedtools-v2.22.1
- BWA-v0.7.12
- igraph-v0.7.1
- `MACS2-v2.1.1.20160309 <https://pypi.python.org/pypi/MACS2/2.1.1.20160309>`_
- `Picard <https://github.com/broadinstitute/picard/releases/tag/2.6.0>`_
- java-v1.8

Optional:

- STAR-v2.5.2b
- `RSEM-1.2.31 <https://github.com/deweylab/RSEM/releases>`_

.. note::
    Software with older versions might not work.


Software installation
=====================

Install stack
-------------

Download the `latest release of
stack <https://github.com/commercialhaskell/stack/releases>`_ for your
platform. For example, if your system is CentOS 6.5, choose "Linux 64-bit,
libgmp4 for CentOS 6.x".

Next, unpack the tarball and move ``stack`` executable to a directory
that is in your system path, e.g., ``/usr/bin``.

::

    tar zxf stack-x.x.x-linux-x64.tar.gz
    mv stack-x.x.x-linux-x64/stack /usr/bin


Download Taiji source and install GHC
--------------------------------------

Download the `latest source code <https://github.com/kaizhang/Taiji/releases>`_ and
unpack it.

::

    tar zxf Taiji-X.X.X.tar.gz

Go into the source code directory and install GHC.

::

    cd Taiji-X.X.X
    stack setup

Once you have a working copy of GHC, you can proceed to install the
dependencies of DBPnet. Under the source code directory, type:

::

    stack build --only-dependencies

and then install Taiji:

::

    stack install


.. _Dependencies:

Dependencies 
-------------

Bowtie2
*****************

**Version:** 2.4.5 or higher, but **don't use the version 2.5.0.**

.. code:: bash

    cd where/you/want/to/install/bowtie2
    wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-source.zip/download/ -O bowtie2-2.4.5-source.zip
    unzip bowtie2-2.4.5-source.zip
    cd bowtie2-2.4.5
    make

Then add the following line to your .bashrc file:

.. code:: bash

    export PATH="${PATH}:/path/to/bowtie2-2.4.5"


.. _dependencies_rnastructure:

RNAstructure
*****************

**About**: `RNAstructure <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-129>`_ is a software package for the prediction and comparison of RNA secondary structures, produced by the `Mathews lab <https://matthewslab.org/>`_. 

**Version**: 6.3.
It is likely that the latest version will work as well.

**Installation**: 
Check the `installation guide <https://rna.urmc.rochester.edu/RNAstructure.html#download>`_ for the latest version.

**Alternative installation**:
We have also provided a .zip file at the root of the repository, which contains the RNAstructure executable and the required files.

Download the zip file and unzip it then make:

.. code:: bash

    cd where/you/want/to/install/RNAstructure
    wget https://github.com/rouskinlab/dreem/raw/main/RNAstructureSource.zip
    unzip RNAstructureSource.zip
    cd RNAstructure
    make all

Then add the following line to your .bashrc file:

.. code:: bash

    export PATH="${PATH}:/path/to/RNAstructure/exe"
    export DATAPATH="/path/to/RNAstructure/data_tables"


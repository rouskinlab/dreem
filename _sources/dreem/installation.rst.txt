=====================
Installation
=====================


Using Docker (recommended)
--------------------------

Installing dependencies can be tricky. Docker handles this for you. Docker is also the easiest way to run DREEM on Windows.

Install Docker
^^^^^^^^^^^^^^

If you haven't already, `install Docker <https://docs.docker.com/get-docker/>`_.
Docker must be running in order to use the docker image. To check if docker is running, run the following command in the terminal::

    docker ps  

Install the Docker image
^^^^^^^^^^^^^^^^^^^^^^^^

There are two options for obtaining and installing the docker image:

**A. By pulling the image from Docker hub (not recommended at this time)**

.. code:: bash

    docker pull rouskinlab/dreem
    mkdir docker-data


        


**B. By cloning the repo from github**

.. code:: bash

    cd where/you/want/dreem
    git clone https://github.com/rouskinlab/dreem.git
    cd dreem
    docker build . -t dreem 
..       
    Run the Docker image
    ^^^^^^^^^^^^^^^^^^^^

    Here's a few examples of how to run the Docker image. See :ref:`API reference<API_reference>` for more.

    .. code:: bash

        # Display help
        docker run dreem --help

        # Run dreem on a single sample
        docker run dreem --fasta /path/to/ref.fasta --fastq1 /path/to/sample_R1.fastq --fastq2 /path/to/sample_R2.fastq --temp-dir /path/to/temp/ --out-dir /path/to/out/ --verbose --log /path/to/dreem_docker_run.log

        # Run dreem on multiple samples
        docker run dreem --fasta /path/to/ref.fasta --fastq1 /path/to/sample1_R1.fastq --fastq2 /path/to/sample1_R2.fastq --fastq1 /path/to/sample2_R1.fastq --fastq2 /path/to/sample2_R2.fastq --temp-dir /path/to/temp/ --out-dir /path/to/out/ --verbose --log /path/to/dreem_docker_run.log



..
    commented that out @Scott because I don't understand why not just giving docker the path to the data folder is not enough.
    **Once done installing and building the image, a container can now be made from that image, in which the analysis can be done:**
        1. First move all the data files that are necessary for the type of analysis that will be done. In most simple cases 
        this will be three files: two fastq files which contains the paired-end reads, and one fasta which contains all the reference sequences.
        
        2. For the docker container to access these, the docker container will mount the dreem/docker-data in the users filesystem, to a folder within the container called var/data/.
        Therefore, for the container to access the desired files, the arguments given to the containers entry point must be given in the form of the path as the container will view it.
        Anything in the docker-data folder will be know to the docker container as var/data/<path/within/docker-data> and therefore the command will be given as:

            ``docker run -iv $(pwd)/docker-data/:/var/data dreem-docker --fasta /var/data/ref.fasta --fastq1 /var/data/sample_R1.fastq --fastq2 /var/data/sample_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log``
        
        3. If run correctly, the docker-data folder should be populated with a folder labeled temp, where intermediary files are stored;
        an out folder, where the analysis results are stored; as well as a log file, which details the steps that occured during this executation of the dreem pipeline.

        The containers access to the files can be a little confusing and unintuitive, so below is an example of the correct arguments to be used for the files available in the repo:

            ``docker run -iv $(pwd)/docker-data/:/var/data dreem-docker --fasta /var/data/3509.fasta --fastq1 /var/data/3509-O-flank_1=bi1-ms2-DB_R1.fastq --fastq2 /var/data/3509-O-flank_1=bi1-ms2-DB_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log``






Using Pypi and virtualenv
-------------------------

Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

If you want to use other modules than ``dreem.draw``, make sure that you have all the dependencies installed. See :ref:`Dependencies` for more information.


Create a virtual environment and install DREEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    python3.10 -m venv dreem-env # you must use python 3.10 for now
    source dreem-env/bin/activate
    pip install dreem


Using Pypi and Conda
-------------------------

Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

If you want to use other modules than ``dreem.draw``, make sure that you have all the dependencies installed. See :ref:`Dependencies` for more information.



Create a virtual environment and install DREEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    conda create -n dreem python=3.10 # you must use python 3.10 for now
    conda activate dreem
    pip install dreem
    conda create -n dreem python=3.10 # you must use python 3.10 for now
    conda activate dreem
    pip install dreem


Using Source
------------------------------------

Best if you want to contribute to the project, or if you want to use the latest version of the code.

Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

Install:

- :ref:`Dependencies` if you want to use other modules than ``dreem.draw``.
- `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ 
- `xcode command line tools <https://www.freecodecamp.org/news/install-xcode-command-line-tools/>`_ (if running on mac OSX).

.. note::

    Unfortunately, Windows OS is not friendly to bioinformatics pipelines because of all the dependencies, so to run on windows it is suggested you run on Docker through `WSL2 <https://docs.docker.com/desktop/windows/wsl/>`_. 





.. code:: bash

   cd path/to/where/you/want/dreem
   git clone https://github.com/rouskinlab/dreem.git
   cd dreem
   conda env create -f dreem/env.yml
   conda activate dreem
   pip install .
   pytest 


.. note::

    The final line ``pytest`` is not required but helpful in ensuring all the correct dependencies have been installed and dreem can access them.


Test installation
-----------------

.. note::

    We still need to implement this feature.

Run:

.. code:: bash

    dreem test



========
Examples
========

Here are some examples to get started with the dreem pipeline. For each example it is assumed the conda enivorment or pyenv is active


---------
CLI
---------
    The following command will run dreem on some testfiles provided in the repo:

    .. code::
        
        dreem --fasta docker-data/3509.fasta \
              --fastq1 docker-data/3509-O-flank_1=bi1-ms2-DB_R1.fastq \
              --fastq2 docker-data/3509-O-flank_1=bi1-ms2-DB_R2.fastq\
               --verbose --autosect --out-dir docker-data/out --temp-dir docker-data/temp
    
    *Command explained:*

        * ``dreem`` - conda's entry point executable for the dreem program
        * ``--fasta`` - path to the reference fasta (reference sequences)
        * ``--fastq1`` - path to 1 of 2 of two files containing sequnece data in the form if "reads" 
        * ``--fastq2`` - path to 1 of 2 of two files containing sequnece data in the form if "reads" 
        * ``--verbose`` - forces dreem to output info on real-time processes
        * ``--out-dir`` - path for which dreem will write outputs of each step of analysis
        * ``--temp-dir`` - path for which dreem will write outputs of each step of analysis ( see documentation of modules I/O for more information on what files each module stores in either the temp or out dir)
    
    Since these test files are directly in the repo this command should be able to be run as is, if run from the repo directory

    .. note::

        The full list of arguments can be found by running ``dreem --help`` or in the :ref:`API reference <API_reference>`.

------------
Python
------------
.. note::
    * Since the paths in this example are given as relative (to the current working directory) 
        in order to run, this python script assume that the user is running it in the repo directory.

    * Running this file also assumes that the user has activated the conda envorment that was installed during the installation steps

The follwing lines can be copy and pasted into a python file and run as a script:

.. code-block:: python
    :linenos:
    
    from dreem.main import run

    if __name__ == '__main__':

        run(
            fasta="docker-data/3509.fasta",
            fastq1="docker-data/3509-O-flank_1=bi1-ms2-DB_R1.fastq",
            fastq2="docker-data/3509-O-flank_1=bi1-ms2-DB_R2.fastq",
            out_dir="docker-data/out/",
            temp_dir="docker-data/temp/",
            autosect=True
            )

To create a simple python file to run user-specified files use the same code but with the paths to your files

.. code-block:: python
    :linenos:   

    from dreem.main import run
    
    if __name__ == '__main__':

        run(
            fasta="docker-data/3509.fasta",
            fastq1="path/to/sample_R1.fastq",
            fastq2="path/to/sample_R2.fastq",
            out_dir="path/to/out/",
            temp_dir="path/to/temp/",
            autosect=True
            )
            
------------
Docker
------------
Docker containers have their own directory and therefore for docker to access files the directory needs to be mounted. 
For simplicity of use, it is recommended the user moves the files to be analyzed, to ``docker-data`` folder in the repo directory.


 .. code::

        docker run -iv $(pwd)/docker-data/:/var/data dreem \
        --fasta /var/data/ref.fasta \
        --fastq1 /var/data/sample_R1.fastq --fastq2 /var/data/sample_R2.fastq \
        --temp-dir /var/data/temp/ --out-dir /var/data/out/ \
        --verbose --log /var/data/dreem_docker_run.log


*Command explained:*
    
    * ``docker`` - executable entry point available as a result of installing docker
    * ``run`` - command for docker telling it to run a container based on the image
    * ``-iv`` - arguments for docker run, **i** indicates to run the container interactively, **v** indicates that docker should mount one of the container directory 
    * ``$(pwd)/docker-data/:/var/data`` - this command creates a mount between one of the containers directory and the hosts directory
        in this case, it mounts the var/data directory in the container to a directory in the current working directory.
    * ``dreem`` - the name of the image from which to build the container
    * ``--fasta`` - path to the reference fasta (reference sequences)
    * ``--fastq1`` - path to 1 of 2 of two files containing sequnece data in the form if "reads" 
    * ``--fastq2`` - path to 1 of 2 of two files containing sequnece data in the form if "reads" 
    * ``--verbose`` - forces dreem to output info on real-time processes
    * ``--out-dir`` - path for which dreem will write outputs of each step of analysis
    * ``--temp-dir`` - path for which dreem will write outputs of each step of analysis ( see documentation of modules I/O for more information on what files each module stores in either the temp or out dir)


Mounting different directories for docker:
------------------------------------------ 

The volume mount logic is as follows:

    ``dir/on/host/:/dir/on/container``

The mount argument could also be changed from:

    ``$(pwd)/docker-data/:/var/data``

into the folder you desire in this format:
    
    ``path/to/data/folder:var/data``

but it is not recommended that you change the /dir/on/container as the /var/data directory is written with the intent of being mounted to a directory on the hosts system
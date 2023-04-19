=====================
Installation
=====================

.. note::

    If you want to use other modules than ``dreem.draw``, make sure that you have all the dependencies installed. See :ref:`Dependencies` for more information.

Using Docker(Work in progress)
------------
    fInstalling the docker in pretty easy. So assuming have docker installed on your local computer run the following sequence of commands:

        ``cd dreem 
        docker build -t dreem .``

    Running a container from the image is a little bit more complicated and will be explained in detail below in the notes. 
    However for quick use follow these steps:
        1. First move all the data files that are necessary for the type of analysis that will be done. In most simple cases 
        this will be three files: two fastq files which contains the paired-end reads, and one fasta which contains all the reference sequences.
        
        2. For the docker container to access these, the docker container will mount the dreem/docker-data in the users filesystem, to a folder within the container called var/data/.
        Therefore, for the container to access the desired files, the arguments given to the containers entry point must be given in the form of the path as the container will view it.
        Anything in the docker-data folder will be know to the docker container as var/data/<path/within/docker-data> and therefore the command will be given as:

            docker run -iv $(pwd)/docker-data/:/var/data dreem --fasta /var/data/ref.fasta --fastq1 /var/data/sample_R1.fastq --fastq2 /var/data/sample_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log
        
        3.If run correctly, the docker-data folder should be populated with a folder labeled temp, where intermediary files are stored;
        an out folder, where the analysis results are stored; as well as a log file, which details the steps that occured during this executation of the dreem pipeline.

        The containers access to the files can be a little confusing and unintuitive, so below is an example of the correct arguments to be used for the files available in the repo:

            docker run -iv $(pwd)/docker-data/:/var/data dreem --fasta /var/data/3509.fasta --fastq1 /var/data/3509-O-flank_1=bi1-ms2-DB_R1.fastq --fastq2 /var/data/3509-O-flank_1=bi1-ms2-DB_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log




.. note::
    BUILD notes
    This first set of commands is pretty simple but it's important to understand what's happening with these commands so that we don't let docker become mystified as it so often is.
    The "cd dreem" ensures we're in the same directory as the Dockerfile. 

    In the second line, "docker" is the main commmand to which the user can provide a number of different functions. 
    In the case of this line the fucntion the user will call is the "build" function. The build function takes a directory
    as an argument which is in this case is the current working diretory i.e "." and then builds a image from the instructions 
    and stores in somewhere in the docker programs dockers data--better explaination here.
    
    The "-t" flag allows the user to give the build command a tag to call the image it builds in this case the user should call 
    the image "dreem" for simpplicty of use, but the user can change this name.

    If this build successfully runs, it will build an image which can then be used to run what is referred to as a Docker container. 
    Docker containers are sometimes referred to as a virtual machine but sadly that is not quite the case. But the main point that 
    needs to be understood here is that containers use the image to create a new isolated envorment to run the programs. 

    RUN notes
    Running a container from the dreem image can be a little confusing so different componants of the command are broken down below.


    At the time of these instructions being written the dreem/docker-data/ folder contains three sample files for testing the docker image:
            dreem/docker-data/3509-O-flank_1=bi1-ms2-DB_R1.fastq, dreem/docker-data/3509-O-flank_1=bi1-ms2-DB_R2.fastq and dreem/docker-data/3509.fasta
    These test files will be used to break down the reccommended run command:
        docker run -iv $(pwd)/docker-data/:/var/data dreem --fasta /var/data/3509.fasta --fastq1 /var/data/3509-O-flank_1=bi1-ms2-DB_R1.fastq --fastq2 /var/data/3509-O-flank_1=bi1-ms2-DB_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log

    The first thing that stands out is the "-iv" flag. This combined flag refers to two flags: 
    the -i flag, short for  --interactive, which signals the container will be run in interative session and the -v flag, short for --volume flag, which siginals the container to mount an outside folder to one in the containers filesystem.
    The information on the use of the --interative flag is not as straight forward as the --volume flag so until the information here is improved it will remain a black-box to most users.
    The --volume flag is used for mount folders outside of the container's file system. In this example the argument given to the -v flag is "$(pwd)/docker-data/:/var/data"
    This assumes the command is being made from the repo directory where the Dockerfile is located. the use of $(pwd) ensure that the absoulute path is given for the docker-data directory. The text following the colon: "/var/data" refers to directory that is built into the image in the docker file. 
    The part in front of the colon can be changed to something else but it is not reccommended for first time docker users.

    The rest of the arugments flags like fasta and --fastq1 are part of the dreem pipelines main arguments.




Using Pypi and pyenv (Work in progress)
---------------------------------------

::

    python3.11 -m venv venv
    source venv/bin/activate
    pip install dreem
    git clone https://github.com/yvesmartindestaillades/dreem
    pip install -r dreem/requirements.txt
    rm -fr dreem

.. note::

    Update dreem using ``pip install dreem -U``    

Using Conda
-----------

::

    conda install -c yvesmartindestaillades dreem
    [TODO]


Using Source (developers only)
------------------------------------

::

   cd path/to/where/you/want/dreem
   git clone https://github.com/yvesmartindestaillades/dreem
   cd dreem
   python3 -m venv venv
   source bin/activate
   pip install -r requirements.txt
   pip install .



=====================
Installation
=====================

.. note::

    If you want to use other modules than ``dreem.draw``, make sure that you have all the dependencies installed. See :ref:`Dependencies` for more information.

Using Docker
------------

[Instructions for installing Docker]

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



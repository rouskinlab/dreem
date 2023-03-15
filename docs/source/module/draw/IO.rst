
I/O files
++++++++++++++++++++++++

**Input**

.. code:: text

    sample_1.json
    sample_2.json
    ...

**Output**

.. code:: text

    sample_1/
      |- construct_1/
          |- section_1/
              |- mutation_fraction.html
              |- mutation_fraction_identity.html
              |- base_coverage.html
           ...
       ...
    sample_2/
       ...

.. note::

    Use the ``--flat`` option to output the files in a flat structure:

    .. code:: text
    
        sample_1/
          |- construct_1__section_1__mutation_fraction.html
          |- construct_1__section_1__mutation_fraction_identity.html
          |- construct_1__section_1__base_coverage.html
            ...
        sample_2/
            ...



Get started
=================

Import your data into a study object
------------------------------------

.. code::

    from dreem.draw import study
    import json

    # Define the list of samples you want to load
    my_samples = ['my_sample_1', 'my_sample_2', 'my_sample_3']

    # Load the data as a list of dictionaries
    data = []
    for sample in my_samples:
        with open('{}.json'.format(sample), 'r') as f:
            data.append(json.load(f))

    # Create a study object
    my_study = study.Study(data)

    # Print the list of available samples
    print(my_study.get_samples())


Plot the coverage across one of your samples
--------------------------------------------

.. code::

    study.

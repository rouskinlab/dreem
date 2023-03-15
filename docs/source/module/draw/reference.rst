
API Reference
++++++++++++++++++++++++

.. code:: text

    # Commands

    --input / -i
        Path to a dreem output format file. Can be specified multiple times.

    --flat / -f
        Flatten the output folder structure. This names your files [reference]__[section]__[plot_name].html
    
    --out_dir / -o
        The output folder to write the plots to. Defaults to ./draw

    --reference / -r
        One or several references to plot. If not specified, all references are plotted

    --section / -s
        One or several sections to plot. If not specified, all sections are plotted
    
    # Plot names

    --mutation_fraction
        Plots mutation_fraction plot. See Plots/gallery.

    --mutation_fraction_identity
        Plots mutation_fraction_identity plot. See Plots/gallery.

    --base_coverage
        Plots base_coverage plot. See Plots/gallery.

    --mutations_in_barcodes
        Plots mutations_in_barcodes plot. See Plots/gallery.

    --mutations_per_read_per_sample
        Plots mutations_per_read_per_sample plot. See Plots/gallery.


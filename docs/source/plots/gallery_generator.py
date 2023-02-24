import os, sys
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(__file__,'../../../..')))
import dreem

gallery_path = os.path.join(os.path.dirname(__file__), 'gallery.rst')

def beautify_title(title):
    title = title.replace('_', ' ')
    title = title[0].upper() + title[1:]
    return title

def strip_extension(filename):
    return '.'.join(filename.split('.')[:-1])


def write_plot(plot):
    name = strip_extension(plot)
    return f"""
{beautify_title(name)}
{"-"*len(name)}
                
.. raw:: html
    :file: plots_figs/{plot}
    
.. dropdown:: :fa:`eye,mr-1` **DOCSTRING**: {name}

    .. autofunction:: dreem.draw.study.Study.{name}
    

    """


def generate_rst():
    
    with open(gallery_path, 'w') as f:
        f.write("""\nGallery\n=========\n\n\n""")
    
    plots =  os.listdir(os.path.join(os.path.dirname(__file__), 'plots_figs'))
    plots.sort()
    
    for plot in plots:
        
        if not plot.endswith('.html'):
            continue
        
        with open(gallery_path, 'a') as f:
            f.write(write_plot(plot))
            


def generate_html():

    data = dreem.draw.load_dataset()

    study = dreem.draw.Study()
    study.df = data
    sample, reference, section, family = study.df.iloc[0][['sample', 'reference', 'section', 'family']]

    path_figs = os.path.join(os.path.dirname(__file__), 'plots_figs')
    
    if not os.path.exists(path_figs):
        os.mkdir(path_figs)
    
    for file in os.listdir(path_figs):
        if file.endswith('.html'):
            os.remove(os.path.join(path_figs, file))

    ################################################################################################
    # Generate HTML plots and save them in the docs/source/plots/plots_figs folder
    ################################################################################################
    
    study.mutation_fraction(
        sample = sample,
        reference = reference,
        section='ROI',
        to_html = os.path.join(path_figs, 'mutation_fraction.html')
    )
  
    study.mutations_in_barcodes(
        to_html = os.path.join(path_figs, 'mutations_in_barcodes.html')
    )
    
    study.deltaG_vs_mut_rates(
        sample = sample,
        section = 'ROI',
        base_type = ['A','C'],
        to_html = os.path.join(path_figs, 'deltaG_vs_mut_rates.html')
    )
    
    study.num_aligned_reads_per_reference_frequency_distribution(
        sample = sample,
        section = 'full',
        to_html = os.path.join(path_figs, 'num_aligned_reads_per_reference_frequency_distribution.html')
    )
    
    study.mutations_per_read_per_sample(
        sample = sample,
        section = 'full',
        to_html = os.path.join(path_figs, 'mutations_per_read_per_sample.html')
    )
    
    ################################################################################################
        
    
    
def main():
    generate_html()
    generate_rst()

            
if __name__ == '__main__':
    main()

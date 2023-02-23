import os

gallery_path = os.path.join(os.path.dirname(__file__), 'gallery.rst')

def beautify_title(title):
    title = title.replace('_', ' ')
    title = title[0].upper() + title[1:]
    return title

def strip_extension(filename):
    return filename[:-4]

with open(gallery_path, 'w') as f:
    f.write("""\nGallery\n=========\n\n\n""")
for plot in os.listdir(os.path.join(os.path.dirname(__file__), 'plots_figs')):
    if not plot.endswith('.png'):
        continue
    
    name = strip_extension(plot)
    with open(gallery_path, 'a') as f:
        f.write(f"""
.. dropdown:: :fa:`eye,mr-1` {beautify_title(name)} 

    .. autofunction:: dreem.draw.study.Study.{name}
    
.. image:: plots_figs/{plot}
    :width: 700
    :align: center
    :alt: {beautify_title(name)}
""")
        


import os

MODULES = [
    'aggregate',
    'align',
    'demultiplex',
    'draw',
    'vector',
    'cluster',
]

pretty_name = {
    'aggregate': 'Aggregation',
    'align': 'Alignment',
    'demultiplex': 'Demultiplexing',
    'draw': 'Drawing',
    'vector': 'Vectorization',
    'cluster': 'Clustering',
}

def copy_module(path, module):
    os.system(f'rm -fr {os.path.join(path, module)}')
    os.system(f'cp -fr {os.path.join(path, "template_module")} {os.path.join(path, module)}')

if __name__ == '__main__':
    for module in MODULES:
        copy_module(os.path.abspath(os.path.dirname(__file__)), module)
        # open main and write module in 2nd line
        with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), module, 'main.rst'), 'r') as f:
            lines = f.readlines()
        with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), module, 'main.rst'), 'w') as f:
            f.write(lines[0])
            f.write(f'{pretty_name[module]} \n')
            for line in lines[2:]:
                f.write(line)
            
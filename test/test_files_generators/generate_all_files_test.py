
import json
import os,sys
sys.path.append(os.getcwd())
import dreem.util as util
import dreem
current_dir = os.path.dirname(os.path.realpath(__file__))
notebooks = [f for f in os.listdir(current_dir) if f.endswith('.ipynb')]
for notebook in notebooks:
    print(notebook)
    with open(os.path.join(current_dir, notebook), 'r') as f:
        cells = json.load(f)['cells']
    for cell in cells:
        if cell['cell_type'] != 'markdown':
            cell['source'] = [line for line in cell['source'] if not line.startswith('%')]
            lines = ''.join(cell['source'])
            print(lines)
            exec(lines)

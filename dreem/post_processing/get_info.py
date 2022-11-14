from util import Path
import yaml

def get_attributes(file, type='txt'):
    path = Path()

    def open_fun(f,type):
        if type == 'txt':
            return f.read()
        if type == 'yml':
            return yaml.safe_load(f)

    if file == 'samples.csv':
        with open(path.sample_attribute_path, 'r') as f:
            attributes = open_fun(f,type)
        f.close()

    if file == 'library.csv':
        with open(path.library_attributes_path, 'r') as f:
            attributes = open_fun(f,type)
        f.close()

    return attributes

def echo_attributes(file):
    attributes = get_attributes(file)
    print(attributes)

def echo_attributes_samples():
    echo_attributes('samples.csv')

def echo_attributes_library():
    echo_attributes('library.csv')
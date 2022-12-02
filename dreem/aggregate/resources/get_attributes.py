import yaml, os

def read_library_attributes():
    with open(os.path.join( os.getcwd(), 'dreem/aggregate/resources/library_attributes.yml')) as f:
        library_attributes = yaml.safe_load(f)
    return library_attributes

def read_sample_attributes():
    with open(os.path.join( os.getcwd(),'dreem/aggregate/resources/sample_attributes.yml')) as f:
        sample_attributes = yaml.safe_load(f)
    return sample_attributes
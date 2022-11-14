
from resources import __file__ as resources_file
import yaml, os, subprocess
import string
import random
import shlex

dreem_ppai_file = __file__  

class Path(object):
    def __init__(self) -> None:    
        self.resources_file = '/'.join(resources_file.split('/')[:-1])
        self.dreem_ppai_file = '/'.join(dreem_ppai_file.split('/')[:-1])
        self.sample_attribute_path = self.resources_file+'/sample_attributes.yml'
        self.library_attributes_path = self.resources_file+'/library_attributes.yml'
        self.config_template = self.resources_file+'/config-template.yml'

def run_command(cmd):
    output, error_msg = None, None
    try:
        output = subprocess.check_output(
                cmd, shell=True, stderr=subprocess.STDOUT
        ).decode("utf8")
    except subprocess.CalledProcessError as exc:
        error_msg = exc.output.decode("utf8")
    return output, error_msg

def run_unix_cmds(cmds):
    cmds = [shlex.split(x) for x in cmds]
    outputs =[]
    for cmd in cmds:
        outputs.append(subprocess.Popen(cmd,
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.STDOUT)\
                                        .communicate())


def format_path(path):
    if path == '.':
        path = os.path.abspath('')
    if path[-1] != '/':
        path = path+'/'
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def get_random_string(length):
    # With combination of lower and upper case
    result_str = ''.join([random.choice(string.ascii_letters) for i in range(length)])
    # return random string
    return result_str


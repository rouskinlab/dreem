import os 
from dreem.util import run_cmd, clear_folder
import one_mut_per_read
if __name__ == '__main__':
    #run_cmd('pytest {} -v'.format(os.path.join(os.getcwd(),'test/test_files_generators')))
    clear_folder(os.path.join(os.getcwd(),'test','test_files'))
    one_mut_per_read.make_files()
    
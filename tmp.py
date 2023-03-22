import dreem
import json


f1 = '/Users/ymdt/src/dreem/test_files/my_test_sample_R1.fastq'
f2 = '/Users/ymdt/src/dreem/test_files/my_test_sample_R2.fastq'

with open(f1, 'r') as fh1, open(f2, 'r') as fh2:
    while True:
        r1 = fh1.readline()
        r2 = fh2.readline()
        
        if r1 == '':
            print('Done, no errors')
            break
        
        if r1[0] == '@':
            assert r1 == r2, 'Reads do not match: {} != {}'.format(r1, r2)

        else:
            if r2[0] == '@':
                raise ValueError('shift in reads')

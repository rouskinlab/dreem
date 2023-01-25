
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import pytest
import pandas as pd
import dreem.util as util
import os
from dreem import clustering
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import json 
import numpy as np

module = 'clustering'

sample_name = 'test_set_1'
MAX_CLUSTERS = 2

reads_partition = [100000]*2

half_sequence_length = 100
unpaired_bases = int(0.3*half_sequence_length)
shared_bases = [int(u*unpaired_bases) for u in np.linspace(0.0, 0.8, 4)]
mu_unpaired = [0.07]

sample_profile = {}
for sc in shared_bases:
    for mu in mu_unpaired:
        profile_name = 'r{}_ub{}_sb{}_mu{}'.format(reads_partition, unpaired_bases, sc, int(mu*100))
        sample_profile[profile_name] = {
            'n_reads': reads_partition,
            'n_AC': half_sequence_length,
            'n_unpaired': unpaired_bases,
            'n_shared': sc,
            'path_bv': os.path.join(test_files_dir, 'input', module, sample_name, 'r{}_ub{}_sb{}_mu{}/0.orc'.format(reads_partition, unpaired_bases, sc, int(mu*100))),
            'path_json': os.path.join(test_files_dir, 'output', module, sample_name, profile_name),
            'mu_unpaired': [mu]*2
        }
    
module_input = os.path.join(input_dir, module)
module_expected = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

inputs = ['clustering']

@pytest.mark.skip(reason="Too bugged")
def test_make_files():
    os.makedirs(os.path.join(test_files_dir, 'input', module, sample_name), exist_ok=True)
    os.makedirs(os.path.join(test_files_dir, 'expected_output', module, sample_name), exist_ok=True)
    files_generator.generate_files(sample_profile, module, inputs, [], test_files_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    
@pytest.mark.skip(reason="Too bugged")
def test_run():
    clustering.run(
        input_dir = os.path.join(module_input, sample_name),
        out_dir = module_output,
        num_runs= 5,
        min_iter=30,
        max_clusters=MAX_CLUSTERS,
        n_cpus=10
    )

@pytest.mark.skip(reason="Too bugged")
def test_assert_files_exist():
    assert os.path.exists(os.path.join(module_output, 'best_cluster_reads.json'))

@pytest.mark.skip(reason="not implemented")
def test_assess_performance():

    for name, sample in sample_profile.items():
    
        with open(os.path.join(module_output, 'best_cluster_reads.json'), 'r') as f:
            result_dict = json.load(f)
            n_reads = [ sum([k[1]==str(idx+1) for k in result_dict[name]['K2'].keys()]) for idx in range(len(sample['n_reads']))]

            result_np = pd.DataFrame.from_dict(result_dict[name]['K'+str(len(sample['n_reads']))]).to_numpy().T
            F1 = 1.0

            if len(n_reads) == 2:
                F1_class = []
                for i in range(2):
                    TP = np.count_nonzero(np.argmax(result_np[:n_reads[0]], axis=1)==i)
                    if TP != 0:
                        precision = TP/np.count_nonzero(np.argmax(result_np, axis=1)==i)
                        recall = TP/n_reads[0]
                        F1_class.append(2*precision*recall/(precision+recall))
                F1 = max(F1_class)
    
            if len(n_reads) == 3:
                n_reads_cumulative = [0] + n_reads[:-1]
                F1 = 0
                indices = [0,1,2]
                
                # Compute F1 score of each class for all combinations
                for c in range(len(n_reads)):
                    F1_class = []
                    start = sum(n_reads_cumulative[:c+1])
                    end = start + n_reads[c]
                    for i in indices:
                        TP = np.count_nonzero(np.argmax(result_np[start:end], axis=1)==i)
                        if TP != 0:
                            precision = TP/np.count_nonzero(np.argmax(result_np, axis=1)==i)
                            recall = TP/n_reads[c]
                            F1_class.append(2*precision*recall/(precision+recall))
                        else:
                            F1_class.append(0.0)
                    indices.pop(np.argmax(F1_class))
                    F1 += np.max(F1_class)
                F1 /= 3

            print('F1 score: ', F1)
            sample['F1'] = F1
            assert F1 > 0.5

                

if __name__ == '__main__':
    test_make_files()
    test_run()
    test_assert_files_exist()
    test_assess_performance()
import dreem
import json

dreem.draw.run(
        inpt = (json.load(open('/Users/ymdt/src/dreem/test_output/my_test_sample.json','r')),),
        out_dir = '/Users/ymdt/src/dreem/test_output',
        mutation_fraction = True,
        mutation_fraction_identity = True,
        base_coverage = True,
        mutations_in_barcodes = False,
        mutations_per_read_per_sample = True,
    )
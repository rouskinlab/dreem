import dreem
import json

s = dreem.draw.Study(
    data = json.load(open('/Users/ymdt/src/dreem/test_output/my_python_sample.json','r'))
)
bowtie2 --vers
a = s.mutation_fraction(
    sample=s.get_samples()[0],
    reference = 'reference_1',
    section='1-50',
    cluster='pop_avg',
    to_html = '/Users/ymdt/src/dreem/test_output/my_python_sample.html'
)
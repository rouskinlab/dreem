
import yaml, sys, os, random
import pandas as pd
import numpy as np
import subprocess
import json
from dreem.aggregate import poisson

test_files_dir = os.getcwd()+'/test_output' #os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..','test_output'))
input_dir = os.path.join(test_files_dir,'input')
prediction_dir = os.path.join(test_files_dir,'predicted_output')
output_dir = os.path.join(test_files_dir,'output')

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def generate_barcodes(barcode_length, n, min_barcode_distance):
    barcodes = []
    while(len(barcodes) < n):
        barcode = ''.join(random.choice('ATCG') for _ in range(barcode_length))
        if all(hamming_distance(barcode, b) > min_barcode_distance for b in barcodes):
            barcodes.append(barcode)
    return barcodes

def invert_sequence(seq):
    return ''.join([{'A':'T','T':'A','C':'G','G':'C'}[s] for s in seq])[::-1]

def next_base(base):
    return {'A':'T','T':'C','C':'G','G':'A',0:1}[base]

def create_sequence(length, bases=['A','T','C','G']):
    return ''.join([random.choice(bases) for _ in range(length)])

def copy_prediction_as_results(path_predicted, path_output):
    assert os.path.exists(path_predicted), 'The predicted output folder doesn\'t exist'
    os.makedirs(path_output, exist_ok=True)
    os.system('cp -r {} {}'.format(path_predicted, path_output))
    
def print_fasta_line(f, id, seq):
    f.write('>{}\n{}\n'.format(id, seq))

def print_fastq_line(f, id, seq, qual):
    f.write('@{}\n{}\n+\n{}\n'.format(id, seq, qual))    

def generate_fastq_files(path, sample_profile, construct=None, max_barcode_muts=None):
    """Write a fastq file with the given parameters
    
    Arguments:
        path {str} -- where to write the fastq files
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'reads' -- sequence to use for each construct [list]
                'number_of_reads' -- number of reads to generate for each construct [list]
                'mutations' -- number of mutations to introduce in each read [list]
                'deletions' -- number of deletions to introduce in each read [list]
                'insertions' -- number of insertions to introduce in each read [list]
        construct {str} -- if specified, only generate fastq files for this construct
    """
    if construct is not None:
        f_prefix = construct
    else:
        f_prefix = path.split('/')[-1]
    fastq1_name, fastq2_name = os.path.join(path, f_prefix + '_R1.fastq'), os.path.join(path, f_prefix + '_R2.fastq')
    with open(fastq1_name, 'w') as f1, open(fastq2_name, 'w') as f2:
        # write placeholder reads
        for c, v in sample_profile.items():
            for i in range(v['number_of_reads']):
                if max_barcode_muts is not None:
                    bs, be = v['barcode_start'], v['barcode_start'] + len(v['barcodes'])
                    if len([m for m in v['mutations'][i] if m >= bs and m < be]) > max_barcode_muts:
                        continue
                sequence = v['reads'][i]
                cigar = make_cigar(len(sequence),v['mutations'][i], v['deletions'][i] , v['insertions'][i])
                
                print_fastq_line(f1, '{}:{}:{}'.format(c, i, cigar), sequence, 'F'*len(sequence))
                print_fastq_line(f2, '{}:{}:{}'.format(c, i, cigar), invert_sequence(sequence), 'F'*len(sequence))
    
def generate_fasta_file(filename, sample_profile):
    """Write a fasta file with the given parameters
    
    Arguments:
        filename {str} -- where to write the fasta file
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'reads' -- sequence to use for each construct [list]
    """
    with open(filename, 'w') as f:
        for c, v in sample_profile.items():
            print_fasta_line(f, c, v['reference'])
            
def generate_library_file(filename, sample_profile):
    """Write a library file with the given parameters

    Arguments:
        filename {str} -- where to write the library file
        sample_profile {dict} -- dictionary with the following keys
            constructs {list} -- list of construct names [list]
                barcodes {list} -- list of barcodes to use for each construct [list]
                barcode_start {int} -- where to start the barcode in the read [int]
    """

    df = pd.DataFrame(sample_profile).T.reset_index().rename(columns={'index':'construct'})
    content_cols = ['construct']
    if 'barcodes' in df.columns:
        df.rename(columns={'barcodes':'barcode'}, inplace=True)
        content_cols += ['barcode_start','barcode']
    if 'sections' in df.columns:
        df.rename(columns={'sections':'section'}, inplace=True)
        content_cols += ['section','section_start','section_end']
        df = df.explode(['section','section_start','section_end'])[content_cols]
    df['some_random_attribute'] = 'some_random_value'
    df.to_csv(filename, index=False)

def generate_demultiplexed_fastq_files(folder, sample_profile):
    """Demultiplex a fastq file based on a library file

    Arguments:
        folder {str} -- where to write the demultiplexed fastq files
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'reads' -- sequence to use for each construct [list]
                'number_of_reads' -- number of reads to generate for each construct [list]
                'mutations' -- number of mutations to introduce in each read [list]
                'deletions' -- number of deletions to introduce in each read [list]
                'insertions' -- number of insertions to introduce in each read [list]
    """

    for c, v in sample_profile.items():
        generate_fastq_files(folder, {c:v}, c, max_barcode_muts=1)


def print_sam_header(f, construct, len_sequence):
    f.write('@HD VN:1.0 SO:unsorted\n' + '@SQ SN:{} LN:{}\n'.format(construct, len_sequence))

def print_sam_lines(f, read_name, sequence, construct, cigar):
    f.write('{}\t3\t{}\t{}\t{}\t{}\t*\t0\t{}\t{}\n'.format(read_name, construct, 1, 255, cigar, len(sequence), sequence))
    f.write('{}\t19\t{}\t{}\t{}\t{}\t*\t0\t{}\t{}\n'.format(read_name, construct, 1, 255, cigar, len(sequence), invert_sequence(sequence)))


def make_cigar(len_sequence, mutations, insertions, deletions):
    """Create a cigar string for a read with the given mutations and indels
    The sequence is assumed to be 1-based.
    The cigar string is 1-based.
    The mutations and indels are 0-based.

    The mutations and indels are assumed to be sorted. This is checked at the beginning of the code.

    Arguments:
        len_sequence {int} -- length of the sequence [int]
        mutations {list} -- list of mutations [list]
        insertions {list} -- list of insertions [list]
        deletions {list} -- list of deletions [list]
    """

    # check that the mutations and indels are sorted
    assert sorted(mutations) == mutations
    assert sorted(insertions) == insertions
    assert sorted(deletions) == deletions

    for i in range(len(mutations)-1):
        assert mutations[i]+1 != mutations[i+1], 'Consecutive mutations are not allowed'
    
    for i in range(len(insertions)-1):
        assert insertions[i]+1 != insertions[i+1], 'Consecutive insertions are not allowed'
    
    for i in range(len(deletions)-1):
        assert deletions[i]+1 != deletions[i+1], 'Consecutive deletions are not allowed'

    # initialize the cigar string
    cigar = ''

    # keep track of the current position in the sequence
    current_position = 1

    # keep track of the current position in the cigar string
    current_cigar_position = 1
    
    # loop through the mutations, insertions, and deletions to create the cigar string. If there's no mutation, insertion, or deletion at a position, then the position is a match.
    while current_position <= len_sequence:
        # if there's a mutation at the current position, then add a mutation to the cigar string
        if current_position-1 in mutations:
            cigar += '{}M1X'.format(current_cigar_position)
            current_cigar_position = 0
            current_position += 1
        # if there's an insertion at the current position, then add an insertion to the cigar string
        elif current_position-1 in insertions:
            cigar += '{}M1I'.format(current_cigar_position)
            current_cigar_position = 0
            current_position += 1
        # if there's a deletion at the current position, then add a deletion to the cigar string
        elif current_position-1 in deletions:
            cigar += '{}M1D'.format(current_cigar_position)
            current_cigar_position = 0
            current_position += 1
        # if there's no mutation, insertion, or deletion at the current position, then add a match to the cigar string
        else:
            current_cigar_position += 1
            current_position += 1
    
    # add the last match to the cigar string
    cigar += '{}M'.format(current_cigar_position)

    return cigar


def generate_sam_files(folder, sample_profile):
    """Generate sam files for the given parameters. 
    The sam files are written from the input parameters, not from an external software. They do not use fastq files.

    Arguments:
        folder {str} -- where to write the sam files
        sample_profile {dict} -- dictionary of parameters for each construct
            - constructs {list} -- list of construct names [list]
                - number_of_reads {list} -- number of reads to generate for each construct [list]
                - mutations {list} -- number of mutations to introduce in each read for each construct [list(list)]
                - reads {list} -- sequence to use for each construct [list]
                - insertions {list} -- number of insertions to introduce in each read for each construct [list(list)]
                - deletions {list} -- number of deletions to introduce in each read for each construct [list(list)]
                - barcode_start {int} -- where to start the barcode in the read [int]
                - barcodes {list} -- list of barcodes to use for each construct [list]            
    """

    for c, v in sample_profile.items():
        with open(os.path.join(folder, '{}.sam'.format(c)), 'w') as f:
            print_sam_header(f, c, len(v['reads']))
            for i in range(v['number_of_reads']):
                cigar = make_cigar(len( v['reads'][i]), v['mutations'][i], v['insertions'][i], v['deletions'][i])
                print_sam_lines(f, '{}:{}'.format(c, str(i).zfill(len(str(v['number_of_reads'])))), v['reads'][i], c, cigar)


def check_sample_profile(reads, mutations, insertions, deletions):
    """ Check that the sample profile is valid.
    Make sure that:
    - the mutations, insertions, and deletions are sorted
    - the mutations, insertions, and deletions are not consecutive
    - the mutations, insertions, and deletions are not at the same position
    - the mutations, insertions, and deletions are within the sequence length    
    """
    for idx, s in enumerate(reads):
        [mutations[idx][i].sort() for i in range(len(mutations[idx]))]
        [insertions[idx][i].sort() for i in range(len(insertions[idx]))]
        [deletions[idx][i].sort() for i in range(len(deletions[idx]))]
        for j in range(len(mutations[idx])):
            for i in range(len(mutations[idx][j])):
                assert mutations[idx][j][i] < len(s[j]), 'Mutation is out of sequence length'
            for i in range(len(insertions[idx][j])):
                assert insertions[idx][j][i] < len(s[j]), 'Insertion is out of sequence length'
            for i in range(len(deletions[idx][j])):
                assert deletions[idx][j][i] < len(s[j]), 'Deletion is out of sequence length'
        

def make_sample_profile(constructs, reads, number_of_reads, mutations, insertions, deletions, barcodes=None, barcode_start=None, sections=None, section_start=None, section_end=None):
    """Make a dictionary of the sample profile.
    The dictionary is used to make sure that the sample profile is valid. 
    The lists are grouped by construct, and the number of reads is grouped by construct.
    """
    check_sample_profile(reads, mutations, insertions, deletions)
    if barcodes is not None:
        assert len(constructs) == len(barcodes)
        for b in barcodes:
            assert len(b) == len(barcodes[0])
    
    sample_profile = {}
               
    for idx, c in enumerate(constructs):
        sample_profile[c] = {}
        sample_profile[c]['reads'] = reads[idx]            
        if barcodes is not None:
            for j in range(len(sample_profile[c]['reads'])): 
                sample_profile[c]['reads'][j] = sample_profile[c]['reads'][j][:barcode_start] + barcodes[idx] + sample_profile[c]['reads'][j][barcode_start+len(barcodes[idx]):]
        sample_profile[c]['reference'] = reads[idx][0] 
        sample_profile[c]['number_of_reads'] = number_of_reads[idx]
        sample_profile[c]['mutations'] = mutations[idx]
        sample_profile[c]['insertions'] = insertions[idx]
        sample_profile[c]['deletions'] = deletions[idx]
        if sections is not None:
            sample_profile[c]['sections'] = sections[idx]
            sample_profile[c]['section_start'] = section_start[idx]
            sample_profile[c]['section_end'] = section_end[idx]
            # assert that the sections are within the sequence length
            for j in range(len(sample_profile[c]['sections'])):
                assert sample_profile[c]['section_start'][j] < len(sample_profile[c]['reads'][j]), 'Section start is out of sequence length'
                assert sample_profile[c]['section_end'][j] < len(sample_profile[c]['reads'][j]), 'Section end is out of sequence length'
        if barcodes is not None:
            sample_profile[c]['barcodes'] = barcodes[idx]
            sample_profile[c]['barcode_start'] = barcode_start
        for i in range(sample_profile[c]['number_of_reads']):
            sequence = sample_profile[c]['reads'][i]
            for j in sample_profile[c]['mutations'][i]:
                sequence = sequence[:j] + next_base(sequence[j]) + sequence[j+1:]
            for j in sample_profile[c]['insertions'][i]:
                j += len([k for k in sample_profile[c]['insertions'][i] if k < j])
                sequence = sequence[:j] + next_base(sequence[j]) + sequence[j:]
            for j in sample_profile[c]['deletions'][i]:
                j += len([k for k in sample_profile[c]['insertions'][i] if k < j])
                j -= len([k for k in sample_profile[c]['deletions'][i] if k < j])
                sequence = sequence[:j] + sequence[j+1:]
            sample_profile[c]['reads'][i] = sequence
    return sample_profile






def generate_bitvector_files(folder, sample_profile, library):
    """Generate a bitvector file for the sample profile. 
    
    The bitvector is an orc file generated using pandas. A bitvector is a binary file that contains the reads for a certain section of a certain construct. The sections of the construct are defined by the library. 
    For each section of each construct, we are interested in the residues inside the [section_start, section_end] 0-indexed interval. 
    The bitvector columns are the residues, and the rows are the reads. The columns are under the format '0-indexed position'+'residue'. Ex: '0A' means that the first residue is an A.

    Parameters
    ----------
    folder : str
        The folder to write the bitvector file to.
    sample_profile : dict
        The sample dictionary.
        Attributes are:
        - constructs
            - reads
            - number_of_reads
            - mutations
            - insertions
            - deletions
    library : str
        The library to use. Columns are 'construct', 'section', 'section_start','section_end'.

    """

    library = pd.read_csv(library)

    # get through the constructs
    for construct in sample_profile:

    # create a folder under the construct's name
        construct_folder = os.path.join(folder, construct)
        if not os.path.exists(construct_folder):
            os.makedirs(construct_folder)
    # get through the sections for that construct
        for section in library[library['construct'] == construct]['section']:
            section_start = library[(library['construct'] == construct) & (library['section'] == section)]['section_start'].values[0]
            section_end = library[(library['construct'] == construct) & (library['section'] == section)]['section_end'].values[0]
            sequence = sample_profile[construct]['reference'][section_start:section_end]
            print(sample_profile[construct]['reference'], section_start, section_end, sequence)
            assert len(sequence) == section_end - section_start, 'Section length is not equal to the length of the sequence'
            columns = [str(i)+sequence[i] for i in range(section_end - section_start)]

            # create a dataframe with the columns
            df = pd.DataFrame(columns=columns)
            for i in range(sample_profile[construct]['number_of_reads']):
                bv = [1]*len(columns)
                for j in range(len(columns)):
                    # if the residue has a deletion, use update_bv_byte
                    if j+section_start in sample_profile[construct]['deletions'][i]:
                        bv[j] = update_bv_byte(bv[j], 'deletion')
                    if j+section_start in sample_profile[construct]['insertions'][i]:
                        bv[j] = update_bv_byte(bv[j], 'insertion_3')
                        bv[j+1] = update_bv_byte(bv[j+1], 'insertion_5')
                    if j+section_start in sample_profile[construct]['mutations'][i]:
                        base = next_base(sample_profile[construct]['reference'][j+section_start])
                        bv[j] = update_bv_byte(bv[j], 'substitution_'+base)
                df.loc[i] = ''.join([str(b) for b in bv])
            df.to_orc(os.path.join(construct_folder, section+'.orc'), index=False)

def update_bv_byte(byte, position):
    """Add a bit to a byte and return the new byte for bitvectoring.
    00000001: match
    00000010: deletion
    00000100: ≥1 base is inserted immediately 3' of this position
    00001000: ≥1 base is inserted immediately 5' of this position
    00010000: substitution to A
    00100000: substitution to C
    01000000: substitution to G
    10000000: substitution to T
    """
    if position == 'match':
        return byte + 2**0
    elif position == 'deletion':
        return byte + 2**1
    elif position == 'insertion_3':
        return byte + 2**2
    elif position == 'insertion_5':
        return byte + 2**3
    elif position == 'substitution_A':
        return byte + 2**4
    elif position == 'substitution_C':
        return byte + 2**5
    elif position == 'substitution_G':
        return byte + 2**6
    elif position == 'substitution_T':
        return byte + 2**7
    else:
        raise ValueError('Position not recognized')            

def generate_clustering_file(file, sample_profile):
    out = {}
    with open (file, 'w') as f:
        json.dump(out, f)

def generate_samples_csv_file(samples_csv_name):
    sample = samples_csv_name.split('/')[-2].split('.')[0]
    df = pd.DataFrame({
        'sample': [sample],
        'user': ['Napoleon'],
        'date': ['1769-08-15'],
        'exp_env': ['in_vitro'],
        'temperature_k': [310],
        'buffer': ['MgCl2'],
        'inc_time_tot_secs': [100]
    }, index=[0])
    df.to_csv(samples_csv_name, index=False)


def count_bases(positions, ss, se):
    out = np.zeros(se-ss, dtype=int)
    for p in positions:
        if p >= ss and p < se:
            out[p-ss] += 1
    return out

def count_mut_indel(mut_indel, ss, se):
    return  np.array([count_bases(mut_indel[p], ss, se) for p in range(len(mut_indel))]).sum(axis=0).tolist()

def count_mut_mod(ref, muts, base, ss, se):
    return np.array([count_bases([b for b in muts[p] if ref[b] == next_base(base)], ss, se) for p in range(len(muts))]).sum(axis=0).tolist()

from dreem.aggregate.rnastructure import RNAstructure

def generate_output_files(file, sample_profile, library, samples, clusters = None, rnastructure_config = None):
    if clusters is None:
        library = pd.read_csv(library)
        samples = pd.read_csv(samples)
        out = samples.to_dict('records')[0]
        out['construct'] = {}
        for construct, v in sample_profile.items():
            out['construct'][construct] = {}
            out['construct'][construct]['num_reads'] = sample_profile[construct]['number_of_reads']
            out['construct'][construct]['num_aligned'] = sample_profile[construct]['number_of_reads']
            out['construct'][construct]['barcode'] = v['barcodes']
            out['construct'][construct]['barcode_start'] = v['barcode_start']
            out['construct'][construct]['some_random_attribute'] = library['some_random_attribute'].values[0]
            out['construct'][construct]['sequence'] = v['reference']
            for s, ss, se in zip(v['sections'], v['section_start'], v['section_end']):
                # DREEM
                out['construct'][construct][s] = {}
                out['construct'][construct][s]['section_start'] = ss
                out['construct'][construct][s]['section_end'] = se
                out['construct'][construct][s]['sequence'] = v['reference'][ss:se]
                out['construct'][construct][s]['num_of_mutations'] = [len([a for a in v['mutations'][b] if (a>=ss) and (a<se)]) for b in range(len(v['mutations']))]
                out['construct'][construct][s]['mut_bases'] =  count_mut_indel(v['mutations'], ss, se)
                out['construct'][construct][s]['del_bases'] =  count_mut_indel(v['deletions'], ss, se)
                out['construct'][construct][s]['ins_bases'] =  count_mut_indel(v['insertions'], ss, se)
                out['construct'][construct][s]['cov_bases'] = [v['number_of_reads']]*(se-ss)
                out['construct'][construct][s]['info_bases'] = [v['number_of_reads']]*(se-ss)
                out['construct'][construct][s]['mut_rates'] = np.array( np.array(out['construct'][construct][s]['mut_bases'])/np.array(out['construct'][construct][s]['cov_bases'])).tolist()
                for base in ['A', 'C', 'G', 'T']:
                    out['construct'][construct][s]['mod_bases_{}'.format(base)] = count_mut_mod(v['reference'], v['mutations'], base, ss, se)
                out['construct'][construct][s]['worst_cov_bases'] = v['number_of_reads']
                out['construct'][construct][s]['skips_short_reads'] = 0
                out['construct'][construct][s]['skips_too_many_muts'] = 0
                out['construct'][construct][s]['skips_low_mapq'] = 0

                # RNAstructure
                if rnastructure_config is not None:
                    rnastructure_config['temp_folder'] = os.path.join(os.path.dirname(file), 'temp')
                    os.makedirs(rnastructure_config['temp_folder'], exist_ok=True)
                    rna = RNAstructure(rnastructure_config)
                    mp = pd.Series({
                        'sequence': out['construct'][construct][s]['sequence'],
                        'construct': construct,
                        'section': s,
                        'info_bases': out['construct'][construct][s]['info_bases'],
                        'mut_bases': out['construct'][construct][s]['mut_bases'],
                        'temperature_k': out['temperature_k']
                    })
                    rna_pred = rna.run(mp, out['sample'])
                    for k, va in rna_pred.items():
                        out['construct'][construct][s][k] = va
                        
                    # Poisson
                    poisson_pred = poisson.compute_conf_interval(mp.info_bases, mp.mut_bases)
                    for k, va in poisson_pred.items():
                        out['construct'][construct][s][k] = va
                    os.system('rm -fr {}'.format(rnastructure_config['temp_folder']))

    else:       
        raise NotImplementedError('Clustering not implemented yet')

    with open(file, 'w') as f:
        json.dump(out, f, indent=4)


def generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name, rnastructure_config=None):
    """Generate the files for a sample profile.

    Parameters
    ----------
    sample_profile : dict
        The sample dictionary.
        Attributes are:
        - constructs
            - reads
            - number_of_reads
            - mutations
            - insertions
            - deletions
    module : str
        The module name.
    inputs : list
        The inputs files.
    outputs : list
        The outputs files.
    test_files_dir : str
        The test files directory.
    sample_name : str
        The sample name.
    rnastructure_config : str
        The RNAstructure config file.
    """
    # make input and output folders
    input_folder = os.path.join(test_files_dir, 'input', module, sample_name)
    output_folder = os.path.join(test_files_dir, 'predicted_output', module, sample_name)
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if ('bitvector' in inputs or 'bitvector' in outputs) and 'library' not in inputs:
        inputs.insert(0,'library')
    if 'library' in inputs:
        if inputs.index('library') > 0:
            inputs.insert(0, inputs.pop(inputs.index('library')))
    for inpt in inputs:
        generate_file_factory(input_folder, inpt, sample_profile)
    for output in outputs:
        generate_file_factory(input_folder, output, sample_profile,output_folder, rnastructure_config)

def generate_file_factory(input_path, file_type, sample_profile, output_path=None, rnastructure_config=None):
    if output_path is None:
        path = input_path
    else:
        path = output_path
    library = os.path.join(input_path, 'library.csv')
    samples = os.path.join(input_path, 'samples.csv')
    if file_type == 'library':
        generate_library_file(library, sample_profile)
    elif file_type == 'fasta':
        generate_fasta_file(os.path.join(path, 'reference.fasta'), sample_profile)
    elif file_type == 'fastq':
        generate_fastq_files(path, sample_profile)
    elif file_type == 'bitvector':
        generate_bitvector_files(path, sample_profile, library)
    elif file_type == 'clustering':
        generate_clustering_file(os.path.join(path, 'clustering.json'), sample_profile)
    elif file_type == 'samples_csv':
        generate_samples_csv_file(samples)
    elif file_type == 'demultiplexed_fastq':
        generate_demultiplexed_fastq_files(path, sample_profile)
    elif file_type == 'sam':
        generate_sam_files(path, sample_profile)
    elif file_type == 'output':
        generate_output_files(path+'.json', sample_profile, library, samples, None, rnastructure_config) # TODO: add clustering
    else:
        raise ValueError('File type not recognized: "{}"'.format(file_type))

def assert_files_exist(sample_profile, module, files, in_out_pred_dir, sample_name):
    """Assert that the files for a sample profile exist.

    Parameters
    ----------
    sample_profile : dict
        The sample dictionary.
        Attributes are:
        - constructs
            - reads
            - number_of_reads
            - mutations
            - insertions
            - deletions
    module : str
        The module name.
    files : list
        The files.
    in_out_pred_dir : str
        The input, output or predicted output directory.
    sample_name : str
        The sample name.
    """
    # make input and output folders
    folder = os.path.join(in_out_pred_dir, module, sample_name)
    assert os.path.exists(folder), 'Folder of {} files does not exist: {}'.format(folder.split('/')[-3], folder)
    
    for file in files:
        assert_file_factory(folder, file, sample_profile, sample_name)
        
def assert_file_factory(path, file_type, sample_profile, sample_name):
    if file_type == 'library':
        assert os.path.exists(os.path.join(path, 'library.csv')), 'Library file does not exist: {}'.format(os.path.join(path, 'library.csv'))
    elif file_type == 'fasta':
        assert os.path.exists(os.path.join(path, 'reference.fasta')), 'Fasta file does not exist: {}'.format(os.path.join(path, 'reference.fasta'))
    elif file_type == 'fastq':
            for read in [1,2]:
                assert os.path.exists(os.path.join(path, '{}_R{}.fastq'.format(sample_name, read))), 'Fastq file does not exist: {}'.format(os.path.join(path, '{}_R{}.fastq'.format(sample_name, read)))
    elif file_type == 'bitvector':
        for construct in sample_profile:
            for section in sample_profile[construct]['sections']:
                assert os.path.exists(os.path.join(path, '{}/{}.orc'.format(construct, section))), 'Bitvector file does not exist: {}'.format(os.path.join(path, '{}/{}.orc'.format(construct, section)))
    elif file_type == 'clustering':
        assert os.path.exists(os.path.join(path, 'clustering.json')), 'Clustering file does not exist: {}'.format(os.path.join(path, 'clustering.json'))
    elif file_type == 'samples_csv':
        assert os.path.exists(os.path.join(path, 'samples.csv')), 'Samples csv file does not exist: {}'.format(os.path.join(path, 'samples.csv'))
    elif file_type == 'demultiplexed_fastq':
        for construct in sample_profile:
            for read in [1,2]:
                assert os.path.exists(os.path.join(path, '{}_R{}.fastq'.format(construct, read))), 'Demultiplexed fastq file does not exist: {}'.format(os.path.join(path, '{}_R{}.fastq'.format(construct, read)))
    elif file_type == 'sam':
        for construct in sample_profile:
            assert os.path.exists(os.path.join(path,'{}.sam'.format(construct))), 'Sam file does not exist: {}'.format(os.path.join(path, '{}.sam'.format(construct)))

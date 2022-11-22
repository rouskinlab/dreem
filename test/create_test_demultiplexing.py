
import pandas as pd
import random
import os


def generate_demultiplexing_test_samples(directory= 'test/demultiplexing', num_reads=1000, num_constructs=10, barcode_length=8, min_barcode_distance=3, barcode_start=10, length_sequence=20, sample_name='test_reference'):
    try:
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

        def print_fastq_line(f, id, seq, qual):
            f.write('@{}\n{}\n+\n{}\n'.format(id, seq, qual))    

        def print_fasta_line(f, id, seq):
            f.write('>{}\n{}\n'.format(id, seq))



        # create two fake fastq files
        barcode_start, barcode_end = barcode_start, barcode_start + barcode_length
        barcodes = generate_barcodes(barcode_length, num_constructs, min_barcode_distance)
        sequence = ''.join(random.choice('ATCG') for _ in range(length_sequence))
        constructs = barcodes

        print('Generating {} reads with {} constructs'.format(num_reads, num_constructs))
        print('creating fasta file')
        # fasta

        df = pd.DataFrame({'construct': constructs, 
                            'barcode_start': [barcode_start for i in range(num_constructs)],
                            'barcode_end': [barcode_end for i in range(num_constructs)],  
                            'barcode_sequence': barcodes}).to_csv(os.path.join(directory, 'library.csv'), index=False)

        print('creating fastq1')
        # fastq1
        with open(directory + '/{}_R1.fastq'.format(sample_name), 'w') as f:
            for barcode in barcodes:
                this_sequence = sequence[:barcode_start] + barcode + sequence[barcode_end:]
                for r in range(num_reads):
                    print_fastq_line(f, 'READ:{}:{}'.format(barcode, r), this_sequence, 'F'*length_sequence)
            f.close()

        print('creating demultiplexed fastq1')
        # demultiplexed fastq1
        for construct, barcode in zip(constructs, barcodes):
            with open(directory + '/{}/{}_R1.fastq'.format(sample_name,construct), 'w') as f:
                this_sequence = sequence[:barcode_start] + barcode + sequence[barcode_end:]
                for r in range(num_reads):
                    print_fastq_line(f, 'READ:{}:{}'.format(barcode, r), this_sequence, 'F'*length_sequence)
                f.close()

        # fastq2
        print('creating fastq2')
        barcode_start, barcode_end = length_sequence - barcode_end, length_sequence - barcode_start
        barcodes = [invert_sequence(barcode) for barcode in barcodes]
        sequence = invert_sequence(sequence)

        with open(directory + '/{}_R2.fastq'.format(sample_name), 'w') as f:
            for barcode in barcodes:
                this_sequence = sequence[:barcode_start] + barcode + sequence[barcode_end:]
                for r in range(num_reads):
                    print_fastq_line(f, 'READ:{}:{}'.format(barcode, r), this_sequence, 'F'*length_sequence)
            f.close()

        print('creating demultiplexed fastq2')
        # demultiplexed fastq2
        for construct, barcode in zip(constructs, barcodes):
            with open(directory + '/{}/{}_R2.fastq'.format(sample_name,construct), 'w') as f:
                this_sequence = sequence[:barcode_start] + barcode + sequence[barcode_end:]
                for r in range(num_reads):
                    print_fastq_line(f, 'READ:{}:{}'.format(barcode, r), this_sequence, 'F'*length_sequence)
                f.close()
        return True
    except Exception as e:
        print(e)
        return False
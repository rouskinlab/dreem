import yaml, sys, os
from dreem import util
from dreem import demultiplexing, alignment, vectoring, clustering, aggregation, drawer

def run(**args):
    """Run DREEM.

     Args:
         args (_type_): _description_
    """
    def verbose_print(x):
        print(x) if args['verbose'] else None

    # make output and temp folders
    for folder in ['output', 'temp']:
        os.makedirs(os.path.join(args['out_dir'], folder), exist_ok=True)
        
    # sort fast pairs
    args['fastq1'], args['fastq2'], args['samples'] = util.sort_fastq_pairs(args['fastq1'], args['fastq2'])
    
    # Run DREEM
    verbose_print("""

    ========================================

                RUNNING   DREEM

    ========================================

    """)
    if args['demultiplexing']:
        verbose_print('\ndemultiplexing \n------------------')
        
        # Run demultiplexing
        demultiplexing.run(**{**args, **{'out_dir': os.path.join(args['out_dir'], 'output', 'demultiplexing')}})
        verbose_print('demultiplexing done')
    
    ## Alignment: 
    verbose_print('\nalignment \n----------------')
    fastq1, fastq2, samples = [], [], []
    if args['demultiplexing']:
        for sample in os.listdir(os.path.join(args['out_dir'], 'output', 'demultiplexing')):
            path = os.path.join(os.path.join(args['out_dir'], 'output', 'demultiplexing'), sample)
            fastq1 = fastq1 + [os.path.join(path, f) for f in os.listdir(path) if f.endswith('_R1.fastq')]
            fastq2 = fastq2 + [os.path.join(path, f) for f in os.listdir(path) if f.endswith('_R2.fastq')]
            samples.append(sample)
    else:
        # samples is the name of the sample, which is the name of the fastq file without the extension
        fastq1, fastq2, samples = args['fastq1'], args['fastq2'],  args['samples']

    for f1, f2, sample in zip(fastq1, fastq2, samples):
        verbose_print('Aligning this fastq pair: ', '\n   ',f1, '\n   ',f2)
        alignment.run(**{**{k:v for k,v in args.items()},\
                      **{'out_dir': os.path.join(args['out_dir'], 'output', 'alignment'),
                         'sample': sample, 
                         'fastq1': f1, 
                         'fastq2': f2}})

    ## Vectoring
    verbose_print('\nvectoring \n------------------')
    for sample in args['samples']:
        vectoring.run(**{**args, 
                    **{'out_dir': os.path.join(args['out_dir'], 'output', 'vectoring', sample),
                       'input_dir': os.path.join(args['out_dir'], 'output', 'alignment', sample),}})
    
    ## Clustering
    if args['clustering']:
        verbose_print('\nclustering \n------------------')
        for sample in args['samples']:
            clustering.run(**{**args, 'input_dir': os.path.join(args['out_dir'], 'output', 'vectoring', sample), 
                                        'out_dir': os.path.join(args['out_dir'], 'output', 'clustering', sample)})

    ## Aggregate
    verbose_print('\naggregating \n------------------')
    for sample in args['samples']:
        args['clustering'] = os.path.join(args['out_dir'], 'output', 'clustering', sample+'.json') if args['clustering'] else None
        aggregation.run(**{**args, 'input_dir': os.path.join(args['out_dir'], 'output', 'vectoring'), 'sample': sample})

    verbose_print('Done!')
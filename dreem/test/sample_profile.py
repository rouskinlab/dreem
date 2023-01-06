from dreem.test import files_generator 

number_of_constructs = 1
number_of_reads = [10]*number_of_constructs
mutations = [ [[]]+[[23]]+[[35]]+[[]]*4+[[37]]+[[32]]+[[33,36]] for n in range(number_of_constructs) ] # 0-based
insertions = [ [[]]*3+[[11]]+[[8, 23]]+[[]]*2+[[15]]+[[]]*2 for n in range(number_of_constructs) ] # 0-based
deletions = [ [[]]*5+[[2]]+[[4, 6]]+[[]]+[[8]]+[[]] for n in range(number_of_constructs) ] # 0-based
no_info = [ [[]]*2+[[2]]+[[4, 6]]+[[]]+[[3]]+[[]]*5 for n in range(number_of_constructs) ] # 0-based

#insertions = [ [[]]*10 for n in range(number_of_constructs) ] # 0-based
#no_info = [ [[]]*10 for n in range(number_of_constructs) ] # 0-based
#deletions = [ [[]]*10 for n in range(number_of_constructs) ] # 0-based
#mutations = [ [[]]*10 for n in range(number_of_constructs) ] # 0-based

#number_of_reads = [2]*number_of_constructs
#mutations = [[[10,14], [12]]  for n in range(number_of_constructs) ] # 0-based
#insertions = [ [[], []] for n in range(number_of_constructs) ] # 0-based
#deletions = [ [[4], [ 4,7]]  for n in range(number_of_constructs) ] # 0-based
no_info = [[[]]*number_of_reads[n] for n in range(number_of_constructs) ] # 0-based

length = [50, 150]
sequences = [[files_generator.create_sequence(length[k])]*number_of_reads[k] for k in range(number_of_constructs)]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 30
len_barcode = 10
barcodes = files_generator.generate_barcodes(len_barcode, number_of_constructs, 3)
sections_start = [[0, 25],[0, 25, 50, 75]] # 0-based
sections_end = [[25, 50],[25, 50, 75, 99]] # 0-based
sections = [['{}-{}'.format(ss+1, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)] # 0-based
sample_profile = files_generator.make_sample_profile(constructs, sequences, number_of_reads, mutations, insertions, deletions, no_info, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)

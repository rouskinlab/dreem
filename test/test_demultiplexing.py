import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import pytest
import dreem
import pandas as pd 

test_files = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_files')
test_files_sample = os.path.join(test_files, 'my_sample_to_demultiplex')
demultiplexed_fastq_dir = os.path.join(test_files_sample,'my_sample_to_demultiplex')
library = os.path.join(test_files_sample, 'library.csv')
fasta =  os.path.join(test_files_sample, 'sample.fasta')
fastq1 = os.path.join(test_files_sample, 'sample_R1.fastq')
fastq2 = os.path.join(test_files_sample, 'sample_R2.fastq')
out_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_output','demultiplexing')
temp_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_temp','demultiplexing')

def test_demultiplexing_python():
    dreem.demultiplex.run(
        library=library,
        fastq1=(fastq1,),
        fastq2=(fastq2,),
        clipped=1,
        index_tolerance=5,
        mismatch_tolerence=1,
        out_dir=os.path.join(out_dir, 'python')+"/",
        temp_dir=os.path.join(temp_dir, 'python')+"/",
        fasta=fasta,
    )
    """dreem.demultiplex.run(
        library="/Users/scottgrote/Documents/absolutely_final_repo/dreem/test_files/my_sample_to_demultiplex/library2.csv",
        fastq1=(fastq1,),
        fastq2=(fastq2,),
        clipped=1,
        index_tolerance=5,
        mismatch_tolerence=1,
        out_dir=os.path.join(out_dir, 'python')+"/",
        temp_dir=os.path.join(temp_dir, 'python')+"/",
        fasta=fasta,
    )"""
print("dumb")
def test_demultiplexing_cli():
    os.system(f"dreem demultiplex --library {library} --fasta {fasta} --fastq1 {fastq1} --fastq2 {fastq2} --out-dir {os.path.join(out_dir, 'cli')}/ --temp-dir {os.path.join(temp_dir, 'cli')}/ --index-tolerance {5} --mismatch-tolerence {1} --clipped {1}")
    
@pytest.mark.parametrize('sample', ['cli', 'python'])
def test_results(sample):
    #directory = os.path.join(out_dir, sample, 'my_sample_to_demultiplex')
    df=pd.read_csv(os.path.join(out_dir,"python","demultiplex_info.csv"))
    for i in df.index:
        if(df.at[i,"fq1_construct_name_"]=="3042-O-flank_1=hp1-DB"):
            assert df.at[i,"fq1_final_count_"]== 10
            assert df.at[i,"fq1_lost_to_multigrep_"]==0
            assert df.at[i,"fq1_init_barcode"]==6
            assert df.at[i,"fq1_barcode_clipped_1"]==1
            assert df.at[i,"fq1_init_rev_barcode"]==3
            assert df.at[i,"fq1_unfiltered"]==10
            assert df.at[i,"fq1_sec_filter"]==7
            assert df.at[i,"fq1_rev_sec_filter"]==3
            assert df.at[i,"fq1_reads_lost_to_filter"]==0
            assert df.at[i,"fq1_pre_union"]==10
            

            assert df.at[i,"fq2_final_count_"]== 10
            assert df.at[i,"fq2_lost_to_multigrep_"]==0
            assert df.at[i,"fq2_init_barcode"]==2
            assert df.at[i,"fq2_barcode_clipped_1"]==1
            assert df.at[i,"fq2_init_rev_barcode"]==7
            assert df.at[i,"fq2_unfiltered"]==10
            assert df.at[i,"fq2_sec_filter"]==3
            assert df.at[i,"fq2_rev_sec_filter"]==7
            assert df.at[i,"fq2_reads_lost_to_filter"]==0
            assert df.at[i,"fq2_pre_union"]==10

        if(df.at[i,"fq1_construct_name_"]=="3043-CC-flank_1=hp1-DB"):
            assert df.at[i,"fq1_final_count_"]== 10
            assert df.at[i,"fq1_lost_to_multigrep_"]==-1
            assert df.at[i,"fq1_init_barcode"]==5
            assert df.at[i,"fq1_barcode_clipped_1"]==0
            assert df.at[i,"fq1_init_rev_barcode"]==5
            assert df.at[i,"fq1_unfiltered"]==10
            assert df.at[i,"fq1_sec_filter"]==5
            assert df.at[i,"fq1_rev_sec_filter"]==4
            assert df.at[i,"fq1_reads_lost_to_filter"]==1
            assert df.at[i,"fq1_pre_union"]==9
            

            assert df.at[i,"fq2_final_count_"]== 10
            assert df.at[i,"fq2_lost_to_multigrep_"]==0
            assert df.at[i,"fq2_init_barcode"]==5
            assert df.at[i,"fq2_barcode_clipped_1"]==0
            assert df.at[i,"fq2_init_rev_barcode"]==5
            assert df.at[i,"fq2_unfiltered"]==10
            assert df.at[i,"fq2_sec_filter"]==5
            assert df.at[i,"fq2_rev_sec_filter"]==5
            assert df.at[i,"fq2_reads_lost_to_filter"]==0
            assert df.at[i,"fq2_pre_union"]==10

        if(df.at[i,"fq1_construct_name_"]=="3414-CC-flank_1=lp6-DB"):
            assert df.at[i,"fq1_final_count_"]== 10
            assert df.at[i,"fq1_lost_to_multigrep_"]==0
            assert df.at[i,"fq1_init_barcode"]==7
            assert df.at[i,"fq1_barcode_clipped_1"]==0
            assert df.at[i,"fq1_init_rev_barcode"]==3
            assert df.at[i,"fq1_unfiltered"]==10
            assert df.at[i,"fq1_sec_filter"]==7
            assert df.at[i,"fq1_rev_sec_filter"]==3
            assert df.at[i,"fq1_reads_lost_to_filter"]==0
            assert df.at[i,"fq1_pre_union"]==10
            

            assert df.at[i,"fq2_final_count_"]== 10
            assert df.at[i,"fq2_lost_to_multigrep_"]==0
            assert df.at[i,"fq2_init_barcode"]==3
            assert df.at[i,"fq2_barcode_clipped_1"]==0
            assert df.at[i,"fq2_init_rev_barcode"]==7
            assert df.at[i,"fq2_unfiltered"]==10
            assert df.at[i,"fq2_sec_filter"]==3
            assert df.at[i,"fq2_rev_sec_filter"]==7
            assert df.at[i,"fq2_reads_lost_to_filter"]==0
            assert df.at[i,"fq2_pre_union"]==10

            
    assert df.at[0,"percent_reads_unused"]==100.0
    assert df.at[0,"orginal_fq_total_read_count"]==30
    #print("passed")
    
if __name__ == '__main__':
    test_demultiplexing_python()
    #test_results('python')
    test_demultiplexing_cli()
    #test_results('cli')
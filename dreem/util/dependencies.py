import os
import subprocess

def run_cmd(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output.decode('utf-8').strip()

#TODO: #21 make this dependent on the env.yml file

def check_fastqc_exists(version = 'v0.12.1'):
    fastqc_installed_version = run_cmd(['fastqc', '--version'])
    assert fastqc_installed_version.startswith('FastQC'), 'fastqc is not installed'
    assert fastqc_installed_version == 'FastQC '+version, 'fastqc version is not correct: {} != {}'.format(fastqc_installed_version, version)
    
def check_cutadapt_exists(version = '4.1'):
    cutadapt_installed_version = run_cmd(['cutadapt', '--version'])
    assert not 'cutadapt: not found' in cutadapt_installed_version, 'cutadapt is not installed'
    assert cutadapt_installed_version == version, 'cutadapt version is not correct: {} != {}'.format(cutadapt_installed_version, version)
    
def check_bowtie2_exists(version = '2.4.5'):
    bowtie2_installed_version = run_cmd(['bowtie2', '--version']).split('\n')[0]
    assert bowtie2_installed_version.split(' ')[1] == 'version', 'bowtie2 is not installed'
    assert bowtie2_installed_version.endswith('bowtie2-align-s version '+version), 'bowtie2 version is not correct: {} != {}'.format(bowtie2_installed_version, version)
    
def check_samtools_exists(version = '1.17'):
    samtools_installed_version = run_cmd(['samtools', '--version'])
    assert samtools_installed_version.startswith('samtools'), 'samtools is not installed'
    assert samtools_installed_version.split('\n')[0].split(' ')[1] == version, 'samtools version is not correct: {} != {}'.format(samtools_installed_version, version)
    
def check_rnastructure_exists(path, version = '6.4'):
    if path and path != None:
        cmd = [os.path.join(path, 'Fold'), '--version']
    else:
        cmd = ['Fold', '--version']
    rnastructure_installed_version = run_cmd(cmd)
    assert rnastructure_installed_version.startswith('Fold: Version '), 'RNAstructure is not found. If the "Fold" command runs in your terminal, try setting the "rnastructure_path" parameter in the config file to the path of the RNAstructure installation (ex: rnastructure_path="/Users/username/RNAstructure/exe").'
    assert rnastructure_installed_version.split(' ')[2] == version, 'RNAstructure version is not correct: {} != {}'.format(rnastructure_installed_version, version)

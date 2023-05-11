from dreem.cluster import run
# from dreem.vector.profile import VectorReader

if __name__ == '__main__':
    result = run(
        mv_file=['/Users/mfa/prj_tmp/dreem-test/output/vectoring/3kb-ASOs-ctrl2/SARS2/13369-13597/report.json'],
        out_dir='/Users/mfa/prj_tmp/dreem-test/output/',
        num_runs=2,
        max_clusters=2,
        signal_thresh=0.005,
        exclude_gu=False,
        include_del=False,
        min_iter=5,
        max_iter=5,
        convergence_cutoff=1e-1)
    
    print(result)

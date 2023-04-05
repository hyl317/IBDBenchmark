# run ancIBD given an input vcf file

import argparse
#import sys
#sys.path.append('/mnt/archgen/users/yilei/tools/ancIBD/package')
from ancIBD.run import hapBLOCK_chroms
from ancIBD.IO.prepare_h5 import vcf_to_1240K_hdf
from pathlib import Path
import h5py
import itertools
import subprocess
import os

ch = 3
covs = ['cov3', 'cov2', 'cov1', 'cov3over4']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ancIBD.')
    parser.add_argument('-r', action="store", dest="r", type=int, required=True,
                        help="which replicate to run")
    args = parser.parse_args()

    r = args.r

    for cov1, cov2 in itertools.combinations(covs, 2):
        vcf1 = f'{cov1}/batch{r}/batch{r}.all.chr3.vcf.gz'
        vcf2 = f'{cov2}/batch{r}/batch{r}.all.chr3.vcf.gz'
        if not os.path.exists(f'{vcf1}.csi'):
            subprocess.run(f'bcftools index {vcf1}'.split())
        if not os.path.exists(f'{vcf2}.csi'):
            subprocess.run(f'bcftools index {vcf2}'.split())
        path2vcf = f'./mixedCoverage/batch{r}/{cov1}_{cov2}.all.chr3.vcf.gz'
        subprocess.run(f'bcftools merge -Oz -o {path2vcf} --force-samples {vcf1} {vcf2}'.split())
        vcf_to_1240K_hdf(in_vcf_path = path2vcf,
                 path_vcf = f"./mixedCoverage/batch{r}/{cov1}_{cov2}.1240k.chr{ch}.vcf",
                 path_h5 = f"./mixedCoverage/batch{r}/{cov1}_{cov2}.chr{ch}.h5",
                 marker_path = f"/mnt/archgen/users/yilei/bin/ancIBD_data/filters/snps_bcftools_ch{ch}.csv",
                 map_path = "/mnt/archgen/users/yilei/bin/ancIBD_data/afs/v51.1_1240k.snp",
                 af_path = f"/mnt/archgen/users/yilei/bin/ancIBD_data/afs/v51.1_1240k_AF_ch{ch}.tsv",
                 col_sample_af = "",
                 buffer_size=20000, chunk_width=8, chunk_length=20000,
                 ch=ch)


        with h5py.File(f"mixedCoverage/batch{r}/{cov1}_{cov2}.chr{ch}.h5", 'r') as f:
            iids = f['samples'][:].astype('str')
        iids1 = [iid for iid in iids if not iid.startswith('2:')]
        iids2 = [iid for iid in iids if iid.startswith('2:')]
        assert(len(iids1) + len(iids2) == len(iids))
        run_iids_all = list(itertools.product(iids1, iids2))
        run_iids = []
        for pair in run_iids_all:
            id1, id2 = pair
            if id1.startswith('2:'):
                id1 = id1[2:]
            elif id2.startswith('2:'):
                id2 = id2[2:]
            if id1 != id2:
                run_iids.append(pair)
        print(f'{run_iids}')

        df_ibd = hapBLOCK_chroms(folder_in=f"./mixedCoverage/batch{r}/{cov1}_{cov2}.chr",
                             iids=iids, run_iids=run_iids,
                             ch=ch, folder_out=f'./mixedCoverage/batch{r}',
                             output=False, prefix_out=f"{cov1}_{cov2}", logfile=False,
                             l_model='hdf5', e_model='haploid_gl', h_model='FiveStateScaled', t_model='standard',
                             ibd_in=1, ibd_out=10, ibd_jump=400,
                             min_cm=6, cutoff_post=0.99, max_gap=0.0075,
                             processes=1)

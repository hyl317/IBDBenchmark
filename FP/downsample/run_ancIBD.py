# run ancIBD given an input vcf file

import argparse
#import sys
#sys.path.append('/mnt/archgen/users/yilei/tools/ancIBD/package')
from ancIBD.run import hapBLOCK_chroms
from ancIBD.IO.prepare_h5 import vcf_to_1240K_hdf
from pathlib import Path
import h5py

ch = 3

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ancIBD.')
    parser.add_argument('--vcf', action="store", dest="vcf", type=str, required=True,
                        help="path to the imputed vcf file")
    parser.add_argument('--prefix', action="store", dest="prefix", type=str, required=True,
                        help="prefix of output .h5 file")
    args = parser.parse_args()

    path2vcf = args.vcf
    prefix = args.prefix
    path = Path(path2vcf)
    pDir = path.parent.absolute()

    vcf_to_1240K_hdf(in_vcf_path = path2vcf,
                 path_vcf = f"{pDir}/{prefix}.1240k.chr{ch}.vcf",
                 path_h5 = f"{pDir}/{prefix}.chr{ch}.h5",
                 marker_path = f"/mnt/archgen/users/yilei/bin/ancIBD_data/filters/snps_bcftools_ch{ch}.csv",
                 map_path = "/mnt/archgen/users/yilei/bin/ancIBD_data/afs/v51.1_1240k.snp",
                 af_path = f"/mnt/archgen/users/yilei/bin/ancIBD_data/afs/v51.1_1240k_AF_ch{ch}.tsv",
                 col_sample_af = "",
                 buffer_size=20000, chunk_width=8, chunk_length=20000,
                 ch=ch)

    with h5py.File(f"{pDir}/{prefix}.chr{ch}.h5", 'r') as f:
        iids = f['samples'][:].astype('str')
    # iids = ['I4893', 'I4596', 'I1583', 'I2978', 'I5838', 'I1507', 'I2861', 'I2520', 'I3758', 'I5077', 'I0708', 'I5233', 'I3123']
    #  iids = ['I4893', 'I4596.chr3.bam', 'I1583.chr3.bam', 'I2978.chr3.bam',
    #     'I5838.chr3.bam', 'I1507.chr3.bam', 'I2861.chr3.bam',
    #     'I2520.chr3.bam', 'I3758.chr3.bam', 'I5077.chr3.bam',
    #     'I0708.chr3.bam', 'I5233.chr3.bam', 'I3123']
    df_ibd = hapBLOCK_chroms(folder_in=f"{pDir}/{prefix}.chr",
                             iids=iids, run_iids=[],
                             ch=ch, folder_out=f'{pDir}',
                             output=False, prefix_out='', logfile=False,
                             l_model='hdf5', e_model='haploid_gl', h_model='FiveStateScaled', t_model='standard',
                             ibd_in=1, ibd_out=10, ibd_jump=400,
                             min_cm=6, cutoff_post=0.99, max_gap=0.0075,
                             processes=1)

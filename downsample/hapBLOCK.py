import argparse
import sys
sys.path.append('/mnt/archgen/users/yilei/tools/hapBLOCK/python3/')
from run import hapBLOCK_all
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('--min_cm1', action="store", dest="min_cm1", type=float, required=True,
                        help="path to the vcf file")
    parser.add_argument('--post', action="store", dest="post", type=float, required=False, default=0.99,
                        help="posterior prob cutoff")
    parser.add_argument('--gap', action="store", dest="gap", type=float, required=False, default=0.0075,
                        help="maximum gap to merge")
    parser.add_argument('--out', action="store", dest="out", type=str, required=False, default='ancIBD',
                        help="output directory")
    parser.add_argument('--infolder', action="store", dest="infolder", type=str, required=False, default='processed_1240k',
                        help="input hdf5 directory")
    args = parser.parse_args()

    hapBLOCK_all(folder_in=f"./{args.infolder}/ch", iids=['I2105', 'I3950', 'I5273', 'I5279'],
                   chs=range(1,23), folder_out=f"./{args.out}", output=False, prefix_out="", logfile=False,
                   l_model="h5", IBD2=False, p_col="variants/AF_ALL", 
                   ibd_in=1, ibd_out=10, ibd_jump=500, min_cm1=args.min_cm1, min_cm2=2,
                   cutoff_post=args.post, max_gap=args.gap, snp_cm=0, maf=0.0, save=3)

    
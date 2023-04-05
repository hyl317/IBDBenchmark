import argparse
import pickle
import pandas as pd
import os

import sys
sys.path.append('/mnt/archgen/users/yilei/tools/ancIBD/notebook/simulate/python')
from multi_mosaic import multi_run_lengths
from modify_h5 import ModifyHDF5Genotypes
from ancIBD.run import hapBLOCK_chroms
import h5py
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate genotypes for ancIBD.')
    parser.add_argument('-l', action="store", dest="l", type=int, required=True,
                        help="length of IBD segments to simulate")
    parser.add_argument('-n', action="store", dest="n", type=int, required=False, default=250,
                        help="number of pairs to simulate")
    args = parser.parse_args()

    # first, simulate groundtruth genotypes
    n_blocks = 1 if args.l > 0 else 0
    # multi_run_lengths(base_path="./simGeno/", pop_list=["TSI", ],
    #                   lengths=[args.l], n_blocks=n_blocks, n=args.n, ch=3, chunk_length=0.0025, IBD2=False)
    
    # second, add errors    
    in_path = f"./simGeno/ch3_{args.l}cm/sim_ch3.h5" # The default dataset
    shuffle_rates = [0.326, 0.485, 0.517, 0.508, 0.384, 0.213]
    for err, shuffle_rate in zip(['cov2', 'cov1', 'cov3over4', 'cov1over2', 'cov1over4', 'cov1over10'], shuffle_rates):
        out_path = f"./simGeno/ch3_{args.l}cm/sim_{err}_ch3.h5"
        # empiricalError = pickle.load(open(f'/mnt/archgen/users/yilei/IBDsim/realMosaicSim_wgs/imputeAccuMap.{err}', 'rb'))
        # m = ModifyHDF5Genotypes(original_path=in_path, save_path=out_path)
        # m.downsample_gt(frac=1.0, ad=False, gp=True, 
        #             mult_alt=False, gt_type="int8", error=0.01,
        #             shuffle_cm=shuffle_rate, cty=0.99, simulated_error=empiricalError, compression="gzip")

        with h5py.File(out_path, 'r') as f:
            iids = f['samples'][:].astype('str')
        
        run_iids = []
        for i in np.arange(len(iids), step=2):
            run_iids.append([iids[i], iids[i+1]])
        print(f'run ancIBD on pairs: {run_iids}')

        hapBLOCK_chroms(folder_in=f"./simGeno/ch3_{args.l}cm/sim_{err}_ch", iids = iids, run_iids = run_iids, 
                ch=3, folder_out=f'./calledIBD_strict/{err}/ch3_{args.l}cm/', output=False, prefix_out='', logfile=False,
                             l_model='hdf5', e_model='haploid_gl', h_model='FiveStateScaled', t_model='standard',
                             ibd_in=1, ibd_out=25, ibd_jump=400, p_col='variants/AF_ALL',
                             min_cm=2, cutoff_post=0.995, max_gap=0.0075,
                             processes=1)
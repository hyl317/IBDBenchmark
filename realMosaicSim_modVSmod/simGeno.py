import argparse
import pickle
import pandas as pd
import os

import sys
sys.path.append('/mnt/archgen/users/yilei/tools/ancIBD/notebook/simulate/python')
from multi_mosaic import multi_run_lengths
from modify_h5 import ModifyHDF5Genotypes
from ancIBD.run import hapBLOCK_chroms
import numpy as np
import h5py
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate genotypes for ancIBD.')
    parser.add_argument('-l', action="store", dest="l", type=int, required=True,
                        help="length of IBD segments to simulate")
    parser.add_argument('-n', action="store", dest="n", type=int, required=False, default=250,
                        help="number of pairs to simulate")
    args = parser.parse_args()

    # first, simulate groundtruth genotypes    
    in_path = f"/mnt/archgen/users/yilei/IBDsim/realMosaicSim_wgs/simGeno/ch3_{args.l}cm/sim_ch3.h5" # The default dataset
    if not os.path.exists(in_path):
        multi_run_lengths(base_path="./simGeno/", pop_list=["TSI", ],
                      lengths=[args.l], n_blocks=1, n=args.n, ch=3, chunk_length=0.0025, IBD2=False)
    else:
        print(f'ground truth genotypes already exist, skip first step...', flush=True)
    
        
    # I set a particularly relaxed IBD1 length cutoff (4cM) here solely for the purpose of benchmarking
    # not recommended in practice!
    with h5py.File(in_path, 'r') as f:
        iids = f['samples'][:].astype('str')
        
    run_iids = []
    for i in np.arange(len(iids), step=2):
        run_iids.append([iids[i], iids[i+1]])
    print(f'run ancIBD on pairs: {run_iids}')

    hapBLOCK_chroms(folder_in=f"/mnt/archgen/users/yilei/IBDsim/realMosaicSim_wgs/simGeno/ch3_{args.l}cm/sim_ch", \
            iids = iids, run_iids = run_iids, 
            ch=3, folder_out=f'./calledIBD/ch3_{args.l}cm/', output=True, prefix_out='', logfile=False,
            l_model='hdf5', e_model='haploid_gl', h_model='FiveStateScaled', t_model='standard',
            ibd_in=1, ibd_out=10, ibd_jump=0.75, p_col='variants/AF_ALL',
            min_cm=2, cutoff_post=0.99, max_gap=0.0075, processes=1)
        

    
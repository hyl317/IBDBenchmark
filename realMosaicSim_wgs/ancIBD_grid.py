import argparse
import pickle
import pandas as pd
import os
from tqdm import tqdm
import sys
from ancIBD.run import hapBLOCK_chroms
import h5py
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate genotypes for ancIBD.')
    parser.add_argument('-l', action="store", dest="l", type=int, required=True,
                        help="length of IBD segments simulated")
    parser.add_argument('-n', action="store", dest="n", type=int, required=False, default=250,
                        help="number of pairs simulated")
    parser.add_argument('--cov', action="store", dest="cov", type=str, required=True,
                        help="a string representing which coverage to run")
    args = parser.parse_args()

    err = args.cov
    out_path = f"./simGeno/ch3_{args.l}cm/sim_{err}_ch3.h5"
    with h5py.File(out_path, 'r') as f:
        iids = f['samples'][:].astype('str')
        
    run_iids = []
    for i in np.arange(len(iids), step=2):
        run_iids.append([iids[i], iids[i+1]])
    print(f'run ancIBD on {len(run_iids)} pairs')

    for i, ibd_in in tqdm(enumerate([0.1, 1, 10])):
        for j, ibd_out in tqdm(enumerate([1, 10, 50])):
            for k, ibd_jump in tqdm(enumerate([100, 250, 500, 1000])):
                for m, gap in tqdm(enumerate([0.005, 0.0075, 0.01])):
                    for n, cutoff in tqdm(enumerate([0.99, 0.995, 0.999])):
                        hapBLOCK_chroms(folder_in=f"./simGeno/ch3_{args.l}cm/sim_{err}_ch", iids = iids, run_iids = run_iids, 
                            ch=3, folder_out=f'./grid/{err}/ch3_{args.l}cm/', output=False, \
                            prefix_out=f'ibd_in{i}_ibd_out{j}_ibd_jump{k}_gap{m}_cutoff_post{n}', logfile=False,
                             l_model='hdf5', e_model='haploid_gl', h_model='FiveStateScaled', t_model='standard',
                             ibd_in=ibd_in, ibd_out=ibd_out, ibd_jump=ibd_jump, p_col='variants/AF_ALL',
                             min_cm=6, cutoff_post=cutoff, max_gap=gap,
                             processes=1)
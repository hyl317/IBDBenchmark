from plot_utility import get_individual_idx, calc_FP_all
import numpy as np
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import time
import h5py
import itertools


path2HDF5 = '/mnt/archgen/users/yilei/IBDsim/groundtruth_IBD_glimpse_transversion/MAC0/hdf5'
iids = ['I2105', 'I3950', 'I5279', 'I5273']
bins = [(5,6), (6,8), (8,12), (12, np.inf)]
covs = ['cov2', 'cov1', 'cov3over4', 'cov1over2', 'cov1over4', 'cov1over10']
cov2str = {'cov2':'2x', 'cov1':'1x', 'cov3over4':'0.75x', 'cov1over2':'0.5x', 'cov1over4':'0.25x', 'cov1over10':'0.1x'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('-i', action="store", dest="i", type=int, required=True,
                        help="which length bin to plot")
    args = parser.parse_args()

    bin = bins[args.i - 1]
    print(f'bins to plot: {bin}')
    fig, axs = plt.subplots(6, 2, sharex=True, sharey=True, figsize=(6, 10.5), tight_layout=True)

    # pre-compute necessary quantity to speed up later computation
    print(f'start pre-computation...', flush=True)
    t1 = time.time()
    bp_vec = {}
    homo_index = defaultdict(lambda:{})
    oppoHomo_index = defaultdict(lambda:{})
    for ch in np.arange(1,23):
        f = h5py.File(f'{path2HDF5}/ch{ch}.h5', 'r')
        filterByAF = np.where(np.array(f['variants/RAF'][:,0])>0.2)[0]
        bp_vec[ch] = np.array(f['variants/POS'])[filterByAF]
        for id1, id2 in itertools.combinations(iids, 2):
            index1, index2 = get_individual_idx(f, id1), get_individual_idx(f, id2)
            gt1 = f['calldata/GT'][:, index1, :]
            gt2 = f['calldata/GT'][:, index2, :]
            gt1 = np.sum(gt1[filterByAF], axis=1)
            gt2 = np.sum(gt2[filterByAF], axis=1)
            assert(len(gt1) == len(gt2))
            assert(len(gt1) == len(bp_vec[ch]))
            homo_index[(min(id1, id2), max(id1, id2))][ch] = np.where(np.logical_and(np.logical_or(gt1==0, gt1==2), np.logical_or(gt2==0, gt2==2)))[0]
            oppoHomo_index[(min(id1, id2), max(id1, id2))][ch] = np.where(np.logical_or(np.logical_and(gt1==0, gt2==2), np.logical_and(gt1==2, gt2==0)))[0]
    print(f'pre-computation finished, takes {time.time()-t1}s', flush=True)




    #### plot WGS BAM results
    for i, cov in enumerate(covs):
        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/wgs/{cov}', bin[0], bin[1], bp_vec, homo_index, oppoHomo_index, threshold=0.5)
        axs[i,0].scatter(fraction_covered, oppoHomoRates, color='blue', alpha=0.3, label='ancIBD', s=5.0)

        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/wgs/{cov}', bin[0], bin[1], bp_vec, homo_index, oppoHomo_index, mode='IBIS_1240k', threshold=0.5)
        axs[i,0].scatter(fraction_covered, oppoHomoRates, color='orange', alpha=0.3, label='IBIS_1240k', s=5.0)

        if i == 0:
            axs[i,0].legend(loc='upper right', fontsize='large')
    
    #### plot 1240k BAM results
    for i, cov in enumerate(covs[:4]):
        i_ = i
        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/1240k/{cov}', bin[0], bin[1], bp_vec, homo_index, oppoHomo_index, threshold=0.5)
        axs[i_,1].scatter(fraction_covered, oppoHomoRates, color='blue', alpha=0.3, label='ancIBD', s=5.0)

        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/1240k/{cov}', bin[0], bin[1], bp_vec, homo_index, oppoHomo_index, mode='IBIS_1240k', threshold=0.5)
        axs[i_,1].scatter(fraction_covered, oppoHomoRates, color='orange', alpha=0.3, label='IBIS_1240k', s=5.0)
    axs[4,1].set_visible(False)
    axs[5,1].set_visible(False)
    #axs[6,1].set_visible(False)

    fig.text(0.2, 1.025, f'Length Bin: {bin[0]} - {bin[1]} cM', fontsize=22)
    fig.text(0.25, 1.01, 'wgs', ha='center', va='center', fontsize=18)
    fig.text(0.75, 1.01, '1240k', ha='center', va='center', fontsize=18)

    # add text to indicate coverages
    for i, cov in enumerate(covs[::-1]):
        fig.text(-0.005, 1/12 + i/6, cov2str[cov], ha='center', va='center', rotation='vertical', fontsize=16)

    fig.text(0.5, -0.01, 'Positive Predictive Value(PPV)', ha='center', va='center', fontsize=22)
    fig.text(-0.075, 0.5, 'Rate of Opposing Homozygotes', ha='center', va='center', rotation='vertical', fontsize=22)
    plt.savefig(f'./sidebyside/{bin[0]}_{bin[1]}.AF20.png', dpi=300, bbox_inches="tight")
    plt.savefig(f'./sidebyside/{bin[0]}_{bin[1]}.AF20.pdf', dpi=300, bbox_inches="tight")
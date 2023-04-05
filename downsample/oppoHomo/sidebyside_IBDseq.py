from plot_utility import calc_FP_all
import numpy as np
import matplotlib.pyplot as plt
import argparse

bins = [(5,6),(6,8), (8,12), (12, np.inf)]
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
    #### plot WGS BAM results
    for i, cov in enumerate(covs):
        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/wgs/{cov}', bin[0], bin[1], threshold=0.5)
        axs[i,0].scatter(fraction_covered, oppoHomoRates, color='blue', alpha=0.3, label='ancIBD', s=5.0)

        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/wgs/{cov}', bin[0], bin[1], mode='IBDseq', threshold=0.5)
        axs[i,0].scatter(fraction_covered, oppoHomoRates, color='orange', alpha=0.3, label='IBDseq', s=5.0)

        if i == 0:
            axs[i,0].legend(loc='upper right', fontsize='large')
    
    #### plot 1240k BAM results
    for i, cov in enumerate(covs[:4]):
        i_ = i
        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/1240k/{cov}', bin[0], bin[1], threshold=0.5)
        axs[i_,1].scatter(fraction_covered, oppoHomoRates, color='blue', alpha=0.3, label='ancIBD', s=5.0)

        density, fraction_covered, oppoHomoRates = calc_FP_all(f'../callIBD/1240k/{cov}', bin[0], bin[1], mode='IBDseq', threshold=0.5)
        axs[i_,1].scatter(fraction_covered, oppoHomoRates, color='orange', alpha=0.3, label='IBDseq', s=5.0)
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
    plt.savefig(f'./sidebyside/{bin[0]}_{bin[1]}_IBDseq.png', dpi=300, bbox_inches="tight")
    plt.savefig(f'./sidebyside/{bin[0]}_{bin[1]}_IBDseq.pdf', dpi=300, bbox_inches="tight")
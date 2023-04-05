import argparse
import sys
sys.path.append('/mnt/archgen/users/yilei/tools/hapBLOCK/python3/')
from run import hapBLOCK_all
import numpy as np
import os
from collections import defaultdict
import itertools
import argparse
import shutil


path2GroundTruth = '/mnt/archgen/users/yilei/IBDsim/groundtruth_IBD_minDP20_wgs/merged.wgs.4cM.seg'
iids = ['I2105', 'I3950', 'I5279', 'I5273']


def read_ibis_output():
    results = defaultdict(lambda: defaultdict(lambda: []))
    with open(path2GroundTruth) as f:
        for line in f:
            id1, id2, ch, _, _, _, start_cm, end_cm, *_ = line.strip().split()
            ch, start_cm, end_cm = int(ch), float(start_cm), float(end_cm)
            results[(min(id1, id2), max(id1, id2))][ch].append((start_cm, end_cm))
    return results

def intersect_intervals(truth, inferred):
    # truth: a list of true ROH blocks
    # inferred: a list of inferred ROH blocks
    # return the length of intersection of these two lists of ROH blocks
    i = 0 
    j = 0
    n = len(truth)
    m = len(inferred)
    accu = 0
    while (i < n and j < m):
        l = max(truth[i][0], inferred[j][0])
        r = min(truth[i][1], inferred[j][1])
        if r > l:
            accu += r - l
 
        if truth[i][1] < inferred[j][1]:
            i += 1
        else:
            j += 1
    return accu

def grabSegments_hapBLOCK(basepath, pair, mode='hapBLOCK'):
    inferred = defaultdict(lambda: [])
    with open(os.path.join(basepath, mode, f'{pair[0]}_{pair[1]}', f'{pair[0]}_{pair[1]}_ibd.tsv')) as f:
        line = f.readline() # skip header line
        line = f.readline()
        while line:
            _, _, startM, endM, _, _, ch, *_ = line.strip().split()
            ch, start_cm, end_cm = int(ch), 100*float(startM), 100*float(endM)
            inferred[ch].append((start_cm, end_cm))
            line = f.readline()
    return inferred

def calc_ppv_sensitivty_one_pair(basepath, pair, min_l, max_l, ground_truth, mode='hapBLOCK'):
    inferred = None
    if mode.startswith('hapBLOCK'):
        inferred = grabSegments_hapBLOCK(basepath, pair, mode)
    else:
        return RuntimeError('Unsupported Mode')
    
    # calculate precision, the fraction of the total length of all inferred IBD segments (within a given length bin) 
    # that overlap any true segment of any size
    ppv_numerator = 0
    ppv_denominator = 0
    for ch, ibds in inferred.items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        ppv_denominator += sum([ibd[1]-ibd[0] for ibd in ibds_within_bin])
        ppv_numerator += intersect_intervals(ground_truth[pair][ch], ibds_within_bin)
    
    # calculate sensitivity (aka recall), the fraction of the total length of all true IBD segments (in a length bin) 
    # that an algorithm calls as identical by descent, considering all called segments of any size.
    sensitivity_numerator = 0
    sensitivity_denominator = 0
    for ch, ibds in ground_truth[pair].items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        sensitivity_denominator += sum([ibd[1]-ibd[0] for ibd in ibds_within_bin])
        sensitivity_numerator += intersect_intervals(ibds_within_bin, inferred[ch])
    return ppv_numerator, ppv_denominator, sensitivity_numerator, sensitivity_denominator


def calc_ppv_sensitivty_one_batch(basepath, min_l, max_l, groundtruth, mode='hapBLOCK'):
    ppv_numerator_accu = 0
    ppv_denominator_accu = 0
    sensitivity_numerator_accu = 0 
    sensitivity_denominator_accu = 0
    for id1, id2 in itertools.combinations(iids, 2):
        pair = (min(id1, id2), max(id1, id2))
        ppv_numerator, ppv_denominator, sensitivity_numerator, sensitivity_denominator = \
                        calc_ppv_sensitivty_one_pair(basepath, pair, min_l, max_l, groundtruth, mode=mode)
        ppv_numerator_accu += ppv_numerator
        ppv_denominator_accu += ppv_denominator
        sensitivity_numerator_accu += sensitivity_numerator
        sensitivity_denominator_accu += sensitivity_denominator
    ppv = ppv_numerator_accu/ppv_denominator_accu if ppv_denominator_accu != 0 else np.nan
    sensitivity = sensitivity_numerator_accu/sensitivity_denominator_accu if sensitivity_denominator_accu != 0 else np.nan
    return ppv, sensitivity

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('--min_cm1', action="store", dest="min_cm1", type=float, required=True,
                        help="path to the vcf file")
    parser.add_argument('--post', action="store", dest="post", type=float, required=False, default=0.99,
                        help="posterior prob cutoff")
    parser.add_argument('--gap', action="store", dest="gap", type=float, required=False, default=0.0075,
                        help="maximum gap to merge")
    args = parser.parse_args()

    groundtruth = read_ibis_output()
    bins = [(5,6), (6,8), (8,10), (10,13), (13,15), (15, np.inf)]
    ###########################  grid search for optimal parameter value ######################################
    # a total of 980 runs
    # with open(f'./grid_summary.tsv', 'w') as f:
    #     f.write(f'ibd_in\tibd_out\tibd_jump\tgap\tprecision_5-6cM\tprecision_6-8cM\tprecision_8-10cM\tprecision_10-13cM\tprecision_13-15cM\tprecision_>15cM\trecall_5-6cM\trecall_6-8cM\trecall_8-10cM\trecall_10-13cM\trecall_13-15cM\trecall_>15cM\n')
    #     for i, ibd_in in enumerate([0.8, 0.9, 1.0, 1.1, 1.2]): #5
    #         for ibd_out in [8, 10, 12, 15]: #4
    #             for ibd_jump in [400, 450, 500, 550, 600, 650, 700]: #7
    #                 for j, gap in enumerate([0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075, 0.008]): #7
    #                     hapBLOCK_all(folder_in="./processed_1240k/ch", iids=['I2105', 'I3950', 'I5273', 'I5279'],
    #                         chs=range(1,23), folder_out=f"./grid/hapBLOCK", output=False, prefix_out="", logfile=False,
    #                         l_model="hdf5", IBD2=False, p_col="variants/AF_ALL", 
    #                         ibd_in=ibd_in, ibd_out=ibd_out, ibd_jump=ibd_jump, min_cm1=args.min_cm1, min_cm2=2,
    #                         cutoff_post=args.post, max_gap=gap, snp_cm=0, save=0)

    #                     precisions = np.full((len(bins),), np.nan)
    #                     recalls = np.full((len(bins),), np.nan)
    #                     for k, bin in enumerate(bins):
    #                         precision, recall = calc_ppv_sensitivty_one_batch('./grid', bin[0], bin[1], groundtruth, mode='hapBLOCK')
    #                         precisions[k] = precision
    #                         recalls[k] = recall
    #                     precisions = "\t".join(map(lambda x: str(round(x,3)), precisions))
    #                     recalls = "\t".join(map(lambda x: str(round(x,3)), recalls))
    #                     f.write(f'{ibd_in}\t{ibd_out}\t{ibd_jump}\t{gap}\t{precisions}\t{recalls}\n')
    #                     shutil.rmtree('./grid/hapBLOCK')

    with open(f'./grid_summary.round2.tsv', 'w') as f:
        f.write(f'ibd_in\tibd_out\tibd_jump\tgap\tpost_cutoff\tsnp_cm\tprecision_5-6cM\tprecision_6-8cM\tprecision_8-10cM\tprecision_10-13cM\tprecision_13-15cM\tprecision_>15cM\trecall_5-6cM\trecall_6-8cM\trecall_8-10cM\trecall_10-13cM\trecall_13-15cM\trecall_>15cM\n')
        for i, ibd_in in enumerate([0.1, 1.0, 10]): #3
            for ibd_out in [5, 10, 20]: #3
                for ibd_jump in [100, 500, 1000]: #3
                    for j, gap in enumerate([0, 0.005, 0.01]): #3
                        for post in [0.99, 0.9925, 0.995]: #3
                            for snp_cm in [180, 200, 220]: #3
                                hapBLOCK_all(folder_in="./processed_1240k/ch", iids=['I2105', 'I3950', 'I5273', 'I5279'],
                                    chs=range(1,23), folder_out=f"./grid/hapBLOCK", output=False, prefix_out="", logfile=False,
                                    l_model="hdf5", IBD2=False, p_col="variants/AF_ALL", 
                                    ibd_in=ibd_in, ibd_out=ibd_out, ibd_jump=ibd_jump, min_cm1=args.min_cm1, min_cm2=2,
                                    cutoff_post=post, max_gap=gap, snp_cm=snp_cm, save=0)

                                precisions = np.full((len(bins),), np.nan)
                                recalls = np.full((len(bins),), np.nan)
                                for k, bin in enumerate(bins):
                                    precision, recall = calc_ppv_sensitivty_one_batch('./grid', bin[0], bin[1], groundtruth, mode='hapBLOCK')
                                    precisions[k] = precision
                                    recalls[k] = recall
                                precisions = "\t".join(map(lambda x: str(round(x,3)), precisions))
                                recalls = "\t".join(map(lambda x: str(round(x,3)), recalls))
                                f.write(f'{ibd_in}\t{ibd_out}\t{ibd_jump}\t{gap}\t{post}\t{snp_cm}\t{precisions}\t{recalls}\n')
                                shutil.rmtree('./grid/hapBLOCK')
                    
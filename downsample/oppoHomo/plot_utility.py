path2GroundTruth = '/mnt/archgen/users/yilei/IBDsim/groundtruth_IBD_glimpse_transversion/MAC0/merged.glimpse.GP99.transversion_only.seg'
path2HDF5 = '/mnt/archgen/users/yilei/IBDsim/groundtruth_IBD_glimpse_transversion/MAC0/hdf5'
iids = ['I2105', 'I3950', 'I5279', 'I5273']
from collections import defaultdict
import itertools
import os
import numpy as np
import h5py
from tqdm import tqdm
from collections import defaultdict
import time

# build a dictionary of 1240k targets for ibis
def build_dict_1240k_targets():
    ret = {}
    for ch in range(1,23):
        listofbp = []
        with open(f'/mnt/archgen/users/yilei/IBDsim/snps/ch{ch}.1240k.snps') as f:
            for line in f:
                _, _, _, bp, _, _ = line.strip().split()
                listofbp.append(int(bp))
        ret[ch] = np.array(listofbp)
    return ret

dict_1240k = build_dict_1240k_targets()

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

def grabSegments_ancIBD(basepath, batch, pair, mode='ancIBD'):
    inferred = defaultdict(lambda: [])
    with open(os.path.join(basepath, f'batch{batch}', mode, f'{pair[0]}_{pair[1]}', f'{pair[0]}_{pair[1]}_ibd.tsv')) as f:
        line = f.readline() # skip header line
        line = f.readline()
        while line:
            _, _, ch, start, end, nSNP, StartM, EndM, _, startBP, endBP = line.strip().split()
            ch, start, end, start_cm, end_cm = int(ch), int(start), int(end), 100*float(StartM), 100*float(EndM)
            inferred[ch].append((start_cm, end_cm, int(startBP), int(endBP)))
            line = f.readline()
    return inferred

def grabSegments_ibis(basepath, batch, pair, mode='IBIS_1240k'):
    inferred = defaultdict(lambda: [])
    with open(os.path.join(basepath, f'batch{batch}', mode, f'batch{batch}.seg')) as f:
        for line in f:
            id1, id2, ch, startbp, endbp, _, startCM, endCM, *_ = line.strip().split()
            id1, id2 = min(id1, id2), max(id1, id2)
            if id1 == pair[0] and id2 == pair[1]:
                ch, startCM, endCM = int(ch), float(startCM), float(endCM)
                startbp, endbp = int(startbp), int(endbp)
                i = np.searchsorted(dict_1240k[ch], startbp)
                j = np.searchsorted(dict_1240k[ch], endbp)
                inferred[ch].append((startCM, endCM, startbp, endbp))
    return inferred

def readHapMap(path2Map):
    # assume the first row is header, so we ignore it
    bps = []
    cMs = []
    with open(path2Map) as f:
        f.readline()
        line = f.readline()
        while line:
            _, bp, _, cM = line.strip().split()
            bps.append(int(bp))
            cMs.append(float(cM))
            line = f.readline()
    return np.array(bps), np.array(cMs)


def bp2Morgan(bp, bps, cMs):
    # bps: a list of basepair position
    # cMs: a list of geneticMap position in cM corresponding to bps
    assert(len(bps) == len(cMs))
    i = np.searchsorted(bps, bp, side='left')

    if i == 0:
        return cMs[0]*(bp/bps[0])
    elif i == len(bps):
        return cMs[-1]
    elif bps[i] == bp:
        return cMs[i]
    else:
        left_bp, right_bp = bps[i-1], bps[i]
        left_cM, right_cM = cMs[i-1], cMs[i]
        return left_cM + (right_cM - left_cM)*(bp - left_bp)/(right_bp - left_bp)


def grabSegments_IBDseq(basepath, batch, pair, bpsMap, cMsMap, mode='IBDseq'):
    inferred = defaultdict(lambda:[])
    with open(os.path.join(basepath, f'batch{batch}', mode, f'batch{batch}.IBDseq.ibd')) as f:
        for line in f:
            id1, _, id2, _, ch, startBP, endBP, _ = line.strip().split()
            id1, id2 = min(id1, id2), max(id1, id2)
            if id1 == pair[0] and id2 == pair[1]:
                ch, startBP, endBP = int(ch), int(startBP), int(endBP)
                startCM, endCM = bp2Morgan(startBP, bpsMap[ch], cMsMap[ch]), bp2Morgan(endBP, bpsMap[ch], cMsMap[ch])
                inferred[ch].append((startCM, endCM, startBP, endBP))
    return inferred


def get_individual_idx(f, iid="", f_col="samples"):
    """Return index of individual iid"""
    samples = f[f_col].asstr()[:]
    idx = (samples == iid)
    if np.sum(idx)!=1:
        raise RuntimeWarning(f"{np.sum(idx)} entries found for {iid}")
    assert(np.sum(idx)>0) # Sanity Check
    idx=np.where(idx)[0][0]
    return idx  

def calc_snpDensity_one_pair(basepath, batch, pair, min_l, max_l, ground_truth, bp_vec, homo_index, oppoHomo_index, mode='ancIBD', threshold=0.5):
    inferred = None
    if mode.startswith('ancIBD'):
        inferred = grabSegments_ancIBD(basepath, batch, pair, mode=mode)
    elif mode.startswith('IBIS'):
        inferred = grabSegments_ibis(basepath, batch, pair, mode=mode)
    elif mode.startswith('IBDseq'):        
        bpsMap = {}
        cMsMap = {}
        for ch in np.arange(1, 23):
            bps, cms = readHapMap(f"/mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{ch}.txt")
            bpsMap[ch] = bps
            cMsMap[ch] = cms
        inferred = grabSegments_IBDseq(basepath, batch, pair, bpsMap, cMsMap, mode=mode)
    else:
        return RuntimeError('Unsupported Mode')
    
    # calculate precision, the fraction of the total length of all inferred IBD segments (within a given length bin) 
    # that overlap any true segment of any size
    density = []
    fraction_covered = []
    oppoHomoRates = []
    for ch, ibds in inferred.items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        for ibd in ibds_within_bin:
            covered_frac = intersect_intervals(ground_truth[pair][ch], [(ibd[0], ibd[1])])/(ibd[1]-ibd[0])
            fraction_covered.append(covered_frac)
            i, j = np.searchsorted(bp_vec[ch], ibd[2]), np.searchsorted(bp_vec[ch], ibd[3])
            j -= 1
            numHomo = np.sum(np.logical_and(homo_index[ch] >= i, homo_index[ch] <= j))
            numOppoHomo = np.sum(np.logical_and(oppoHomo_index[ch] >= i, oppoHomo_index[ch] <= j))
            oppoHomorate = numOppoHomo/numHomo if numHomo != 0 else 0.0
            oppoHomoRates.append(oppoHomorate)
            if covered_frac > 0.8 and oppoHomorate > 0.05:
                # suspicious segments
                print(f'{pair}, ch{ch}, {ibd} has PPV {covered_frac} with oppoHomo {oppoHomorate}')
                
    return density, fraction_covered, oppoHomoRates

def calc_FP_one_batch(basepath, batch, min_l, max_l, groundtruth, bp_vec, homo_index, oppoHomo_index, mode='ancIBD', threshold=0.5):
    # return the number of false positive segments (within a given length bin) in this batch
    density_all = []
    fraction_covered_all = []
    oppoHomo_rates_all = []
    for id1, id2 in itertools.combinations(iids, 2):
        pair = (min(id1, id2), max(id1, id2))
        density, fraction_covered, oppoHomoRates = calc_snpDensity_one_pair(basepath, batch, pair, min_l, max_l, groundtruth, bp_vec, homo_index[pair], oppoHomo_index[pair], mode=mode, threshold=threshold)
        density_all.extend(density)
        fraction_covered_all.extend(fraction_covered)
        oppoHomo_rates_all.extend(oppoHomoRates)
    return density_all, fraction_covered_all, oppoHomo_rates_all


from scipy import stats

def calc_FP_all(basepath, min_l, max_l, bp_vec, homo_index, oppoHomo_index, mode='ancIBD', threshold=0.5):
    groundtruth = read_ibis_output()
    density_all = []
    fraction_covered_all = []
    oppoHomo_all = []

    for b in tqdm(range(1, 51)):
        density, fraction_covered, oppoHomo = calc_FP_one_batch(basepath, b, min_l, max_l, groundtruth, bp_vec, homo_index, oppoHomo_index, mode=mode, threshold=threshold)
        density_all.extend(density)
        fraction_covered_all.extend(fraction_covered)
        oppoHomo_all.extend(oppoHomo)

    return density_all, fraction_covered_all, oppoHomo_all
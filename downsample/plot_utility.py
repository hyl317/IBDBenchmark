iids = ['I2105', 'I3950', 'I5279', 'I5273']
from collections import defaultdict
import itertools
import os
import numpy as np
from collections import defaultdict
import gzip

def read_ibis_output(path2groundtruth):
    results = defaultdict(lambda: defaultdict(lambda: []))
    with open(path2groundtruth) as f:
        for line in f:
            if line.startswith('#'):
                continue
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

def intersect_intervals_relaxed(set1, set2, threshold=0.5, verbose=False):
    # return the number of segments in set2 that are covered by segments in set1 over 50% (or whatever percentage specified by the threshold) of its length
    count = 0
    for seg in set2:
        overlap = intersect_intervals(set1, [seg])
        if overlap/(seg[1] - seg[0]) > threshold:
            count += 1
        else:
            if verbose:
                print(f'Segment {round(seg[0],3)}-{round(seg[1],3)} is not detected')
    return count

def readMaskTrack(path2mask):
    dict = defaultdict(lambda:[])
    with open(path2mask) as f:
        for line in f:
            ch, startBP, endBP, startCM, endCM = line.strip().split()
            ch, startBP, endBP, startCM, endCM = int(ch), int(startBP), int(endBP), float(startCM), float(endCM)
            dict[ch].append((startCM, endCM))
    return dict

def maskSegments(unmasked_segments, masks):
    masked_segments = defaultdict(lambda:[])
    for ch, segs in unmasked_segments.items():
        for seg in segs:
            start_cm, end_cm = seg[0], seg[1]
            for mask in masks[ch]:
                mask_start, mask_end = mask[0], mask[1]
                if end_cm <= mask_start or start_cm >= mask_end:
                    continue
                elif start_cm >= mask_start: # implicitly, start_cm < mask_end, guaranteed by the above condition
                    start_cm = mask_end
                    break
                elif end_cm <= mask_end: # implicitly, end_cm > mask_start
                    end_cm = mask_start
                    break
                else:
                    assert(start_cm < mask_start)
                    assert(end_cm > mask_end)
                    start_cm = (start_cm, mask_end)
                    end_cm = (mask_start, end_cm)
                    break
            if isinstance(start_cm, tuple):
                masked_segments[ch].append((start_cm[0], end_cm[0]))
                masked_segments[ch].append((start_cm[1], end_cm[1]))
            else:
                if end_cm > start_cm:
                    masked_segments[ch].append((start_cm, end_cm))
    return masked_segments


def grabSegments_ancIBD(basepath, batch, pair, mode='ancIBD', masks=None, snpcm=220):
    inferred = defaultdict(lambda: [])
    with open(os.path.join(basepath, f'batch{batch}', mode, f'{pair[0]}_{pair[1]}', f'{pair[0]}_{pair[1]}_ibd.tsv')) as f:
        line = f.readline() # skip header line
        line = f.readline()
        while line:
            _, _, ch, _, _, nSNP, StartM, EndM, *_ = line.strip().split()
            ch, start_cm, end_cm, nSNP = int(ch), 100*float(StartM), 100*float(EndM), int(nSNP)
            if nSNP/(end_cm - start_cm) < snpcm:
                line = f.readline()
                continue
            inferred[ch].append((start_cm, end_cm))
            line = f.readline()
    if masks:
        inferred = maskSegments(inferred, masks)
    return inferred

def grabSegments_ibis(basepath, batch, pair, mode='IBIS', masks=None):
    inferred = defaultdict(lambda: [])
    with open(os.path.join(basepath, f'batch{batch}', mode, f'batch{batch}.seg')) as f:
        for line in f:
            id1, id2, ch, _, _, _, startCM, endCM, *_ = line.strip().split()
            id1, id2 = min(id1, id2), max(id1, id2)
            if id1 == pair[0] and id2 == pair[1]:
                ch, startCM, endCM = int(ch), float(startCM), float(endCM)
                inferred[ch].append((startCM, endCM))
    if masks:
        inferred = maskSegments(inferred, masks)
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


def grabSegments_hapIBD(basepath, batch, pair, bpsMap, cMsMap, mode='hapIBD', masks=None):
    inferred = defaultdict(lambda: [])
    with gzip.open(os.path.join(basepath, f'batch{batch}', mode, f'batch{batch}.ibd.gz'), mode='rt') as f:
        for line in f:
            id1, _, id2, _, ch, startBP, endBP, length = line.strip().split()
            id1, id2 = min(id1, id2), max(id1, id2)
            if id1 == pair[0] and id2 == pair[1]:
                ch, startBP, endBP, length = int(ch), int(startBP), int(endBP), float(length)
                startCM = bp2Morgan(startBP, bpsMap[ch], cMsMap[ch])
                endCM = bp2Morgan(endBP, bpsMap[ch], cMsMap[ch])
                inferred[ch].append((startCM, endCM))
    if masks:
        inferred = maskSegments(inferred, masks)
    return inferred

def grabSegments_IBDseq(basepath, batch, pair, bpsMap, cMsMap, mode='IBDseq', masks=None):
    inferred = defaultdict(lambda: [])
    with open(os.path.join(basepath, f'batch{batch}', mode, f'batch{batch}.IBDseq.ibd')) as f:
        for line in f:
            id1, _, id2, _, ch, startBP, endBP, _ = line.strip().split()
            id1, id2 = min(id1, id2), max(id1, id2)
            if id1 == pair[0] and id2 == pair[1]:
                ch, startBP, endBP = int(ch), int(startBP), int(endBP)
                startCM = bp2Morgan(startBP, bpsMap[ch], cMsMap[ch])
                endCM = bp2Morgan(endBP, bpsMap[ch], cMsMap[ch])
                if endCM - startCM > 5.0:
                    inferred[ch].append((startCM, endCM))
    if masks:
        inferred = maskSegments(inferred, masks)
    return inferred

def calc_ppv_sensitivty_one_pair(basepath, batch, pair, min_l, max_l, ground_truth, mode='ancIBD', bpsMap=None, cMsMap=None, masks=None, snpcm=220):
    ###############
    # the definition of precision and recall used in this function is taken from the IBIS paper
    ###############
    inferred = None
    if mode.startswith('ancIBD'):
        inferred = grabSegments_ancIBD(basepath, batch, pair, mode=mode, masks=masks, snpcm=snpcm)
    elif mode.startswith('IBIS'):
        inferred = grabSegments_ibis(basepath, batch, pair, mode=mode, masks=masks)
    elif mode.startswith('hapIBD'): 
        inferred = grabSegments_hapIBD(basepath, batch, pair, bpsMap, cMsMap, mode=mode, masks=masks)
    elif mode.startswith('IBDseq'):
        inferred = grabSegments_IBDseq(basepath, batch, pair, bpsMap, cMsMap, mode=mode, masks=masks)
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


def calc_ppv_sensitivty_one_pair_relaxed(basepath, batch, pair, min_l, max_l, ground_truth, mode='ancIBD', bpsMap=None, cMsMap=None, masks=None, snpcm=220):
    ###############
    # this function uses a more relaxed def of precision and recall
    ###############
    inferred = None
    if mode.startswith('ancIBD'):
        inferred = grabSegments_ancIBD(basepath, batch, pair, mode=mode, masks=masks, snpcm=snpcm)
    elif mode.startswith('IBIS'):
        inferred = grabSegments_ibis(basepath, batch, pair, mode=mode, masks=masks)
    elif mode.startswith('hapIBD'):
        inferred = grabSegments_hapIBD(basepath, batch, pair, bpsMap, cMsMap, mode=mode, masks=masks)
    elif mode.startswith('IBDseq'):
        inferred = grabSegments_IBDseq(basepath, batch, pair, bpsMap, cMsMap, mode=mode, masks=masks)
    else:
        return RuntimeError('Unsupported Mode')
    
    ppv_numerator = 0
    ppv_denominator = 0
    for ch, ibds in inferred.items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        ppv_denominator += len(ibds_within_bin)
        ppv_numerator += intersect_intervals_relaxed(ground_truth[pair][ch], ibds_within_bin)
    
    sensitivity_numerator = 0
    sensitivity_denominator = 0
    for ch, ibds in ground_truth[pair].items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        sensitivity_denominator += len(ibds_within_bin)
        sensitivity_numerator += intersect_intervals_relaxed(inferred[ch], ibds_within_bin)
    return ppv_numerator, ppv_denominator, sensitivity_numerator, sensitivity_denominator


def calc_ppv_sensitivty_one_batch(basepath, batch, min_l, max_l, groundtruth, mode='ancIBD', bpsMap=None, cMsMap=None, masks=None, snpcm=220):
    ppv_numerator_accu = 0
    ppv_denominator_accu = 0
    sensitivity_numerator_accu = 0 
    sensitivity_denominator_accu = 0
    for id1, id2 in itertools.combinations(iids, 2):
        pair = (min(id1, id2), max(id1, id2))
        ppv_numerator, ppv_denominator, sensitivity_numerator, sensitivity_denominator = \
                        calc_ppv_sensitivty_one_pair_relaxed(basepath, batch, pair, min_l, max_l, groundtruth, mode=mode, \
                        bpsMap=bpsMap, cMsMap=cMsMap, masks=masks, snpcm=snpcm)
        ppv_numerator_accu += ppv_numerator
        ppv_denominator_accu += ppv_denominator
        sensitivity_numerator_accu += sensitivity_numerator
        sensitivity_denominator_accu += sensitivity_denominator
    ppv = ppv_numerator_accu/ppv_denominator_accu if ppv_denominator_accu != 0 else 1.0
    sensitivity = sensitivity_numerator_accu/sensitivity_denominator_accu if sensitivity_denominator_accu != 0 else 1.0
    return ppv, sensitivity


from scipy import stats

def calc_ppv_sensitivity_all(path2groundtruth, basepath, min_l, max_l, mode='ancIBD', masks=None, snpcm=220):
    groundtruth = read_ibis_output(path2groundtruth)
    ppv_ = []
    sensitivity_ = []

    if mode.startswith('hapIBD') or mode.startswith('IBDseq'):
        # read maps
        bpsMap = {}
        cMsMap = {}
        for ch in np.arange(1, 23):
            bps, cms = readHapMap(f"/mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{ch}.txt")
            bpsMap[ch] = bps
            cMsMap[ch] = cms
    else:
        bpsMap = None
        cMsMap = None

    for b in range(1, 51):
        ppv, sensitivity = calc_ppv_sensitivty_one_batch(basepath, b, min_l, max_l, groundtruth, mode=mode, bpsMap=bpsMap, cMsMap=cMsMap, masks=masks, snpcm=snpcm)
        ppv_.append(ppv)
        sensitivity_.append(sensitivity)
    return np.nanmean(ppv_), stats.sem(ppv_, nan_policy='omit'), np.nanmean(sensitivity_), stats.sem(sensitivity_, nan_policy='omit')

def get_ppv_sensitivity_allCov(covs, min_l, max_l, path2groundtruth, datatype='wgs', mode='ancIBD', masks=None, snpcm=220):
    ppv_ = np.zeros(len(covs))
    ppv_se_ = np.zeros(len(covs))
    sensitivity_ = np.zeros(len(covs))
    sensitivity_se_ = np.zeros(len(covs))
    for k, cov in enumerate(covs):
        ppv, ppv_se, sensitivity, sensitivity_se = calc_ppv_sensitivity_all(path2groundtruth, \
            f'./callIBD/{datatype}/{cov}', min_l, max_l, mode=mode, masks=masks, snpcm=snpcm)
        ppv_[k] = ppv
        ppv_se_[k] = ppv_se
        sensitivity_[k] = sensitivity
        sensitivity_se_[k] = sensitivity_se
    return ppv_, ppv_se_, sensitivity_, sensitivity_se_
    
####################################################################################################################################
################################################### Examine Length Bias ############################################################

def find_overlap(target_segment, inferred_segments, threshold=0.5):
    # return the segment in the set of inferred segments that overlaps with the target segment
    # if no such segment exist, return None
    start0, end0 = target_segment
    length = end0 - start0
    for segment in inferred_segments:
        start, end = segment
        overlap_start = max(start, start0)
        overlap_end = min(end, end0)
        if overlap_start < overlap_end and (overlap_end - overlap_start) >= threshold*length:
            return segment
    return None

def calc_lengthBias_one_pair(basepath, batch, pair, min_l, max_l, ground_truth, mode='ancIBD', masks=None):
    inferred = None
    if mode.startswith('ancIBD'):
        inferred = grabSegments_ancIBD(basepath, batch, pair, mode, masks)
    elif mode.startswith('IBIS'):
        inferred = grabSegments_ibis(basepath, batch, pair, mode, masks)
    else:
        return RuntimeError('Unsupported Mode')

    length_diff = []
    for ch, ibds in ground_truth[pair].items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        for ibd in ibds_within_bin:
            overlapping_segment = find_overlap(ibd, inferred[ch])
            if overlapping_segment != None:
                length_diff.append((overlapping_segment[1] - overlapping_segment[0]) - (ibd[1] - ibd[0]))
                #print(f'found an overlap segment of {ibd}: {overlapping_segment}, {pair}, batch{batch}')
    return length_diff

def calc_lengthBias_one_batch(basepath, batch, min_l, max_l, groundtruth, mode='ancIBD', masks=None):
    length_diff = []
    for id1, id2 in itertools.combinations(iids, 2):
        pair = (min(id1, id2), max(id1, id2))
        length_diff.extend(calc_lengthBias_one_pair(basepath, batch, pair, min_l, max_l, groundtruth, mode=mode, masks=masks))
    return length_diff

def calc_lengthBias_all(path2groundtruth, basepath, min_l, max_l, mode='ancIBD', masks=None):
    groundtruth = read_ibis_output(path2groundtruth)
    length_diff = []
    for b in range(1, 51):
        length_diff.extend(calc_lengthBias_one_batch(basepath, b, min_l, max_l, groundtruth, mode=mode, masks=masks))
    return length_diff

########################################## Examine Power ####################################################
def calc_power_one_pair(basepath, batch, pair, min_l, max_l, ground_truth, mode='ancIBD', masks=None):
    """
    For a given length bin, return the number of segments that are detected, and the total number of segments fall into this length bin
    """
    inferred = None
    if mode.startswith('ancIBD'):
        inferred = grabSegments_ancIBD(basepath, batch, pair, mode, masks)
    elif mode.startswith('IBIS'):
        inferred = grabSegments_ibis(basepath, batch, pair, mode, masks)
    else:
        return RuntimeError('Unsupported Mode')

    denominator = 0
    numerator = 0
    for ch, ibds in ground_truth[pair].items():
        ibds_within_bin = [ibd for ibd in ibds if ibd[1] - ibd[0] >= min_l and ibd[1] - ibd[0] <= max_l]
        denominator += len(ibds_within_bin)
        for ibd in ibds_within_bin:
            overlapping_segment = find_overlap(ibd, inferred[ch])
            if overlapping_segment != None:
                numerator += 1
    return numerator, denominator

def calc_power_one_batch(basepath, batch, min_l, max_l, groundtruth, mode='ancIBD', masks=None):
    numerator_sum = 0
    denominator_sum = 0
    for id1, id2 in itertools.combinations(iids, 2):
        pair = (min(id1, id2), max(id1, id2))
        numerator, denominator = calc_power_one_pair(basepath, batch, pair, min_l, max_l, groundtruth, mode=mode, masks=masks)
        numerator_sum += numerator
        denominator_sum += denominator
    return numerator_sum, denominator_sum

def calc_power_all(path2groundtruth, basepath, min_l, max_l, mode='ancIBD', masks=None):
    groundtruth = read_ibis_output(path2groundtruth)
    numerator_sum = 0
    denominator_sum = 0
    for b in range(1, 51):
        numerator, denominator = calc_power_one_batch(basepath, b, min_l, max_l, groundtruth, mode=mode, masks=masks)
        numerator_sum += numerator
        denominator_sum += denominator
    return numerator_sum/denominator_sum
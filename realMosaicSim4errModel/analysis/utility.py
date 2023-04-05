import pandas as pd
import numpy as np
from collections import defaultdict

# some utility function for analyzing ancIBD's performance on simulated data

def readGroundtruthIBD(path2file):
    """
    Read the ibd_info.csv file, and return a dictionary such that
    key is the sample pair
    value is a tuple representing the groundtruth IBD: (startM, endM)
    """
    with open(path2file) as f:
        record = {}
        f.readline() # ignore headerline
        line = f.readline()
        while line:
            beginM, endM, iid1, iid2, *_ = line.strip().split()
            beginM, endM = float(beginM), float(endM)
            iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
            record[(iid1, iid2)] = (beginM, endM)
            line = f.readline()
    return record

def readInferredIBD(path2file):
    # Note that different from groundtruth, in inferred IBD, a pair of samples may share more than one IBD segments
    record = defaultdict(lambda:[])
    with open(path2file) as f:
        f.readline() # ignore headerline
        line = f.readline()
        while line:
            Start, End, StartM, EndM, length, lengthM, ch, iid1, iid2 = line.strip().split('\t')
            line = f.readline()
            iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
            StartM,EndM = float(StartM), float(EndM)
            record[(iid1, iid2)].append((StartM, EndM))
    return record

def overlapPercentage_precision(inferred, truth):
    # return the percentage of inferred segment that overlaps with the true segment
    s_inf, e_inf = inferred
    s_tru, e_tru = truth
    if e_inf <= s_tru or s_inf >= e_tru:
        return 0.0
    else:
        return (min(e_inf, e_tru) - max(s_inf, s_tru))/(e_inf - s_inf)

def overlapPercentage_power(inferred, truth):
    # return the percentage of the true segment that overlaps with the inferred segment
    s_inf, e_inf = inferred
    s_tru, e_tru = truth
    if e_inf <= s_tru or s_inf >= e_tru:
        return 0.0
    else:
        return (min(e_inf, e_tru) - max(s_inf, s_tru))/(e_tru - s_tru)

def lengthDiffOneBatch(inferredCSV, simulatedCSV, verbose=False):
    # return a list of length difference between the inferred length and the simulated true length (inferred - simulated)
    gtIBD = readGroundtruthIBD(simulatedCSV)
    diffs = []
    df = pd.read_csv(inferredCSV, sep='\t')
    df = df.reset_index()
    for index, row in df.iterrows():
        iid1, iid2 = row['iid1'], row['iid2']
        iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
        startM_inferred, endM_inferred = row['StartM'], row['EndM']
        if overlapPercentage_power((startM_inferred, endM_inferred), gtIBD[(iid1, iid2)]) > 0.8:
            s_tru, e_tru = gtIBD[(iid1, iid2)]
            diff = (endM_inferred - startM_inferred) - (e_tru - s_tru)
            diffs.append(diff)
            if diff > 5.5/100:
                print(f'true segment: {gtIBD[(iid1, iid2)]}, inferred segment: {startM_inferred} - {endM_inferred}')
    diffs = 100*np.array(diffs) # convert to cM
    diffs_pos = diffs[diffs >= 0]
    diffs_neg = diffs[diffs < 0]
    if verbose:
        print(f'percentage of overshoot: {len(diffs_pos)/len(diffs)}')
        print(f'mean of overshoot: {np.mean(diffs_pos)}')
        print(f'mean of undershoot: {np.mean(diffs_neg)}')
    return diffs, diffs_pos, diffs_neg 

def FPOneBatch(inferredCSV, simulatedCSV, cutoff=0.25):
    gtIBD = readGroundtruthIBD(simulatedCSV)
    fps = []
    df = pd.read_csv(inferredCSV, sep=',')
    df = df.reset_index()
    for index, row in df.iterrows():
        iid1, iid2 = row['iid1'], row['iid2']
        iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
        startM_inferred, endM_inferred = row['StartM'], row['EndM']
        if overlapPercentage_precision((startM_inferred, endM_inferred), gtIBD[(iid1, iid2)]) < cutoff:
            fps.append(endM_inferred - startM_inferred)
    return fps


def PowerOneBatch(inferredCSV, simulatedCSV, threshold=0.1):
    """
    Return the percentage of segments in simulatedCSV that are detected in the inferred CSV 
    """
    inferred = readInferredIBD(inferredCSV)
    simulated = readGroundtruthIBD(simulatedCSV)
    count = 0
    for pair, simulated_seg in simulated.items():
        inferred_segs = inferred[pair]
        for inferred_seg in inferred_segs:
            if overlapPercentage_power(inferred_seg, simulated_seg) > threshold:
                count += 1
                break
    return count/len(simulated)

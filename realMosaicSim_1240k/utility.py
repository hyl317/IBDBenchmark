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
            iid1,iid2,ch,Start,End,length,StartM,EndM,lengthM,StartBP,EndBP = line.strip().split(',')
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
    # return the percentage of true segment that overlaps with an inferred segment
    s_inf, e_inf = inferred
    s_tru, e_tru = truth
    if e_inf <= s_tru or s_inf >= e_tru:
        return 0.0
    else:
        return (min(e_inf, e_tru) - max(s_inf, s_tru))/(e_tru - s_tru)

def getOverlap(seg1, seg2):
    # get the length of overlapping regions of the two given segments
    s1, e1 = seg1
    s2, e2 = seg2
    if e1 <= s2 or s1 >= e2:
        return 0
    else:
        return min(e1, e2) - max(s1, s2)


def precisionOneBatch(inferredCSV, simulatedCSV, minL=-np.inf, maxL=np.inf):
    minL, maxL = minL/100, maxL/100
    gtIBD = readGroundtruthIBD(simulatedCSV)
    df = pd.read_csv(inferredCSV, sep=',')
    df = df.reset_index()
    df = df[np.logical_and(df['lengthM'] >= minL, df['lengthM'] <= maxL)]
    # iterate over each of the inferred segments
    tot = 0.0
    overlaps = 0.0
    for index, row in df.iterrows():
        iid1, iid2 = row['iid1'], row['iid2']
        iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
        startM_inferred, endM_inferred = row['StartM'], row['EndM']
        overlap = getOverlap((startM_inferred, endM_inferred), gtIBD[(iid1, iid2)])
        tot += endM_inferred - startM_inferred
        overlaps += overlap
    return overlaps/tot

def powerOneBatch(inferredCSV, simulatedCSV):
    inferredIBD = readInferredIBD(inferredCSV)
    df = pd.read_csv(simulatedCSV, sep='\t')
    df = df.reset_index()
    # iterate over each of the true segments
    tot = 0.0
    overlaps = 0.0
    for index, row in df.iterrows():
        iid1, iid2 = row['iid1'], row['iid2']
        iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
        startM_sim, endM_sim = row['IBD_Begin'], row['IBD_End']
        for seg in inferredIBD[(iid1, iid2)]:
            overlaps += getOverlap((startM_sim, endM_sim), seg)
        tot += endM_sim - startM_sim
    return overlaps/tot

##################################################################################################

def lengthCalledOneBatch(inferredCSV, simulatedCSV, threshold=0.5):
    # return a list of length difference between the inferred length and the simulated true length (inferred - simulated)
    gtIBD = readGroundtruthIBD(simulatedCSV)
    lengthCalled = []
    df = pd.read_csv(inferredCSV, sep='\t')
    df = df.reset_index()
    count = 0
    for index, row in df.iterrows():
        iid1, iid2 = row['iid1'], row['iid2']
        iid1, iid2 = min(iid1, iid2), max(iid1, iid2)
        startM_inferred, endM_inferred = row['StartM'], row['EndM']
        precision = overlapPercentage_precision((startM_inferred, endM_inferred), gtIBD[(iid1, iid2)])
        if precision > threshold:
            lengthCalled.append((endM_inferred - startM_inferred)*100)
            if overlapPercentage_power((startM_inferred, endM_inferred), gtIBD[(iid1, iid2)]) > 0.8:
                count += 1
    return lengthCalled, count/len(gtIBD)

def lengthCallAll(bl_lens, cov, basepath, threshold=0.5):
    all = []
    powers = []
    for l in bl_lens:
        path2file = f'{basepath}/{cov}/ch3_{l}cm/ch3.tsv' if len(cov) > 0 else f'{basepath}/ch3_{l}cm/ch3.tsv'
        ret, power = lengthCalledOneBatch(path2file, \
            f'/mnt/archgen/users/yilei/IBDsim/realMosaicSim_wgs/simGeno/ch3_{l}cm/ibd_info.csv')
        all.append(ret)
        powers.append(power)
    return all, powers


###################################################################################################

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

from collections import defaultdict
import pandas as pd
import numpy as np
from numba import njit
import json
from math import exp, ceil
import os 
import gzip
def get_mut_rates(mut_full, window_size, windows, obs_count, chr):
    '''
    window: window index
    a chromsome
    e.g. 0 1000000 1.5
         1000000 2000000 2.5

    '''
    mut = []
    mut_full = mut_full[mut_full['chrom'] == chr]
    assert len(mut_full)!=0, f"mutation rate missing for chromsome {chr}"
    current_window = 0
    inter = 0   #interval of mut rates
    while current_window <= len(windows):
        while inter < len(mut_full):  
            '''
            window 100
            100 * 1000 - 101 * 1000
            '''


            ''' go to current mut interval'''
            if ((mut_full.iloc[inter]['start'] <= (windows[current_window] * window_size)) and (mut_full.iloc[inter]['end'] >= (windows[current_window] + 1) * window_size)):
                if (mut_full.iloc[inter]['mut_rate'] == 0) and (obs_count[current_window] > 0):   
                    '''
                    might cause numerical problem. so far all 0 mut_rates are manully modifyed
                    '''
                    mut.extend([0.000125] * obs_count[current_window])
                else:   
                    mut.extend([mut_full.iloc[inter]['mut_rate']] * obs_count[current_window])
                current_window += 1
                break
            inter += 1
        if ((inter >= len(mut_full)) or (current_window >= len(windows))):
            break
    return mut
    '''
    seems to load correctly. checked.
    '''
def load_fasta(fasta_file):
    '''
    Read a fasta file with a single chromosome in and return the sequence as a string
    '''
    fasta_sequence = ''
    with open(fasta_file) as data:
        for line in data:
            if not line.startswith('>'):
                fasta_sequence += line.strip().upper()

    return fasta_sequence
                
def get_weights(bedfile, window_size, mut_bed): # return weights from the first window to the last window
    first_window = 0
    last_window = 0
    window_size = int(window_size)
    weights = defaultdict(lambda: defaultdict(float))
    with open(f"{bedfile}") as data:
        print(f"loading mask file {bedfile}")
        for line in data:
            chr, start, end = line.strip().split() 
            start = int(start)
            end = int(end)
            start_window = ceil((start+0.5) / window_size) - 1
            end_window = ceil((end - 0.5) / window_size) - 1
            if start_window == end_window: # same window   start window 1, end window 1
                weights[chr][start_window] += (end - start) / window_size
            else:
                weights[chr][start_window] += (window_size*(start_window+1) - start) / window_size
            
                if end_window > start_window + 1:       # fill in windows in the middle   start_window 1, end_window 4, so window 2 window3 will be filled as 1
                    for window_tofill in range(start_window + 1, end_window):
                        weights[chr][window_tofill] += float(1)
                    weights[chr][end_window] += (end - end_window * window_size) / window_size
                    
                else:    # e.g. start window 1, end window 2
                    weights[chr][start_window + 1] += (end - end_window * window_size) / window_size
        if mut_bed is None:
            print("no mut file provided, assuming constant mutation rates (as 1 )")
            return weights
        else:
            print(f"loading mut file {mut_bed}")
            mut = pd.read_csv(mut_bed, names = ['chr', 'start', 'end', 'rate'], dtype = {'chr':str, 'start':int, 'end':int, 'rate':float}, sep='\s+') 
            for chr in list(weights.keys()):
                mut_chr = mut[mut['chr'] == chr].reset_index(drop = True)
                mut_loop = mut_chr.iterrows()
                mut_row = next(mut_loop)[1]
                for window in list(weights[chr].keys()):
                    while mut_row.end < window * window_size:  # 100 000 - 101 000
                        mut_row = next(mut_loop)[1]
                    weights[chr][window] *= mut_row.rate
            return weights

def sort_chrom(value):
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)

def load_observations_gt(gt_file, mask_file, window_size, max_variants, data_type,mut_bed, phased = False):  # return number of variants
    chr = []
    chr_index = []  #chrs
    windows = []
    weights = []  # weights    [1,0,1,0,2,...] from the first bin with data to last bin with data
    obs = []    # [0,0,0,1,0,0...][0,1,1,0,...] from the first bin with data to last bin with data
    call_index = []    #[0,1,3,5,...]   bins with data
    # assuming full length 100M
    window_all = []
    # observation is always sorted by chr.
    call = get_weights(mask_file, window_size, mut_bed)

    chr_unsorted = list(call.keys())
    if not chr_unsorted[0].startswith("chr"):
        chr_sorted = sorted(chr_unsorted, key=sort_chrom)
    else:
        chr_sorted = sorted(chr_unsorted, key=lambda x: int(x[3:]))
    print(f"all chromosomes loading after sorting: {chr_sorted}")
    chr_total_sort = []
    for chr in chr_sorted:
        chr_total_sort.append(chr)
        first_w = list(call[chr].keys())[0]  # the first window callable
        last_w = list(call[chr].keys())[-1]  # the last window callable
        weights_ = np.zeros(last_w - first_w + 1)
        call_index_ = []
        for i in list(call[chr].keys()):   # loop through all windows with data
            i = int(i)
            call_index_.append(i - first_w)
        for w in range(first_w, last_w+1):
            weights_[w-first_w] = call[chr][w]
        weights.append(weights_)
        call_index.append(call_index_)
    snp = defaultdict(lambda:defaultdict(lambda:defaultdict(int)))
    #print(f"weights loaded for all chrs {chr_total_sort}")
    with open(gt_file, 'r') as data:
        if phased == True:
            print("loading phased data")
            for line in data:
                chr, pos, _, phase = line.strip().split()[0:4]
                window = ceil(int(pos) / window_size) - 1
                snp[chr][int(phase)][window]+=1
        else:
            print("loading unphased data")
            for line in data:
                #chr, pos = line.strip().split()
                chr, pos = line.strip().split()[0:2]
                window = ceil(int(pos) / window_size) - 1
                snp[chr][0][window]+=1            
    for chr in list(snp.keys()):   # has to start from (chr)1
        
        #print(chr)
        chr_index.append(chr)
        first_w = sorted(list(call[chr].keys()))[0]  # assuming all are chr1
        last_w = sorted(list(call[chr].keys()))[-1]
        window_all.append(list(range(first_w, last_w+1)))
        if phased == True:
            p = 2
        else:
            p = 1
        snps = np.zeros((p, last_w - first_w + 1))
        for p_ in range(p):
            for window in list(snp[chr][p_].keys()):
                snps[p_][window-first_w] = snp[chr][p_][window]   # snp from the first bin with data.
                if snps[p_][window-first_w] > max_variants:
                    snps[p_][window-first_w] = max_variants
                    print(f"window {window}, chr{chr} has more than {max_variants} variants, set to {max_variants}, file {gt_file}")
            #assert weights[int(chr)-1][window-first_w] >0 , f"weights for window {window} is 0 but there are derived in this window."

        obs.append(snps)
    return chr_index, weights, obs, call_index, window_all

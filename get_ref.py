import pandas as pd
import numpy as np
def make_ref1(vcf, ref):  # 1:  all segragating sites between den3 and 3 neas
    v = pd.read_csv(vcf, sep='\t')
    v['afr_alt'] = 0
    for i in range(0, 200):
        v['afr_alt'] += v[f'afr_{i}']
    v['afr_ref'] = 200 - v['afr_alt']   #all segragating sites between den3 and 3 neas
    v['freq_afr'] = v['afr_alt']/200
    v['nea_alt'] = v['nea_0'] + v['nea_1'] + v['nea_2']+v['nea_3'] + v['nea_out_0'] + v['nea_out_1']
    v['nea_ref'] = 6 - v['nea_alt']
    v['freq_nea'] = v['nea_alt']/6
    v['d3_alt'] = v['den3_0']+v['den3_1']
    v['d3_ref'] = 2 - v['d3_alt']
    v['freq_d3'] = v['d3_alt']/2
    v['max_diff_d3'] = v[['freq_afr','freq_nea','freq_d3']].max(axis=1) - v[['freq_afr','freq_nea','freq_d3']].min(axis=1)
    """
    hapmap = pd.read_csv(hap_map, sep = '\t', header = None,
    names = ['chr', 'snp', 'gen_pos', 'physical_pos'])
    loop_hapmap = hapmap.iterrows()
    row = next(loop_hapmap)[1]
    row_1 = next(loop_hapmap)[1]
    pos_v = []
    for i in range(0, len(pos)):
        while (row_1.physical_pos <= pos[i]):
            row, row_1 = row_1, next(loop_hapmap)[1]
        pos_v.append(row.gen_pos + (pos[i] - row.physical_pos)*(row_1.gen_pos - row.gen_pos)/(row_1.physical_pos - row.physical_pos))
    v['map'] = pos_v
    """
    pos = v.pos.to_list()
    v['map'] = pos
    v['map'] = v['map'] * 1.25e-6
    v['ref'] = "A"
    v['alt'] = "G"
    v['anc_ref'] = 1
    v['anc_alt'] = 0
    vd3 = v[v['max_diff_d3'] > 0.01].copy()
    vd3 = vd3[['chrom','pos','ref','alt','map','nea_ref','nea_alt','d3_ref','d3_alt','afr_ref','afr_alt', 'anc_ref','anc_alt']]
    vd3 = vd3[vd3['pos']>0]
    vd3.to_csv(ref, index = False)
def make_ref2(vcf, ref):  # 2:  all segragating sites between Den3, Den25 and 3 Neas
    v = pd.read_csv(vcf, sep='\t')
    v['afr_alt'] = 0
    for i in range(0, 200):
        v['afr_alt'] += v[f'afr_{i}']
    v['afr_ref'] = 200 - v['afr_alt']   #all segragating sites between den25 den3 3nea
    v['freq_afr'] = v['afr_alt']/200
    v['nea_alt'] = v['nea_0'] + v['nea_1'] + v['nea_2']+v['nea_3'] + v['nea_out_0'] + v['nea_out_1']
    v['nea_ref'] = 6 - v['nea_alt']
    v['freq_nea'] = v['nea_alt']/6
    v['den_alt'] = v['den3_0']+v['den3_1']+v['den25_0']+v['den25_1']
    v['den_ref'] = 4 - v['den_alt']
    v['freq_den'] = v['den_alt']/4
    v['max_diff_den'] = v[['freq_afr','freq_nea','freq_den']].max(axis=1) - v[['freq_afr','freq_nea','freq_den']].min(axis=1)
    """
    hapmap = pd.read_csv(hap_map, sep = '\t', header = None,
    names = ['chr', 'snp', 'gen_pos', 'physical_pos'])
    loop_hapmap = hapmap.iterrows()
    row = next(loop_hapmap)[1]
    row_1 = next(loop_hapmap)[1]
    pos_v = []
    for i in range(0, len(pos)):
        while (row_1.physical_pos <= pos[i]):
            row, row_1 = row_1, next(loop_hapmap)[1]
        pos_v.append(row.gen_pos + (pos[i] - row.physical_pos)*(row_1.gen_pos - row.gen_pos)/(row_1.physical_pos - row.physical_pos))
    v['map'] = pos_v
    """
    pos = v.pos.to_list()
    v['map'] = pos
    v['map'] = v['map'] * 1.25e-6
    v['ref'] = "A"
    v['alt'] = "G"
    v['anc_ref'] = 1
    v['anc_alt'] = 0
    vden = v[v['max_diff_den'] > 0.01].copy()
    vden = vden[['chrom','pos','ref','alt','map','nea_ref','nea_alt','den_ref','den_alt','afr_ref','afr_alt', 'anc_ref','anc_alt']]
    vden = vden[vden['pos']>0]
    vden.to_csv(ref, index = False)

def make_ref3(vcf, ref):   #3:  all segragating sites between Den3, and Altai, requiring near fixation in afrs. reflecting Archaic admixture
    v = pd.read_csv(vcf, sep='\t')
    v['afr_alt'] = 0
    for i in range(0, 200):
        v['afr_alt'] += v[f'afr_{i}']
    v['afr_ref'] = 200 - v['afr_alt']   #actually reflect the design of archaic admixture, so nearly fixed can mean that some polymorphic pos are not detected,
    v['freq_afr'] = v['afr_alt']/200
    v['nea_alt'] = v['nea_0'] + v['nea_1'] + v['nea_2']+v['nea_3'] + v['nea_out_0'] + v['nea_out_1']
    v['nea_ref'] = 6 - v['nea_alt']
    v['freq_nea'] = v['nea_alt']/6
    v['d3_alt'] = v['den3_0']+v['den3_1']
    v['d3_ref'] = 2 - v['d3_alt']
    v['freq_d3'] = v['d3_alt']/2
    v['freq_altai'] = v['nea_out_0'] / 2 + v['nea_out_1'] / 2
    v['max_diff_aa'] = v[['freq_afr','freq_altai','freq_d3']].max(axis=1) - v[['freq_afr','freq_altai','freq_d3']].min(axis=1) # must be ~0, =0.5, ~1
    """
    hapmap = pd.read_csv(hap_map, sep = '\t', header = None,
    names = ['chr', 'snp', 'gen_pos', 'physical_pos'])
    loop_hapmap = hapmap.iterrows()
    row = next(loop_hapmap)[1]
    row_1 = next(loop_hapmap)[1]
    pos_v = []
    for i in range(0, len(pos)):
        while (row_1.physical_pos <= pos[i]):
            row, row_1 = row_1, next(loop_hapmap)[1]
        pos_v.append(row.gen_pos + (pos[i] - row.physical_pos)*(row_1.gen_pos - row.gen_pos)/(row_1.physical_pos - row.physical_pos))
    v['map'] = pos_v
    """
    pos = v.pos.to_list()
    v['map'] = pos
    v['map'] = v['map'] * 1.25e-6
    v['ref'] = "A"
    v['alt'] = "G"
    v['anc_ref'] = 1
    v['anc_alt'] = 0
    # 1 / 200 = 0.005
    v_aa = v[((v['freq_afr'] > 0.994) | (v['freq_afr'] < 0.006)) & (v['max_diff_aa'] > 0.4)].copy()
    v_aa = v_aa[['chrom','pos','ref','alt','map','nea_ref','nea_alt','d3_ref','d3_alt','afr_ref','afr_alt', 'anc_ref','anc_alt']]
    v_aa = v_aa[v_aa['pos']>0]
    v_aa.to_csv(ref, index = False)

def make_ref4(vcf, ref): #4:  all segragating sites between Den3 and 3 Neas. Archaic admixture like.  but not used in the end for sums 
    v = pd.read_csv(vcf, sep='\t')
    v['afr_alt'] = 0
    for i in range(0, 200):
        v['afr_alt'] += v[f'afr_{i}']
    v['afr_ref'] = 200 - v['afr_alt']   #all segragating sites between den3 and 3 neas, but afr should be nearly fixed
    v['freq_afr'] = v['afr_alt']/200
    v = v[(v['freq_afr'] > 0.95) | (v['freq_afr'] < 0.05)]
    v['nea_alt'] = v['nea_0'] + v['nea_1'] + v['nea_2']+v['nea_3'] + v['nea_out_0'] + v['nea_out_1']
    v['nea_ref'] = 6 - v['nea_alt']
    v['freq_nea'] = v['nea_alt']/6
    v['d3_alt'] = v['den3_0']+v['den3_1']
    v['d3_ref'] = 2 - v['d3_alt']
    v['freq_d3'] = v['d3_alt']/2
    v['max_diff_d3'] = v[['freq_afr','freq_nea','freq_d3']].max(axis=1) - v[['freq_afr','freq_nea','freq_d3']].min(axis=1)
    """
    hapmap = pd.read_csv(hap_map, sep = '\t', header = None,
    names = ['chr', 'snp', 'gen_pos', 'physical_pos'])
    loop_hapmap = hapmap.iterrows()
    row = next(loop_hapmap)[1]
    row_1 = next(loop_hapmap)[1]
    pos_v = []
    for i in range(0, len(pos)):
        while (row_1.physical_pos <= pos[i]):
            row, row_1 = row_1, next(loop_hapmap)[1]
        pos_v.append(row.gen_pos + (pos[i] - row.physical_pos)*(row_1.gen_pos - row.gen_pos)/(row_1.physical_pos - row.physical_pos))
    v['map'] = pos_v
    """
    
    pos = v.pos.to_list()
    v['map'] = pos
    v['map'] = v['map'] * 1.25e-6
    v['ref'] = "A"
    v['alt'] = "G"
    v['anc_ref'] = 1
    v['anc_alt'] = 0
    vd3 = v[v['max_diff_d3'] >= 0.05].copy()
    vd3 = vd3[['chrom','pos','ref','alt','map','nea_ref','nea_alt','d3_ref','d3_alt','afr_ref','afr_alt', 'anc_ref','anc_alt']]
    vd3 = vd3[vd3['pos']>0]
    vd3.to_csv(ref, index = False)
def make_ref5(vcf, ref): # 5.  all segragating sites between Den3, Den25 and 3 Neas. Archaic admixture like. 
    v = pd.read_csv(vcf, sep='\t')
    v['afr_alt'] = 0
    for i in range(0, 200):
        v['afr_alt'] += v[f'afr_{i}']
    v['afr_ref'] = 200 - v['afr_alt']   #all segragating sites between den25 den3 3nea, but afr should be nearly fixed
    v['freq_afr'] = v['afr_alt']/200
    v = v[(v['freq_afr'] > 0.95) | (v['freq_afr'] < 0.05)]
    v['nea_alt'] = v['nea_0'] + v['nea_1'] + v['nea_2']+v['nea_3'] + v['nea_out_0'] + v['nea_out_1']
    v['nea_ref'] = 6 - v['nea_alt']
    v['freq_nea'] = v['nea_alt']/6
    v['den_alt'] = v['den3_0']+v['den3_1']+v['den25_0']+v['den25_1']
    v['den_ref'] = 4 - v['den_alt']
    v['freq_den'] = v['den_alt']/4
    v['max_diff_den'] = v[['freq_afr','freq_nea','freq_den']].max(axis=1) - v[['freq_afr','freq_nea','freq_den']].min(axis=1)
    """
    hapmap = pd.read_csv(hap_map, sep = '\t', header = None,
    names = ['chr', 'snp', 'gen_pos', 'physical_pos'])
    loop_hapmap = hapmap.iterrows()
    row = next(loop_hapmap)[1]
    row_1 = next(loop_hapmap)[1]
    pos_v = []
    for i in range(0, len(pos)):
        while (row_1.physical_pos <= pos[i]):
            row, row_1 = row_1, next(loop_hapmap)[1]
        pos_v.append(row.gen_pos + (pos[i] - row.physical_pos)*(row_1.gen_pos - row.gen_pos)/(row_1.physical_pos - row.physical_pos))
    v['map'] = pos_v
    """
    pos = v.pos.to_list()
    v['map'] = pos
    v['map'] = v['map'] * 1.25e-6
    v['ref'] = "A"
    v['alt'] = "G"
    v['anc_ref'] = 1
    v['anc_alt'] = 0
    vden = v[v['max_diff_den'] >= 0.05].copy()
    vden = vden[['chrom','pos','ref','alt','map','nea_ref','nea_alt','den_ref','den_alt','afr_ref','afr_alt', 'anc_ref','anc_alt']]
    vden = vden[vden['pos']>0]
    vden.to_csv(ref, index = False)
def make_ref6(vcf, ref):  # 6.  a fake Archaic panel.
    v = pd.read_csv(vcf, sep='\t')
    v['afr_alt'] = 0
    for i in range(0, 200):
        v['afr_alt'] += v[f'afr_{i}']
    v['afr_ref'] = 200 - v['afr_alt']   
    v['freq_afr'] = v['afr_alt']/200
    v = v[(v['freq_afr'] > 0.95) | (v['freq_afr'] < 0.05)]
    v['nea_alt'] = v['nea_0'] + v['nea_1'] + v['nea_2']+v['nea_3'] + v['nea_out_0'] + v['nea_out_1']
    v['nea_ref'] = 6 - v['nea_alt']
    v['freq_nea'] = v['nea_alt']/6
    v['den_alt'] = v['den3_0']+v['den3_1']+v['den25_0']+v['den25_1']
    v['den_ref'] = 4 - v['den_alt']
    v['freq_den'] = v['den_alt']/4
    v['max_diff_den'] = v[['freq_afr','freq_nea','freq_den']].max(axis=1) - v[['freq_afr','freq_nea','freq_den']].min(axis=1)
    v['N1'] = v['nea_0'] + v['nea_1']
    v['N2'] = v['nea_2'] + v['nea_3']
    v['N3'] = v['nea_out_0'] + v['nea_out_1']
    v['D1'] = v['den3_0'] + v['den3_1']
    v['D2'] = v['den25_0'] + v['den25_1']
    Archaic_max = np.max([v['N1'],v['N2'],v['N3'],v['D1'],v['D2']], axis = 0)
    Archaic_min = np.min([v['N1'],v['N2'],v['N3'],v['D1'],v['D2']], axis = 0)

    v['Archaic_alt'] = np.where(v['afr_alt'] < 100, Archaic_max, Archaic_min)
    v['Archaic_ref'] = 2 - v['Archaic_alt']
    """
    hapmap = pd.read_csv(hap_map, sep = '\t', header = None,
    names = ['chr', 'snp', 'gen_pos', 'physical_pos'])
    loop_hapmap = hapmap.iterrows()
    row = next(loop_hapmap)[1]
    row_1 = next(loop_hapmap)[1]
    pos_v = []
    for i in range(0, len(pos)):
        while (row_1.physical_pos <= pos[i]):
            row, row_1 = row_1, next(loop_hapmap)[1]
        pos_v.append(row.gen_pos + (pos[i] - row.physical_pos)*(row_1.gen_pos - row.gen_pos)/(row_1.physical_pos - row.physical_pos))
    v['map'] = pos_v
    """
    pos = v.pos.to_list()
    v['map'] = pos
    v['map'] = v['map'] * 1.25e-6
    v['ref'] = "A"
    v['alt'] = "G"
    v['anc_ref'] = 1
    v['anc_alt'] = 0
    vden = v[v['max_diff_den'] >= 0.05].copy()
    vden = vden[['chrom','pos','ref','alt','map','afr_ref','afr_alt', 'anc_ref','anc_alt','Archaic_ref','Archaic_alt']]
    vden = vden[vden['pos']>0]
    vden.to_csv(ref, index = False)
def make_infile(pop, ind, ref, vcf, out):
    hap1 = ind*2
    hap2 = ind*2+1
    vcf = pd.read_csv(vcf, sep='\t')  
    ref = pd.read_csv(ref)
    pos = ref.pos.to_list()
    vcf = vcf[vcf['pos'].isin(pos)]
    vcf = vcf[['chrom', 'pos', f'{pop}_{hap1}', f'{pop}_{hap2}']]
    vcf['talt'] = vcf[f'{pop}_{hap1}'] + vcf[f'{pop}_{hap2}']
    vcf['tref'] = 2 - vcf['talt']
    vcf = vcf[['chrom', 'pos', 'tref', 'talt']]
    vcf.to_csv(out, index = False)

def make_infile_hmmix(pop, ind, vcf, out):
    snp = pd.read_csv(vcf, sep = '\t')
    snp['afr_total'] = 0
    for i in range(0, 400):
        snp['afr_total'] += snp[f'afr_{i}']
    snp = snp[snp['afr_total'] == 0]
    snp = snp[snp['pos'] > 0]
    hap1 = ind*2
    hap2 = ind*2 + 1
    hap_1 = f"{pop}_{hap1}"
    hap_2 = f"{pop}_{hap2}"
    print(hap_1, hap_2)
    snp['target'] = snp[f'{hap_1}'] + snp[f'{hap_2}']
    snp = snp[snp['target'] > 0]
    
    with open(out, 'w') as f:
        for index, row_vcf in snp.iterrows():
            print(row_vcf.chrom, row_vcf.pos, sep = '\t', file = f)  

def make_infile_hmmix12(pop, ind, vcf, out):
    snp = pd.read_csv(vcf, sep = '\t')
    snp['afr_total'] = 0
    for i in range(0, 400):
        snp['afr_total'] += snp[f'afr_{i}']
    snp = snp[snp['afr_total'] == 0]
    snp = snp[snp['pos'] > 0]
    hap1 = ind*2
    hap2 = ind*2 + 1
    hap_1 = f"{pop}_{hap1}"
    hap_2 = f"{pop}_{hap2}"
    print(hap_1, hap_2)
    snp1 = snp[snp[f'{hap_1}'] > 0]
    snp2 = snp[snp[f'{hap_2}'] > 0] 
    
    with open(out, 'a') as f:
        for index, row_vcf in snp1.iterrows():
            print(row_vcf.chrom, row_vcf.pos, sep = '\t', file = f)  
    with open(out, 'a') as f:
        for index, row_vcf in snp2.iterrows():
            print(row_vcf.chrom+22, row_vcf.pos, sep = '\t', file = f)  
            

def anno(callfile, infile, vcf, out,chr, window_size = 1000, rec_rate = 1.25e-8, rec = None, min_match = 0.05, min_diff = 0.05, hap = False):
    c = pd.read_csv(callfile)
    chr = int(chr)
    chr_full = [chr]
    if hap == True:
        chr_full = [chr, chr+22]
    for chr_hap in chr_full:    #chr_hap can > 22
        c = c[c['chrom'] == chr_hap].sort_values('start').reset_index(drop=True)
        c['start'] = c['start'] * window_size
        c['end'] = c['end'] * window_size + window_size
        if len(c) == 0:
            return
        else:
            with open(out, 'a') as f:
                with open(out+".pos", 'a') as f1:
                    vcf = pd.read_csv(vcf, sep = '\s+')
                    vcf = vcf[['pos', 'nea_out_0', 'nea_out_1', 'nea_0', 'nea_1', 'nea_2', 'nea_3',
                            'den3_0', 'den3_1', 'den25_0', 'den25_1']]
                    in_obs = pd.read_csv(infile, sep = '\t', names = ['chrom', 'pos'])
                    in_obs = in_obs[in_obs['chrom'] == chr_hap].sort_values('pos').reset_index(drop=True)
                    if not rec is None:
                        rec_m = pd.read_csv(rec, sep = '\s+', names=['chrom', 'pos', 'gen_pos'])
                        rec_m = rec_m[rec_m['chrom'] == chr_hap].sort_values('pos').reset_index(drop=True)
                        print(rec_m)
                        loop_rec = rec_m.iterrows()
                        row_rec = next(loop_rec)[1]
                        row_rec_ = next(loop_rec)[1]
                    loop_c = c.iterrows()
                    while True:
                        try:
                            row = next(loop_c)[1]
                            if not rec is None:
                                while row_rec_.pos < row.start:
                                    row_rec, row_rec_ = row_rec_, next(loop_rec)[1]
                                start_gen = row_rec.gen_pos + (row.start - row_rec.pos)*(row_rec_.gen_pos - row_rec.gen_pos)/(row_rec_.pos - row_rec.pos) 
                                while row_rec_.pos < row.end:
                                    row_rec, row_rec_ = row_rec_, next(loop_rec)[1]
                                end_gen = row_rec.gen_pos + (row.end - row_rec.pos)*(row_rec_.gen_pos - row_rec.gen_pos)/(row_rec_.pos - row_rec.pos)
                            else:
                                start_gen = row.start * rec_rate * 100 #in cm
                                end_gen = row.end * rec_rate * 100
                            len_gen = end_gen - start_gen
                            in_obs_ = in_obs[(in_obs['pos'] > row.start) & (in_obs['pos'] <= row.end)]
                            if len(in_obs_) == 0:
                                print(row.chrom, row.start, row.end, row.end-row.start, len_gen, "0/0", "0/0", "0/0", "0/0", "0/0", "NA", "NA", 
                                    sep = '\t', file = f)
                            else:
                                nea_out = 0
                                nea1 = 0
                                nea2 = 0
                                den3 = 0
                                den25 = 0
                                LEN = len(in_obs_)
                                for pos_ in in_obs_.pos.tolist():
                                    nea_out_ = 0
                                    nea1_ = 0
                                    nea2_ = 0
                                    den3_ = 0
                                    den25_ = 0
                                    for a in ['nea_out_0', 'nea_out_1']:
                                        if vcf[vcf['pos'] == pos_][a].values[0] > 0:
                                            nea_out_ = 1
                                    nea_out += nea_out_
                                    for a in ['nea_0', 'nea_1']:
                                        if vcf[vcf['pos'] == pos_][a].values[0] > 0:
                                            nea1_ = 1
                                    nea1 += nea1_
                                    for a in ['nea_2', 'nea_3']:
                                        if vcf[vcf['pos'] == pos_][a].values[0] > 0:
                                            nea2_ = 1
                                    nea2 += nea2_
                                    for a in ['den3_0', 'den3_1']:
                                        if vcf[vcf['pos'] == pos_][a].values[0] > 0:
                                            den3_ = 1
                                    den3 += den3_
                                    for a in ['den25_0', 'den25_1']:
                                        if vcf[vcf['pos'] == pos_][a].values[0] > 0:
                                            den25_ = 1
                                    den25 += den25_
                                    if (nea_out_ + nea1_ + nea2_ > 0) and (den3_ + den25_ == 0):
                                        print(int(row.chrom),pos_,"Nea", sep = '\t', file = f1)
                                    elif (nea_out_ + nea1_ + nea2_ == 0) and (den3_ + den25_ > 0):
                                        print(int(row.chrom),pos_,"Den", sep = '\t', file = f1)
                                    elif (nea_out_ + nea1_ + nea2_ > 0) and (den3_ + den25_ > 0):
                                        print(int(row.chrom),pos_,"Shared", sep = '\t', file = f1)
                            samples = ['nea_out', 'nea_1', 'nea_2', 'den3', 'den25']
                            M = [round(nea_out/LEN, 3), 
                                    round(nea1/LEN, 3), 
                                    round(nea2/LEN, 3), 
                                    round(den3/LEN, 3), 
                                    round(den25/LEN, 3)]
                            L1 = np.argmax(np.array(M))
                            if M[L1] < min_match:  # little match
                                Source1 = "NA"
                                Source2 = "NA"
                            else:
                                M1 = M[L1]
                                Source1 = samples[L1]
                                M.pop(L1)
                                samples.pop(L1)
                                L2 = np.argmax(np.array(M))
                                if M[L2] < min_match:
                                    Source2 = "NA"
                                else:
                                    Source2 = samples[L2]
                            NEA_max = np.max([nea_out/LEN, nea1/LEN, nea2/LEN])
                            DEN_max = np.max([den3/LEN, den25/LEN])
                            if abs(NEA_max - DEN_max) < min_diff:
                                Source1 = "Ambiguous"
                                Source2 = "Ambiguous"
                            print(round(row.chrom), round(row.start), round(row.end), round(row.end-row.start), round(len_gen,6),
                                f"{nea_out}/{LEN}", f"{nea1}/{LEN}", f"{nea2}/{LEN}", f"{den3}/{LEN}", f"{den25}/{LEN}",
                                    Source1, Source2, sep = '\t', file = f)
                            print(round(row.chrom), round(row.start), round(row.end), round(row.end-row.start), round(len_gen,6),
                                f"{nea_out}/{LEN}", f"{nea1}/{LEN}", f"{nea2}/{LEN}", f"{den3}/{LEN}", f"{den25}/{LEN}",
                                    Source1, Source2, sep = '\t')
                            
                        except StopIteration:
                            break

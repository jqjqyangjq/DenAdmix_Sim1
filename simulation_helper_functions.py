import msprime
import gzip
import pandas as pd
from scipy.stats import binom, poisson
import numpy as np
from scipy import stats
import tskit
import numpy as np
import collections
from collections import defaultdict


def write_all_snps(file_name, ts, sample_names, chang_dict=None, chrom='1'):
    """
    This function writes genotype matrix to a gzipped file with user defined pop names.
    The population names must be given in order of msprime sample id (ascending order)
    """
    """Generates the individual sample names from the user definded pop names"""
    samples = []
    for p in enumerate(sample_names):
        n_pop = len(ts.get_samples(p[0]))
        samples.extend([f"{p[1]}_{i}" for i in range(n_pop)])
    print(samples)
    if chang_dict is not None:
        samples = [chang_dict.get(item,item)  for item in samples]
    print(samples)
    """Write genotype matrix to file."""
    with gzip.open(file_name, 'wt') as f:
        f.write("chrom\tpos\t")
        f.write("\t".join(samples))
        f.write("\n")
        for v in ts.variants():
            gt = '\t'.join(str(gt) for gt in v.genotypes)
            f.write(f"{chrom}\t{int(v.site.position)}\t{gt}\n")    
    
def write_eigenstrat(file_prefix, ts, sample_names, ploidy=2, fixed_ancestral = None, fixed_derived = None, polymorphic = None, target_meta=None, rec_rate=None, rec_map=None, chrom='1'):
    """
    This function takes a tree sequence and writes the three eigenstrat files .geno .snp and .ind
    The population names must be given in order of msprime sample id (ascending order)
    """
    with open('{}.ind'.format(file_prefix),'w') as outind:
        samples = []
        for p in enumerate(sample_names):
            if len(ts.get_samples(p[0])) >= int(ploidy):
                n_pop = ts.get_samples(p[0])[::int(ploidy)]
            else:
                n_pop = ts.get_samples(p[0])
            for i in n_pop:
                if target_meta is not None and i in np.array(target_meta.sim_index):
                    Ind = int((i - n_pop[0])/int(ploidy))
                    samples.extend([f"Ancient_{p[1]}_{Ind}_sampled_{int(target_meta[target_meta['sim_index']==i]['sample_target'])}\tU\tAncient_{p[1]}\n"])
                else:
                    Ind = int((i - n_pop[0])/int(ploidy))
                    samples.extend([f"{p[1]}_{Ind}\tU\t{p[1]}\n"])
        ind = ''.join(str(x) for x in samples)
        outind.write(ind)
        outind.close()
    if all(i is None for i in [fixed_ancestral,fixed_derived,polymorphic]):
        print('No ascertainment')

    if fixed_ancestral is not None:
        sample_ID = [x[0] for x in enumerate(sample_names) if x[1] in fixed_ancestral]
        fixed_ancestral_IDs = np.concatenate([ts.get_samples(i) for i in sample_ID])
        for i in fixed_ancestral:
            print('ascertaining {} to be fixed ancestral'.format(i))

    if fixed_derived is not None:
        sample_ID = [x[0] for x in enumerate(sample_names) if x[1] in fixed_derived]
        fixed_derived_IDs = np.concatenate([ts.get_samples(i) for i in sample_ID])
        for i in fixed_derived:
            print('ascertaining {} to be fixed derived'.format(i))

    if polymorphic is not None:
        sample_ID = [x[0] for x in enumerate(sample_names) if x[1] in polymorphic]
        polymorphic_IDs = np.concatenate([ts.get_samples(i) for i in sample_ID])
        for i in polymorphic:
            print('ascertaining {} to be fixed derived or polymorphic'.format(i))

    if all(i is None for i in [rec_map,rec_rate]):
        print("Using default recombination rate of 1e-8")
    elif rec_map is None:
        print("Using constant recombination rate of {}".format(rec_rate))
    else:
        print("Using recombination map")
    with gzip.open('{}.geno.gz'.format(file_prefix),'wt') as outgeno , gzip.open('{}.snp.gz'.format(file_prefix),'wt') as snpfile:
        for v in ts.variants():
            if fixed_ancestral is not None and np.sum(v.genotypes[fixed_ancestral_IDs]) != 0:
                pass
            elif fixed_derived is not None and np.sum(v.genotypes[fixed_derived_IDs]) != len(v.genotypes[fixed_derived_IDs]):
                pass
            elif polymorphic is not None and np.sum(v.genotypes[polymorphic_IDs]) == 0:
                pass
            else:
                pos = int(v.site.position)
                if v.site.id == 0:
                    pos_genetic = 0
                else:
                    if all(i is None for i in [rec_map,rec_rate]):

                        pos_genetic = pos * 1e-8
                    elif rec_map is None:
                        pos_genetic = pos * rec_rate
                    else:
                        pos_genetic = np.interp(pos, rec_map.pos, rec_map.map)
                snpfile.write('rs{0}\t{1}\t{2}\t{3}\n'.format(pos,str(chrom),pos_genetic,pos))
                gt = np.add.reduceat(v.genotypes, np.arange(len(v.genotypes))[::int(ploidy)])
                geno = ''.join(str(x) for x in gt)
                outgeno.write('{0}\n'.format(geno))
        snpfile.close()
        outgeno.close()

def write_snps_to_eigenstrat(file_prefix, snps, sample_names, rec_rate=None, rec_map=None, target=False):
    """
    This function takes a snps table created by the write_all_snps function and writes the three eigenstrat
    files .geno .snp and .ind.
    The population names must be given in order of msprime sample id (ascending order)
    """
    if target == True:
        print("reshaping for target")
        target_cols = [col for col in snps if col.startswith('target')]
        snps=snps[[col for col in snps if col not in target_cols] + target_cols]

    with open('{}.geno'.format(file_prefix),'w') as genofile:
        gt_h=snps.iloc[:,2:snps.shape[1]]
        gt_d = pd.DataFrame(np.add.reduceat(gt_h.values, np.arange(len(gt_h.columns))[::2], axis=1).astype(int))
        for gt_d_row in gt_d.values:
            geno = ''.join(str(x) for x in gt_d_row)
            genofile.write('{0}\n'.format(geno))

    with open('{}.snp'.format(file_prefix),'w') as snpfile:
        D_snp= dict()
        D_snp["rs"] = ["rs" + str(i) for i in snps.pos]
        D_snp["chrom"] = snps.chrom
        if all(i is None for i in [rec_map,rec_rate]):
            D_snp['genetic_pos'] = snps.pos * 1e-8
        elif rec_map is None:
            D_snp['genetic_pos'] = snps.pos * rec_rate
        else:
            D_snp['genetic_pos'] = np.interp(snps.pos, rec_map.pos, rec_map.map)
        D_snp["pos"] = snps.pos
        snp = pd.DataFrame.from_dict(D_snp)
        snp.loc[:0,'genetic_pos'] = 0
        for snp_row in snp.values:
            snp_str = '\t'.join(str(x) for x in snp_row)
            snpfile.write('{0}\n'.format(snp_str))

    with open('{}.ind'.format(file_prefix),'w') as indfile:
        samples = []
        for p in sample_names:
            pop_x_cols = [col for col in snps.columns if col.startswith(f"{p}")]
            n_pop = int(len(pop_x_cols))
            if n_pop < 2:
               n_pop = n_pop
            else:
                n_pop = int(len(pop_x_cols)/2)
            samples.extend([f"{p}_{i}\tU\t{p}\n" for i in range(n_pop)])
        ind = ''.join(str(x) for x in samples)
        indfile.write(ind)
        indfile.close()


def get_introgressed_segments(ts,mixed_pop_ID,introgressing_pop_ID,split_time,file_name=None,chrom='1'):
    """
    This function takes a tree sequence object (while running the msprim sim, record_migrations must have been set to True!!!),
    and returns all (merged) segments introgressed from the introgressing_pop found in sampled individuals of the mixed_pop.
    If a file_name is given it writes a file with individual id, chrom, start, end of intriogressed segment.
    """
    Testpopulation = ts.get_samples(mixed_pop_ID)

    de_seg = {i: [] for i in Testpopulation}
    print("new function used")
    """
    This is a helper function and merges neighboring segments in an individuals genome to one.
    """
    def combine_segs(segs, get_segs = False):
        merged = np.empty([0, 2])
        if len(segs) == 0:
            if get_segs:
                return([])
            else:
                return(0)
        sorted_segs = segs[np.argsort(segs[:, 0]), :]
        for higher in sorted_segs:
            if len(merged) == 0:
                merged = np.vstack([merged, higher])
            else:
                lower = merged[-1, :]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1, :] = (lower[0], upper_bound)
                else:
                    merged = np.vstack([merged, higher])
        if get_segs:
            return(merged)
        else:
            return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)
    """
    Looping through the tree sequence tables of migrants. If a migrant is from the
    introgressing_pop found in the mixed_pop, the sub interval of the tree is taken.
    From this interval all leaves from gthe migrant node (i.e. descendedts of the migrant)
    and are appended to our introgressed seg object if they are in any sampled ind of the mixed_pop.
    """
    mig_events = []
    for mr in ts.migrations():
        if mr.source == mixed_pop_ID and mr.dest == introgressing_pop_ID:
            if mr.time == split_time: # migraton caused by merge
                pass
            mig_events.append(mr)
    mig_events.sort(key = lambda x: x.left)
    mig_events.sort(key = lambda x: x.right)
    for mr in mig_events:
            for tree in ts.trees(leaf_lists=True):
                if mr.left > tree.get_interval()[0]:
                    continue
                if mr.right <= tree.get_interval()[0]:
                    break
                for l in tree.leaves(mr.node):
                    if l in Testpopulation:
                        de_seg[l].append(tree.get_interval())

    """
    Merging all length neighboring intervals of introgressed segments using the combine_segs
    helper function.
    """
    true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]
    """
    Writing the introgressed segs of a individual (haploid genome id) with chrom, and start and stop in bp to a file.
    """
    if file_name is not None:
        with open('{}'.format(file_name),'w') as out_true:
            for haplotype, archaic_segments in enumerate(true_de_segs):
                for archaic_segment in archaic_segments:
                    out_true.write('{}\tchr{}\t{}\t{}\n'.format(haplotype,str(chrom), int(archaic_segment[0]), int(archaic_segment[1])))
    else:
        return(true_de_segs)

def admixfrog_sample_in(sample_reads,min_len,bin_len):
    max_len= sample_reads.len
    sample_reads.loc[sample_reads['deam'] == 0, 'deam'] = "deam"
    sample_reads.loc[sample_reads['deam'] == -1, 'deam'] = "nodeam"
    i = 1
    while (min_len + bin_len*i) <= max(max_len + bin_len):
        sample_reads.loc[(sample_reads['len'] >= (min_len + bin_len*(i-1))) & (sample_reads['len'] <= (min_len + bin_len*(i))), 'len'] = i-1
        i=i+1
    admixfrog_sample_in = sample_reads[['chrom', 'pos']].copy()
    admixfrog_sample_in['lib'] = [p1 + '_' + p2 + '_' + p3 for p1, p2, p3 in zip(sample_reads['lib'], sample_reads['len'].astype(str), sample_reads['deam'])]
    admixfrog_sample_in['tref'] = sample_reads['tref']
    admixfrog_sample_in['talt'] = sample_reads['talt']
    admixfrog_sample_in['tdeam'] = np.repeat(0, admixfrog_sample_in.shape[0], axis=0)
    admixfrog_sample_in['tother'] = np.repeat(0, admixfrog_sample_in.shape[0], axis=0)
    admixfrog_sample_in=admixfrog_sample_in.groupby(['chrom','pos','lib']).sum().reset_index()
    return admixfrog_sample_in

def get_all_pop_allele_counts(snps, sample_names, rec_rate=None, rec_map=None, chrom='1'):
    """
    This function creates a panda table with allele counts per SNP for each population from a
    snps file table (created by write_all_snps function)
    The population names must be given in order of msprime sample id (ascending order)
    """

    """"Reading in data from csv file containing the variants (snp file)"""
    if 'pos' not in snps:
        snps['pos'] = snps.index.astype(int)
    snps.pos = snps.pos.astype(int)
    snps2 = snps.drop_duplicates(subset=['pos'])


    D = dict()
    D['chrom'] = chrom
    D['pos'] = snps2.pos.astype(int)
    D['map'] = snps2.pos * 1e-6
    D['anc'] = 'A'
    D['der'] = 'G'

    for p in sample_names:
        pop_x_cols = [col for col in snps2.columns if col.startswith(f"{p}")]
        n_pop_x = len(pop_x_cols)
        """
        Determining the count of derived (1) and ancestral (0) alleles per population:
        the number of haploid individuals minus the number of derived alleles determines
        the number of ancestral alleles per population
        """
        D[f"{p}_der"] = np.sum(snps2[pop_x_cols], 1)
        D[f"{p}_anc"] = n_pop_x - D[f"{p}_der"]


    """
    Putting it into a data frame and check if number of input and poutput snps is the same
    and dropping potential duplicates
    """
    allele_table = pd.DataFrame.from_dict(D)
    assert allele_table.shape[0] == snps2.shape[0]
    assert allele_table.drop_duplicates(subset=['chrom', 'pos']).shape[0] == allele_table.shape[0]
    return(allele_table)

def ascertain_all_allele_counts(allele_table, fixed_der=None,fixed_anc=None,polymorphic=None):
    """
    This function ascertaines the allele table object (generated by the get_all_pop_allele_counts function)
    by user defined populations and outputs the filtered table.
    The populations ascertainment can be done using different categories:
    either be fixed for the derived/ancestral allele or polymorphic or both.
    Multiple populations can be defined per category.
    Optionally a snp file can be provided and is filtered in the same way and outputted.
    """

    """First the variants are filtered out which are not fixed"""
    if all(v is None for v in (fixed_der,fixed_anc)):
        filter_1 = (allele_table.pos == allele_table.pos)
    else:
        if fixed_der == None:
            fixed_anc_der = [f"{i}_der" for i in fixed_anc]
            fixed_anc_anc = [f"{i}_anc" for i in fixed_anc]
            print([f"Ascertaining for fixed ancestral sites in {i}" for i in fixed_anc])
            filter_1 = (allele_table[fixed_anc_der].sum(axis=1) == 0)
        if fixed_anc == None:
            fixed_der_der = [f"{i}_der" for i in fixed_der]
            fixed_der_anc = [f"{i}_anc" for i in fixed_der]
            print([f"Ascertaining for fixed derived sites in {i}" for i in fixed_der])
            filter_1 = (allele_table[fixed_der_anc].sum(axis=1) == 0)
        elif all(v is not None for v in (fixed_der,fixed_anc)):
            fixed_anc_der = [f"{i}_der" for i in fixed_anc]
            fixed_anc_anc = [f"{i}_anc" for i in fixed_anc]
            fixed_der_der = [f"{i}_der" for i in fixed_der]
            fixed_der_anc = [f"{i}_anc" for i in fixed_der]
            print([f"Ascertaining for fixed derived sites in {i}" for i in fixed_der])
            print([f"Ascertaining for fixed ancestral sites in {i}" for i in fixed_anc])
            filter_1 = (allele_table[fixed_anc_der].sum(axis=1) + allele_table[fixed_der_anc].sum(axis=1) == 0) | (allele_table[fixed_der_der].sum(axis=1) + allele_table[fixed_anc_anc].sum(axis=1) == 0)
    allele_table_asc_1 = allele_table[filter_1]

    """Second, filtering for polymorphic variants"""
    if polymorphic == None:
        filter_2 = (allele_table_asc_1.pos == allele_table_asc_1.pos)
    else:
        polymorphic_der = [f"{i}_der" for i in polymorphic]
        polymorphic_anc = [f"{i}_anc" for i in polymorphic]
        filter_2 = (allele_table_asc_1[polymorphic_der].sum(axis=1) !=0) & (allele_table_asc_1[polymorphic_anc].sum(axis=1) !=0)
        print([f"Ascertaining for polymorphic sites in {i}" for i in polymorphic])
    allele_table_asc_2 = allele_table_asc_1[filter_2]
    filter_all=allele_table_asc_2.index
    if allele_table_asc_2.shape[0] == allele_table.shape[0]:
        print("Data not ascertained")

    return(allele_table_asc_2,filter_all)

def coverage_simulation(ids, allele_table ,snp_table, coverage, contaminant_pop, contamination_percent, libs):
    """
    This function takes 1 or more diploid individuals specified by their haploid ids from the snp file,
    the allele table file and the snp file to simulates coverage for each SNP and the composition of
    endogenous and contaminant read per SNP observed. The user defines the source of contamination by
    its pop lable. The coverage, contamination_percent and libs (string) parameters must be given as arrays.
    """
    if isinstance(len(ids)/2, float) is not True:
        raise ValueError("Number of ids must be two or multiples of two")
    S = []
    for cov, cont, lib in zip(coverage, contamination_percent, libs):
        """Print input parameters to stdout for user to check"""
        print(f'Lib:{lib}\tCov:{cov}\tcont:{contaminant_pop}\tcont_%:{cont}', end="\t")
        print(f'Lib:{lib}\tCov:{cov}\tcont:{contaminant_pop}\tcont_%:{cont}', end="\t")

        """Determining the true allele frequency:
        Counts the msprime simulated alternative/derived alleles of a single deploid sample.
        """
        data = allele_table[['chrom', 'pos']].copy()
        data['true_der'] = np.sum(snp_table[ids],1)
        data['true_anc'] = 2 - data['true_der']
        data['lib'] = lib

        """Sampling the number of reads per SNP:
        The number of endogenous and contaminant reads are sampled for a poisson distribution with
        lambda (mean read count) being the user defined coverage (cov).
        Number of endogenous reads (cov_real) are simulated with lambda =  cov * 1 - proportion of contamination.
        Number of contaminat reads (cov_cont) are simulated with lambda =  cov * proportion of contamination.
        """
        cov_real = poisson.rvs(cov * (1. - cont), size=data.shape[0])
        cov_cont = poisson.rvs(cov * cont, size=data.shape[0])

        """Sampling the alternative/reference allele per read:
        The frequency of the derived allele (p) from the true data is used to sample
        the amount of derived reads from all endogenous reads (cov_real) per SNP using a binomial distribution
        with mean = p per read.
        The frequency of the contaminant derived allele (p_cont) taken from the contaminant population allele frequency
        is used to sample the amount of derived reads from all contaminant reads (cov_cont) per SNP using a binomial distribution
        with mean = p_cont per read.
        """
        p = data['true_der'] / (data['true_anc'] + data['true_der'])
        p_cont = allele_table[f"{contaminant_pop}_der"] / (allele_table[f"{contaminant_pop}_anc"] + allele_table[f"{contaminant_pop}_der"])
        data['ralt'] = binom.rvs(cov_real, p)
        data['rref'] = cov_real - data['ralt']
        data['calt'] = binom.rvs(cov_cont, p_cont)
        data['cref'] = cov_cont - data['calt']

        """Observed allele counts:
        Are simply the sum of contaminant + endogenous reads having the derived or ancestral allele
        """
        data['talt'] = data.ralt + data.calt
        data['tref'] = data.rref + data.cref


        print(f"alt:\t{lib}\t{np.mean(data['ralt']):.3f}\t{np.mean(data['calt']):.3f}\t{np.mean(data['talt']):.3f}")
        print(f"ref:\t{lib}\t{np.mean(data['rref']):.3f}\t{np.mean(data['cref']):.3f}\t{np.mean(data['tref']):.3f}")

        """Excludes all SNPs with 0 reads"""
        data = data[data.tref+data.talt>0]
        S.append(data)
    """Write the coverage file to csv"""
    data = pd.concat(S).sort_values(['chrom', 'pos', 'lib'])

    return(data)

def admixfrog_input_gt(ids, allele_table ,snp_table):
    """
    This function takes 1 or more diploid individuals specified by their haploid ids from the snp file,
    the allele table file and the snp file to simulates coverage for each SNP and the composition of
    endogenous and contaminant read per SNP observed. The user defines the source of contamination by
    its pop lable. The coverage, contamination_percent and libs (string) parameters must be given as arrays.
    """
    if isinstance(len(ids)/2, float) is not True:
        raise ValueError("Number of ids must be two or multiples of two")
    S = []

    data = allele_table[['chrom', 'pos']].copy()
    data['talt'] = np.sum(snp_table[ids],1)
    data['tref'] = 2 - data['talt']

    """Excludes all SNPs with 0 reads"""
    data = data[data.tref+data.talt>0]
    S.append(data)
    """Write the coverage file to csv"""
    data = pd.concat(S).sort_values(['chrom', 'pos'])

    return(data)

def deamination_and_read_len_simulation(data,endo_read_params,cont_read_params,read_length_cutoff):
    """
    This function takes the output of the coverage_simulation function and simulates read length
    and deamination for the endogenpous and contaminant reads. The nuber of libraries is determined
    from the input data. The endo_read_params (array of 3), cont_read_params (array of 3) and read_length_cutoff (integer)
    parameters must be provided as arrays either length one or samle length as libraries in the input data.
    """
    print(endo_read_params)
    print(cont_read_params)
    print(read_length_cutoff)
    S_reads = []
    """If the number of libraries and endo_read_params is not the same the same parameters will be used for each simulation"""
    for l in enumerate(np.unique(data.lib)):
        if (len(np.unique(data.lib)) == len(endo_read_params)):
            erp = endo_read_params[l[0]]
            crp = cont_read_params[l[0]]
            rlc = read_length_cutoff[l[0]]
        else:
            erp = endo_read_params[0]
            crp = cont_read_params[0]
            rlc = read_length_cutoff[0]

        """Reading in the parameters:
        Each set of reads (endogenous and contaminant) have 3 parameters:
        p_d_x = probability of diamination per read
        shape_l_x = shape parameter for the gamma distributions the read length are sampled from
        scale_l_x = scale parameter for the gamma distributions the read length are sampled from
        """
        p_d_endo , shape_l_endo , scale_l_endo = erp
        p_d_cont , shape_l_cont , scale_l_cont = crp

        print(f"Simulating endogenous reads of {l[1]} with prob. deam.:{p_d_endo}, min cutoff: {rlc} and mean length: {rlc[0]+shape_l_endo*scale_l_endo}")
        print(f"Simulating contaminant reads of {l[1]} with prob. deam.:{p_d_cont}, min cutoff: {rlc} and mean length: {rlc[0]+shape_l_cont*scale_l_cont}")
        data_x = data[data.lib == l[1]]


        """Demination simulation:
        A binomial distribution is used to determin if a read is deaminated or not.
        This is done seperatly for endogenous and contaminant reads,
        according to the user defined probability of being deaminated
        """
        deam_endo=binom.rvs(1, p_d_endo, size= (sum(data_x.rref)+sum(data_x.ralt)))-1
        deam_cont=binom.rvs(1, p_d_cont, size= (sum(data_x.cref)+sum(data_x.calt)))-1

        """Read length simulation:
        A gamma distribution is used to simulate read length. The read length is always the user defined
        read length cutoff + the outcome of the sampling from a gamme distribution using different parameters
        for endogenous and contaminat reads.
        """
        len_endo=(np.random.gamma(shape=shape_l_endo,scale=scale_l_endo,size=(sum(data_x.rref)+
                    sum(data_x.ralt)))+rlc).astype(int)
        len_cont=(np.random.gamma(shape=shape_l_cont,scale=scale_l_cont,size=(sum(data_x.cref)+
                    sum(data_x.calt)))+rlc).astype(int)

        """Putting it together in a dataframe:
        Enodgenous reads. Each read gets an entry in the data frame.
        """
        sim_reads_endo = pd.DataFrame({'chrom': np.repeat(data_x.chrom, (data_x.rref+data_x.ralt))})
        sim_reads_endo['pos'] = np.concatenate((np.repeat(data_x.pos, (data_x.rref)), np.repeat(data_x.pos, (data_x.ralt))), axis=None)
        sim_reads_endo['tref'] = np.concatenate((np.repeat(1, sum(data_x.rref)), np.repeat(0, sum(data_x.ralt))), axis=None)
        sim_reads_endo['talt'] = np.concatenate((np.repeat(0, sum(data_x.rref)), np.repeat(1, sum(data_x.ralt))), axis=None)
        sim_reads_endo['lib'] = np.repeat(data_x.lib, (data_x.rref+data_x.ralt))
        sim_reads_endo['len'] = len_endo
        sim_reads_endo['deam'] = deam_endo

        """Putting it together in a dataframe:
        Contaminant reads. Each read gets an entry in the data frame.
        """
        sim_reads_cont = pd.DataFrame({'chrom': np.repeat(data_x.chrom, (data_x.cref+data_x.calt))})
        sim_reads_cont['pos'] = np.concatenate((np.repeat(data_x.pos, (data_x.cref)), np.repeat(data_x.pos, (data_x.calt))), axis=None)
        sim_reads_cont['tref'] = np.concatenate((np.repeat(1, sum(data_x.cref)), np.repeat(0, sum(data_x.calt))), axis=None)
        sim_reads_cont['talt'] = np.concatenate((np.repeat(0, sum(data_x.cref)), np.repeat(1, sum(data_x.calt))), axis=None)
        sim_reads_cont['lib'] = np.repeat(data_x.lib, (data_x.cref+data_x.calt))
        sim_reads_cont['len'] = len_cont
        sim_reads_cont['deam'] = deam_cont

        """Putting it all together in a dataframe:
        Including dmgsites all 0 (I forgott what that is Ben P. knows)
        """
        sim_reads_result = pd.concat([sim_reads_endo,sim_reads_cont])
        sim_reads_result['dmgsite']=np.repeat(0,(sum(data_x.tref)+sum(data_x.talt)))
        sim_reads_result=sim_reads_result.sort_values(by=['pos'],ignore_index=True)
        S_reads.append(sim_reads_result)
    sim_reads_result_all = pd.concat(S_reads).sort_values(['chrom', 'pos', 'lib'])

    return sim_reads_result_all

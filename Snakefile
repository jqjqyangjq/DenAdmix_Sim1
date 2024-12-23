from curses import window
from inspect import isdatadescriptor
import os
from pickle import TRUE
from random import random, sample
from sre_parse import expand_template
import sys
from textwrap import indent
from turtle import pos
from wsgiref.validate import InputWrapper
import simulation_helper_functions as ts_help
import sim_NEA_DEN_EMH as sim_demo
import msprime
import tskit
import pandas as pd
import numpy as np
from collections import defaultdict, namedtuple, Counter
from os import path
import gzip
from scipy.stats import binom, poisson, bernoulli
import random
from get_ref import make_ref1,make_ref2,make_ref3, make_ref4,make_ref5,make_ref6, make_infile, make_infile_hmmix,make_infile_hmmix12, anno

"""
from make_ref:
1:  all segragating sites between Den3 and 3 Neas
2:  all segragating sites between Den3, Den25 and 3 Neas
3:  all segragating sites between Den3, and Altai, requiring near fixation in afrs. reflecting Archaic admixture
4:  all segragating sites between Den3 and 3 Neas. Archaic admixture like.  but not used in the end for sums 
5.  all segragating sites between Den3, Den25 and 3 Neas. Archaic admixture like. 
6.  a fake Archaic panel.
"""
ref_state = {
    "1":"afr nea d3",
    "2":"afr nea den",
    "3":"afr nea d3",
    "4":"afr nea d3",
    "5":"afr nea den",
    "6":"afr Archaic"
}

ref_prefix = {
    "1":"d3",
    "2":"den",
    "3":"d_AA",
    "4":"d3_AA",
    "5":"den_AA"
}
pop=["N", "S"]
chr = list(range(1,23))
window_size = [1000, 5000]
genome_cov = [0, 0.5, 1]
split_time = [90, 150, 210, 270, 330, 390]
frequency = [0.01]
individual = [0,2,4]
contamination = [0.02]
#coverage = [0.2,0.5,1,5,10,30]
coverage = [0.5,1,2,5,10]
penalty = [0.2, 0.5]
frog_out = ["bin.xz", "cont.xz", "log", "res2.xz", "res.xz", "snp.xz", "pars.yaml"]
plot_out = ["tpfpfsrc.png", "sensi.png"]
chrom_len = 80000000
#mut_rate_c = [mut_rate, mut_rate, mut_rate]
#cov_ll = [0.1]
divergence = [90, 150, 210, 270, 330, 390]
missing_mask = "mask_missing.chr1_20" #chr1 mask from map35_100 + 1kg accessibility


rule run_sim:
    priority: 1
    input:
        #rec=get_rec_true
        rec_file = "/mnt/diversity/jiaqi/hmm/hmm_extend_sim/map/100M_r_map_chr1.txt"
    priority : 1
    output:
        ts="sim/D_{divergence}/chrom_{chrom}_ts.trees",
        target_sample_list="sim/D_{divergence}/chrom_{chrom}_list.txt",
        #mut_rates_file = "sim/D_{divergence}/chrom_{chrom}_mut_rates.txt"
    run:
        random_seed_sim = int(wildcards.chrom) + 419
        r_rate_map = msprime.RateMap.read_hapmap(input.rec_file)
        sample_list = {
              "nea_out": {"t_sample": [120], "ind": 1, "group": "NEA"},
              "nea":  {"t_sample": [60,80], "ind": 1, "group": "NEA"},
              "Intro_nea": {"t_sample": [40], "ind": 1, "group": "NEA"},
              "den3":  {"t_sample": [80], "ind": 1, "group": "DEN"},
              "den25":  {"t_sample": [200], "ind": 1, "group": "DEN"},
              "ea": {"t_sample": [40, 40, 40, 40, 40, 20, 20, 20, 20, 20, 0, 0, 0, 0, 0], "ind": 1, "group": "EA"},
              "afr": {"t_sample": [0], "ind": 200, "group": "AFR"},
              "Intro_den_S": {"t_sample": [40], "ind": 1, "group": "DEN"}
          }
        samples = sim_demo.define_samples(sample_list)
        Demography= sim_demo.demo_archaic_introgression(N_p = 0.025, D_p = 0.0025, D_s=int(wildcards.divergence))
        recomb_rate = 1.25e-8
        
        seq_len = chrom_len
        print("Simulation started")
        simulation_result_tree = msprime.sim_ancestry(
        random_seed=random_seed_sim,
        samples=samples,
        sequence_length=seq_len,
        #recombination_rate = r_rate_map, varing
        recombination_rate=1.25e-8, # constant,
        #recombination_rate=recomb_rate,
        demography=Demography,
        record_migrations=True,
        ploidy =  2)   ##################

        population = []
        individual = []
        group = []
        for i, pop in enumerate(sample_list):
            for n in range(len(sample_list[pop]["t_sample"])):
                population.extend([f"{pop}" for i in range(sample_list[pop]["ind"])])
                individual.extend([f"{pop}_{ind}" for ind in range(sample_list[pop]["ind"])])
                group.extend(["{}".format(sample_list[pop]["group"]) for i in range(sample_list[pop]["ind"])])
        target_df = pd.DataFrame({'group': group, 'population': population, 'individual': individual})

        # Add mutations on tree sequence
        simulation_result = msprime.sim_mutations(simulation_result_tree, rate=1.25e-8,model=msprime.BinaryMutationModel(),random_seed=random_seed_sim)
        print("Simulation finished")
        target_df.to_csv('{}'.format(output.target_sample_list), float_format="%.5f", index=False)
        #simulation_result.dump('{}'.format(output[0]))
        simulation_result.dump('{}'.format(output.ts))
        print(samples)

rule get_true_segments:
    input:
        ts="sim/D_{divergence}/chrom_{chrom}_ts.trees"
    output:
        TrueSeg_den_S="sim/D_{divergence}/chrom_{chrom}_den_S.bed",
        TrueSeg_nea_S="sim/D_{divergence}/chrom_{chrom}_nea_S.bed"
    run:
        ts=tskit.load('{}'.format(input.ts))
        ts_help.get_introgressed_segments(ts,5,7,600000//29,file_name=output.TrueSeg_den_S,chrom=wildcards.chrom)
        ts_help.get_introgressed_segments(ts,5,2,600000//29,file_name=output.TrueSeg_nea_S,chrom=wildcards.chrom)

rule merge_true_segments_S:
    input:
        seg_den = expand("sim/D_{{divergence}}/chrom_{chrom}_den_S.bed", chrom = chr),
        seg_nea = expand("sim/D_{{divergence}}/chrom_{chrom}_nea_S.bed", chrom = chr)
    output:
        seg_all_chr = "sim/D_{divergence}/chrom_true_introgressed_all.bed.xz"
    run:
        seg_merged = pd.DataFrame(columns=["hap", "chr", "start", "end", "source"])
        for i in range(len(input.seg_den)):
            den = pd.read_csv('{}'.format(input.seg_den[i]), sep='\t', names = ["hap", "chr","start", "end"])
            den["source"] = "den_S"
            seg_merged = pd.concat([seg_merged, den])
        for i in range(len(input.seg_nea)):
            nea = pd.read_csv('{}'.format(input.seg_nea[i]), sep='\t', names = ["hap", "chr","start", "end"])
            nea["source"] = "nea"
            seg_merged = pd.concat([seg_merged, nea])
        seg_merged.to_csv('{}'.format(output.seg_all_chr),index=False)
rule get_snp_table:
    input:
        ts="sim/D_{divergence}/chrom_{chrom}_ts.trees",
        target_sample_list="sim/D_{divergence}/chrom_{chrom}_list.txt"
    output:
        data="sim/D_{divergence}/chrom_{chrom}.snps.tsv.gz"
    run:
        file_name = '{}'.format(output.data)
        ts=tskit.load('{}'.format(input.ts))
        #target_sample_list = pd.read_csv('{}'.format(input.target_sample_list), sep=",")
        # Choose sample names (in correct order)
        sample_names = ['nea_out',
                        'nea',
                        'Intro_nea',
                        "den3",
                        "den25",
                        "ea",
                        "afr",
                        "Intro_den_S"]
        chrom = wildcards.chrom
        ts_help.write_all_snps(
                        file_name=file_name,
                        ts=ts,
                        sample_names=sample_names,
                        chrom=chrom)

rule write_all_snp_table:
    input:
        expand("sim/D_{divergence}/chrom_{chrom}.snps.tsv.gz", chrom = chr, divergence = divergence)

rule output_ref:
    input:
        vcf = "sim/D_{divergence}/chrom_{chrom}.snps.tsv.gz",
    output:
        ref_file = "sim/D_{divergence}/chrom_{chrom}_ref{ref_type}.xz"
    wildcard_constraints:
        ref_type = "|".join(['1','2','3','4','5', '6'])
    run:
        if wildcards.ref_type == "1":
            make_ref1(vcf = input.vcf, ref = output.ref_file)
        elif wildcards.ref_type == "2":
            make_ref2(vcf = input.vcf, ref = output.ref_file)
        elif wildcards.ref_type == "3":
            make_ref3(vcf = input.vcf, ref = output.ref_file)
        elif wildcards.ref_type == "4":
            make_ref4(vcf = input.vcf, ref = output.ref_file)
        elif wildcards.ref_type == "5":
            make_ref5(vcf = input.vcf, ref = output.ref_file)
        elif wildcards.ref_type == "6":
            make_ref6(vcf = input.vcf, ref = output.ref_file)

rule output_ref_missing:
    input:
        ref_file = "sim/D_{divergence}/chrom_{chrom}_ref{ref_type}.xz",
        mask = missing_mask
    output:
        ref_file = "sim/D_{divergence}/chrom_{chrom}_ref{ref_type}.missing.xz"
    wildcard_constraints:
        ref_type = "|".join(['1','2','3','4','5','6'])
    run:
        s = "xzcat {input.ref_file} | sed -n '1p' | xz > {output.ref_file}"
        shell(s)
        s = "bedtools intersect -a <(xzcat {input.ref_file} | sed 1d | column -t -s , | awk -v OFS='\t' '{{print $1,$2-1,$2,$0}}' ) -b {input.mask} | cut -f 4- | sed -E 's/[[:space:]]+/,/g' | xz >> {output.ref_file}"
        shell(s)

#ref1  Den3 included, 3 neas included
#ref2  Den3 Den25 included, 3 neas included
#ref3  Den3 and Altai included, working as 1740K.sites
rule all_ref:
    input:
        expand("sim/D_{{divergence}}/chrom_{chrom}_ref{{ref_type}}.xz", chrom = chr)
    output:
        ref1 = "sim/D_{divergence}/chrAll_ref{ref_type}.xz"
    wildcard_constraints:
        ref_type = "|".join(['1','2','3','4','5','6'])
    run:
        s = "xzcat sim/D_{wildcards.divergence}/chrom_1_ref{wildcards.ref_type}.xz  > sim/D_{wildcards.divergence}/chrAll_ref{wildcards.ref_type}"
        shell(s)
        s = " for i in {{2..22}}; do xzcat  sim/D_{wildcards.divergence}/chrom_${{i}}_ref{wildcards.ref_type}.xz| sed 1d >> sim/D_{wildcards.divergence}/chrAll_ref{wildcards.ref_type}; done"
        shell(s)
        s = "xz sim/D_{wildcards.divergence}/chrAll_ref{wildcards.ref_type}"
        shell(s)
        
rule all_ref_missing:
    input:
        expand("sim/D_{{divergence}}/chrom_{chrom}_ref{{ref_type}}.missing.xz", chrom = chr)
    output:
        ref1 = "sim/D_{divergence}/chrAll_ref{ref_type}.missing.xz"
    wildcard_constraints:
        ref_type = "|".join(['1','2','3','4','5','6'])
    run:
        s = "xzcat sim/D_{wildcards.divergence}/chrom_1_ref{wildcards.ref_type}.missing.xz  > sim/D_{wildcards.divergence}/chrAll_ref{wildcards.ref_type}.missing"
        shell(s)
        s = " for i in {{2..22}}; do xzcat  sim/D_{wildcards.divergence}/chrom_${{i}}_ref{wildcards.ref_type}.missing.xz| sed 1d >> sim/D_{wildcards.divergence}/chrAll_ref{wildcards.ref_type}.missing; done"
        shell(s)
        s = "xz sim/D_{wildcards.divergence}/chrAll_ref{wildcards.ref_type}.missing"
        shell(s)
    
rule concat_all_ref:
    input:
        expand("sim/D_{divergence}/chrAll_ref{ref_type}.xz", chrom = chr,ref_type = ['1', '2', '3','4','5'], divergence = divergence),
        expand("sim/D_{divergence}/chrAll_ref{ref_type}.missing.xz", chrom = chr,ref_type = ['1', '2', '3','4','5'], divergence = divergence),
        

rule infile:
    input:
        vcf = "sim/D_{divergence}/chrom_{chrom}.snps.tsv.gz",
        ref1 = "sim/D_{divergence}/chrom_{chrom}_ref{ref_type}.xz"
    output:
        vcf1 = "sim/D_{divergence}/infile/{pop}.ind_{ind}.chr_{chrom}.in{ref_type}.xz"
    run:
        make_infile(pop=wildcards.pop, ind=int(wildcards.ind), ref=input.ref1, vcf=input.vcf, out=output.vcf1)
        
rule infile_missing:
    input:
        vcf1 = "sim/D_{divergence}/infile/{pop}.ind_{ind}.chr_{chrom}.in{ref_type}.xz",
        mask = missing_mask
    output:
        vcf1 = "sim/D_{divergence}/infile_missing/{pop}.ind_{ind}.chr_{chrom}.in{ref_type}.xz"
    run:
        s = "xzcat {input.vcf1} | sed -n '1p' | xz > {output.vcf1}"
        shell(s)
        s = "bedtools intersect -a <(xzcat {input.vcf1} | sed 1d | column -t -s , | awk -v OFS='\t' '{{print $1,$2-1,$2,$0}}') -b {input.mask} | cut -f 4- | sed -E 's/[[:space:]]+/,/g' | xz >> {output.vcf1}"
        shell(s)
rule concat_infile:
    input:
        expand("sim/D_{{divergence}}/infile/{{pop}}.ind_{{ind}}.chr_{chrom}.in{{ref_type}}.xz",  chrom = chr),
    output:
        infile1= "sim/D_{divergence}/infile/{pop}.ind_{ind}.chrAll.in{ref_type}.xz"
    run:
        s = "xzcat sim/D_{wildcards.divergence}/infile/{wildcards.pop}.ind_{wildcards.ind}.chr_1.in{wildcards.ref_type}.xz > sim/D_{wildcards.divergence}/infile/{wildcards.pop}.ind_{wildcards.ind}.chrAll.in{wildcards.ref_type}; "
        s += " for i in {{2..22}}; do echo \"$i\"; xzcat  sim/D_{wildcards.divergence}/infile/{wildcards.pop}.ind_{wildcards.ind}.chr_$i.in{wildcards.ref_type}.xz| sed 1d >>  sim/D_{wildcards.divergence}/infile/{wildcards.pop}.ind_{wildcards.ind}.chrAll.in{wildcards.ref_type}; done"
        shell(s)
        s1 = "xz sim/D_{wildcards.divergence}/infile/{wildcards.pop}.ind_{wildcards.ind}.chrAll.in{wildcards.ref_type}"
        shell(s1)
rule concat_infile_missing:
    input:
        expand("sim/D_{{divergence}}/infile_missing/{{pop}}.ind_{{ind}}.chr_{chrom}.in{{ref_type}}.xz",  chrom = chr),
    output:
        infile1= "sim/D_{divergence}/infile_missing/{pop}.ind_{ind}.chrAll.in{ref_type}.xz"
    run:
        s = "xzcat sim/D_{wildcards.divergence}/infile_missing/{wildcards.pop}.ind_{wildcards.ind}.chr_1.in{wildcards.ref_type}.xz > sim/D_{wildcards.divergence}/infile_missing/{wildcards.pop}.ind_{wildcards.ind}.chrAll.in{wildcards.ref_type}; "
        s += " for i in {{2..22}}; do echo \"$i\"; xzcat  sim/D_{wildcards.divergence}/infile_missing/{wildcards.pop}.ind_{wildcards.ind}.chr_$i.in{wildcards.ref_type}.xz| sed 1d >>  sim/D_{wildcards.divergence}/infile_missing/{wildcards.pop}.ind_{wildcards.ind}.chrAll.in{wildcards.ref_type}; done"
        shell(s)
        s1 = "xz sim/D_{wildcards.divergence}/infile_missing/{wildcards.pop}.ind_{wildcards.ind}.chrAll.in{wildcards.ref_type}"
        shell(s1)        
rule admixfrog_call:
    input:
        infile = "sim/D_{divergence}/infile/{pop}.ind_{ind}.chrAll.in{ref_type}.xz",
        ref = "sim/D_{divergence}/chrAll_ref{ref_type}.xz"
    output:
        bin = "sim/D_{divergence}/admixfrog_call/{ref_type}/{pop}.ind_{ind}.rle.xz"
    run:
        s = "admixfrog --infile {input.infile} --ref {input.ref} -o sim/D_{wildcards.divergence}/admixfrog_call/"
        s += wildcards.ref_type
        s += "/{wildcards.pop}.ind_{wildcards.ind} --ll-tol 1e-3 --bin-size 5000 --ancestral anc --max-iter 200 --filter-pos 50 --filter-map 0"
        s += " --states "
        s += ref_state[wildcards.ref_type] 
        s += " --run-penalty 0.2 --n-post-replicates 200 --dont-est-contamination --gt-mode"
        shell(s)


rule admixfrog_call_missing:
    input:
        infile = "sim/D_{divergence}/infile_missing/{pop}.ind_{ind}.chrAll.in{ref_type}.xz",
        ref = "sim/D_{divergence}/chrAll_ref{ref_type}.missing.xz",
    output:
        bin = "sim/D_{divergence}/admixfrog_call_missing/{ref_type}/{pop}.ind_{ind}.rle.xz"
    run:
        s = "admixfrog --infile {input.infile} --ref {input.ref} -o sim/D_{wildcards.divergence}/admixfrog_call_missing/"
        s += wildcards.ref_type
        s += "/{wildcards.pop}.ind_{wildcards.ind} --ll-tol 1e-3 --bin-size 5000 --ancestral anc --max-iter 200 --filter-pos 50 --filter-map 0"
        s += " --states "
        s += ref_state[wildcards.ref_type] 
        s += " --run-penalty 0.2 --n-post-replicates 200 --dont-est-contamination --gt-mode"
        shell(s)
        
rule admixfrog_call_all:
    input:
        nead3 = expand("sim/D_{divergence}/admixfrog_call/{asc}/{pop}.ind_{ind}.rle.xz", pop = ["ea"], ind=list(range(0,15)), divergence = divergence, asc = ['4','5','6'])
        
rule all_admixfrog_infile:
    input:
        in1_missing = expand("sim/D_{divergence}/infile_missing/{pop}.ind_{ind}.chrAll.in{ref_type}.xz", pop = ["ea",], ind=list(range(0,15)), ref_type = ['1','2','3','4','5'], divergence = divergence),
        in1 = expand("sim/D_{divergence}/infile/{pop}.ind_{ind}.chrAll.in{ref_type}.xz", pop = ["ea"], ind=list(range(0,15)), ref_type = ['1','2','3','4','5'], divergence = divergence),
        

rule all_admixfrog_call:
    input:
        call = expand("sim/D_{divergence}/admixfrog_call/{ref_type}/{pop}.ind_{ind}.rle.xz", pop = ["ea"], ind=list(range(0,15)), ref_type = ['1','2','3','4','5'], divergence = divergence),
        call_missing = expand("sim/D_{divergence}/admixfrog_call_missing/{ref_type}/{pop}.ind_{ind}.rle.xz", pop = ["ea"], ind=list(range(0,15)), ref_type = ['1','2','3','4','5'], divergence = divergence),
        
        


# rule admixfrog_bed:
#     input:
#         rle_nea_den = "sim/D_{divergence}/admixfrog_call/{pop}.ind_{ind}.rle.xz",
#         rle_archaic = "sim/D_{divergence}/admixfrog_call/{pop}.ind_{ind}.Archaic.rle.xz",
#         rle_den = "sim/D_{divergence}/admixfrog_call/Den.{pop}.ind_{ind}.rle.xz"
#     output:
#         bed_nea_den = "sim/D_{divergence}/admixfrog_call/{pop}.ind_{ind}.nea_den.bed",
#         bed_archaic = "sim/D_{divergence}/admixfrog_call/{pop}.ind_{ind}.Archaic.bed",
#         bed_den = "sim/D_{divergence}/admixfrog_call/Den.{pop}.ind_{ind}.den.bed"
#     run:
#         s = "xzcat {input.rle_nea_den} | grep -v 'afr' | grep 'state' | column -t -s , | awk -v OFS='\t' '{{print $1,$8,$11,$5}}' > {output.bed_nea_den} "
#         shell(s) 
#         ss = "xzcat {input.rle_archaic} | grep -v 'afr' | grep 'state' | column -t -s , | awk -v OFS='\t' '{{print $1,$8,$11,$5}}' > {output.bed_archaic} "
#         shell(ss)
#         sss = "xzcat {input.rle_den} | grep -v 'afr' | grep 'state' | column -t -s , | awk -v OFS='\t' '{{print $1,$8,$11,$5}}' > {output.bed_den} "
#         shell(sss)

# rule all_bed:
#     input:
#         bed_nea_den = expand("sim/D_{divergence}/admixfrog_call/{pop}.ind_{ind}.nea_den.bed", pop = ["ea"], ind=list(range(0,45))),
#         bed_archaic = expand("sim/D_{divergence}/admixfrog_call/{pop}.ind_{ind}.Archaic.bed", pop = ['ea'], ind=list(range(0,45))),
#         bed_den = expand("sim/D_{divergence}/admixfrog_call/Den.{pop}.ind_{ind}.den.bed", pop = ['ea'], ind=list(range(0,45)))

rule make_hmmix_input:
    input:
        data="sim/D_{divergence}/chrom_{chr}.snps.tsv.gz"
    output:
        hmmix_input = temp("sim/D_{divergence}/hmmix_infile/{pop}_ind{ind}.chr{chr}.in")
    run:
        make_infile_hmmix(pop=wildcards.pop, 
        ind = int(wildcards.ind),
        vcf = f"{input.data}", 
        out = f"{output.hmmix_input}")
rule make_hmmix_input_hap:
    input:
        data="sim/D_{divergence}/chrom_{chr}.snps.tsv.gz"
    output:
        hmmix_input = temp("sim/D_{divergence}/hmmix_infile/{pop}_ind{ind}.chr{chr}.in.hap")
    run:
        make_infile_hmmix12(pop=wildcards.pop, 
        ind = int(wildcards.ind),
        vcf = f"{input.data}", 
        out = f"{output.hmmix_input}")        
        
rule hmmix_input_concat:
    input:
        hmmix_input = expand("sim/D_{{divergence}}/hmmix_infile/{{pop}}_ind{{ind}}.chr{chr}.in", chr = chr)
    output:
        hmmix_input = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in"
    run:
        s = "for i in {{1..22}}; do cat sim/D_{wildcards.divergence}/hmmix_infile/{wildcards.pop}_ind{wildcards.ind}.chr$i.in >>{output.hmmix_input}; done"
        shell(s)

rule hmmix_input_concat_missing:
    input:
        hmmix_input = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in"
    output:
        hmmix_input = "sim/D_{divergence}/hmmix_infile_missing/{pop}.ind{ind}.in"
    run:
        s = "bedtools intersect -a <(awk -v OFS='\t' '{{print $1,$2-1,$2}}' {input.hmmix_input}) -b "
        s += missing_mask
        s += " | awk -v OFS='\t' '{{print $1,$3}}' > {output.hmmix_input}"
        shell(s)
    
rule hmmix_input_concat_hap:
    input:
        hmmix_input = expand("sim/D_{{divergence}}/hmmix_infile/{{pop}}_ind{{ind}}.chr{chr}.in.hap", chr = chr)
    output:
        hmmix_input = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in.hap"
    run:
        s = "for i in {{1..22}}; do cat sim/D_{wildcards.divergence}/hmmix_infile/{wildcards.pop}_ind{wildcards.ind}.chr$i.in.hap >>{output.hmmix_input}; done"
        shell(s)
        
rule all_hmmix_input:
    input:
        nomissing = expand("sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in", pop = ["ea"], ind=list(range(0,15)), chr = chr, divergence = divergence),
        missing = expand("sim/D_{divergence}/hmmix_infile_missing/{pop}.ind{ind}.in", pop = ["ea"], ind=list(range(0,15)), chr = chr, divergence = divergence)



rule hmmix_call:
    input:
        infile = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in"
    output:
        param = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.posterior.called"
    wildcard_constraints:
        ind = "|".join([str(i) for i in list(range(0,100))])
    run:
        s = "hmm/main.py gt_mode -data_type modern -count_file {input.infile} -mask_file mask_no_missing.chr1_20 "
        s += " -out sim/D_{wildcards.divergence}/hmmix_outfile/{wildcards.pop}.ind{wildcards.ind}"
        shell(s)
rule hmmix_call_missing:
    input:
        infile = "sim/D_{divergence}/hmmix_infile_missing/{pop}.ind{ind}.in"
    output:
        param = "sim/D_{divergence}/hmmix_outfile_missing/{pop}.ind{ind}.posterior.called"
    wildcard_constraints:
        ind = "|".join([str(i) for i in list(range(0,100))])
    run:
        s = "hmm/main.py gt_mode -data_type modern -count_file {input.infile} -mask_file "
        s+= missing_mask
        s += " -out sim/D_{wildcards.divergence}/hmmix_outfile_missing/{wildcards.pop}.ind{wildcards.ind}"
        shell(s)
rule hmmix_call_hap:
    input:
        infile = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in.hap"
    output:
        param = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.hap.posterior.called"
    wildcard_constraints:
        ind = "|".join([str(i) for i in list(range(0,100))])
    run:
        s = "hmm/main.py gt_mode -data_type modern -count_file {input.infile} -mask_file mask_no_missing.chr1_20.hap "
        s += " -out sim/D_{wildcards.divergence}/hmmix_outfile/{wildcards.pop}.ind{wildcards.ind}.hap"
        shell(s)

rule hmmix_anno:
    input:
        called = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.posterior.called",
        infile = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in",
        rec = "constant_rec.chr1_20"
    wildcard_constraints:
        ind = "|".join([str(i) for i in list(range(0,100))])
    output:
        anno_file = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.posterior.called.anno"
    run:
        with open(output.anno_file, 'w') as f:
            print('chrom','start','end','length', 'map_len', 'nea_out', 'nea_1', 'nea_2', 'den3', 'den25', 'source1', 'source2',
                  sep = '\t', file = f)
        for i in chr:
            anno(input.called, input.infile, f"sim/D_{wildcards.divergence}/chrom_{i}.snps.tsv.gz", output.anno_file, i, rec = input.rec)
            
rule hmmix_anno_missing:
    input:
        called = "sim/D_{divergence}/hmmix_outfile_missing/{pop}.ind{ind}.posterior.called",
        infile = "sim/D_{divergence}/hmmix_infile_missing/{pop}.ind{ind}.in",
        rec = "constant_rec.chr1_20"
    wildcard_constraints:
        ind = "|".join([str(i) for i in list(range(0,100))])
    output:
        anno_file = "sim/D_{divergence}/hmmix_outfile_missing/{pop}.ind{ind}.posterior.called.anno"
    run:
        with open(output.anno_file, 'w') as f:
            print('chrom','start','end','length', 'map_len', 'nea_out', 'nea_1', 'nea_2', 'den3', 'den25', 'source1', 'source2',
                  sep = '\t', file = f)
        for i in chr:
            anno(input.called, input.infile, f"sim/D_{wildcards.divergence}/chrom_{i}.snps.tsv.gz", output.anno_file, i, rec = input.rec)

rule hmmix_anno_hap:
    input:
        called = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.hap.posterior.called",
        infile = "sim/D_{divergence}/hmmix_infile/{pop}.ind{ind}.in.hap",
        rec = "constant_rec.chr1_20"
    wildcard_constraints:
        ind = "|".join([str(i) for i in list(range(0,100))])
    output:
        anno_file = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.hap.posterior.called.anno"
    run:
        with open(output.anno_file, 'w') as f:
            print('chrom','start','end','length', 'map_len', 'nea_out', 'nea_1', 'nea_2', 'den3', 'den25', 'source1', 'source2',
                  sep = '\t', file = f)
        for i in chr:
            anno(input.called, input.infile, f"sim/D_{wildcards.divergence}/chrom_{i}.snps.tsv.gz", output.anno_file, i, rec = input.rec, hap = True)
rule hmmix_anno_sep:
    input:
        anno_file = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.posterior.called.anno"
    output:
        hmm_den_match = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.posterior.called.anno.den",
        hmm_nea_match = "sim/D_{divergence}/hmmix_outfile/{pop}.ind{ind}.posterior.called.anno.nea"
    run:
        a = pd.read_csv(input.anno_file,sep = '\s+',na_values='', keep_default_na=False)
        a_den = a[(a['source1'] == "den25") | (a['source1'] == "den3")].copy()
        a_nea = a[(a['source1'] == "nea_out") | (a['source1'] == "nea_1") | (a['source1'] == "nea_2")].copy()
        a_den['MinMatch'] = a_den.apply(lambda x: min(int(x['den3'].split("/")[0]), int(x['den25'].split("/")[0])), axis=1)    # might need to get the maximum
        
        a_den.to_csv(output.hmm_den_match, sep = '\t', index = False)
        a_nea['MinMatch'] = a_nea.apply(lambda x: min(int(x['nea_out'].split("/")[0]), int(x['nea_1'].split("/")[0]), int(x['nea_2'].split("/")[0])), axis=1)# same here
        a_nea.to_csv(output.hmm_nea_match, sep = '\t', index = False)
rule hmmix_anno_sep_missing:
    input:
        anno_file = "sim/D_{divergence}/hmmix_outfile_missing/{pop}.ind{ind}.posterior.called.anno"
    output:
        hmm_den_match = "sim/D_{divergence}/hmmix_outfile_missing/{pop}.ind{ind}.posterior.called.anno.den",
        hmm_nea_match = "sim/D_{divergence}/hmmix_outfile_missing/{pop}.ind{ind}.posterior.called.anno.nea"
    run:
        a = pd.read_csv(input.anno_file,sep = '\s+',na_values='', keep_default_na=False)
        a_den = a[(a['source1'] == "den25") | (a['source1'] == "den3")].copy()
        a_nea = a[(a['source1'] == "nea_out") | (a['source1'] == "nea_1") | (a['source1'] == "nea_2")].copy()
        a_den['MinMatch'] = a_den.apply(lambda x: max(int(x['den3'].split("/")[0]), int(x['den25'].split("/")[0])), axis=1)
        
        a_den.to_csv(output.hmm_den_match, sep = '\t', index = False)
        a_nea['MinMatch'] = a_nea.apply(lambda x: max(int(x['nea_out'].split("/")[0]), int(x['nea_1'].split("/")[0]), int(x['nea_2'].split("/")[0])), axis=1)
        a_nea.to_csv(output.hmm_nea_match, sep = '\t', index = False)
        
#sum stats:
#1. sensitivity
#2. accuracy
#3. time to introgression
#4. length of introgression
#5. divergence of haplotypes

#1.total called
#2. how many total called hit the true introgressed
#3. how many total called hit the true introgressed, but hit the wrong source. misassingment
#4. how many total called hit the true introgressed, no matter of the source
#5. total length of simulation of the target source
#6. total length of simulation of the all sources  . test senesitivity and accuracy when do not care about the source
#7. total length of simulation of the target source, after filtered by length applied
#8. minimum map_len filter
#9. minimum match filter

rule sensi_:
    input:
        hmm_den_match = "sim/D_{divergence}/hmmix_outfile/ea.ind{ind}.posterior.called.anno.den",
        hmm_nea_match = "sim/D_{divergence}/hmmix_outfile/ea.ind{ind}.posterior.called.anno.nea",
        frog_call_AA = "sim/D_{divergence}/admixfrog_call/3/ea.ind_{ind}.rle.xz",
        frog_call_d3 = "sim/D_{divergence}/admixfrog_call/1/ea.ind_{ind}.rle.xz",
        frog_call_den = "sim/D_{divergence}/admixfrog_call/2/ea.ind_{ind}.rle.xz",
        frog_call_d3_AA = "sim/D_{divergence}/admixfrog_call/4/ea.ind_{ind}.rle.xz",
        frog_call_den_AA = "sim/D_{divergence}/admixfrog_call/5/ea.ind_{ind}.rle.xz",
        true_call = "sim/D_{divergence}/chrom_true_introgressed_all.bed.xz"
    output:
        sums = "sim/D_{divergence}/sum_stats/ea.ind{ind}.sums", # cut offs, 0.01cm, 0.02cm, 0.05cm. sensi, acc, mean_len
        #dis = "sim/D_{divergence}/sum_stats/{pop}.ind{ind}.all_len",
        out1 = "sim/D_{divergence}/sum_stats/ea.ind{ind}.frog_call",
        true_temp1 = "sim/D_{divergence}/sum_stats/ea.ind{ind}.true.nea.bed",
        true_temp2 = "sim/D_{divergence}/sum_stats/ea.ind{ind}.true.den.bed",
        true_temp3 = "sim/D_{divergence}/sum_stats/ea.ind{ind}.true.all.bed"
    run:
        s = "xzcat {input.frog_call_AA} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"AA\",i}}' >>{output.out1} "
        shell(s)
        s = "xzcat {input.frog_call_d3} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"d3\",i}}' >>{output.out1} "
        shell(s)
        s = "xzcat {input.frog_call_den} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"den\",i}}' >>{output.out1}"
        shell(s)
        s = "xzcat {input.frog_call_d3_AA} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"d3_AA\",i}}' >>{output.out1}"
        shell(s)
        s = "xzcat {input.frog_call_den_AA} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"den_AA\",i}}' >>{output.out1}"
        shell(s)
        td = "xzcat {input.true_call} | sed 1d | column -t -s , | sed 's/chr//g' | "
        td += "awk -v ind={wildcards.ind} -v OFS='\t' '((($1==(ind*2))||($1==(ind*2+1)))&&($5~\"nea\")){{print $2,$3,$4}}' |"
        td += "bedtools sort | bedtools merge > {output.true_temp1}"
        shell(td)
        tn = "xzcat {input.true_call} | sed 1d | column -t -s , | sed 's/chr//g' | "
        tn += "awk -v ind={wildcards.ind} -v OFS='\t' '((($1==(ind*2))||($1==(ind*2+1)))&&($5~\"den\")){{print $2,$3,$4}}' |"
        tn += "bedtools sort | bedtools merge > {output.true_temp2}"
        shell(tn)
        ta = "cat {output.true_temp1} {output.true_temp2} | bedtools sort | bedtools merge > {output.true_temp3}"
        shell(ta)
        for g_len in [0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1]:
            print(g_len)
            for min_m in [2,3,4,5,6]:
                """
                'called','hit','hit_wrong','hit_all_source','hit_none_source','simulation',
                                                'simulation_all_source','simulation_filter_by_len','min_map','min_match','ind','source','method'
                                                """
                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}" 
                s1 += " '($13>=m && $5>=(l-0.001)){{sum+=$3-$2}}END{{print sum}}' <(sed 1d {input.hmm_nea_match}) | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp3} -A| bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '{min_m}') <(echo '{wildcards.ind}')  <(echo 'nea') <(echo 'hmm')"
                s1 += ">> {output.sums}"
                shell(s1)
                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}" 
                s1 += " '($13>=m && $5>=(l-0.001)){{sum+=$3-$2}}END{{print sum}}' <(sed 1d {input.hmm_den_match}) | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp3} -A| bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)*800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '{min_m}') <(echo '{wildcards.ind}') <(echo 'den') <(echo 'hmm')"
                s1 += ">> {output.sums}"
                shell(s1)

                
                
            for ty in ['AA','d3','den','d3_AA','den_AA']:                    
                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v t="
                s1 += f"\"{ty}\" "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{sum+=$3-$2}}END{{print sum}}' {output.out1} | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} -A| bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '0') <(echo '{wildcards.ind}')  <(echo 'den') <(echo 'frog_{ty}')"
                s1 += ">> {output.sums}"
                shell(s1)

                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v t="
                s1 += f"\"{ty}\" "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{sum+=$3-$2}}END{{print sum}}' {output.out1} | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} -A| bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '0') <(echo '{wildcards.ind}')  <(echo 'nea') <(echo 'frog_{ty}') "
                s1 += ">> {output.sums}"
                shell(s1)

rule sensi_missing:
    input:
        hmm_den_match = "sim/D_{divergence}/hmmix_outfile_missing/ea.ind{ind}.posterior.called.anno.den",
        hmm_nea_match = "sim/D_{divergence}/hmmix_outfile_missing/ea.ind{ind}.posterior.called.anno.nea",
        frog_call_AA = "sim/D_{divergence}/admixfrog_call_missing/3/ea.ind_{ind}.rle.xz",
        frog_call_d3 = "sim/D_{divergence}/admixfrog_call_missing/1/ea.ind_{ind}.rle.xz",
        frog_call_den = "sim/D_{divergence}/admixfrog_call_missing/2/ea.ind_{ind}.rle.xz",
        frog_call_d3_AA = "sim/D_{divergence}/admixfrog_call_missing/4/ea.ind_{ind}.rle.xz",
        frog_call_den_AA = "sim/D_{divergence}/admixfrog_call_missing/5/ea.ind_{ind}.rle.xz",
        true_call = "sim/D_{divergence}/chrom_true_introgressed_all.bed.xz"
    output:
        sums = "sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.sums",
        out1 = "sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.frog_call",
        true_temp1 = "sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.true.nea.bed",
        true_temp2 = "sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.true.den.bed",
        true_temp3 = "sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.true.all.bed"
    run:
        s = "xzcat {input.frog_call_AA} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"AA\",i}}' >>{output.out1} "
        shell(s)
        s = "xzcat {input.frog_call_d3} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"d3\",i}}' >>{output.out1} "
        shell(s)
        s = "xzcat {input.frog_call_den} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"den\",i}}' >>{output.out1}"
        shell(s)
        s = "xzcat {input.frog_call_d3_AA} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"d3_AA\",i}}' >>{output.out1}"
        shell(s)
        s = "xzcat {input.frog_call_den_AA} | column -t -s , | awk -v i={wildcards.ind} -v OFS='\t' '($6==\"state\"&&$5!=\"afr\"){{print $1,$8,$11,$14,$5,\"den_AA\",i}}' >>{output.out1}"
        shell(s)
        
        td = "xzcat {input.true_call} | sed 1d | column -t -s , | sed 's/chr//g' | "
        td += "awk -v ind={wildcards.ind} -v OFS='\t' '((($1==(ind*2))||($1==(ind*2+1)))&&($5~\"nea\")){{print $2,$3,$4}}' |"
        td += "bedtools sort | bedtools merge > {output.true_temp1}"
        shell(td)
        tn = "xzcat {input.true_call} | sed 1d | column -t -s , | sed 's/chr//g' | "
        tn += "awk -v ind={wildcards.ind} -v OFS='\t' '((($1==(ind*2))||($1==(ind*2+1)))&&($5~\"den\")){{print $2,$3,$4}}' |"
        tn += "bedtools sort | bedtools merge > {output.true_temp2}"
        shell(tn)
        ta = "cat {output.true_temp1} {output.true_temp2} | bedtools sort | bedtools merge > {output.true_temp3}"
        shell(ta)
        for g_len in [0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1]:
            for min_m in [2,3,4,5,6]:
                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}" 
                s1 += " '($13>=m && $5>=(l-0.001)){{sum+=$3-$2}}END{{print sum}}' <(sed 1d {input.hmm_nea_match}) | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_nea_match})) -b {output.true_temp3} -A | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '{min_m}') <(echo '{wildcards.ind}')  <(echo 'nea') <(echo 'hmm')"
                s1 += ">> {output.sums}"
                shell(s1)
                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}" 
                s1 += " '($13>=m && $5>=(l-0.001)){{sum+=$3-$2}}END{{print sum}}' <(sed 1d {input.hmm_den_match}) | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v m="
                s1 += f"{min_m}"
                s1 += " '(($13>=m) && ($5>=(l-0.001))){{print $1,$2,$3}}' <(sed 1d {input.hmm_den_match})) -b {output.true_temp3} -A | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "                
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '{min_m}') <(echo '{wildcards.ind}') <(echo 'den') <(echo 'hmm')"
                s1 += ">> {output.sums}"
                shell(s1)

                
                
            for ty in ['AA','d3','den','d3_AA','den_AA']:                    
                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v t="
                s1 += f"\"{ty}\" "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{sum+=$3-$2}}END{{print sum}}' {output.out1} | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"d\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} -A| bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "         
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp2} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '0')  <(echo '{wildcards.ind}')  <(echo 'den') <(echo 'frog_{ty}')"
                s1 += ">> {output.sums}"
                shell(s1)

                s1 = "paste <(awk -v l="
                s1 += f"{g_len} -v t="
                s1 += f"\"{ty}\" "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{sum+=$3-$2}}END{{print sum}}' {output.out1} | sed 's/$/.0/')"
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp1} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp2} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(bedtools intersect -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} | bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(bedtools subtract -a <(awk -v OFS='\t' -v l="
                s1 += f"{g_len} -v t="
                s1 += f"{ty} "
                s1 += " '($4>=(l-0.001) && $5~\"nea\" && $6==t){{print $1,$2,$3}}' {output.out1}) -b {output.true_temp3} -A| bedtools sort | bedtools merge | awk '{{sum+=$3-$2}}END{{print sum}}' | sed 's/$/.0/') "              
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/') "
                s1 += " <(awk '{{sum+=$3-$2}}END{{print sum}}' {output.true_temp3} | sed 's/$/.0/') "
                s1 += f" <(awk -v l={g_len} "
                s1 += " '(($3-$2)>=(l-0.001)* 800000){{sum+=$3-$2}}END{{print sum}}' {output.true_temp1} | sed 's/$/.0/' ) "
                s1 += f"<(echo '{g_len}') <(echo '0') <(echo '{wildcards.ind}')  <(echo 'nea') <(echo 'frog_{ty}') "
                s1 += ">> {output.sums}"
                shell(s1)
rule merge_call:
    input:
        expand("sim/D_{{divergence}}/hmmix_outfile_missing/ea.ind{ind}.posterior.called.anno.{archaic}", archaic = ['nea','den'], ind = list(range(0,15))),
        expand("sim/D_{{divergence}}/hmmix_outfile/ea.ind{ind}.posterior.called.anno.{archaic}", archaic = ['nea','den'], ind = list(range(0,15))),
        frog_call1 = expand("sim/D_{{divergence}}/sum_stats/ea.ind{ind}.frog_call", ind = list(range(0,15))),
        frog_call2 = expand("sim/D_{{divergence}}/sum_stats_missing/ea.ind{ind}.frog_call", ind = list(range(0,15))),
        hmmix_call1 = expand("sim/D_{{divergence}}/sum_stats/ea.ind{ind}.sums", ind = list(range(0,15))),
        hmmix_call2 = expand("sim/D_{{divergence}}/sum_stats_missing/ea.ind{ind}.sums", ind = list(range(0,15)))
    output:
        frog_call1 = "sim/D_{divergence}/sum_stats/ea.frog_call.{phase}",
        frog_call2 = "sim/D_{divergence}/sum_stats_missing/ea.frog_call.{phase}",
        hmm_call1 = "sim/D_{divergence}/sum_stats/ea.posterior.called.anno.{phase}",
        hmm_call2 = "sim/D_{divergence}/sum_stats_missing/ea.posterior.called.anno.{phase}",
        true_call1 = "sim/D_{divergence}/sum_stats/ea.true.bed.{phase}",
        true_call2 = "sim/D_{divergence}/sum_stats_missing/ea.true.bed.{phase}",
        
    run:
        if wildcards.phase == "0":
            l = "for i in {{0..4}};"
        elif wildcards.phase == "1":
            l = "for i in {{5..9}};"
        elif wildcards.phase == "2":
            l = "for i in {{10..14}};"
        else:
            print("wrong wildcards phase")
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/sum_stats/ea.ind$i.frog_call >> {output.frog_call1}; done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/hmmix_outfile/ea.ind$i.posterior.called.anno.nea | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"nea\",ind}}' >>{output.hmm_call1};done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/hmmix_outfile/ea.ind$i.posterior.called.anno.den | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"den\",ind}}' >>{output.hmm_call1};done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/sum_stats/ea.ind$i.true.den.bed | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"den\",ind}}' >>{output.true_call1};done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/sum_stats/ea.ind$i.true.nea.bed | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"nea\",ind}}' >>{output.true_call1};done"
        shell(s)
### missing
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/sum_stats_missing/ea.ind$i.frog_call >> {output.frog_call2}; done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/hmmix_outfile_missing/ea.ind$i.posterior.called.anno.nea | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"nea\",ind}}' >>{output.hmm_call2};done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/hmmix_outfile_missing/ea.ind$i.posterior.called.anno.den | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"den\",ind}}' >>{output.hmm_call2};done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/sum_stats_missing/ea.ind$i.true.den.bed | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"den\",ind}}' >>{output.true_call2};done"
        shell(s)
        s = ""
        s += l
        s += "do cat sim/D_{wildcards.divergence}/sum_stats_missing/ea.ind$i.true.nea.bed | sed 1d | awk -v ind=$i -v OFS='\t' '{{print $0,\"nea\",ind}}' >>{output.true_call2};done"
        shell(s)
        
divergence = [90,150,210,270,330,390]
rule sensi_all:
    input:
        expand("sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.frog_call", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.sums", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.true.nea.bed", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats_missing/ea.ind{ind}.true.den.bed", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats_missing/ea.true.bed.{phase}", phase = [0,1,2], divergence = divergence),
        expand("sim/D_{divergence}/sum_stats_missing/ea.frog_call.{phase}", phase = [0,1,2], divergence = divergence),
        expand("sim/D_{divergence}/sum_stats_missing/ea.posterior.called.anno.{phase}", phase = [0,1,2], divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.ind{ind}.frog_call", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.ind{ind}.sums", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.ind{ind}.true.nea.bed", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.ind{ind}.true.den.bed", ind=list(range(0,15)), divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.true.bed.{phase}", phase = [0,1,2], divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.frog_call.{phase}", phase = [0,1,2], divergence = divergence),
        expand("sim/D_{divergence}/sum_stats/ea.posterior.called.anno.{phase}", phase = [0,1,2], divergence = divergence),
#!/usr/bin/env python3
import argparse
import numpy as np
from helper_f import load_observations_gt
from anno import anno, anno_strict
from hmm import write_HMM_to_file, read_HMM_parameters_from_file
import pandas as pd
from call import Call_LS, get_runs
from hmm_gt import TrainModel_GT
import time
def main():
    parser = argparse.ArgumentParser(add_help = True)

    subparser = parser.add_subparsers(dest = 'mode')

    train_gt_subparser = subparser.add_parser('gt_mode', help='Train HMM')
    train_gt_subparser.add_argument("-data_type", metavar='',help="[required] data type: modern genotypes or ancient genotyes", type=str, default = "modern")

    train_gt_subparser.add_argument("-not_est_transition",  action='store_true', help = "estimate transition parameter or not", default= False)
    train_gt_subparser.add_argument("-transition_1", metavar='',  help = "transition_param_1", default= None, type=float)
    train_gt_subparser.add_argument("-transition_2", metavar='',  help = "transition_param_2", default= None, type=float)
    train_gt_subparser.add_argument("-not_est_starting_prob",  action='store_true', help = "estimate starting probabilities or not", default= False)
    train_gt_subparser.add_argument("-starting_prob",  metavar='', help = "starting probabilities for archaic state", default= None, type = float)
    train_gt_subparser.add_argument("-count_file", metavar='', 
                                 help="[required] count of derived alleles", type=str, required = True)
    train_gt_subparser.add_argument("-max_variants_per_window", metavar='',
                                 help="maximum number of variants per window allowed", type=int, default = 50)
    train_gt_subparser.add_argument("-mask_file", metavar='',
                                 help="file with integrated mask", type=str, default=None)
    train_gt_subparser.add_argument("-mut_file", metavar='',
                                 help="file with mutation rates (default is mutation rate is uniform)", type=str, default=None)
    train_gt_subparser.add_argument("-param", metavar='',
                                 help="markov parameters file (default is human/neanderthal like parameters)", type=str, default= None)
    train_gt_subparser.add_argument("-out", metavar='',
                                 help="outputfile (default is a file named trained.json)", default = 'trained')
    train_gt_subparser.add_argument("-window_size", metavar='',
                                 help="size of bins (default is 1000 bp)", type=int, default = 1000)
    train_gt_subparser.add_argument("-iteration", metavar='', help = "max iteration for EM", type = int, default = 1000)
    train_gt_subparser.add_argument("-penalty", metavar='', help = "penalty score used for calling fragments", type = float, default = 0.2)
    train_gt_subparser.add_argument("-phased", action='store_true', help = "whether phased or unphased", default= False)
    call_subparser = subparser.add_parser('call', help='call fragments based on posterior')
    call_subparser.add_argument("-posterior", metavar='',
                                 help="output of posterior from fwd-bwd")
    #call_subparser.add_argument("-param", metavar='',
    #                             help="markov parameters file (default is human/neanderthal like parameters)", type=str, default= None)
    #call_subparser.add_argument("-called", metavar='',
    #                                help="output of called fragments", default = "called.txt")
    #call_subparser.add_argument("-window_size", metavar='', help="window_size", default= 1000, type = int)
    call_subparser.add_argument("-penalty", metavar='', help='penalty score used', default=0.2)
    call_subparser.add_argument("-call_type", metavar='', help='Penalty, or Average', default= None)
    call_subparser.add_argument("-call_out", metavar='', help='output file', default= None)
    anno_subparser = subparser.add_parser('anno', help='call fragments based on posterior')
    anno_subparser.add_argument("-called", metavar='',
                                 help="called file", default=None, type = str)
    anno_subparser.add_argument("-anno", metavar='',
                                 help="annotated file", type=str, default= None)
    anno_subparser.add_argument("-sample", metavar='',
                                    help="archaic individuals used for annotation", default = "")
    anno_subparser.add_argument("-group1", metavar='',
                                    help="archaic individuals by group used for annotation, e.g. Neanderthals", default = None)
    anno_subparser.add_argument("-group2", metavar='',
                                    help="archaic individuals by group used for annotation, e.g. Denisovans", default = None)
    anno_subparser.add_argument("-vcf", metavar='',
                                    help="vcf used for comparison", default = "")
    anno_subparser.add_argument("-map_file", metavar='', help = "map file used for annotation. one map with all chrs", default = None)
    anno_subparser.add_argument("-map", metavar='', help = "map used for annotation", default = "AA_Map")
    anno_subparser.add_argument("-window_size", metavar='', help = "window size", type = int, default = 1000)
    anno_subparser.add_argument("-count_file", metavar='', help = "observation", type = str, default = None)
    anno_subparser.add_argument("-phased", action='store_true', help = "whether phased or unphased", default= False)
    anno_subparser.add_argument("-strict", action='store_true', help = "whether anno and refine in a strict way", default= False)
    args = parser.parse_args()

    if args.mode == 'gt_mode':   # same as Larits. But starts with the first callable position.
        if not hasattr(args, 'count_file'):
            parser.print_help()
            return

        print("gt mode")
        hmm_parameters = read_HMM_parameters_from_file(args.param)
        print(args.param)
        #print(f"maximum number of variants per window allowed is {args.max_variants_per_window}")
        t1 = time.time()

        print(f"loading input data {args.count_file} with {args.data_type} mask {args.mask_file}, with window size {args.window_size}, phasing info: {args.phased}")
        print(f"maximum number of variants per window allowed is {args.max_variants_per_window}")
        chr_index, weights, obs, call_index, w = load_observations_gt(args.count_file, 
        args.mask_file, args.window_size, args.max_variants_per_window, args.data_type, args.mut_file,args.phased)
        nnn = len(chr_index)
        nnnn = len(weights)
        nnnnn = len(call_index)
        print(f"number of chr (with phasing info {args.phased}) : {nnn},{nnnn},{nnnnn}")
        #check for zero:
        t2 = time.time()
        print(f"loading time: {t2 - t1}")
        print(f"finished loading {args.count_file}")
        print(f"finished loading {args.mask_file}")
        print('-' * 40)
        if args.transition_1 is not None:
            hmm_parameters.transitions[0] = [1-args.transition_1, args.transition_1]
        if args.transition_2 is not None:
            hmm_parameters.transitions[1] = [args.transition_2, 1-args.transition_2]
        if args.not_est_starting_prob:
            print("not estimating starting probabilities")
            hmm_parameters.starting_probabilities = np.array([1-args.starting_prob, args.starting_prob])
        print(hmm_parameters)
        print('> Output is',args.out) 
        print('> Window size is',args.window_size, 'bp') 
        print('-' * 40)
        with open(args.out+".log", 'w') as log_file:
            log_file.write(f"gt mode\n")
            log_file.write(f"input type {args.data_type}\n")
            log_file.write(f"input_file {args.count_file}\n")
            log_file.write(f"{hmm_parameters}\n")
            log_file.write(f"maximum number of variants per window allowed is {args.max_variants_per_window}\n")
            log_file.write(f"output file prefix: {args.out}\n")
            log_file.write(f"window size: {args.window_size}\n")
            log_file.write(f"mut rate file: {args.mut_file}\n")
            log_file.write(f"posteiors file: {args.out}.posterior\n")
        hmm_parameters = TrainModel_GT(obs = obs,
                                weights = weights,
                                chr_index = chr_index,
                                w= w,
                                window_size = args.window_size,
                                call_index= call_index,
                                post_file= args.out+".posterior",
                                log_file = args.out+".log",
                                hmm_parameters = hmm_parameters,
                                maxiterations=args.iteration,
                                not_est_trans = args.not_est_transition,
                                not_est_starting_prob= args.not_est_starting_prob,
                                phased = args.phased)
        write_HMM_to_file(hmm_parameters, args.out+".param")
        get_runs(args.out+".posterior", args.penalty)     
        #Call(args.out+".posterior", hmm_parameters, args.out+".called", args.window_size)
        print("call finished")
    if args.mode == "call":
        print("calling fragments using posterior")
        if args.call_type == "Penalty":
            get_runs(args.posterior, args.penalty)
            print("call finished")
        elif args.call_type == "Average":
            print("Using averagel posteriors" + f"output is {args.call_out}")
            Call_LS(args.posterior, args.call_out)
        else:
            print("please specify the call type, either 'Penalty' or 'Average'")
    if args.mode == "anno":
        print(f"annotating fragments using observation file {args.count_file} and called file {args.called}")
        print(f"using archaic vcf {args.vcf}")
        print(f"using archaic individuals {args.sample}")
        print(f"sample names {args.sample}")
        print(f"out file {args.anno}.anno")
        print(f"annotating and refine called segs (at pos level) in strict way: {args.strict}")
        if args.strict is False:
            anno(gt_file = args.count_file, 
                called = args.called,
                vcf = args.vcf,
                out = args.anno,
                group1 = args.group1,
                group2 = args.group2,
                samples_input = args.sample,
                map_file = args.map_file,
                map = args.map,
                window_size = args.window_size,
                phased = args.phased)
        else:
            anno_strict(gt_file = args.count_file, 
                called = args.called,
                vcf = args.vcf,
                out = args.anno,
                group1 = args.group1,
                group2 = args.group2,
                samples_input = args.sample,
                map_file = args.map_file,
                map = args.map,
                window_size = args.window_size,
                phased = args.phased)
if __name__ == "__main__":
    main()

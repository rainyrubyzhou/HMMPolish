#################################################################################$$
# Driver Program for graph construction and viterbi search in graph. 
#################################################################################$$
import sys
import time
import os
import random
import argparse
import polish_main
import graph_prepare

#example running command: 
#python driver.py --read ../data/SRR13951181.fa --seed ../data/canu_ass.fa  --hmm ../data/7_profile.hmm

if __name__ == '__main__':
    #parse args
    parser = argparse.ArgumentParser()
    #parser.add_argument('-r', type=str, required=True, help = "Reads file for graph construction (in fasta format).")
    parser.add_argument('--read', type=str, required =True, help = "Reads file for graph construction (in fasta format).")
    
    #parser.add_argument('-s', type=int, default = 500,  required=False, help = "Seed file (in fasta format).")
    parser.add_argument('--seed', type=str, required = True,  help = "Seed file (in fasta format).")
    #parser.add_argument('-m', type=int, default = 500,  required=False, help = "Indicated Profile HMMs (in .hmm format).")
    parser.add_argument('--hmm', type=str, required = True, help = "Indicated Profile HMMs (in .hmm format).")
    parser.add_argument('-o', type=str, required = False, help = "Output filename.")
    parser.add_argument('--wei', type=float, default = 0.9,  required = False, help = "Weight of viterbi score in the recursive function. (Default = 0.9)")
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    args = parser.parse_args()
    print("Input files are:\n", " reads: ",  args.read, "\n  seed: ", args.seed, "\n  hmm:", args.hmm)
    seq_file = args.read
    seed_file = args.seed
    hmm_file = args.hmm    
    hmm_coef = args.wei
    output_file = args.o

    polish_main.prepare(seq_file, output_file, seed_file, hmm_file, hmm_coef, args.verbose)

from Bio import SeqIO
from numpy.lib import utils
import graph_prepare
import graph_main
import graph_viterbi
#import hmmdagcon_pacbio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
import time
import argparse
import re
import math
import aminoacid as aa
import aligngraph as ag
import numpy as np
import align
import utils
import os


"""
seq_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_145/sample_01.fa"
output_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_145/fp_out_supp2.fa"
seed_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_145/7622.fa"

seq_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_094/SRR13951184.fa"
output_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_094/fp_out_supp2.fa"
seed_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_094/7502.fa"

seq_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_097/SRR13951181.fa"
output_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_097/fp_out_ass.fa"
seed_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_097/canu_ass.fa"

seq_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_097/SRR13951181.fa"
output_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_097/fp_out_raw.fa"
seed_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_097/7553.fa"
"""
seq_file = "/home/runzhouyu2/work/frame-pro/real_data/HIV_HeLa/mapped.fa"
output_file = "/home/runzhouyu2/work/frame-pro/real_data/HIV_HeLa/fp_out.fa"
seed_file = "/home/runzhouyu2/work/frame-pro/real_data/HIV_HeLa/canu_ass.fa"

#hmm_file = "/home/runzhouyu2/work/frame-pro/real_data/noro_BMH19_145/hmm/all.hmm" 
hmm_file = "/home/runzhouyu2/work/frame-pro/HIV/9_cds/9_cds.hmm" 

dag_coef = 0.1
hmm_coef = 0.9 

ori_aln_graph = None
#def ori_graph():
path, filename = graph_prepare.find_dir(output_file)
consensus_file = path + "/" + filename.split(".fasta")[0] + "_consensus.fasta"
output_handle = open(output_file, "w")
#consensus_handle = open(consensus_file, "w")  
network_ext = 3
consensus_path, m5_dict, con_start_dict = graph_prepare.produce_concensus(seed_file,seq_file)
p_seq = graph_prepare.dna_transalate(consensus_path)
hmmer_domtbl = graph_prepare.hmmer_consensus(hmm_file, p_seq)
print("start read hmm", time.asctime( time.localtime(time.time()) ))
hmms = graph_viterbi.read_hmm_file(hmm_file)
print("finish read hmm", time.asctime( time.localtime(time.time()) ))

seq_dict = SeqIO.index(seq_file, "fasta")
consnesus_dict = SeqIO.index(consensus_path, "fasta")
target_dict = graph_main.read_hmmer_result(hmmer_domtbl, seq_dict, consnesus_dict, hmms)
for seq_index in target_dict:
    for hmm in target_dict[seq_index]:
        print("Successfully loaded HMMs:",hmm)
        target_hit = target_dict[seq_index][hmm]
        aligns = graph_main.split_alignment(target_hit, seq_dict[target_hit.id], m5_dict[target_hit.id], con_start_dict[seq_index][0], con_start_dict[seq_index][1],network_ext)
        aln_graph = graph_main.construct_network(aligns)
        ori_aln_graph = aln_graph
        # run with penalty for stop codon
        output = graph_viterbi.viterbi(hmms[hmm], ori_aln_graph, dag_coef, hmm_coef, psc = 1)
        print(dag_coef, hmm_coef, hmm)
        # print output["alignment"]
        seq_out = Seq(output["alignment"])
        record = SeqRecord(seq_out)
        record.id = aligns.backbone_seq.id + "_" + hmm + "_dag_" + str(dag_coef) + "_hmm_" + str(hmm_coef)
        record.description = ""
        SeqIO.write(record, output_handle, "fasta")
        
output_handle.close()


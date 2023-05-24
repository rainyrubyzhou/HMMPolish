from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from numpy.lib import utils
import graph_prepare
import graph_main
import graph_viterbi
import time
import subprocess

#example input
#seq_file = "SRR13951181.fa"
#output_file = "fp_out_raw.fa"
#seed_file = "7553.fa"
#hmm_file = "9_cds.hmm" 

def prepare(seq_file, output_file, seed_file, hmm_file, hmm_coef, verbose):
    hmm_coef = 1
    dag_coef = 1 - hmm_coef
    ori_aln_graph = None
    #def ori_graph():
 
    try: 
        #merge_cmd = "cat " + seed_file + " >>"  + seq_file
        #p1 = subprocess.run(merge_cmd, shell = True, check = True)
        merge_cmd2 = "cat " + seed_file + " " +  seq_file + "> merge.fa"
        p = subprocess.run(merge_cmd2, shell = True, check = True)

        #remove dup
        #rm_dup_cmd = "cat "+ seq_file + " | seqkit rmdup -n -o " + seq_file
        rm_dup_cmd = "cat merge.fa | seqkit rmdup -n -o " + seq_file
        #print(rm_dup_cmd)
        p = subprocess.run(rm_dup_cmd, shell = True, check = True)
        
        rm_merge = "rm merge.fa"
        p = subprocess.run(rm_merge, shell = True, check = True)

    except subprocess.CalledProcessError:
        print("Please input correct file paths for reads and seed files.")
        return

    path, filename = graph_prepare.find_dir(seq_file)
    consensus_file = path + "/" + filename.split(".fasta")[0] + "_consensus.fasta"
    if output_file == None:
        output_file = path + "/" + "polished.fa"

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
            print("\nSuccessfully loaded HMMs:",hmm)
            target_hit = target_dict[seq_index][hmm]
            aligns = graph_main.split_alignment(target_hit, seq_dict[target_hit.id], m5_dict[target_hit.id], con_start_dict[seq_index][0], con_start_dict[seq_index][1],network_ext)
            aln_graph = graph_main.construct_network(aligns)
            ori_aln_graph = aln_graph

            print("Score weight:", round(hmm_coef, 2), round(dag_coef, 2))
            output = graph_viterbi.viterbi(hmms[hmm], ori_aln_graph, dag_coef, hmm_coef, 1)
            # print output["alignment"]
            """ output nt seq"""
            seq_out = Seq(output["alignment"])
            record = SeqRecord(seq_out)
            record.id = aligns.backbone_seq.id + "_" + hmm + "_dag_" + str(round(dag_coef, 2)) + "_hmm_" + str(round(hmm_coef, 2))
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")
                
    output_handle.close()
    """verbose version: 
    1. output path weights
    2. 6-frame translation and hmmscan"""
    if verbose:
        out_p_seq = graph_prepare.dna_transalate(output_file)
        path_root = output_file.split(".fa")[0]
        pep_out = path_root + "_pep.fa"
        dir_path, old_name = graph_prepare.find_dir(out_p_seq)
        domtbl_path = dir_path+"/" + old_name.split(".fasta")[0] + "_domtbl.out"
        hmmer_cmd = "hmmscan --domtblout " + domtbl_path + " -E 1000 " + hmm_file + " " + pep_out + " >log"
        p2 = subprocess.run(hmmer_cmd, shell = True, check = True)

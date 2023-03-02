__author__ = 'Nan'

from Bio import SeqIO
import subprocess
import process
import aligngraph as ag
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def find_dir(path):
    """
    Given a file path, determine the directory
    :param path:
    :return:
    """
    file_part = path.split("/")
    dir_path = "/".join(file_part[0:len(file_part)-1])
    file_name = file_part[-1]
    return dir_path, file_name



def classify_seq(file_path, threshold = 6000):
    handle = open(file_path)
    dir_path, old_name = find_dir(file_path)
    seed_output = dir_path + "/" + old_name.split(".fasta")[0] + "_seed.fasta"
    short_output = dir_path + "/" + old_name.split(".fasta")[0] + "_short.fasta"
    seed_handle = open(seed_output, "w")
    short_handle = open(short_output, "w")
    records = SeqIO.parse(handle,"fasta")
    for record in records:

        if len(record)>threshold:
            SeqIO.write(record, seed_handle, "fasta")
        else:
            SeqIO.write(record, short_handle, "fasta")

    seed_handle.close()
    short_handle.close()

    return seed_output, short_output


def blasr_compare(seed_path, short_path):
    dir_path, old_name = find_dir(seed_path)
    out_path = dir_path + "/" + old_name.split("_seed.fasta")[0] + "_blasr.m5"
    
    blasr_cmd = "blasr " + seed_path +" " + short_path+" --bestn 200 -m 5 --out " + out_path
    # print blasr_cmd
    p = subprocess.run(blasr_cmd, shell = True, check = True)
    """
    p = subprocess.Popen(blasr_cmd , shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.communicate()
    print(out)
    """
    return out_path


def pair_m5_seq(seq, m5_align_list):
    """
    Given sequence and m5 result, form dag_align for construct network
    :param seq: back bone sequence for network
    :param m5_align_list:
    :return: dag_align:
    """
    dag_align = process.DagAlign(seq)
    for m5_align in m5_align_list:
        dag_align.add_alignment(m5_align)

    return dag_align


def dna_transalate(fasta_path):
    """

    :param fasta_path:
    :return: p_fasta:
    """
    path_root = fasta_path.split(".fasta")[0]
    p_fasta = path_root + "_pep.fasta"
    dna2pep_path = "/home/runzhouyu2/tools/dna2pep-1.1/dna2pep.py"
    
    dna2pep_cmd = "python " + dna2pep_path +" -r all --fasta " + p_fasta + " " + fasta_path
    p = subprocess.run(dna2pep_cmd, shell = True, check = True)
    """
    p = subprocess.Popen(dna2pep_cmd , shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.communicate()
    """
    return p_fasta
    
def hmmer_consensus(hmm_path, consensus_path):
    dir_path, old_name = find_dir(consensus_path)
    # print old_name
    domtbl_path = dir_path+"/" + old_name.split(".fasta")[0] + "_domtbl.out"
    print(domtbl_path)
    
    #hmmpress_cmd = "hmmpress " + hmm_path 
    #p1 = subprocess.run(hmmpress_cmd, shell = True, check = True)
    
    hmmer_cmd = "hmmscan --domtblout " + domtbl_path + " -E 1000 " + hmm_path + " " + consensus_path
    p2 = subprocess.run(hmmer_cmd, shell = True, check = True)
    """
    p = subprocess.Popen(hmmer_cmd , shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.communicate()
    """
    return domtbl_path

def produce_concensus(seed_file, fasta_file):
    dir_path, old_name = find_dir(fasta_file)
    consensus_path = dir_path + "/" + old_name.split(".fasta")[0] + "_consensus.fasta"
    consensus_handle = open(consensus_path, "w")

    #seed, short = classify_seq(fasta_file, 9300) # will not use the "short" file, all reads will be aligned
    #seed = dir_path + "/EB_canu_ass.fa"
    seed = seed_file
    
    #m5_align = blasr_compare(seed, short) #will include only shorter reads
    m5_align = blasr_compare(seed, fasta_file) #will include all reads 
    seq_dict = SeqIO.index(seed, "fasta")
    m5_handle = open(m5_align)
    m5_dict = process.read_m5(m5_handle)
    seq_con_pos_dict = {}
    for seq_key in seq_dict:
        print(seq_key)
        if seq_key in m5_dict:
            network_align = pair_m5_seq(seq_dict[seq_key], m5_dict[seq_key])
            aln_graph = process.construct_network(network_align)
            s,c, cov, c_start, c_end = aln_graph.generate_consensus()
            #print("nodes:", len(aln_graph.nodes), "\t edges:", len(aln_graph.edges) ,"\nconsensus:", s[0:10])
            seq_con_pos_dict[seq_key] = [c_start, c_end]
            seq = Seq(s)
            record = SeqRecord(seq)
            record.id = seq_key + "_consensus"
            record.description = ""
            SeqIO.write(record, consensus_handle, "fasta")

    return consensus_path, m5_dict, seq_con_pos_dict


if __name__ == "__main__":
    file = "reads.fasta"
    seed, short = classify_seq(file)
    m5_align = blasr_compare(seed, short)
    seq_dict = SeqIO.index(seed, "fasta")
    m5_handle = open(m5_align)
    m5_dict = process.read_m5(m5_handle)
    for seq_key in seq_dict:
        print(seq_key)
        print(m5_dict)
        network_align = pair_m5_seq(seq_dict[seq_key], m5_dict[seq_key])
        aln_graph = process.construct_network(network_align)
        s,c, cov = aln_graph.generate_consensus()
        print(s)


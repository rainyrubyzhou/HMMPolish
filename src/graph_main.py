__author__ = 'Nan'
from Bio import SeqIO
from Bio import SearchIO
import aligngraph


class DomHit(object):
    def __init__(self, HMM_model, bit_score, target_score, hmm_from, hmm_to, ali_from, ali_to, reading_frame, direction,
                 index):
        self.model = HMM_model
        self.bit_score = bit_score
        self.target_score = target_score
        self.hmm_from = hmm_from
        self.hmm_to = hmm_to
        self.ali_from = ali_from
        self.ali_to = ali_to
        self.frame = reading_frame
        self.direction = direction
        self.index = index  # index of same HMM domain
        # self.pep_len = hmm_query_len


class TargetHit(object):
    def __init__(self, seq_id, HMM_model, seq_length, con_len, hmm_len):
        self.id = seq_id
        self.model = HMM_model
        self.hit_start = None
        self.hit_end = None
        self._doms = []
        self.direction = None
        self.length = 0
        self.seq_len = seq_length
        self.con_len = con_len
        self.hmm_len = hmm_len

    def add_dom(self, dom):
        """
        add dom_hit object to the domian list
        :param dom: dom_hit object
        :return: None
        """
        if dom.model == self.model:
            self._doms.append(dom)
            if dom.frame == 1:
                new_start = dom.ali_from * 3
                new_end = dom.ali_to * 3
            elif dom.frame == 2:
                new_start = dom.ali_from * 3 + 1
                new_end = dom.ali_to * 3 + 1
            elif dom.frame == 3:
                new_start = dom.ali_from * 3 + 2
                new_end = dom.ali_to * 3 + 2

            elif dom.frame == -1:
                new_start = self.con_len - dom.ali_to * 3
                new_end = self.con_len - dom.ali_from * 3

            elif dom.frame == -2:
                new_start = self.con_len - dom.ali_to * 3 - 1
                new_end = self.con_len - dom.ali_from * 3 - 1

            elif dom.frame == -3:
                new_start = self.con_len - dom.ali_to * 3 - 2
                new_end = self.con_len - dom.ali_from * 3 - 2

            #python3 won't accept None VS int
            #if new_start < self.hit_start or self.hit_start == None:
            if self.hit_start == None or new_start < self.hit_start:
                self.hit_start = new_start

            #if new_end > self.hit_end or self.hit_end == None:
            if self.hit_end == None or new_end > self.hit_end:
                self.hit_end = new_end

            self.length = self.hit_end - self.hit_start

    def det_direction(self):
        """
        determine target direction by simply compare sum of each direction's bit score
        :return:
        """
        score_pos = 0
        score_rev = 0
        for dom in self._doms:

            if dom.direction == "+":
                score_pos += dom.bit_score
            if dom.direction == "-":
                score_rev += dom.bit_score

        if score_pos > score_rev:
            self.direction = "+"
        if score_rev > score_pos:
            self.direction = "-"


class DagAlign(object):
    def __init__(self, seq):
        self.backbone_seq = seq
        self._alignments = []
        self.direction = None

    def add_alignment(self, align):
        """

        :param align:
        :return:
        """
        self._alignments.append(align)

    def get_alignments(self):
        """

        :return:
        """
        return self._alignments

    def set_direction(self, direction):
        """

        :param direction:
        :return:
        """
        if direction == "+":
            self.direction = "+"
        if direction == "-":
            self.direction = "-"


class M5Align(object):
    def __init__(self, qName, qLength, new_q_start, new_q_end, qStrand, tName, tLength, new_t_start, new_t_end, tStrand,
                 new_qAlignedSeq,
                 new_matchPattern, new_tAlignedSeq):
        self.query_name = qName
        self.query_length = qLength
        self.query_start = new_q_start
        self.query_end = new_q_end
        self.query_strand = qStrand
        self.target_name = tName
        self.target_length = tLength
        self.target_start = new_t_start
        self.target_end = new_t_end
        self.target_strand = tStrand
        self.q_aligned_seq = new_qAlignedSeq
        self.match_pattern = new_matchPattern
        self.t_aligned_seq = new_tAlignedSeq


class ProgramError(Exception):
    def __init__(self, value):
        self.value = value

    @property
    def __str__(self):
        if self.value == 1:
            return "input error"


def read_hmmer_result(handle_domtbl, seq_dict, consensus_dict, hmm_dict):
    """

    :param handle_domtbl:
    :param seq_dict:
    :return:
    """
    domtbl = SearchIO.parse(handle_domtbl, "hmmscan3-domtab")
    target_dict = {}
    for hmm_query in domtbl:
        seq_id = hmm_query.id.split("_consensus_")[0]
        con_id =  hmm_query.id.split("_rframe")[0]
        
        if not (seq_id in target_dict):
            target_dict[seq_id] = {}

        reading_frame = hmm_query.id.split("_rframe")[1]
        if reading_frame[0] == "-":
            direction = "-"
        else:
            direction = "+"

        for hmm_hit in hmm_query:
            if not (hmm_hit.id in target_dict[seq_id]):
                target_dict[seq_id][hmm_hit.id] = TargetHit(seq_id, hmm_hit.id, len(seq_dict[seq_id]), len(consensus_dict[con_id]), hmm_dict[hmm_hit.id].length)

            for hmm_hsp in hmm_hit:
                dom_hit = DomHit(hmm_hit.id, hmm_hsp.bitscore, hmm_hit.bitscore, hmm_hsp.hit_start, hmm_hsp.hit_end,
                                 hmm_hsp.query_start, hmm_hsp.query_end, int(reading_frame), direction,
                                 hmm_hsp.domain_index)
                target_dict[seq_id][hmm_hit.id].add_dom(dom_hit)
                # notice: the start position is the position before first peptide, so 0 is the start position

    return target_dict


def read_m5(m5_handle):
    """

    :param m5_handle:
    :return: m5_alignment: dict of m5_align object
    """
    m5_alignment = {}
    lines = m5_handle.readlines()

    for line in lines:
        line = line.strip()
        line_split = line.split()
        if line_split[0] == "Warning:":
            continue
        else:
            qName = line_split[0].split("/")[0]
            qLength = int(line_split[1])
            qStart = int(line_split[2])
            qEnd = int(line_split[3])
            qStrand = line_split[4]
            tName = line_split[5]
            tLength = int(line_split[6])
            tStart = int(line_split[7])
            tEnd = int(line_split[8])
            tStrand = line_split[9]
            # score = line_split[10]
            # numMatch = line_split[11]
            # numMismatch = line_split[12]
            # numIns = line_split[13]
            # numDel = line_split[14]
            # mapQV = line_split[15]
            qAlignedSeq = line_split[16]
            matchPattern = line_split[17]
            tAlignedSeq = line_split[18]
            m5_align = M5Align(qName, qLength, qStart, qEnd, qStrand, tName, tLength, tStart, tEnd, tStrand,
                               qAlignedSeq, matchPattern, tAlignedSeq)
            if not (m5_align.query_name in m5_alignment):
                m5_alignment[m5_align.query_name] = []

            m5_alignment[m5_align.query_name].append(m5_align)
    return m5_alignment


def split_alignment(target_hit, seq, m5_align_list, con_start = None, con_end = None,cov_scale=3):
    # m5_parse_output = open("/mnt/home/dunan/Job/HMMDAGCON_large_test_Feb_2016/20160319_new_state_worst_case/S1_5503_program_new.m5", "w")
    
    con_length = target_hit.con_len
    
    # to make sure the right region of the seq
  
    if con_start == None:
        con_start = 0
    if con_end == None:
        con_end = len(seq)    
    
    
    seq_length =  con_end -  con_start
    
    length = int(float(target_hit.length)/con_length * seq_length)

    if cov_scale < 1:
        raise ProgramError(1)
    else:
        leng_to_tend = int(length*(float(cov_scale) - 1) / 2)
        if leng_to_tend > 2 * target_hit.hmm_len:
            leng_to_tend = 2 * target_hit.hmm_len
        network_start = con_start + int(float(target_hit.hit_start)/con_length * seq_length) - leng_to_tend
        
        network_end = con_start + int(float(target_hit.hit_end)/con_length * seq_length) + leng_to_tend

        if network_start < 0:
            network_start = 0
        if network_end > len(seq):
            network_end = len(seq)
        print(network_start, network_end)
        target_hit.det_direction()
        seq_cut = seq[network_start:network_end]
        seq_cut.id = seq.id + "_" + str(network_start) + "_" +str(network_end)
        dag_align = DagAlign(seq_cut)
        dag_align.set_direction(target_hit.direction)

        for m5_align in m5_align_list:
            qName = m5_align.query_name
            qLength = m5_align.query_length
            qStart = m5_align.query_start
            qEnd = m5_align.query_end
            qStrand = m5_align.query_strand
            tName = m5_align.target_name
            tLength = m5_align.target_length
            tStart = m5_align.target_start
            tEnd = m5_align.target_end
            tStrand = m5_align.target_strand
            
            score = "0"
            numMatch = "0"
            numMismatch = "0"
            numIns = "0"
            numDel = "0"
            mapQV = "0"
            
            qAlignedSeq = m5_align.q_aligned_seq
            matchPattern = m5_align.match_pattern
            tAlignedSeq = m5_align.t_aligned_seq

            # rename query sequence name
            qName = qName + "_" + str(network_start) + "_" + str(network_end)
            
            #seq_start = network_start + 1
            seq_start = network_start
            seq_end = network_end - 1

            if (seq_start >= qEnd) or (seq_end <= qStart):
                # q is not in the range
                continue

            else:
                # if q begin later than seq we want, start from 0, other with start from diffrence
                if seq_start < qStart:
                    q_start_pos = 0
                else:
                    q_start_pos = seq_start - qStart

                # if q end earlier than seq we want, end from 0, other with end from diffrence
                if seq_end > qEnd:
                    q_end_pos_rev = 0
                else:
                    q_end_pos_rev = qEnd - seq_end

                # find the real start position in given sequence with "-"
                index_query = 0
                start_modify = 0
                for i, charac_i in enumerate(qAlignedSeq):
                    if index_query == q_start_pos:
                        real_start = i
                        break
                    elif charac_i == "-":
                        continue
                    else:
                        index_query += 1

                # find the real end position in given sequence with "-"
                index_query = 0
                for j, charac_j in reversed(list(enumerate(qAlignedSeq))):
                    if index_query == q_end_pos_rev:
                        real_end = j + 1
                        break
                    elif charac_j == "-":
                        continue
                    else:
                        index_query += 1

                # find new qStart and qEnd

                new_q_start = qStart + q_start_pos - seq_start

                new_q_end = qEnd - seq_start - q_end_pos_rev

                # find new tStart and tEnd
                if tStrand == "+":
                    t_start_pos = 0
                    for i2, charac_i2 in enumerate(tAlignedSeq):
                        if i2 == real_start:
                            break
                        elif charac_i2 == "-":
                            continue
                        else:
                            t_start_pos += 1

                    t_end_pos_rev = 0
                    for j2, charac_j2 in reversed(list(enumerate(tAlignedSeq))):
                        if j2 == real_end - 1:
                            break
                        elif charac_j2 == "-":
                            continue
                        else:
                            t_end_pos_rev += 1


                elif tStrand == "-":
                    t_start_pos = 0
                    for i2, charac_i2 in enumerate(tAlignedSeq):

                        if i2 == real_start:

                            break
                        elif charac_i2 == "-":
                            continue
                        else:
                            t_start_pos += 1

                    t_end_pos_rev = 0
                    for j2, charac_j2 in reversed(list(enumerate(tAlignedSeq))):

                        if j2 == real_end - 1:
                            break
                        elif charac_j2 != "-":
                            t_end_pos_rev += 1

                new_t_start = tStart + t_start_pos
                new_t_end = tEnd - t_end_pos_rev

                new_q_length = seq_end - seq_start
                new_t_length = new_t_end - new_t_start

                alignment = M5Align(qName, new_q_length, new_q_start, new_q_end, qStrand, tName, new_t_length,
                                    new_t_start, new_t_end, tStrand, qAlignedSeq[real_start:real_end],
                                    matchPattern[real_start:real_end],
                                    tAlignedSeq[real_start:real_end])

                dag_align.add_alignment(alignment)

                """
                output = [qName, str(seq_end - seq_start), str(new_q_start),
                              str(new_q_end), qStrand, tName,
                              str(new_t_end - new_t_start), str(new_t_start),
                              str(new_t_end), tStrand, score, numMatch, numMismatch,
                              numIns, numDel, mapQV, qAlignedSeq[real_start:real_end],
                              matchPattern[real_start:real_end],
                              tAlignedSeq[real_start:real_end]]
                
                """

    return dag_align


def construct_network(dag_align):
    """

    :param dag_align:
    :return:
    """
    alignments = dag_align.get_alignments()
    aln_graph = aligngraph.AlnGraph(dag_align.backbone_seq.seq)
    back_len = len(dag_align.backbone_seq)

    for align in alignments:
        tAln, qAln = aligngraph.convert_mismatches(align.t_aligned_seq, align.q_aligned_seq)
        
        qAln_strip = qAln.rstrip("-")
        
        if len(qAln_strip) < len(qAln):
            tail = len(qAln) - len(qAln_strip)
            tAln = tAln[:-tail]
            qAln = qAln_strip
        
        aln = ((0, align.target_end - align.target_start, tAln), (align.query_start, align.query_end, qAln))
        aln_graph.add_alignment(aln)

    aln_graph.merge_nodes()

    if dag_align.direction == "-":
        aligngraph.reverse_graph(aln_graph)

    return aln_graph
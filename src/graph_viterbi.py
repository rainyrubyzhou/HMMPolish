import re
import math
import aminoacid as aa
import aligngraph as ag
import numpy as np

path_compensate = 0
backbone_penalty = 10.00
background_path = -1.38629
AA_index = "ACDEFGHIKLMNPQRSTVWY"


class HmmModel:
    def __init__(self):
        self.length = 0
        self.name = None
        self.accession = ""
        self.match_emission = []
        self.insert_emission = []
        self.transition = []
        self.reference = []

    def set_length(self, length):
        """
        :param length: hmm model length
        """
        self.length = length
        self.match_emission = [[None for i in range(21)] for j in range(length + 1)]
        self.insert_emission = [[None for i in range(21)] for j in range(length + 1)]
        self.transition = [[None for i in range(7)] for j in range(length + 1)]
        self.reference = [None for i in range(length + 1)]

    def read_hmm(self, file_handle):
        """
        :param file_handle: handle of hmmer3 hmm model file
        :return is_end, to check if there is remain hmm model in the file handle
        """
        lines = []
        hmm_pos = 1
        is_end = False
        
        while True:
            line = file_handle.readline()
            if line == "":
                is_end = True
                break

            if line.strip() == '//':
                is_end = False
                break
            else:
                lines.append(line)

        # hmm_model = HmmModel()
        if is_end is False: 
            for i, line in enumerate(lines):
                line = line.strip()
                line_split = re.split(" +", line)
                if line_split[0] == "NAME":
                    self.name = line_split[1]
    
                elif line_split[0] == "ACC":
                    self.accession = line_split[1]
    
                elif line_split[0] == "LENG":
                    self.set_length(int(line_split[1]))
                elif len(line_split) > 1 and line_split[0] == "COMPO":
                    """
                    when found compo, record the initial insert_emission and transition, also mark the line position
                    """
    
                    hmm_pos = 0
                    # pos = 0, the match_emission is background emission
                    print(" hmm model length: ", len(self.match_emission), len(self.match_emission[0]))
                    for j in range(20):
                        line_split1 = re.split(" +", lines[i].strip("\n").strip("\r"))
                        if line_split1[j + 2] == "*":
                            self.match_emission[hmm_pos][j] = float("-inf")
                        else:
                            self.match_emission[hmm_pos][j] = -float(line_split1[j + 2])
                    for j in range(20):
                        line_split2 = re.split(" +", lines[i + 1].strip("\n").strip("\r"))
                        if line_split2[j + 1] == "*":
                            self.insert_emission[hmm_pos][j] = float("-inf")
                        else:
                            self.insert_emission[hmm_pos][j] = -float(line_split2[j + 1])
                    for j in range(7):
                        line_split3 = re.split(" +", lines[i + 2].strip("\n").strip("\r"))
                        if line_split3[j + 1] == "*":
                            self.transition[hmm_pos][j] = float("-inf")
                        else:
                            self.transition[hmm_pos][j] = -float(line_split3[j + 1])
                    hmm_pos = 1
                elif len(line_split) > 1 and line_split[0] == str(hmm_pos):
    
                    for j in range(20):
                        line_split1 = re.split(" +", lines[i].strip("\n").strip("\r"))
                        if line_split1[j + 2] == "*":
                            self.match_emission[hmm_pos][j] = float("-inf")
                        else:
                            self.match_emission[hmm_pos][j] = -float(line_split1[j + 2])
                    for j in range(20):
                        line_split2 = re.split(" +", lines[i + 1].strip("\n").strip("\r"))
                        if line_split2[j + 1] == "*":
                            self.insert_emission[hmm_pos][j] = float("-inf")
                        else:
                            self.insert_emission[hmm_pos][j] = -float(line_split2[j + 1])
                    for j in range(7):
                        line_split3 = re.split(" +", lines[i + 2].strip("\n").strip("\r"))
                        if line_split3[j + 1] == "*":
                            self.transition[hmm_pos][j] = float("-inf")
                        else:
                            self.transition[hmm_pos][j] = -float(line_split3[j + 1])
                    hmm_pos += 1
            # give value for the "*"
            self.match_emission[0][20] = 0.00000
           
            for i in range(self.length + 1):
                
                if i != 0:
                    #No.21 elements are None now in the list , cannot get the max value
                    self.reference[i] = AA_index[self.match_emission[i].index(max(filter(None.__ne__,self.match_emission[i])))]
                    self.match_emission[i][20] = -3.06027 
                self.insert_emission[i][20] = -3.06027
        return is_end


def read_hmm_file(file_path):
    """
    :param file_path: hmm model file path
    :return: hmm_dict: a dict of all hmm models in the file
    """
    handle = open(file_path)
    is_end = False
    hmm_dict = {}
    while True:
        hmm = HmmModel()
        is_end = hmm.read_hmm(handle)
        if hmm.name == None:
            break
        else:
            hmm_dict[hmm.name] = hmm
    handle.close()
    return hmm_dict


def find_ancestor(aln_node):
    """
    :param graph: the graph
    :param node: the id of the node
    :return node_list: 2d list include all node and its parent and grandparent combination, also may include the
    ancestor node if exist.
    """
    node_list = []
    for in_edge in aln_node._in_edges:
        parent = in_edge.in_node

        for in_edge_parent in parent._in_edges:
            grandparent = in_edge_parent.in_node

            if grandparent._in_edges:
                for in_edge_grandparent in grandparent._in_edges:
                    ancestor = in_edge_grandparent.in_node
                    node = [parent, grandparent, ancestor]
                    node_list.append(node)
            else:
                node = [parent, grandparent]
                node_list.append(node)

    return node_list


def viterbi(hmm_model, aln_graph, c_dag=0, c_hmm = 1, psc = 0.8):
    """
    :param hmm_model: hmm model
    :param sorted_nodes: a list of sorted nodes in graph
    """
    sorted_nodes = aln_graph.get_sorted_nodes()
    aln_graph.filter_graph()

    m_to_m = 0
    i_to_m = 3
    m_to_i = 1
    m_to_d = 2
    i_to_i = 4
    d_to_m = 5
    d_to_d = 6

    hmm_len = hmm_model.length
    dna_len = len(sorted_nodes)
    print(dna_len)

    scores = [[[float("-inf") for x in range(4)] for x in range(hmm_len + 1)] for x in
              range(dna_len + 1)]  # score matrix for DP, score[i][0]:m score[i][1]:i score[i][2]:d score[i][3]:x
    state_x = [[float("-inf") for x in range(2)] for x in range(dna_len + 1)]
    flags = [[[{"path": "", "state": ""} for x in range(4)] for x in range(hmm_len + 1)] for x in
             range(dna_len + 1)]  # flag matrix for DP
    flag_x = [[{"path": "", "state": ""} for x in range(2)]for x in range(dna_len + 1)]
    begin_to_match = [0 for x in range(hmm_len + 1)]
    """
    end_flag[i][0] stores the value of j which gives largest END score;
    end_flag[i][1] indicates whether E state leads to largest score.
    """
    end_flag = [[0 for x in range(2)] for x in range(dna_len + 1)]

    # background = hmm_model.match_emission[0]
    

    state_x[0][0] = 0.0
    state_x[0][1] = 0.0

    max_end_score = float("-Inf")
    max_end_score_pos = 0
    codon = ["" for x in range(3)]
    # aa_emission = ""
    # aa_index = 0
    check_node_order = {}
    
    # test how many possible path for one node
    len_pre = []

    for i, node_i in enumerate(sorted_nodes):
        check_node_order[node_i.ID] = i
        if i > 0:
            # state_x_max = [float("-inf"), float("-inf")]
            for in_edge in node_i._in_edges:
                in_node = in_edge.in_node
                in_node_index = check_node_order[in_node.ID]
                if in_node.backbone_node.coverage < 5:
                    backbone_node_coverage = in_node.backbone_node.coverage + 1
                else:
                    backbone_node_coverage = in_node.backbone_node.coverage
                if node_i.is_backbone is True and node_i.weight == 1:
                    path_score = - backbone_penalty  # ln(0.5) - ln(0.45)
                else:
                    path_score = in_edge.count - backbone_node_coverage * 0.5

                new_x_score = state_x[in_node_index][0] + c_dag * path_score
                if new_x_score>state_x[i][0]:
                    state_x[i][0] = new_x_score
                    flag_x[i][0] = {"path": [in_node], "state": 4}
                new_x_score = state_x[in_node_index][1] + c_dag * path_score

                if new_x_score>state_x[i][1]:
                    state_x[i][1] = new_x_score
                    flag_x[i][1] = {"path": [in_node], "state": 4}

                for j in range(1, hmm_len + 1):
                    new_x_score = c_dag * path_score + scores[in_node_index][j][0]
                    if new_x_score > state_x[i][1]:
                        state_x[i][1] = new_x_score
                        flag_x[i][1] = {"path": [in_node], "state": 0}
                        end_flag[in_node_index][0] = j

        if (i >= 3) and (i < dna_len - 1):
            predecessors = find_ancestor(node_i)
            if len(predecessors) == 0:
                len_pre.append(1)
            else:
                len_pre.append(len(predecessors))
            for j in range(1, hmm_len + 1):
                score_max = [float("-inf"), float("-inf"), float("-inf"), float("-inf")]
                score_max_em = [float("-inf"), float("-inf"), float("-inf"), float("-inf")]




                for predecessor in predecessors:

                    if len(predecessor) > 2:
                        codon[0] = predecessor[1].base
                        codon[1] = predecessor[0].base
                        codon[2] = node_i.base

                        aa_emission = aa.translation(codon[0], codon[1],
                                                     codon[2])

                        aa_index = aa.aminoacidindex(aa_emission)

                        # print "restart", score_temp
                        # only when the find ancestor return 3 node, we need to consider transition

                        last_score_node = predecessor[2]
                        path_node = [predecessor[2], predecessor[1], predecessor[0], node_i]

                        path_score = 0
                        extreme_path = False

                        for i1, node in enumerate(path_node):
                            for in_edge in node._in_edges:
                                in_node = in_edge.in_node

                                if i1 != 0 and in_node == path_node[i1 - 1]:

                                    if in_node.backbone_node.coverage < 5:
                                        backbone_node_coverage = in_node.backbone_node.coverage + 1
                                    else:
                                        backbone_node_coverage = in_node.backbone_node.coverage

                                        
                                    if node.is_backbone is True and node.weight == 1:
                                        path_score -= backbone_penalty  # ln(0.5) - ln(0.45)
                                    else:
                                        path_score += in_edge.count - backbone_node_coverage * 0.5
                                        
                                        # path_score += round(math.log(in_edge.count / (backbone_node_coverage + 1.0)),5) - background_path

                        penalty_sc =  1
                        score_em = -1
                        score_temp = -1
                        
                        if aa.aminoacidindex == 20:
                            penalty_sc = psc #penalty for stop codon
                            #score_em = -1
                            #score_temp = -1
                        
                        score_em = penalty_sc * (c_hmm * (hmm_model.match_emission[j][aa_index] - hmm_model.match_emission[0][
                            aa_index]) + c_dag * path_score + state_x[check_node_order[last_score_node.ID]][0])
                        score_temp = penalty_sc * (c_dag * path_score + state_x[check_node_order[last_score_node.ID]][0])

                        if score_temp > score_max[0]:
                            score_max[0] = score_temp
                            score_max_em[0] = score_em
                            flags[i][j][0] = {"path": [predecessor[2],predecessor[1], predecessor[0]], "state": 4, "aa": aa_emission}

                            # from match to match

                        score_em = penalty_sc *(c_hmm * (
                                hmm_model.transition[j - 1][m_to_m] +
                                hmm_model.match_emission[j][aa_index] - hmm_model.match_emission[0][aa_index]) \
                                       + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j - 1][0])
                        score_temp = penalty_sc *(c_hmm * (hmm_model.transition[j - 1][m_to_m]) \
                                         + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j - 1][0])

                        if score_temp > score_max[0]:
                            score_max[0] = score_temp
                            score_max_em[0] = score_em
                            flags[i][j][0] = {"path": [predecessor[2], predecessor[1], predecessor[0]], "state": 0,
                                                  "aa": aa_emission}

                            # from insertion to match
                        score_em = penalty_sc *(c_hmm * (hmm_model.transition[j - 1][i_to_m] +
                                                hmm_model.match_emission[j][aa_index] - hmm_model.match_emission[0][
                                                    aa_index]) \
                                       + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j - 1][1])
                        score_temp = penalty_sc *(c_hmm * (hmm_model.transition[j - 1][i_to_m]) \
                                         + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j - 1][1])
                        if score_temp > score_max[0]:
                            score_max[0] = score_temp
                            score_max_em[0] = score_em
                            flags[i][j][0] = {"path": [predecessor[2], predecessor[1], predecessor[0]], "state": 1,
                                                  "aa": aa_emission}

                        score_em = penalty_sc *(c_hmm * (hmm_model.transition[j - 1][d_to_m] +
                                                hmm_model.match_emission[j][aa_index] - hmm_model.match_emission[0][
                                                    aa_index]) \
                                       + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j - 1][2])
                        score_temp = penalty_sc *(c_hmm * (hmm_model.transition[j - 1][d_to_m]) \
                                         + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j - 1][2])
                        if score_temp > score_max[0]:
                            score_max[0] = score_temp
                            score_max_em[0] = score_em
                            flags[i][j][0] = {"path": [predecessor[2], predecessor[1], predecessor[0]], "state": 2,
                                                  "aa": aa_emission}

                            """
                            insertion state
                            """
                            # from match
                        score_em = penalty_sc *(c_hmm * (
                                hmm_model.transition[j][m_to_i] + hmm_model.insert_emission[j][aa_index] -
                                hmm_model.match_emission[0][aa_index]) + c_dag * path_score + \
                                       scores[check_node_order[last_score_node.ID]][j][0])
                        score_temp = penalty_sc *(c_hmm * (hmm_model.transition[j][m_to_i]) + c_dag * path_score + \
                                         scores[check_node_order[last_score_node.ID]][j][0])

                        if score_temp > score_max[1]:
                            score_max[1] = score_temp
                            # score_max_em[1] = score_em
                            score_max_em[1] = score_temp
                            flags[i][j][1] = {"path": [predecessor[2], predecessor[1], predecessor[0]], "state": 0,
                                                  "aa": aa_emission}

                            # from insertion
                        score_em = penalty_sc *(c_hmm * (hmm_model.transition[j][
                                                    i_to_i] + hmm_model.insert_emission[j][aa_index] -
                                                hmm_model.match_emission[0][
                                                    aa_index]) \
                                       + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j][1])

                        score_temp = penalty_sc *(c_hmm * (hmm_model.transition[j][i_to_i]) \
                                         + c_dag * path_score + scores[check_node_order[last_score_node.ID]][j][1])

                        if score_temp > score_max[1]:
                            score_max[1] = score_temp
                                # score_max_em[1] = score_em
                            score_max_em[1] = score_temp
                            flags[i][j][1] = {"path": [predecessor[2], predecessor[1], predecessor[0]], "state": 1,
                                                  "aa": aa_emission}

                            """
                            deletion state
                            """
                            # from match
                        score_temp = scores[i][j - 1][0] + c_hmm * hmm_model.transition[j - 1][m_to_d]
                        if score_temp > score_max[2]:
                            score_max[2] = score_temp
                            score_max_em[2] = score_temp
                            flags[i][j][2] = {"path": [], "state": 0, "aa": "-"}
                            # from deletion
                        score_temp = scores[i][j - 1][2] + c_hmm * hmm_model.transition[j - 1][d_to_d]
                        if score_temp > score_max[2]:
                            score_max[2] = score_temp
                            score_max_em[2] = score_temp
                            flags[i][j][2] = {"path": [], "state": 2, "aa": "-"}

                scores[i][j] = score_max_em

    score = state_x[dna_len - 1][1]
    print(score)
        # / math.log10(2.0)
    # print flags[3791][47][1]
    
    
    print("ave possible path", np.mean(len_pre))
    
    """
    traceback
    """
    # traceback = {}
    i_end = dna_len - 1
    i_track = i_end
    score_track = []
    j_end = end_flag[i_end][0]
    j_track = 0
    alignment = ""
    aa_alignment = ""
    state = 4
    state_alignment = ""
    ref_alignment = ""
    i_start = 0
    j_start = 0

    pass_hmm = False

    while (i_track >= 0):
        state_alignment = str(state) + state_alignment
        if state ==  4:
            if pass_hmm is False:
                path = flag_x[i_track][1]["path"]
                last_state = flag_x[i_track][1]["state"]
            else:
                path = flag_x[i_track][0]["path"]
                last_state = flag_x[i_track][0]["state"]
            if sorted_nodes[i_track].base == "B":
                break
            elif sorted_nodes[i_track].base != "E":
               alignment = sorted_nodes[i_track].base + alignment
            # new i_track
            i_track = check_node_order[path[0].ID]
            if last_state == 0:
                j_track = end_flag[i_track][0]

        else:
            aa_emit = flags[i_track][j_track][state]["aa"]
            path = flags[i_track][j_track][state]["path"]
            score_track.append(scores[i_track][j_track][state])
            last_state = flags[i_track][j_track][state]["state"]

            if state == 2:
                # deletion case
                aa_alignment = aa_emit + aa_alignment
                ref_alignment = hmm_model.reference[j_track] + ref_alignment
                j_track -= 1
            else:
                alignment = sorted_nodes[i_track].base + alignment
                alignment = path[2].base + alignment
                alignment = path[1].base + alignment
                i_track = check_node_order[path[0].ID]
                if state == 0:
                    if last_state == 4:
                        pass_hmm = True
                        i_start = check_node_order[path[1].ID]
                        j_start = j_track
                    # match state, j-1 traceback, otherwise means insertion, no need -1 for j
                    aa_alignment = aa_emit + aa_alignment
                    ref_alignment = hmm_model.reference[j_track] + ref_alignment
                    j_track -= 1
                else:
                    aa_alignment = aa_emit.lower() + aa_alignment
                    ref_alignment = "." + ref_alignment



        state = last_state

    """
    while (i_track >= 0) and (j_track >= 0):
        # print dag_graph.node[i_track]

        aa_emit = flags[i_track][j_track][state]["aa"]
        state_alignment = str(state) + state_alignment
        path = flags[i_track][j_track][state]["path"]

        score_track.append(scores[i_track][j_track][state])
        last_state = flags[i_track][j_track][state]["state"]
        if state == 2:
            # deletion case
            aa_alignment = aa_emit + aa_alignment
            ref_alignment = hmm_model.reference[j_track] + ref_alignment
            j_track -= 1

        else:
            alignment = sorted_nodes[i_track].base + alignment
            if len(path) == 3:
                alignment = path[2].base + alignment
                alignment = path[1].base + alignment
                i_track = check_node_order[path[0].ID]
                if state == 0:
                    # match state, j-1 traceback, otherwise means insertion, no need -1 for j
                    aa_alignment = aa_emit + aa_alignment
                    ref_alignment = hmm_model.reference[j_track] + ref_alignment
                    j_track -= 1
                else:
                    aa_alignment = aa_emit.lower() + aa_alignment
                    ref_alignment = "." + ref_alignment

            else:
                # only two node in path, means found beginning
                aa_alignment = aa_emit + aa_alignment
                ref_alignment = hmm_model.reference[j_track] + ref_alignment
                alignment = path[1].base + alignment
                alignment = path[0].base + alignment
                i_start = check_node_order[path[0].ID]
                j_start = j_track
                break

        state = last_state
        """
    traceback = {"alignment": alignment, "state": state_alignment, "score": score, "hmm_start": j_start,
                 "hmm_end": j_end, "seq_start": i_start, "seq_end": i_end, "aa": aa_alignment,
                 "score_track": score_track, "ref": ref_alignment}

    return traceback


def alignment_score(hmm, alignment):
    """

    :param hmm: hmm model
    :param alignment: protein sequence alignment with hmm
    :return:
    """
    score = 0
    m_to_m = 0
    i_to_m = 3
    m_to_i = 1
    m_to_d = 2
    i_to_i = 4
    d_to_m = 5
    d_to_d = 6
    seq = alignment["aa"]
    state = alignment["state"]  # 0 for match 1 for insetion 2 for deletion
    start = alignment["hmm_start"]  # the first amino acid should be indexed as 1
    # end = alignment["hmm_end"]
    hmm_pos = start

    for i in range(len(state)):

        if i == 0:
            # begin to match
            aa_index = aa.aminoacidindex(seq[i])
            score += hmm.match_emission[hmm_pos][aa_index] - hmm.match_emission[0][aa_index]
        else:
            if int(state[i]) == 0:
                # match state
                hmm_pos += 1
                aa_index = aa.aminoacidindex(seq[i])
                if int(state[i - 1]) == 0:
                    score += hmm.transition[hmm_pos - 1][m_to_m] + hmm.match_emission[hmm_pos][aa_index] - \
                             hmm.match_emission[0][aa_index]

                elif int(state[i - 1]) == 1:
                    score += hmm.transition[hmm_pos - 1][i_to_m] + hmm.match_emission[hmm_pos][aa_index] - \
                             hmm.match_emission[0][aa_index]

                elif int(state[i - 1]) == 2:
                    score += hmm.transition[hmm_pos - 1][d_to_m] + hmm.match_emission[hmm_pos][aa_index] - \
                             hmm.match_emission[0][aa_index]

            elif int(state[i]) == 1:
                # insertion
                aa_index = aa.aminoacidindex(seq[i])
                if int(state[i - 1]) == 0:
                    score += hmm.transition[hmm_pos][m_to_i] + hmm.insert_emission[hmm_pos][aa_index] - \
                             hmm.match_emission[0][aa_index]

                elif int(state[i - 1]) == 1:
                    score += hmm.transition[hmm_pos][i_to_i] + hmm.insert_emission[hmm_pos][aa_index] - \
                             hmm.match_emission[0][aa_index]

            elif int(state[i]) == 2:
                hmm_pos += 1
                if int(state[i - 1]) == 0:
                    score += hmm.transition[hmm_pos - 1][m_to_d]

                elif int(state[i - 1]) == 2:
                    score += hmm.transition[hmm_pos - 1][d_to_d]

    return score

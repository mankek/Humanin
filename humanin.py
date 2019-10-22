import pandas as pd


def add_aa(dict_in, set_in):
    for place in dict_in.keys():
        for i in set_in:
            if i not in dict_in[place].keys():
                dict_in[place][i] = 0
    return dict_in


def sort_positions(text_in, text_dict, set_in):
    for place, acid in enumerate(text_in):
        text_dict[place].setdefault(acid, 0)
        text_dict[place][acid] += 1
        set_in.add(acid)
    return text_dict, set_in


def read_proteins():
    aa_set = set()
    pos_dict = dict()
    for i in range(0, 24):
        pos_dict[i] = dict()
    with open("humanin.txt", "r") as file_in:
        for i in file_in.readlines():
            pos_dict, aa_set = sort_positions(i.rstrip(), pos_dict, aa_set)
    pos_dict = add_aa(pos_dict, aa_set)
    return pos_dict, aa_set


pos_dict_out, aa_set_out = read_proteins()
# print(pos_dict_out)


def get_counts(pos_dict, aa_set):
    counts_dict = dict()
    for aa in aa_set:
        counts_dict.setdefault(aa, [])
    for place in pos_dict.keys():
        for aa in pos_dict[place].keys():
            counts_dict[aa].append(pos_dict[place][aa])
    counts_df = pd.DataFrame(counts_dict)
    probs_df = counts_df.applymap(lambda x: x/16 if x != 0 else x)
    consen_dict = {"Sequence": probs_df.idxmax(1), "Probabilities": probs_df.max(1)}
    consen_df = pd.DataFrame(consen_dict)
    consensus = "".join(counts_df.idxmax(1).values)
    return probs_df, consen_df, consensus


probability_df, consensus_df, consensus_text = get_counts(pos_dict_out, aa_set_out)


import pandas as pd
from subprocess import Popen, PIPE
import os
from optparse import OptionParser


def get_sequence(seq_file_path):
    sequence = ""
    with open(seq_file_path, "r") as file_in:
        for line in file_in.readlines()[1:]:
            sequence = sequence + line.rstrip()
    return sequence


def get_translation(seq_in):
    my_env = os.environ.copy()
    my_env["PATH"] = "/usr/sbin:/sbin:" + my_env["PATH"]
    c = Popen(["python", "emboss_transeq.py", "--email", "mankek@alumni.msoe.edu", "--sequence", seq_in, "--frame", "6", "--codontable", "0"], env=my_env, shell=True, stdout=PIPE, stderr=PIPE)
    proc_out, proc_errs = c.communicate()
    job_id = proc_errs.decode().rstrip().split(" ")[-1]
    return job_id
    # print(proc_errs)
    # trans_file = proc_out.decode().split("\r")[1].split(" ")[-1]
    # print(trans_file)
    # return trans_file


def read_translation(trans_file):
    seq_list = dict()
    frame = "none"
    seq = ""
    with open(trans_file + ".out.txt", "r") as file_in:
        for i in file_in.readlines():
            line = i.rstrip()
            if line.startswith(">"):
                seq_list[frame] = seq
                frame = line.split("_")[-1]
                seq = ""
            else:
                seq = seq + line
    seq_list[frame] = seq
    seq_list.pop("none")
    os.remove(trans_file + ".out.txt")
    os.remove(trans_file + ".sequence.txt")
    return seq_list


def add_aa(dict_in, set_in):
    for place in dict_in.keys():
        for i in set_in:
            if i not in dict_in[place].keys():
                dict_in[place][i] = 0
    return dict_in


def count_positions(text_in, text_dict, set_in):
    for place, acid in enumerate(text_in):
        text_dict[place].setdefault(acid, 0)
        text_dict[place][acid] += 1
        set_in.add(acid)
    return text_dict, set_in


def read_proteins(seq_file):
    aa_set = set()
    pos_dict = dict()
    line_count = 0
    for i in range(0, 25):  # 25: the number of amino acids in the peptide (plus stop codon)
        pos_dict[i] = dict()
    with open(seq_file, "r") as file_in:
        for i in file_in.readlines():
            line_count += 1
            pos_dict, aa_set = count_positions(i.rstrip(), pos_dict, aa_set)
    pos_dict = add_aa(pos_dict, aa_set)
    return pos_dict, aa_set, line_count

# print(pos_dict_out)


def get_counts(pos_dict, aa_set, line_count):
    counts_dict = dict()
    for aa in aa_set:
        counts_dict.setdefault(aa, [])
    for place in pos_dict.keys():
        for aa in pos_dict[place].keys():
            counts_dict[aa].append(pos_dict[place][aa] + 1)
    counts_df = pd.DataFrame(counts_dict)
    probs_df = counts_df.applymap(lambda x: (x/(line_count * 2)) if x != 0 else x)  # 32: twice the number of sequences considered to find a consensus
    consen_dict = {"Sequence": probs_df.idxmax(1), "Probabilities": probs_df.max(1)}
    consen_df = pd.DataFrame(consen_dict)
    consensus = "".join(counts_df.idxmax(1).values)
    return probs_df, consen_df, consensus


def find_prob(prob_df, text_in, line_count):
    aa_list = list(prob_df.columns)
    text_probs = 1
    for index, acid in enumerate(text_in):
        if acid in aa_list:
            prob = prob_df.at[index, acid]
            text_probs = text_probs * prob
        else:
            prob = 1/(line_count * 2)
            text_probs = text_probs * prob
    return text_probs


def find_humanin(sequence_in, prob_df, line_count):
    prob = 0
    likely_hum = ""
    sequence_in.replace("-", "*")
    for i in range(0, len(sequence_in) - 25):
        subseq = sequence_in[i:(i + 25)]
        sub_prob = find_prob(prob_df, subseq, line_count)
        if sub_prob > prob:
            prob = sub_prob
            likely_hum = subseq
    return prob, likely_hum


def find_dna(sequence_in, frame, translation, motif):
    # print(translation)
    motif_ind = translation.find(motif)
    # print(motif_ind)
    if int(frame) > 0:
        seq_start = int(frame) - 1
        sequence_framed = sequence_in[seq_start:]
    else:
        seq_start = int(frame)
        sequence_framed = sequence_in[seq_start::-1]
    # print(sequence_in)
    # print(sequence_framed)
    dna_seq = sequence_framed[(motif_ind*3):((motif_ind*3) + (3*25))]
    return dna_seq


def yield_results(frame_opt, file_in, output_file, prob_opt, prob_df, con_prob, line_count):
    seq = get_sequence(file_in)
    trans = get_translation(seq)
    seqs = read_translation(trans)
    if frame_opt == "All":
        for key, value in seqs.items():
            prob, motif = find_humanin(value, prob_df, line_count)
            if prob_opt:
                rel_prob = str((prob/con_prob) * 100) + "% "
            else:
                rel_prob = ""
            output_file.write(file_in + " - Frame: " + key + ": " + rel_prob + motif + "\n")
            output_file.write("DNA Sequence: " + find_dna(seq, key, seqs[key], motif) + "\n\n")
    elif frame_opt == "Best":
        max_key = ""
        max_prob = 0
        max_motif = ""
        for key, value in seqs.items():
            prob, motif = find_humanin(value, prob_df, line_count)
            if prob > max_prob:
                max_prob = prob
                max_key = key
                max_motif = motif
        if prob_opt:
            rel_prob = str((max_prob/con_prob) * 100) + "% "
        else:
            rel_prob = ""
        output_file.write(file_in.split("\\")[-1] + " - Frame " + max_key + ": " + rel_prob + max_motif + "\n")
        output_file.write("DNA Sequence: " + find_dna(seq, max_key, seqs[max_key], max_motif) + "\n\n")
    else:
        for key, value in seqs.items():
            prob, motif = find_humanin(value, prob_df, line_count)
            if prob_opt:
                rel_prob = str((prob / con_prob) * 100) + "% "
            else:
                rel_prob = ""
            if key == str(frame_opt):
                output_file.write(file_in + " - Frame: " + key + ": " + rel_prob + motif + "\n")
                output_file.write("DNA Sequence: " + find_dna(seq, key, seqs[key], motif) + "\n\n")


def main():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Read data from FASTA file; can also be a directory path if "
                                                            "the -d option is used")
    parser.add_option("-o", "--output", dest="output", help="Prefix of the output file")
    parser.add_option("-m", "--matrix", dest="matrix", default="humanin.txt", help="File containing the motif sequences "
                                                                                   "used to create the probability "
                                                                                   "matrix; each line should hold a "
                                                                                   "separate sequence; all sequences "
                                                                                   "should be the same length. Default file is humanin.txt")
    parser.add_option("-d", "--dir", action="store_true", dest="directory", default=False, help="Read data from FASTA "
                                                                                                "files in this "
                                                                                                "directory")
    parser.add_option("-f", "--frames", dest="frames", default="All", help="Determines for which frames the most "
                                                                           "probable humanin sequence is returned; "
                                                                           "options: 1-6, All, or Best")
    parser.add_option("-p", "--probs", action="store_true", dest="probabilities", default=False, help="Includes relative"
                                                                                                      " probability of "
                                                                                                      "the motif (as "
                                                                                                      "compared to "
                                                                                                      "probability of "
                                                                                                      "consensus "
                                                                                                      "sequence)")

    (options, args) = parser.parse_args()
    pos_dict_out, aa_set_out, line_count = read_proteins(options.matrix)
    probability_df, consensus_df, consensus_text = get_counts(pos_dict_out, aa_set_out, line_count)
    consensus_prob = find_prob(probability_df, consensus_text, line_count)
    # if len(args) < 1:
    #     parser.error("Input option required")
    if options.input:
        if options.directory:
            file_directory = options.input
            print("Directory: " + file_directory)
            with open(options.output + "_humanin_results.txt", "w") as file_out:
                for file in os.listdir(file_directory):
                    file_in = os.path.join(file_directory, file)
                    print("File: " + file_in)
                    yield_results(options.frames, file_in, file_out, options.probabilities, probability_df, consensus_prob, line_count)
        else:
            if os.path.isdir(options.input):
                print("Please use the -d option with a path to a directory")
            else:
                with open(options.output + "_humanin_results.txt", "w") as file_out:
                    file = options.input
                    print("File: " + file)
                    yield_results(options.frames, file, file_out, options.probabilities, probability_df, consensus_prob, line_count)
    else:
        print("No input file specified.\nUse 'python humanin.py -h' for the help message.")


if __name__ == "__main__":
    main()


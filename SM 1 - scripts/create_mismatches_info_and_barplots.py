import pickle
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
import sys

gene_name = sys.argv[1]
fasta_file = sys.argv[2]  # "dd_Smed_v6.fasta"
sample = sys.argv[3]  # binary
log_scale = sys.argv[4]  # binary
file_prefix = ''


def set_file_prefix(gene_name):
    num = gene_name.split('_')[3]
    global file_prefix
    file_prefix = gene_name + "/" + num + "_"


def create_mismatches_info_file(gene_name, fasta_file, sample):
    print("\nStarting creating mismatches info file")
    mm_col_a = "good quality mismatch amount"
    mm_col_p = "good quality mismatch %"
    final_df = pd.DataFrame(columns=['position', 'reference', 'coverage',  mm_col_a, mm_col_p, 'A', 'T', 'C', 'G'])

    if sample:
        print("Using sample file to calculate depth")
        depth_df = pd.read_csv(file_prefix + "sample_depth.tsv", sep='\t', names=['Gene', 'position', 'coverage'])
        d_file = open(file_prefix + "sample_pos_good_quality_mm_dict.pkl", "rb")
        output_file_name = file_prefix + "sample_good_quality_mismatch_info.csv"
    else:
        depth_df = pd.read_csv(file_prefix + "depth.tsv", sep='\t', names=['Gene', 'position', 'coverage'])
        d_file = open(file_prefix + "pos_good_quality_mm_dict.pkl", "rb")
        output_file_name = file_prefix + "good_quality_mismatch_info.csv"

    final_df['position'] = depth_df['position']
    final_df['coverage'] = depth_df['coverage']

    record_dict = SeqIO.index(fasta_file, "fasta")
    record = record_dict[gene_name]
    reference_seq = [nuc for nuc in record.seq]
    final_df['reference'] = reference_seq

    good_dict = pickle.load(d_file)
    mm_count = [len(tup[1]) for tup in list(good_dict.items())]
    final_df[mm_col_a] = mm_count
    final_df[mm_col_p] = final_df[mm_col_a] / final_df['coverage'] * 100

    nucs = ['A', 'T', 'C', 'G']
    pos_nuc_cnt_dict = {key: [] for key in nucs}
    for dict_tup in list(good_dict.items()):
        nuc_list = [tup[0] for tup in dict_tup[1]]  # all nucleotides at the position
        for nuc in nucs:
            pos_nuc_cnt_dict[nuc].append(nuc_list.count(nuc)) # count of each nucleotide

    for nuc in nucs:
        final_df[nuc] = pos_nuc_cnt_dict[nuc]
    for nuc in nucs:
        final_df[nuc + " /mismatch amount (%)"] = final_df[nuc] / final_df[mm_col_a] * 100
    for nuc in nucs:
        final_df[nuc + " /total reads (%)"] = final_df[nuc] / final_df['coverage'] * 100
        # if nuc==reference then the value will be 0, because we are looking only at mismatches data

    final_df = final_df.fillna(0) # replace all possible Nan values with 0 (when coverage or mismatch amount is 0)
    final_df.to_csv(output_file_name, index=False)
    d_file.close()
    print("Done creating mismatches info file")
    print("Results were saved to '{}'".format(output_file_name))


def plot_barplot_divided_by_nucleodites(sample, log_scale):
    print("\nCreating barplot of good quality mismatches per position, divided by nucleotides")
    if sample:
        print("Using sample file for ploting")
        input_file_name = file_prefix + "sample_good_quality_mismatch_info.csv"
        plot_file_name = file_prefix + 'sample_barplot_mm_by_nucleotides.pkl'
        plot_title = "Gene {} - Sample - count of good quality mismatches per position, divided by nucleotides".format(gene_name)
    else:
        input_file_name = file_prefix + "good_quality_mismatch_info.csv"
        plot_file_name = file_prefix + 'barplot_mm_by_nucleotides.pkl'
        plot_title = "Gene {} - count of good quality mismatches per position, divided by nucleotides".format(gene_name)
    mm_df = pd.read_csv(input_file_name)
    gene_len = len(mm_df)
    x = range(1, gene_len+1) # all positions in the gene sequence
    bar_width = 0.65
    fig, ax = plt.subplots()
    ax.bar(x, mm_df['A'], bar_width, label='A', color='blue')
    ax.bar(x, mm_df['T'], bar_width, bottom=mm_df['A'], label='T', color='red')
    ax.bar(x, mm_df['C'], bar_width, bottom=mm_df['T'], label='C', color='lime')
    ax.bar(x, mm_df['G'], bar_width, bottom=mm_df['C'], label='G', color='violet')
    if log_scale:
        plt.yscale('log', base=10)
        plt.ylabel("Amount of reads with mismatches (log10 scale)")
    else:
        plt.ylabel("Amount of reads with mismatches")
    plt.xticks(range(1, len(x)+1, 50))
    plt.xlabel("Position in reference sequence")
    plt.title(plot_title)

    # add to the barplot the coverage for each position:
    plt.plot(mm_df['position'], mm_df['coverage'], lw=0.7, color='gray', label='Coverage') # pos x in the reference is pos x-1
    # because the positions in the reference sequence starts from 1, but in python sequences start from 0.

    # reorder the legend:
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])

    # save the plot as a pkl file:
    plot_file = open(plot_file_name, 'wb')
    pickle.dump(ax, plot_file)
    plot_file.close()
    plt.clf()
    print("The plot was saved as '{}'. It needs to be loaded with pickle for display".format(plot_file_name))


def plot_barplot_precentage(sample):
    print("\nCreating barplot of good quality mismatch precentage per position")
    if sample:
        print("Using sample file for ploting")
        input_file_name = file_prefix + "sample_good_quality_mismatch_info.csv"
        plot_file_name = file_prefix +'sample_barplot_mm_precentage.pkl'
        plot_title = "Gene {} - Sample - Precentage of good quality mismatches per position".format(gene_name)
    else:
        input_file_name = file_prefix + "good_quality_mismatch_info.csv"
        plot_file_name = file_prefix + 'barplot_mm_precentage.pkl'
        plot_title = "Gene {} - Precentage of good quality mismatches per position".format(gene_name)
    mm_df = pd.read_csv(input_file_name)

    # plotting a bar plot
    ax = plt.bar(mm_df['position'], mm_df['good quality mismatch %'], width=0.8, color=['dodgerblue', 'magenta'])
    plt.xticks(range(1, len(mm_df)+1, 50))
    plt.xlabel("Position in reference sequence")
    plt.ylabel("Precentage of reads with good quality mismatches")
    plt.title(plot_title)
    plt.grid(linewidth=0.5)

    # save the plot as a pkl file:
    plot_file = open(plot_file_name, 'wb')
    pickle.dump(ax, plot_file)
    plot_file.close()
    print("The plot was saved as '{}'. It needs to be loaded with pickle for display".format(plot_file_name))


def do_work_for_gene(gene_name, sample, log_scale):
    set_file_prefix(gene_name)
    create_mismatches_info_file(gene_name, fasta_file, sample)
    plot_barplot_divided_by_nucleodites(sample, log_scale)
    plot_barplot_precentage(sample)
    print("\nDone working on gene '{}'. You may use show_plot.py <plot_file.pkl> to see the saved plots".format(gene_name))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("There should be exactly 4 arguments, but {} were passed. Exiting.".format(len(sys.argv)))
        exit()
    do_work_for_gene(gene_name, sample, log_scale)

